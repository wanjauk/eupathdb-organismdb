#!/usr/bin/env Rscript
###############################################################################
##
## EuPathDB OrgDB package generation
##
## This script uses resources from EuPathDB to generate a Bioconductor Orgdb
## annotation package.
##
###############################################################################

suppressMessages(library(yaml))
suppressMessages(library(tools))
suppressMessages(library(rtracklayer))
suppressMessages(library(AnnotationForge))
suppressMessages(library(RSQLite))
source('helper.R')

options(stringsAsFactors=FALSE)

# default settings file
config_file <- "config.yaml"

# if alternate config file specified, use that instead
args <- commandArgs(TRUE)

if (length(args) > 0) {
    config_file <- args[1]
} else {
    message("Defaulting to the configuration in 'config.yaml'.")
}

#
# MAIN
#
# Load settings
settings <- yaml.load_file(config_file)

rversion <- paste0(R.Version()$major,  '.', R.Version()$minor)

build_dir <- file.path(settings$build_dir, rversion)
output_dir <- file.path(settings$output_dir, rversion)

# Create build and output directories
if (!file.exists(build_dir)) {
    dir.create(build_dir, recursive=TRUE)
}
if (!file.exists(output_dir)) {
    dir.create(output_dir, recursive=TRUE)
}
build_basename <- file.path(build_dir,
                           sub('.gff', '', basename(settings$gff)))

# Parse GFF
gff <- import.gff3(settings$gff)
genes <- gff[gff$type == 'gene']

#
# 1. Gene name and description
#
gene_file <- sprintf("%s_gene_info.txt", build_basename)

# Check to see if caches version of metadata exists and load it if it does
if (file.exists(gene_file)) {
    message(paste0("Reading pre-existing gene file: ", gene_file))
    gene_info <- read.delim(gene_file)
} else {
    # Otherwise, parse GFF to get basic gene featurs
    message(paste0("Parsing gene information from: ", settings$gff))
    gene_info <- as.data.frame(elementMetadata(genes))

    # GFF field names changed from TriTrypDB 28 -> 29
    if ('description' %in% colnames(gene_info)) {
        #
        # TriTrypDB 28 and earlier
        #

        # Normalize columns names
        colnames(gene_info) <- toupper(colnames(gene_info))
        colnames(gene_info)[colnames(gene_info) == "ID"] <- "GID"

        ## Move gid to the front of the line.
        gid_index <- grep("GID", colnames(gene_info))
        gene_info <- gene_info[, c(gid_index, (1:ncol(gene_info))[-gid_index])]
        colnames(gene_info) <- paste0("GENE", colnames(gene_info))
        colnames(gene_info)[1] <- "GID"
    }

    # Convert form-encoded description string to human-readable
    gene_info$GENEDESCRIPTION <- gsub("\\+", " ", gene_info$GENEDESCRIPTION)

    num_rows <- nrow(gene_info)

    # Remove any newlines present in GENEALIAS field;
    # as.character inserts newlines for objects with >500 characters.
    # https://stat.ethz.ch/R-manual/R-devel/library/base/html/character.html
    if ("GENEALIAS" %in% colnames(gene_info)) {
        gene_info[["GENEALIAS"]] <- as.character(gene_info[["GENEALIAS"]])
        gene_info[["GENEALIAS"]] <- gsub('\n', '', gene_info[["GENEALIAS"]])
    }

    ## Get rid of character(0) and NA entries
    is.empty <- function(stuff) {
        (length(stuff) == 0) && (typeof(stuff) == "character")
    }
    is.empty.col <- function(x) {
        y <- if (length(x)) {
                 do.call("cbind", lapply(x, "is.empty"))
             } else {
                 matrix(FALSE, length(row.names(x)), 0)
             }
        if (.row_names_info(x) > 0L) {
            rownames(y) <- row.names(x)
        }
        y
    }

    # Check for empty columns and remove them
    for (col in colnames(gene_info)) {
        tmp_col <- gene_info[[col]]
        empty_index <- is.empty.col(tmp_col)
        tmp_col[empty_index] <- NA
        gene_info[[col]] <- tmp_col
        if (sum(!is.na(gene_info[[col]])) == num_rows) {
            gene_info[, !(colnames(gene_info) %in% col)]
        }
        gene_info[[col]] <- as.character(gene_info[[col]])
    }

    # Replace remaining NA entries with the "null" string
    gene_info[is.na(gene_info)] <- "null"
    gene_info[gene_info == "NA"] <- "null"

    # Save results
    write.table(gene_info, gene_file, sep='\t', quote=FALSE, row.names=FALSE)
}

#
# 2. Chromosome information
#
chr_file <- sprintf("%s_chr_info.txt", build_basename)

if (file.exists(chr_file)) {
    chr_info <- read.delim(chr_file)
} else {
    message(paste0("Parsing chromosome information from ", chr_file))
    chr_info <- data.frame(
        'GID' = genes$ID,
        'CHR' = as.character(seqnames(genes))
    )
    write.table(chr_info, chr_file, sep='\t', quote=FALSE, row.names=FALSE)
}

#
# 3. Gene type information
#
gene_type_file <- sprintf("%s_gene_type.txt", build_basename)

if (file.exists(gene_type_file)) {
    gene_types <- read.delim(gene_type_file)
} else {
    message("Retrieving gene types from TriTrypDB web API")
    gene_types <- get_gene_types('tritrypdb', settings$description)
    write.table(gene_types, gene_type_file, sep='\t', quote=FALSE, row.names=FALSE)
}

#
# 4. Gene ontology information
#
go_file <- sprintf("%s_go_table.txt", build_basename)

if (file.exists(go_file)) {
    go_table <- read.delim(go_file)
} else {
    library('GO.db')

    print("Parsing GO annotations...")

    # 2016/07/03 recent versions of TriTrypDB no longer include
    # automatically inferred (IEA) GO term annotations in the txt files;
    # swithing to web service calls to retrieve all annotations.
    #go_table <- parse_go_terms(settings$txt)

    organism <- settings$description
    go_table <- retrieve_go_terms('tritrypdb.org', organism)

    # Check to make sure result is valid
    if (nrow(go_table) == 0) {
        # If error occurs, check to make sure the organism name is in
        # settings$description is formatted the same as on EuPathDB
        # e.g. http://tritrypdb.org/tritrypdb/webservices/GeneQuestions/GenesByGoTerm.wadl
        stop("Failed to retrieve GO terms for the specified organism.")
    }

    # Map from non-primary IDs to primary GO ids;
    # non-primary IDs are filtered out by makeOrgPackage
    problem_rows <- go_table[!go_table$GO %in% keys(GO.db),]
    synonyms <- problem_rows$GO

    # Create a mapping data frame
    synonym_mapping <- data.frame()

    for (syn in synonyms) {
        if (!is.null(GOSYNONYM[[syn]])) {
            synonym_mapping <- rbind(synonym_mapping, c(syn, GOSYNONYM[[syn]]@GOID))
        }
    }

    # replace alternative GO term identifiers
    if (nrow(synonym_mapping) > 0) {
        colnames(synonym_mapping) <- c('synonym', 'primary')
        synonym_mapping <- unique(synonym_mapping)

        go_table$GO[!go_table$GO %in% keys(GO.db)] <- synonym_mapping$primary[match(synonyms, synonym_mapping$synonym)]
        go_table <- unique(go_table[complete.cases(go_table),])
    }

    write.table(go_table, go_file, sep='\t', quote=FALSE, row.names=FALSE)
}

#
# 5. KEGG information
#
kegg_mapping_file  <- sprintf("%s_kegg_mapping.txt", build_basename)
kegg_pathways_file <- sprintf("%s_kegg_pathways.txt", build_basename)

convert_kegg_gene_ids <- function(kegg_ids, kegg_id_mapping) {
    result <- c()
    for (kegg_id in kegg_ids) {
        if (substring(kegg_id, 1, 4) == 'tbr:') {
            # T. brucei
            result <- append(result,
                gsub('tbr:', '', kegg_id))
        } else if (substring(kegg_id, 1, 4) == 'tcr:') {
            # T. cruzi
            result <- append(result,
                gsub('tcr:', 'TcCLB.', kegg_id))
        } else if (substring(kegg_id, 1, 4) == 'tgo:') {
            # T. gondii
            result <- append(result, gsub('tgo:', '', gsub('_', '.', kegg_id)))
        } else if (substring(kegg_id, 1, 9) == 'lbz:LBRM_') {
            # L. braziliensis (lbz:LBRM_01_0080)
            result <- append(result, gsub('LBRM', 'LbrM',
                     gsub("_", "\\.", substring(kegg_id, 5))))
        } else if (substring(kegg_id, 1, 9) == 'lma:LMJF_') {
            # L. major (lma:LMJF_11_0100)
            result <- append(result,
                gsub('LMJF', 'LmjF',
                     gsub("_", "\\.", substring(kegg_id, 5))))
        } else if (substring(kegg_id, 1, 8) == 'lma:LMJF') {
            # L. major (lma:LMJF10_TRNALYS_01)
            parts <- unlist(strsplit(kegg_id, "_"))
            result <- append(result,
                sprintf("LmjF.%s.%s.%s",
                        substring(kegg_id, 9, 10),
                        parts[2], parts[3]))
        } else {
            print(sprintf("Skipping KEGG id: %s", kegg_id))
            result <- append(result, NA)
        }
    }
    return(result)
} ## End convert_kegg_gene_ids

if (!file.exists(kegg_mapping_file)) {
    library(KEGGREST)

    # KEGG Organism abbreviation (e.g. "lma")
    org_abbreviation <- paste0(tolower(substring(settings$genus, 1, 1)),
                              substring(settings$species, 1, 2))

    # Overides for cases where KEGG abbreviation differes from the above
    # pattern.

    # L. braziliensis
    if (org_abbreviation == 'lbr') {
        org_abbreviation <- 'lbz'
    }

    # For some species, it is necessary to map gene ids from KEGG to what is
    # currently used on TriTrypDB.
    #
    # TODO: Generalize if possible
    #
    if (org_abbreviation == 'tbr') {
        # Load GeneAlias file and convert entry in KEGG results
        # to newer GeneDB/TriTrypDB identifiers.
        fp <- file(settings$aliases)
        rows <- strsplit(readLines(fp), "\t")
        close(fp)

        kegg_id_mapping <- list()

        # example alias file entries
        #Tb927.10.2410  TRYP_x-70a06.p2kb545_720  Tb10.70.5290
        #Tb927.9.15520  Tb09.244.2520  Tb09.244.2520:mRNA
        #Tb927.8.5760   Tb08.26E13.490
        #Tb10.v4.0258   Tb10.1120
        #Tb927.11.7240  Tb11.02.5150  Tb11.02.5150:mRNA  Tb11.02.5150:pep
        for (row in rows) {
            # get first and third columns in the alias file
            old_ids <- row[2:length(row)]

            for (old_id in old_ids[grepl('Tb\\d+\\.\\w+\\.\\d+', old_ids)]) {
                kegg_id_mapping[old_id] <- row[1]
            }
        }

    } else if (org_abbreviation == 'lma') {
        # L. major
        #
        # Convert KEGG identifiers to TriTrypDB identifiers
        #
        # Note that this currently skips a few entries with a different
        # format, e.g. "md:lma_M00359", and "bsid:85066"
        #
    } else if (org_abbreviation == 'tcr') {
        fp <- file(settings$aliases)
        rows <- strsplit(readLines(fp), "\t")
        close(fp)

        kegg_id_mapping <- list()

        for (row in rows) {
            # get first and third columns in the alias file
            kegg_id_mapping[row[3]] <- row[1]
        }

        # Example: "tcr:509463.30" -> ""
        ##convert_kegg_gene_ids <- function(kegg_ids) {
        ##    kegg_to_genedb(kegg_ids, kegg_id_mapping)
        ##}
    }

    # data frame to store kegg gene mapping and pathway information
    kegg_mapping <- data.frame()
    kegg_pathways <- data.frame()

    pathways <- unique(keggLink("pathway", org_abbreviation))

    # Iterate over pathways and query genes for each one
    for (pathway in pathways) {
        message(sprintf("Processing genes for KEGG pathway %s", pathway))

        # Get pathway info
        meta <- keggGet(pathway)[[1]]
        pathway_desc  <- ifelse(is.null(meta$DESCRIPTION), '', meta$DESCRIPTION)
        pathway_class <- ifelse(is.null(meta$CLASS), '', meta$CLASS)
        kegg_pathways <- rbind(kegg_pathways,
                               data.frame(pathway=pathway,
                                         name=meta$PATHWAY_MAP,
                                         class=pathway_class,
                                         description=pathway_desc))

        # Get genes in pathway
        kegg_ids <- as.character(keggLink(org_abbreviation, pathway))
        gene_ids <- convert_kegg_gene_ids(kegg_ids)

        # Map old T. brucei gene names
        if (org_abbreviation == 'tbr') {
            old_gene_ids <- gene_ids
            gene_ids <- c()

            for (x in old_gene_ids) {
                if (x %in% names(kegg_id_mapping)) {
                    gene_ids <- append(gene_ids, kegg_id_mapping[[x]])
                } else {
                    gene_ids <- append(gene_ids, x)
                }
            }
        }

        if (!is.null(gene_ids)) {
            kegg_mapping <- unique(rbind(kegg_mapping,
                data.frame(GID=gene_ids, pathway=pathway)))
        }
    }
    # Save KEGG mapping
    write.csv(kegg_mapping, file=kegg_mapping_file, quote=FALSE,
              row.names=FALSE)
    write.table(kegg_pathways, file=kegg_pathways_file, quote=FALSE,
              row.names=FALSE, sep='\t')
} else {
    # Otherwise load saved version
    kegg_mapping <- read.csv(kegg_mapping_file)
    kegg_pathways <- read.delim(kegg_pathways_file)
}

# Drop columns with unrecognized identifiers
kegg_mapping <- kegg_mapping[complete.cases(kegg_mapping),]

# Combined KEGG table
kegg_table <- merge(kegg_mapping, kegg_pathways, by='pathway')
colnames(kegg_table) <- c("KEGG_PATH", "GID", "KEGG_NAME", "KEGG_CLASS",
                         "KEGG_DESCRIPTION")

# reorder so GID comes first
kegg_table <- kegg_table[,c(2, 1, 3, 4, 5)]

# R package versions must be of the form "x.y"
db_version <- paste(settings$db_version, '0', sep='.')

message("Building OrgDB")

# Compile list of arguments for makeOrgPackage call
orgdb_args <- list(
    gene_info  = gene_info,
    chromosome = chr_info,
    go         = go_table,
    type       = gene_types,
    version    = db_version,
    author     = settings$author,
    maintainer = settings$maintainer,
    outputDir  = output_dir,
    tax_id     = settings$tax_id,
    genus      = settings$genus,
    species    = settings$species,
    goTable    = "go"
)

# Add KEGG mapping if one exists
if (nrow(kegg_table) > 0) {
    orgdb_args[['kegg']] <- kegg_table
}
org_result <- do.call('makeOrgPackage', orgdb_args)

