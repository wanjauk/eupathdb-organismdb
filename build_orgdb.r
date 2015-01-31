###############################################################################
#
# TriTrypDB OrgDB package generation
#
# This script uses resources from TriTrypDB to generate a Bioconductor Organism
# annotation package.
# 
###############################################################################
library(yaml)
library(tools)
library(rtracklayer)
library(AnnotationForge)
source('helper.r')

options(stringsAsFactors=FALSE)

#
# MAIN
#
# Load settings
settings = yaml.load_file("config.yaml")

# Create build directory
if (!file.exists(settings$build_dir)) {
    dir.create(settings$build_dir)
}
build_basename = file.path(settings$build_dir,
                            sub('.gff', '', basename(settings$gff)))

# Parse GFF
gff = import.gff3(settings$gff)
genes = gff[gff$type == 'gene']

#
# 1. Gene name and description
#
gene_file = sprintf("%s_gene_info.txt", build_basename)

if (file.exists(gene_file)) {
    gene_info = read.delim(gene_file)
} else {
    print("Parsing gene information...")
    gene_info = as.data.frame(elementMetadata(genes[,c('ID', 'description')]))

    # Convert form-encoded description string to human-readable
    gene_info$description = gsub("\\+", " ", gene_info$description)
    colnames(gene_info) <- c("GID", "GENENAME")

    write.table(gene_info, gene_file, sep='\t', quote=FALSE, row.names=FALSE)
}

#
# 2. Chromosome information
#
chr_file = sprintf("%s_chr_info.txt", build_basename)

if (file.exists(chr_file)) {
    chr_info = read.delim(chr_file)
} else {
    print("Parsing chromosome information...")
    chr_info = data.frame(
        'GID' = genes$ID,
        'CHROMOSOME' = as.character(seqnames(genes))
    )
    write.table(chr_info, chr_file, sep='\t', quote=FALSE, row.names=FALSE)
}

#
# 3. Gene ontology information
#
go_file = sprintf("%s_go_table.txt", build_basename)

if (file.exists(go_file)) {
    go_table = read.delim(go_file)
} else {
    print("Parsing GO annotations...")
    go_table = parse_go_terms(settings$txt)
    write.table(go_table, go_file, sep='\t', quote=FALSE, row.names=FALSE)
}

#
# 4. KEGG information
#
kegg_mapping_file  = sprintf("%s_kegg_mapping.txt", build_basename)
kegg_pathways_file = sprintf("%s_kegg_pathways.txt", build_basename)

if (!file.exists(kegg_mapping_file)) {
    library(KEGGREST)

    # KEGG Organism abbreviation (e.g. "lma")
    org_abbreviation = paste0(tolower(substring(settings$genus, 1, 1)), 
                              substring(settings$species, 1, 2))

    # TODO: Generalize if possible

    # L. major
    if (org_abbreviation == 'lma') {
        #
        # Convert KEGG identifiers to TriTrypDB identifiers
        #
        # Note that this currently skips a few entries with a different
        # format, e.g. "md:lma_M00359", and "bsid:85066"
        #
        convert_kegg_gene_ids = function(kegg_ids) {
            result = c()
            for (kegg_id in kegg_ids) {
                if (substring(kegg_id, 1, 9) == 'lma:LMJF_') {
                    # lma:LMJF_11_0100
                    result = append(result, 
                                    gsub('LMJF', 'LmjF', 
                                    gsub("_", "\\.", substring(kegg_id, 5))))
                } else if (substring(kegg_id, 1, 8) == 'lma:LMJF') {
                    # lma:LMJF10_TRNALYS_01
                    parts = unlist(strsplit(kegg_id, "_"))
                    result = append(result,
                                    sprintf("LmjF.%s.%s.%s",
                                    substring(kegg_id, 9, 10),
                                    parts[2], parts[3]))
                } else {
                    print(sprintf("Skipping KEGG id: %s", kegg_id))
                    result = append(result, NA)
                }
            }
            return(result)
        }
    } else if (org_abbreviation == 'tcr') {
        # Load GeneAlias file and convert entry in KEGG results
        # to newer GeneDB/TriTrypDB identifiers.
        fp = file(settings$aliases)
        rows = strsplit(readLines(fp), "\t")
        close(fp)

        kegg_id_mapping = list()

        for (row in rows) {
            # get first and third columns in the alias file
            kegg_id_mapping[row[3]] = row[1]
        }

        # Example: "tcr:509463.30" -> ""
        convert_kegg_gene_ids = function(kegg_ids) {
            kegg_to_genedb(kegg_ids, kegg_id_mapping)   
        }
    }

    # data frame to store kegg gene mapping and pathway information
    kegg_mapping = data.frame()
    kegg_pathways = data.frame()

    pathways = unique(keggLink("pathway", org_abbreviation))

    # Iterate over pathways and query genes for each one
    for (pathway in pathways) {
        message(sprintf("Processing genes for KEGG pathway %s", pathway))

        # Get pathway info
        meta = keggGet(pathway)[[1]]
        pathway_desc  = ifelse(is.null(meta$DESCRIPTION), '', meta$DESCRIPTION)
        pathway_class = ifelse(is.null(meta$CLASS), '', meta$CLASS)
        kegg_pathways = rbind(kegg_pathways, 
                              data.frame(pathway=pathway,
                                         name=meta$PATHWAY_MAP,
                                         class=pathway_class,
                                         description=pathway_desc))
        
        # Get genes in pathway
        result = keggLink(pathway) 
        kegg_ids = result[,2]
        gene_ids = convert_kegg_gene_ids(kegg_ids)
        kegg_mapping = unique(rbind(kegg_mapping,
                             data.frame(GID=gene_ids, pathway=pathway)))
    }
    # Save KEGG mapping
    write.csv(kegg_mapping, file=kegg_mapping_file, quote=FALSE,
              row.names=FALSE)
    write.table(kegg_pathways, file=kegg_pathways_file, quote=FALSE,
              row.names=FALSE, sep='\t')
} else {
    # Otherwise load saved version
    kegg_mapping = read.csv(kegg_mapping_file)
    kegg_pathways = read.delim(kegg_pathways_file)
}

# Drop columns with unrecognized identifiers
kegg_mapping = kegg_mapping[complete.cases(kegg_mapping),]

# Combined KEGG table
kegg_table = merge(kegg_mapping, kegg_pathways, by='pathway')
colnames(kegg_table) = c("KEGG_PATH", "GID", "KEGG_NAME", "KEGG_CLASS",
                         "KEGG_DESCRIPTION")

# reorder so GID comes first
kegg_table = kegg_table[,c(2, 1, 3, 4, 5)]

# Generate package
makeOrgPackage(
    gene_info  = gene_info,
    chromosome = chr_info,
    go         = go_table,
    kegg       = kegg_table,
    version    = settings$tritrypdb_version,
    author     = settings$author,
    maintainer = settings$maintainer,
    outputDir  = settings$output_dir,
    tax_id     = settings$tax_id,
    genus      = settings$genus,
    species    = settings$species,
    goTable    = "go"
)

