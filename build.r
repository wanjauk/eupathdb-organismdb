###############################################################################
#
# TriTrypDB OrganismDB package generation
#
# This script uses resources from TriTrypDB to generate a Bioconductor Organism
# annotation package.
# 
###############################################################################
library(yaml)
library(tools)
library(rtracklayer)
library(AnnotationForge)

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

# Generate package
makeOrgPackage(
    gene_info  = gene_info,
    chromosome = chr_info,
    go         = go_table,
    version    = settings$tritrypdb_version,
    author     = settings$author,
    maintainer = settings$maintainer,
    outputDir  = settings$output_dir,
    tax_id     = settings$tax_id,
    genus      = settings$genus,
    species    = settings$species,
    goTable    = "go"
)

#'
#' TriTrypDB gene information table GO term parser
#'
#' @author Keith Hughitt
#'
#' @param filepath Location of TriTrypDB gene information table.
#' @param verbose  Whether or not to enable verbose output.
#' @return Returns a dataframe where each line includes a gene/GO terms pair
#'         along with some addition information about the GO term. Note that
#'         because each gene may have multiple GO terms, a single gene ID may
#'         appear on multiple lines.
#'
parse_go_terms = function (filepath) {
    if (file_ext(filepath) == 'gz') {
        fp = gzfile(filepath, open='rb')
    } else {
        fp = file(filepath, open='r')
    }

    # Create empty vector to store dataframe rows
    N = 1e5
    gene_ids = c()
    go_rows = data.frame(GO=rep("", N),
                         ONTOLOGY=rep("", N), GO_TERM_NAME=rep("", N),
                         SOURCE=rep("", N), EVIDENCE=rep("", N),
                         stringsAsFactors=FALSE)

    # Counter to keep track of row number
    i = j = 1

    # Iterate through lines in file
    while (length(x <- readLines(fp, n=1, warn=FALSE)) > 0) {
        # Gene ID
        if(grepl("^Gene ID", x)) {
            gene_id = .get_value(x)
            i = i + 1
        }

        # Gene Ontology terms
        else if (grepl("^GO:", x)) {
            gene_ids[j] = gene_id
            go_rows[j,] = c(head(unlist(strsplit(x, '\t')), 5))
            j = j + 1
        }
    }

    # get rid of unallocated rows
    go_rows = go_rows[1:j-1,]

    # drop unneeded columns
    go_rows = go_rows[,c('GO', 'EVIDENCE')]

    # add gene id column
    go_rows = cbind(GID=gene_ids, go_rows)

    # close file pointer
    close(fp)

    # TODO: Determine source of non-unique rows in the dataframe
    return(unique(go_rows))
}

#
# Parses a key: value string and returns the value
#
.get_value = function(x) {
    return(gsub(" ","", tail(unlist(strsplit(x, ': ')), n=1), fixed=TRUE))
}

# TriTrypDB download URL
# http://tritrypdb.org/common/downloads/Current_Release/LmajorFriedlin/
# http://tritrypdb.org/common/downloads/Current_Release/LmajorFriedlin/gff/data/TriTrypDB-8.1_LmajorFriedlin.gff
# http://tritrypdb.org/common/downloads/Current_Release/LmajorFriedlin/txt/TriTrypDB-8.1_LmajorFriedlinGene.txt

