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


# TriTrypDB download URL
# http://tritrypdb.org/common/downloads/Current_Release/LmajorFriedlin/
# http://tritrypdb.org/common/downloads/Current_Release/LmajorFriedlin/gff/data/TriTrypDB-8.1_LmajorFriedlin.gff
# http://tritrypdb.org/common/downloads/Current_Release/LmajorFriedlin/txt/TriTrypDB-8.1_LmajorFriedlinGene.txt

