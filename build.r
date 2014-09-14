###############################################################################
#
# TriTrypDB OrganismDB package generation
#
# This script uses resources from TriTrypDB to generate a Bioconductor Organism
# annotation package.
# 
###############################################################################
library(yaml)
library(rtracklayer)

# Load settings
settings = yaml.load_file("config.yaml")

# Parse GFF
gff = import.gff3(settings$gff)
genes = gff[gff$type == 'gene']

#
# 1. Gene name and description
#
gene_info = elementMetadata(genes[,c('ID', 'description')])

# Convert form-encoded description string to human-readable
gene_info$description = gsub("\\+", " ", gene_info$description)
colnames(gene_info) <- c("GID", "GENENAME")

#
# 2. Chromosome information
#
chr_info = data.frame(
    'GID' = genes$ID,
    'CHROMOSOME' = as.character(seqnames(genes))
)

# Generate package
makeOrgPackage(
    gene_info  = gene_info,
    chromosome = chr_info,
    version    = settings$tritrypdb_version,
    author     = settings$author,
    maintainer = settings$maintainer,
    outputDir  = '.',
    tax_id     = settings$tax_id,
    genus      = settings$genus,
    species    = settings$species
)

# TriTrypDB download URL
# http://tritrypdb.org/common/downloads/Current_Release/LmajorFriedlin/
# http://tritrypdb.org/common/downloads/Current_Release/LmajorFriedlin/gff/data/TriTrypDB-8.1_LmajorFriedlin.gff
# http://tritrypdb.org/common/downloads/Current_Release/LmajorFriedlin/txt/TriTrypDB-8.1_LmajorFriedlinGene.txt

