###############################################################################
#
# TriTrypDB OrganismDB package generation
#
# This script uses resources from TriTrypDB to generate a Bioconductor Organism
# annotation package.
# 
###############################################################################
library(yaml)

# Load settings
settings = yaml.load_file("config.yaml")

# Generate package
makeOrgPackage(
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

