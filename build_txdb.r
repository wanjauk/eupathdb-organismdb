###############################################################################
#
# TriTrypDB TxDB package generation
#
# This script uses resources from TriTrypDB to generate a Bioconductor
# Transcript annotation package.
# 
###############################################################################
library(yaml)
library(GenomicFeatures)

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

txdb = makeTranscriptDbFromGFF(
    file=settings$gff,
    format='gff3',
    exonRankAttributeName=NA,
    dataSource=sprintf('TriTrypDB %s', settings$tritrypdb_version),
    species=paste(settings$genus, settings$species)
)

# Save transcript database
short_name = paste0(substring(tolower(settings$genus), 1, 1), settings$species)
saveDb(txdb, file=sprintf("%s.sqlite", short_name))

# Build TxDB package
makeTxDbPackage(
    txdb,
    version=settings$tritrypdb_version,
    maintainer=settings$maintainer,
    author=settings$author,
    license='Artistic-2.0'
)

