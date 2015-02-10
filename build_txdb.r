#!/usr/bin/env Rscript-devel
###############################################################################
#
# TriTrypDB TxDB package generation
#
# This script uses resources from TriTrypDB to generate a Bioconductor
# Transcript annotation package.
# 
###############################################################################
library(yaml)
library(rtracklayer)
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

# chromosome info
gff = import.gff3(settings$gff)
ch = gff[gff$type == 'chromosome']

chrom_info = data.frame(
    chrom=ch$ID,
    length=as.numeric(ch$size),
    is_circular=NA
)

# WORK-AROUND 2015/02/29
# Added 'useGenesAsTranscripts' flag to force inclusion of ncRNAs
# https://github.com/elsayed-lab/eupathdb-organismdb/issues/1
txdb = makeTxDbFromGFF(
    file=settings$gff,
    format='gff3',
    chrominfo=chrom_info,
    exonRankAttributeName=NA,
    dataSource=sprintf('TriTrypDB %s', settings$db_version),
    species=paste(settings$genus, settings$species),
    gffTxName=settings$gff_txname
)

# Save transcript database
short_name = paste0(substring(tolower(settings$genus), 1, 1), settings$species)
saveDb(txdb, file=file.path(settings$build_dir, sprintf("%s.sqlite", short_name)))

# Build TxDB package
makeTxDbPackage(
    txdb,
    destDir=settings$output_dir,
    version=settings$db_version,
    maintainer=settings$maintainer,
    author=settings$author,
    license='Artistic-2.0'
)

