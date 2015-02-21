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

build_dir = file.path(settings$build_dir,
                      paste0(R.Version()$major,  '.', R.Version()$minor))

# Create build directory
if (!file.exists(build_dir)) {
    dir.create(build_dir)
}
build_basename = file.path(build_dir,
                            sub('.gff', '', basename(settings$gff)))

# chromosome info
gff = import.gff3(settings$gff)
ch = gff[gff$type %in% c('apicoplast_chromosome', 'chromosome', 'supercontig')]

#genes = gff[gff$type == 'gene']
#gene_ch = unique(as.character(chrom(genes)))

chrom_info = data.frame(
    chrom=ch$ID,
    length=as.numeric(ch$size),
    is_circular=NA
)

# WORK-AROUND 2015/02/29
# Added 'useGenesAsTranscripts' flag to force inclusion of ncRNAs
# https://github.com/elsayed-lab/eupathdb-organismdb/issues/1

# 2015/02/10 Switching back to Bioconductor 3.0 until 3.1 is out in the wild
#txdb = makeTxDbFromGFF(
txdb = makeTranscriptDbFromGFF(
    file=settings$gff,
    format='gff3',
    chrominfo=chrom_info,
    exonRankAttributeName=NA,
    dataSource=sprintf('TriTrypDB %s', settings$db_version),
    species=paste(settings$genus, settings$species),
    useGenesAsTranscripts=TRUE
)

#gffTxName=settings$gff_txname

# Save transcript database
short_name = paste0(substring(tolower(settings$genus), 1, 1), settings$species)
saveDb(txdb, file=file.path(build_dir, sprintf("%s.sqlite", short_name)))

# Build TxDB package
makeTxDbPackage(
    txdb,
    destDir=settings$output_dir,
    version=settings$db_version,
    maintainer=settings$maintainer,
    author=settings$author,
    license='Artistic-2.0'
)

