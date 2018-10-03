#!/usr/bin/env Rscript
###############################################################################
##
## EuPathDB TxDB package generation
##
## This script uses resources from EuPathDB to generate a Bioconductor
## Transcript annotation package.
##
###############################################################################
suppressMessages(library(yaml))
suppressMessages(library(rtracklayer))
suppressMessages(library(Biostrings))
suppressMessages(library(GenomicFeatures))

options(stringsAsFactors=FALSE)

config_file <- "config.yaml"
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

# Create build directory
if (!file.exists(build_dir)) {
    dir.create(build_dir)
}
build_basename <- file.path(build_dir,
                            sub('.gff', '', basename(settings$gff)))

#
# chromosome info
#
# NOTE: In the future, perhaps it would be better instead to first filter the
# source GFF and remove all non-standard entries?
#
# non-standard types:
#
# apicoplast_chromosome (T. gondii ME49)
# random_sequence       (T. brucei TREU927)
# geneontig             (T. brucei TREU927)
#

# load GFF
gff <- import.gff3(settings$gff)

# load FASTA
fasta <- readDNAStringSet(settings$fasta)

# Get chromosome length information from FASTA file
chrom_info <- data.frame(
    chrom=sapply(strsplit(names(fasta), ' '), '[[', 1),
    length=width(fasta),
    is_circular=NA
)

# Bioconductor 3.5
if (BiocInstaller::biocVersion() >= 3.5) {
    txdb <- makeTxDbFromGFF(
        file=settings$gff,
        format='gff3',
        chrominfo=chrom_info,
        dataSource=sprintf('%s %s', settings$db_name, settings$db_version),
        organism=paste(settings$genus, settings$species),
        metadata=c('Resource URL', settings$db_url)
    )
} else {
    # Bioconductor 3.4 (work-around)

    # generate metadata dataframe
    metadata <- GenomicFeatures:::.prepareGFFMetadata(
        settings$gff,
        sprintf('%s %s', settings$db_name, settings$db_version),
        paste(settings$genus, settings$species),
        settings$tax_id)

    # Add required 'Resource URL' field to the metadata
    metadata <- rbind(metadata, c('Resource URL', settings$db_url))

    gr <- import(settings$gff, format='gff3',
                colnames=GenomicFeatures:::GFF3_COLNAMES,
                feature.type=GenomicFeatures:::GFF_FEATURE_TYPES)

    circ_seqs <- GenomicFeatures:::DEFAULT_CIRC_SEQS
    gr <- GenomicFeatures:::.tidy_seqinfo(gr, circ_seqs, NULL) 

    ## work-around 2016/02/08
    txdb <- makeTxDbFromGRanges(gr, metadata=metadata)

    # Save transcript database
    short_name <- paste0(substring(tolower(settings$genus), 1, 1), settings$species)
    saveDb(txdb, file=file.path(build_dir, sprintf("%s.sqlite", short_name)))
}

# R package versions must be of the form "x.y"
db_version <- paste(settings$db_version, '0', sep='.')

# Build TxDB package
result <- tryCatch({
    makeTxDbPackage(
        txdb,
        version=db_version,
        maintainer=settings$maintainer,
        author=settings$author,
        destDir=output_dir,
        license='Artistic-2.0',
        pkgname=settings$txdb_name,
        provider=settings$db_name,
        providerVersion=settings$db_version
    )
}, error=function(e) {
    # above fails in recent versions of Bioconductor due to 'inst/extdata'
    # directly not being present
    db_path <- file.path(output_dir, settings$txdb_name, "inst",
                         "extdata", paste(settings$txdb_name,"sqlite",sep="."))
    dir.create(dirname(db_path), recursive=TRUE)
    saveDb(txdb, file=db_path)
})

