#!/usr/bin/env Rscript
###############################################################################
#
# TriTrypDB OrganismDB package generation
#
# This script uses resources from TriTrypDB to generate a Bioconductor
# OrganismDB annotation package.
#
# Note that in order for this to work properly, TranscriptDB and OrgDB
# databases must first be generated using the other build_xx.r scripts.
# 
###############################################################################
library(yaml)
library(GenomicFeatures)
library(OrganismDbi)

options(stringsAsFactors=FALSE)

#
# MAIN
#
# Load settings
settings = yaml.load_file('config.yaml')

# Create build directory
if (!file.exists(settings$build_dir)) {
    dir.create(settings$build_dir)
}
build_basename = file.path(settings$build_dir,
                            sub('.gff', '', basename(settings$gff)))

# Load organism-specific packages
library(settings$orgdb_name, character.only=TRUE)
library(settings$txdb_name,  character.only=TRUE)

# mapping
graph_data = list(
    join1=c(GO.db='GOID', orgdb='GO'),
    join2=c(orgdb='GID',  txdb='GENEID')
)

names(graph_data$join1) = c('GO.db', settings$orgdb_name)
names(graph_data$join2) = c(settings$orgdb_name, settings$txdb_name)

makeOrganismPackage(
    pkgname=settings$organismdb_name,
    graphData=graph_data,
    organism=paste(settings$genus, settings$species),
    version=settings$db_version,
    maintainer=settings$maintainer,
    author=settings$author,
    destDir=settings$output_dir,
    license='Artistic-2.0'
)
