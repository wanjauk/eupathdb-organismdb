install:	prereq orgdb txdb organismdb

orgdb:
	Rscript ./build_orgdb.R
txdb:
	Rscript ./build_txdb.R
organismdb:
	sh prepare_dbs.sh
	Rscript ./build_organismdb.R
	sh finalize.sh

prereq:
	Rscript -e "yaml_test = try(suppressPackageStartupMessages(library('yaml')));\
if (class(yaml_test) == 'try-error') { source('http://bioconductor.org/biocLite.R'); biocLite(yaml) };\
rtracklayer_test = try(suppressPackageStartupMessages(library('rtracklayer')));\
if (class(rtracklayer_test) == 'try-error') { source('http://bioconductor.org/biocLite.R'); biocLite('rtracklayer') };\
gf_test = try(suppressPackageStartupMessages(library('GenomicFeatures')));\
if (class(gf_test) == 'try-error') { source('http://bioconductor.org/biocLite.R'); biocLite('GenomicFeatures') };\
org_test = try(suppressPackageStartupMessages(library('OrganismDbi')));\
if (class(org_test) == 'try-error') { source('http://bioconductor.org/biocLite.R'); biocLite('OrganismDbi') };"
