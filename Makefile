SPECIES=tcruzi_clbrener
ORG=org.TCCLB.tritryp.db

install:	prereq forge orgdb

forge:
	cp settings/${SPECIES}.yaml ./config.yaml
	./build_orgdb.r
	./build_txdb.r
	./build_organismdb.r
	sh prepare_dbs.sh

orgdb:
	Rscript -e "install.packages("./${ORG}", repos=NULL)"
	Rscript build_organismdb.r

prereq:
	Rscript -e "yaml_test = try(suppressPackageStartupMessages(library('yaml')));\
if (class(yaml_test) == 'try-error') { source('http://bioconductor.org/biocLite.R'); biocLite(yaml) };\
rtracklayer_test = try(suppressPackageStartupMessages(library('rtracklayer')));\
if (class(rtracklayer_test) == 'try-error') { source('http://bioconductor.org/biocLite.R'); biocLite('rtracklayer') };\
gf_test = try(suppressPackageStartupMessages(library('GenomicFeatures')));\
if (class(gf_test) == 'try-error') { source('http://bioconductor.org/biocLite.R'); biocLite('GenomicFeatures') };\
org_test = try(suppressPackageStartupMessages(library('OrganismDbi')));\
if (class(org_test) == 'try-error') { source('http://bioconductor.org/biocLite.R'); biocLite('OrganismDbi') };"
