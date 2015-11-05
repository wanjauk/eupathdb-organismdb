install: orgdb txdb organismdb

orgdb:
	Rscript ./build_orgdb.R
txdb:
	Rscript ./build_txdb.R
organismdb:
	sh prepare_dbs.sh
	Rscript ./build_organismdb.R
	sh finalize.sh
