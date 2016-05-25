install: orgdb txdb organismdb

orgdb:
	Rscript ./build_orgdb.R
txdb:
	Rscript ./build_txdb.R
organismdb:
	./prepare_dbs.sh
	./build_organismdb.R
	./sh finalize.sh
