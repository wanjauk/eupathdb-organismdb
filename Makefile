.PHONY: clean

%:	%.yaml
	@rm -r build && mkdir build
	@echo "Generating orgdb package."
	./orgdb.R $*.yaml
	@echo "Finished generating $* orgdb package."
	@echo "Generating TxDb package."
	./txdb.R $*.yaml
	@echo "Finished generating $* txdb package."
	@echo "Installing TxDb and orgdbs."
	./prepare.sh $*.yaml
	@echo "Generating organismDbi package."
	./organismdb.R $*.yaml
	@echo "Finished generating $* organismdb package, installing it."
	./final.sh $*.yaml
	@echo "Installed packages in R."

clean:
	@rm -r build && mkdir build

install: orgdb txdb organismdb

orgdb:
	./orgdb.R
txdb:
	./txdb.R
organismdb:
	./prepare.sh
	./organismdb.R
	./final.sh
