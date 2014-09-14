TriTrypDB OrganismDB package generator
======================================

Overview
--------

This script uses resources from TriTrypDB to generate a Bioconductor Organism
annotation package. 

Requirements
------------

For this script to work, you must have installed the following R libraries:

- [yaml](http://cran.r-project.org/web/packages/yaml/index.html)
- [rtracklayer](http://www.bioconductor.org/packages/release/bioc/html/rtracklayer.html)

Usage
-----

The use this script, begin by modifying the example configuration file,
`config.example.yaml` and save it as `config.yaml`.

From there, simply run `build.r`:

```sh
$ Rscript build.r
```

The script should generate an orgDB package in the working directory. You can
then use the `install.packages` command to install the annotations database
locally, e.g.:

```r
install.packages("./org.Lmajor.eg.db", repos=NULL)
```

Getting Help
------------

This script has been designed to work with any TriTrypDB organisms for
which both GFF and gene TXT annotations are available. Currently, however, it
has only been tested on L. major Friedlin annotations. If you encounter issues
using it for your organism of interest, feel free to submit fixes via pull
requests, or report the problem as an issue.

For questions, feel free to contact me at [khughitt@umd.edu](khughitt@umd.edu).

