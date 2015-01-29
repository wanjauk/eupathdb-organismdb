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
- [GenomicFeatures](http://www.bioconductor.org/packages/release/bioc/html/GenomicFeatures.html)
- [OrganismDbi](http://www.bioconductor.org/packages/release/bioc/html/OrganismDbi.html)

Usage
-----

The use this script, begin by modifying or creating a new YAML configuration
file specifying the details for the species you wish to process. Examples
configuration files for several species can be found in the `settings` folder.
Next, modify the `build_xx.r` scripts to point to the setting file you wish you
use.

From there, simply run the build scripts in the following order:

```sh
$ Rscript build_orgdb.r
$ Rscript build_txdb.r
$ Rscript build_organismdb.r
$ sh finalize.sh
```

This will generate OrgDb, TranscriptDb, and OrganismDb packages in the current
working directory. The final script performs some post-processing to replace
default names with desired ones, etc.

You can then use the `install.packages` command to install the annotations
database locally, e.g.:

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

See Also
--------

- [AnnotationForge - Making Organism packages](http://www.bioconductor.org/packages/release/bioc/vignettes/AnnotationForge/inst/doc/MakingNewOrganismPackages.html)
- [Bioconductor - Making and Utilizing TxDb Objects](http://www.bioconductor.org/packages/release/bioc/vignettes/GenomicFeatures/inst/doc/GenomicFeatures.pdf)
- [Bioconductor - Making an OrganismDb package](http://www.bioconductor.org/help/workflows/annotation/#Making-an-OrganismDb-package)
- [Bioconductor - AnnotationForge](http://www.bioconductor.org/packages/release/bioc/html/AnnotationForge.html)
- [Bioconductor - GO.db](http://www.bioconductor.org/packages/release/data/annotation/html/GO.db.html)
- [Bioconductor - OrganismDbi](http://www.bioconductor.org/packages/release/bioc/html/OrganismDbi.html)

TODO
----

1. Rename packages (e.g. 'org.LmajorFriedlin.tritryp.db')
2. Add strain information
3. Add KEGG annotations
4. Add UTRs
5. Fix DBSCHEMA information

