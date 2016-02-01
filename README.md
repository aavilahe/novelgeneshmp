## About

This package accompanies the manuscript NOVELGENES_TITLE_HERE. It includes
functions used to clean, filter, and analyze data for associations between
FUnkSFAM presence and subject metadata in the manuscript. See the
[vignette](./novelgeneshmp/vignettes/example-run.md) for an example.

## Install

Download the tarball and install with

```bash

R CMD INSTALL novelgeneshmp_0.0.0.9000.tar.gz

```

Or with [devtools](https://github.com/hadley/devtools)

```r

library(devtools)
install_github('aavilaherrera/novelgeneshmp@dev', subdir = 'novelgeneshmp')

# To build vignettes locally (takes a few extra seconds):
install_github('aavilaherrera/novelgeneshmp@dev',
               subdir = 'novelgeneshmp',
               build_vignettes = TRUE
               )

```

## License

GPLv3. See [LICENSE](./novelgeneshmp/LICENSE).
