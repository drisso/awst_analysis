# Code to reproduce the analyses of the AWST paper

The code in this repository can be used to reproduce the results of Risso and Pagnotta (2020).

The folder `SyntheticExperiment` contains the code to reproduce the results in Figure 2, Figure 3, Supplementary Figures 1-5, and Supplementary Tables 1-2, 4.

The folder `scMixology` contains the code to reproduce the results in Figure 4 and Supplementary Table 3.

The folder `LowerGradeGliomaExperiment` contains the code to reproduce the results in Figure 5 and Supplementary Figures 6-20.

The folder `CordBloodMononuclearCellsExperiment` contains the code to reproduce the results in Figure 6 and Supplementary Figures 21-26.

## Package versions

In order to perfectly reproduce the analyses of the paper, you need version 0.0.4 of the `awst` package, available [here](https://github.com/drisso/awst/releases/tag/v0.0.4). Note that this is a frozen version in which no further development nor bug fixes will be done.

The latest, stable version can be downloaded through the Bioconductor project at https://bioconductor.org/packages/awst.

Note that the biggest change between version 0.0.4 and the Bioconductor version is that in the latter the `awst` function will return a matrix of transformed expressions with genes in row and samples in column. This behavior is more consistent with the Bioconductor ecosystem. If you want to run the code in this repository using the latest version of `awst` you have to transpose the result of `awst` with the `t()` function.

## References

Risso and Pagnotta (2021). Per-sample standardization and asymmetric winsorization lead to accurate clustering of RNA-seq expression profiles. [Bioinformatics](https://doi.org/10.1093/bioinformatics/btab091). [[Preprint at biorXiv](https://doi.org/10.1101/2020.06.04.134916).]
