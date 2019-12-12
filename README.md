# DecompPipeline
*DecompPipeline* provides a comprehensive list of preprocessing functions for performing reference-free deconvolution of complex DNA methylation data. It is an integral part of a recently published protocol [doi:10.1101/853150](https://doi.org/10.1101/853150) for reference-free deconvolution and is independent of the deconvolution tool used. 

![Overview of the reference-free deconvolution tool](pictures/protocol_overview.png)

# Installing DecompPipeline
*DecompPipeline* can be directly installed from GitHub within an R session:
```
install.packages("devtools")
devtools::install_github("CompEpigen/DecompPipeline")
```

# Using Decomp
*DecompPipeline* includes three major steps, all of them are extensively documented. A more detailed introduction into *DecompPipeline* can be found in the package [vignette](vignettes/DecompPipeline.md) and in the [protocol](vignettes/DeconvolutionProtocol.Rmd) .

## 1. CpG filtering
There are dedicated preprocessing steps for both array-based data sets (```prepare_data```) and sequencing-based data sets (```prepare_data_BS```).

## 2. Selecting CpG subsets
To select a subset of CpGs for downstream deconvolution analysis, the function ```prepare_CG_subsets``` can be used.

## 3. Starting MeDeCom
After these preprocessing steps, a deconvolution run can be started using *DecompPipeline* by envoking ```start_decomp_pipeline```.

## Combining the above steps
We also provide a wrapper functions that executes all the above steps in a single command (```start_decomp_pipeline```).
