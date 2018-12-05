# DecompPipeline
Large automated pipeline for running MeDeCom 

# Using Decomp
*DecompPipeline* includes three major steps, all of them are extensively documented. A more detailed introduction into *DecompPipeline* can be found in the package vignette (https://github.com/lutsik/DecompPipeline/blob/master/vignettes/DecompPipeline.md).

## CpG filtering
There are dedicated preprocessing steps for both array-based data sets (```prepare_data```) and sequencing-based data sets (```prepare_data_BS```).

## Selecting CpG subsets
To select a subset of CpGs for downstream MeDeCom analysis, the function ```prepare_CG_subsets``` can be used.

## Starting MeDeCom
After these preprocessing steps, a MeDeCom run can be started using *DecompPipeline* by envoking ```start_decomp_pipeline```.

## Combining the above steps
We also provide a wrapper functions that executes all the above steps in a single command (```start_decomp_pipeline```).
