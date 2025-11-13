# Metabolomics in incident diabetes mellitus

Code used to analyse and obtain the results in the paper: 

Barranco-Altirriba M. et al. Guanine and pregnenolone sulfate are significantly associated with incident type 2 diabetes in two independent populations.

The files should be executed in the following order:

* The `discovery_analysis.R` file contains the code to obtain the statistical results of the discovery stage.

* The `validation_analysis.R` file contains the code to obtain the statistical results of the validation stage.

* The `mwise_annotation.R` file uses mWISE R package to give putative annotations to all peaks in the discovery data.

* The `FELLA_enrichment.R` file performs pathway enrichment using FELLA from the previously obtained mWISE annotations and plots the results.

* The `Script_plot_results.R` file contains the code to obtain the results shown in the paper.

Some of these files read auxiliary functions from the `Functions_discovery.R` and the `Functions_validation.R` files.
