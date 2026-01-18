# Metabolomics in incident diabetes mellitus

Code used to analyse and obtain the results in the paper: 

Barranco-Altirriba M, Granado-Casas M, Yanes O, Capellades J, Junza A, Franch-Nadal J, Vendrell J, Llauradó G, Valdés S, García-Escobar E, Bermúdez-López M, Valdivielso JM, López-Lifante V-M, Herrero-Alonso C, Falguera M, Vilanova MB, Arteaga I, Torán-Monserrat P, Perera-Lluna A, Castelblanco E and Mauricio D (2025) Guanine and pregnenolone sulfate are associated with incident type 2 diabetes in two independent populations. Front. Endocrinol. 16:1706886. doi: 10.3389/fendo.2025.1706886

## Data download

After cloning this repository, run the following commands **from inside the repository directory** to download the dataset:

```
chmod +x download_data.sh
./download_data.sh
```

## Code execution

The files should be executed in the following order:

* The `discovery_analysis.R` file contains the code to obtain the statistical results of the discovery stage.

* The `validation_analysis.R` file contains the code to obtain the statistical results of the validation stage.

* The `mwise_annotation.R` file uses mWISE R package to give putative annotations to all peaks in the discovery data.

* The `FELLA_enrichment.R` file performs pathway enrichment using FELLA from the previously obtained mWISE annotations and plots the results.

* The `Script_plot_results.R` file contains the code to obtain the results shown in the paper.

Some of these files read auxiliary functions from the `Functions_discovery.R` and the `Functions_validation.R` files.
