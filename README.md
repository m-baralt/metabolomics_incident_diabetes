# Metabolomics in incident diabetes mellitus

Code used to analyse and obtain the results in the paper: 

Barranco-Altirriba M, Granado-Casas M, Yanes O, Capellades J, Junza A, Franch-Nadal J, Vendrell J, Llauradó G, Valdés S, García-Escobar E, Bermúdez-López M, Valdivielso JM, López-Lifante V-M, Herrero-Alonso C, Falguera M, Vilanova MB, Arteaga I, Torán-Monserrat P, Perera-Lluna A, Castelblanco E and Mauricio D (2025) Guanine and pregnenolone sulfate are associated with incident type 2 diabetes in two independent populations. Front. Endocrinol. 16:1706886. doi: 10.3389/fendo.2025.1706886

## Data download

Git clone this repository in your local folder using:

```
git clone https://github.com/m-baralt/metabolomics_incident_diabetes.git
```

After cloning, run the following commands **from inside the repository directory** to download the dataset:

```
chmod +x download_data.sh
./download_data.sh
```

## Code execution

The R scripts should be executed in the following order. Make sure to update any paths to your local repository as needed.

* `discovery_analysis.R` Performs the statistical analysis of the discovery stage. Update the repository path in the first lines of the script before running.

* The `validation_analysis.R` file contains the code to obtain the statistical results of the validation stage.

* The `mwise_annotation.R` file uses mWISE R package to give putative annotations to all peaks in the discovery data.

* The `FELLA_enrichment.R` file performs pathway enrichment using FELLA from the previously obtained mWISE annotations and plots the results.

* The `Script_plot_results.R` file contains the code to obtain the results shown in the paper.

Some of these files read auxiliary functions from the `Functions_discovery.R` and the `Functions_validation.R` files.
