# Cheminformatics Based Approach to Exploring and Modeling Trait-Associated Metabolite Profiles
*Jeremy R. Ash, Melaine A. Kuenemann, Daniel Rotroff, Alison Motsinger-Reif, and Denis Fourches*

## Usage
To run any of the r scripts provided, clone the repository and open the **metabochem.Rproj** file in Rstudio.
Run the scripts within the project, so that they can find the relevant paths on your machine.

## Files

* The *analyses* folder contains three rmarkdown files and one r script, which perform the majority of the analyses reported in the paper.
	* **differential_analyses.rmd**: Performs the differential analysis on cancer and healthy patient metabolites.  Outputs the significance results. See already rendered report: **differential_analysis.pdf.**
	* **metab_clus.r**: Filters metabolites by significance and clusters on chemical structure. Outputs clustered metabolite concentration profiles.
	* **metab_classifier_plasma.rmd**: Takes the clustered plasma metabolite concentration profiles as input.  Builds and validates classification models predicting patient cancer status.  See already rendered report: **metab_classifier_plasma.html.**
	* **metab_classifier_serum.rmd**: Takes the clustered serum metabolite concentration profiles as input.  Builds and validates classification models predicting patient cancer status.  See already rendered report: **metab_classifier_serum.html.**
* The *analyses/data folder* contains the data used by each analysis script.
	* **sample_metabolites_training_excol_fix.csv**, **sample_metabolites_test.txt**: raw patient metabolite concentration profiles.
	* **training_set_tq_normalize.rda**, **test_set_tq_normalize.rda**: fully processed and normalized patient metabolite concentration profiles.
	* **sample_factors_train.txt**, **ssample_factors_test.txt**: patient characteristic data, including cancer status.
	* **common_test_training_molecule-v2.sdf**: metabolite chemical structures.
	* **metabolics_fingerprint.csv**: metabolite fingerprints.
	* **ML_data.RDATA**: clustered metabolites concentration profiles.
* The *analyses/data folder/signficance_results* folder contains the signficance_results from the differential analysis.
* The *analyses/figure_scripts* folder contains the scripts to recreate a few of the figures in the paper.  The file names contain the corresponding figure numbers.
