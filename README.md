#    Phylogenetic-modeling-of-evolutionary-trait-jumps

**Brief Description**

This repository contains all the required scripts, and files associated with our preprint:

This git repository includes:

README.md - It contains a brief description about all the required scripts, and files.

funtions_TPCM.R - Custom functions required to reproduce the results of this study.

Preprocess_script.R - This script needs to be executed before running Trait_jump_model_paper.rmd to generate the file "norm_exp_15spe_4tips_new.RData" with intermediate results.

rphylopic_test_fig1_new.R - The script is written to generate figure 1 (phylopic image).

trait_jump_output_tau.csv - Using "levolution", evolutionary jumps in tissue-specifcity were estimated for branches of all our gene trees. This file contains mapped gene tree ids, branch information where trait jump occured, and their corresponding jump probabilities.

jump_model_analysis_tau.R - This script uses "trait_jump_output_tau.csv" file to perform jump model analyses for tissue-specificity.

trait_jump_output_avg_exp.csv - Using "levolution", evolutionary jumps in average gene expression levels  were estimated for branches of all our gene trees. This file contains mapped gene tree ids, branch information where trait jump occured, and their corresponding jump probabilities.

jump_model_analysis_avgexp.R - This script uses "trait_jump_output_avg_exp.csv" file to perform jump model analyses for average gene expression levels.

results_Brain.csv - Using "levolution", evolutionary jumps in brain were estimated for branches of all our gene trees. This file contains mapped gene tree ids, branch information where trait jump occured, and their corresponding jump probabilities.

trait_jump_brain.R - This script uses "results_Brain.csv" file to perform jump model analyses for expression in brain.

results_Heart.csv - Using "levolution", evolutionary jumps in heart were estimated for branches of all our gene trees. This file contains mapped gene tree ids, branch information where trait jump occured, and their corresponding jump probabilities.

trait_jump_heart.R - This script uses "results_Heart.csv" file to perform jump model analyses for expression in heart.

results_Kidney.csv - Using "levolution", evolutionary jumps in kidney were estimated for branches of all our gene trees. This file contains mapped gene tree ids, branch information where trait jump occured, and their corresponding jump probabilities.

trait_jump_kidney.R - This script uses "results_Kidney.csv" file to perform jump model analyses for expression in kidney.

results_Liver.csv - Using "levolution", evolutionary jumps in liver were estimated for branches of all our gene trees. This file contains mapped gene tree ids, branch information where trait jump occured, and their corresponding jump probabilities.

trait_jump_liver.R - This script uses "results_Liver.csv" file to perform jump model analyses for expression in liver.

results_Muscle.csv - Using "levolution", evolutionary jumps in muscle were estimated for branches of all our gene trees. This file contains mapped gene tree ids, branch information where trait jump occured, and their corresponding jump probabilities.

trait_jump_muscle.R - This script uses "results_Muscle.csv" file to perform jump model analyses for expression in muscle.

results_Testis.csv - Using "levolution", evolutionary jumps in testis were estimated for branches of all our gene trees. This file contains mapped gene tree ids, branch information where trait jump occured, and their corresponding jump probabilities.

trait_jump_testis.R - This script uses "results_Testis.csv" file to perform jump model analyses for expression in testis.

map_selectome_results.R - This script uses "selectome_trees_output.tsv" file (output from Selectome database) to detect and map positive selection on the branches of our gene trees. The mapped output is provided in "mapped.selectome.results.tsv".

BGEE_data_process.R- This script uses "devStageIds_BGEE.txt" to collect processed adult stage data from BGEE.v14 for eight species. Such processed data are provided in the "Processed_expression_data" folder. After normalization, the expression data are provided in the "Normalized_expression_data" folder.

formatted_preprocessed.R - This script uses preporcessed data for dog, rabbit and ferret provided by Chen et al. (2019) to format it. To collect data, we first used "sra2sh.pl". The output was then stored in "geo_extract.sh". We ran "geo_extract.sh" to obtain tissue wise expression data for the three species. We ran "formatted_preprocessed.R" script then to process tissue data, which we provided in the "Processed_expression_data" folder. After normalization, the expression data are provided in the "Normalized_expression_data" folder. 

Fish_data_process.R - This script is used to process kallisto outputs of medaka, cavefish, northern pike and tilapia, which we provided in the "Processed_expression_data" folder. Raw RNA-seq data were collected using Phylofish_2_kallisto.pl. After normalization, the expression data are provided in the "Normalized_expression_data" folder.

tau.R - To compute tissue specificity and average gene expression levels on processed expression data for all the fifteen species, we used this script.

TPCM.bib - Bibliography file required to knit "Trait_jump_model_paper.rmd".

molecular-biology-and-evolution.csl - This file required to knit "Trait_jump_model_paper.rmd" to format citations and bibliographies.

Trait_jump_model_paper.rmd - The manuscript text and source code for presenting our analyses using the saved output "norm_exp_15spe_4tips_new.RData"  archived at https://doi.org/10.5281/zenodo.10418895. This Rmarkdown file includes the computed results and the source code for generating figures.
