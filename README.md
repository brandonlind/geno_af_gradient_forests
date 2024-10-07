# A comparison of genomic forecasts based on genotypes versus allele frequencies

Accelerating land use and climate change threaten to disrupt relationships between adaptive variation and environmental optima of many species. Consequently, management must increasingly identify non-local genetic sources for restoration programs. Genomic offset methods, like gradientForests, have shown promise in identifying these sources using genomic data, potentially bypassing the need for traditional, time-consuming transplant experiments. However, previous studies primarily used population-level allele frequencies (AF) for training and population-mean fitness for evaluation, ignoring individual variation within populations. Here, we used simulation data to compare the accuracy of genotype- and AF-based models, each evaluated using both individual and population-mean fitness. With over 810,000 evaluations of such models, we found that the number of loci had little impact on model performance. As expected, population-level evaluation provided an optimistic view of predictive performance for both genomic inputs. While genotype- and AF-based models showed similar qualitative and quantitative aspects, genotype-based models improved predictions in landscapes that differed from strict environmental clines by incorporating additional loci beyond those used by AF-based models. This suggests genotype-based models may enhance offset predictions in environments that are discontinuous and have multiple populations in geographically distinct yet similar environments. We close with recommendations for future use and evaluation of these tools.

# Funding

This research was funded by NSF-2043905 (KEL) and Northeastern University.

# Citation

Lind & Lotterhos (2024) A comparison of genomic forecasts based on genotypes versus allele frequencies.

contact information:

Brandon Lind - lind dot brandon dot m (at) gmail dot com

Katie Lotterhos - k dot lotterhos (at) northeastern dot edu

# Conda environments

Conda environments are the the same as that used in [Lind & Lotterhos (2024)](https://github.com/brandonlind/mvp-offsets). Specifically, we used the mvp_env.yml and r35.yml environments. Package and coding versions are available at the top of the notebooks described below. Data used to train models has been [archived previously](https://doi.org/10.26008/1912/bco-dmo.889769.1).

# Code Descriptions

`runtime_API.py` - this file is imported into many of the notebooks and used to load data, metadata, and arguments for making figures.

Below are the descriptions of notebooks in this repo. Notebooks can be viewed in the repository but are best viewed at https://nbviewer.jupyter.org (hyperlinks below).

### 00_create_datasets
[00_set_up_loci_sets](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/00_create_datasets/00_set_up_loci_sets.ipynb)

create random sets of SNP files for 3 reps from each of 225 simulation seeds

### 01_individual_runs
Train GF models using genotypes, calculate performance at the individual level.

[01_kick_off_individual_GF_runs](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/01_individual_runs/01_kick_off_individual_GF_runs.ipynb)

create genotype runs of GF using the sets of random loci assigned to individual runs

[02_submit_remaining_individual_jobs](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/01_individual_runs/02_submit_remaining_individual_jobs.ipynb)

check on currently submitted jobs for genotype runs and resubmit any jobs that failed

[03_gather_individual_scores](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/01_individual_runs/03_gather_individual_scores.ipynb)

gather performance scores from geno runs into one object

### 02_pooled_runs
Train GF models using allele frequencies, calculate performance of models at the population level

[01_kick_off_pooled_GF_runs](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/02_pooled_runs/01_kick_off_pooled_GF_runs.ipynb)

create allele frequency runs of GF using the sets of random loci assigned to individual runs

[02_submit_remaining_pooled_jobs](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/02_pooled_runs/02_submit_remaining_pooled_jobs.ipynb)

check on currently submitted jobs for AF runs and resubmit any jobs that failed

[03_gather_pooled_scores](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/02_pooled_runs/03_gather_pooled_scores.ipynb)

gather scores from pooled runs for the runtime project

### 03_calculate_geno-ind_performance

Calculate performance at the individual level using genotype models

[00_calc_ind-averaged_performance](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/03_calculate_geno-ind_performance/00_calc_ind-averaged_performance.ipynb)

see how averaging across individuals affects perceived performance.

### 04_main_questions

Answer main questions outlined in the manuscript.

[01_Q1_effect_of_marker_set_size](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/04_main_questions/01_Q1_effect_of_marker_set_size.ipynb)

Answer Q1 of the manuscript: How does the number of markers used as input affect performance across population- and individual-level datasets?

[02_Q2_Q3_effect_of_genetic_source](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/04_main_questions/02_Q2_Q3_effect_of_genetic_source.ipynb)

Answer Q2 and Q3 of the manuscript:
  Q2 How does the format of evaluation data affect performance?
  Q3 How does the format of the training data affect model performance?


[03_Q4_computational_requirements](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/04_main_questions/03_Q4_computational_requirements.ipynb)

Answer Q4 from the manuscript: How does the size of the dataset affect computational time and memory requirements?

### 05_supplemental

[01_pca_loading_corrs.ipynb](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/05_supplement/01_pca_loading_corrs.ipynb)

why do GF runs using 500 markers do about as well as runs using 10k-20k markers? 

[03_all_compare_workflows](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/05_supplement/03_all_compare_workflows.ipynb)

see which datasets differed the most between AF and genotype data

[04_check_overlap_of_loci](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/05_supplement/04_check_overlap_of_loci.ipynb)

check the overlap of loci used by *GO<sub>geno,ind</sub>* and *GO<sub>AF,pop</sub>* models

[05_explore_r2_geno_ind_loci](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/05_supplement/05_explore_r2_geno_ind_loci.ipynb)

explore differences in R2 from loci used by GF_geno models but not GF_AF models

[06_determine_levels_of_failed_replicates](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/05_supplement/06_determine_levels_of_failed_replicates.ipynb)

figure out if there is any commonality regarding the simulation levels for the 269 replicates that died when using 20k loci encoded as individual genotypes

[07_conceptual_fig](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/05_supplement/07_conceptual_fig.ipynb)

create a 4x4 conceptual figure that shows a scatter plot between offset and fitness for seed 1231422. The 4x4 is evaluation level (rows - ind or pop) and input data (columns - AF or geno)

[08_lotterhos_env_figs](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/05_supplement/08_lotterhos_env_figs.ipynb)

recreate the environmental figures from Lotterhos 2023

### 06_calc_af_ind

Calulate performance at the individual level using allele frequency models

[00_af_ind_performance.ipynb](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/06_calc_af_ind/00_af_ind_performance.ipynb)


---

`pythonimports` in notebooks can be found here: https://github.com/brandonlind/pythonimports















