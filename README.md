[![DOI](https://zenodo.org/badge/756374985.svg)](https://doi.org/10.5281/zenodo.13899117)

# A comparison of genomic forecasts based on genotypes versus allele frequencies

Accelerating land use and climate change threaten to disrupt relationships between adaptive variation and environmental optima of many species. Consequently, management must increasingly identify non-local genetic sources for restoration programs. Genomic offset methods, like gradientForests, have shown promise in identifying these sources using genomic data, potentially bypassing the need for traditional, time-consuming transplant experiments. However, previous studies primarily used population-level allele frequencies (AF) for training and population-mean fitness for evaluation, ignoring individual variation within populations. Here, we used simulation data to compare the accuracy of genotype- and AF-based models, factorially evaluated using both individual and population-mean fitness. With over 810,000 evaluations of such models, we found that the number of loci had little impact on model performance. As expected, population-level evaluation provided an optimistic view of predictive performance for both genomic inputs. While genotype- and AF-based models showed similar qualitative and quantitative aspects, genotype-based models improved predictions in landscapes that differed from strict environmental clines by incorporating additional loci beyond those used by AF-based models. This suggests genotype-based models may enhance offset predictions in environments that are discontinuous and have multiple populations in geographically distant yet similar environments. We close with recommendations for future use and evaluation of these tools.

# Funding

This research was funded by NSF-2043905 (KEL) and Northeastern University.

# Citation

Manuscript:
Lind & Lotterhos (2025) A comparison of genomic forecasts based on genotypes versus allele frequencies. <i> The American Naturalist</i>.

Archive:
Lind (2025) GitHub.com/brandonlind/geno_af_gradient_forests. Revision release (v1.0.1). Zenodo.

# contact information:

Brandon Lind - lind dot brandon dot m (at) gmail dot com

Katie Lotterhos - k dot lotterhos (at) northeastern dot edu

# Conda environments

Conda environments are the the same as that used in [Lind & Lotterhos (2024)](https://github.com/brandonlind/mvp-offsets). Specifically, we used the mvp_env.yml and r35.yml environments. Package and coding versions are available at the top of the notebooks described below. Data used to train models has been [archived previously](https://doi.org/10.26008/1912/bco-dmo.889769.1).

# Raw data

- Spatially discrete simulation output files that were formatted for GF training here are archived at [https://www.bco-dmo.org/data-set/889769](https://www.bco-dmo.org/dataset/889769)
- Spatially continuous simulation output files that were formatted for GF training here are archived at [available here](https://marineomics.github.io/RDAtraitPredictionTutorial.html).

# Evaluation results

The main results of our manuscript are available in gzip, tab-delimited format in the `/data` directory. Files are separated by the spatially continuous or spatially discrete subdirectories. Column and entry metadata are described for each file in the [`/data/REAMDE.md`](/data)

- data/spatially_continuous/validation.txt.gz
    - evaluation results for the continuous space simulation. This includes results from both individual-level environmental training data $GF_{geno, ind}(ind-env)$, as well as population-level environmental training data, $GF_{geno, ind}(pop-env)$
    - this file was created and saved in [07_continuous_space_sims/02_validate_continuous_sims_offset.ipynb](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/07_continuous_space_sims/02_validate_continuous_sims_offset.ipynb)
- data/spatially_discrete/af_ind/af_ind_results.txt.gz
    - evaluation results for the spatially discrete evaluation using allele frequencies and individual-level fitnesses, $GF_{AF, ind}$
    - this file was created and saved in [06_calc_af_ind/00_af_ind_performance.ipynb](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/06_calc_af_ind/00_af_ind_performance.ipynb)
- data/spatially_discrete/af_pop/pooled_performance.txt.gz
    - evaluation results for the spatially discrete evaluation using allele frequencies and population-mean fitnesses, $GF_{AF, pop}$
    - this file was created and saved in [02_pooled_runs/03_gather_pooled_scores.ipynb](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/02_pooled_runs/03_gather_pooled_scores.ipynb)
    - within notebooks this file is often loaded, concatenated with `ind-averaged_results.txt`, via the `runtime_API.load_results` function.
- data/spatially_discrete/geno_ind/ind-averaged_results.txt.gz
    - evaluation results for the spatially discrete evaluation using individual genotypes and individual-level fitnesses, $GF_{geno, ind}$
    - this file was created and saved in [03_calculate_geno-ind_performance/00_calc_geno-ind_performance.ipynb](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/03_calculate_geno-ind_performance/00_calc_geno-ind_performance.ipynb)
- data/spatially_discrete/geno_pop/ind_performance.txt.gz
    - evaluation results for the spatially discrete evaluation using individual genotypes and individual-level fitnesses, $GF_{geno, pop}$
    - this file was created and saved in [01_individual_runs/03_gather_individual_scores.ipynb](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/01_individual_runs/03_gather_individual_scores.ipynb)
    - within notebooks this file is often loaded, concatenated with `pooled_performance.txt`, via the `runtime_API.load_results` function.

# Code Descriptions

`runtime_API.py` - this file is imported into many of the notebooks and used to load data, metadata, and arguments for making figures and nesting results within data objects. Functions are described in each functions' docstring found in this file. In some ways, this can be considered a configuration file for the notebook scripts.

Below are the descriptions of notebooks in this repo. Notebooks were used to carry out the formatting of data and main results of this manuscript. Notebooks can be viewed in the repository but are best viewed at https://nbviewer.jupyter.org (hyperlinks below).

### 00_create_datasets

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[00_set_up_loci_sets](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/00_create_datasets/00_set_up_loci_sets.ipynb)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;create random sets of SNP files for 3 reps from each of 225 simulation seeds (only 3 reps of 180 simulation levels are analyzed, these are filtered out in `03_gather_individual_scores`)

### 01_individual_runs
Train and evaluate $GF_{geno, pop}$ models

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[01_kick_off_individual_GF_runs](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/01_individual_runs/01_kick_off_individual_GF_runs.ipynb)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;create genotype runs of GF using the sets of random loci assigned to individual runs

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[02_submit_remaining_individual_jobs](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/01_individual_runs/02_submit_remaining_individual_jobs.ipynb)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;check on currently submitted jobs for genotype runs and resubmit any jobs that failed

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[03_gather_individual_scores](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/01_individual_runs/03_gather_individual_scores.ipynb)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;gather performance scores from geno runs into one object

### 02_pooled_runs
Train and evaluate $GF_{AF, pop}$

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[01_kick_off_pooled_GF_runs](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/02_pooled_runs/01_kick_off_pooled_GF_runs.ipynb)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;create allele frequency runs of GF using the sets of random loci assigned to individual runs

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[02_submit_remaining_pooled_jobs](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/02_pooled_runs/02_submit_remaining_pooled_jobs.ipynb)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;check on currently submitted jobs for AF runs and resubmit any jobs that failed

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[03_gather_pooled_scores](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/02_pooled_runs/03_gather_pooled_scores.ipynb)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;gather scores from pooled runs for the runtime project

### 03_calculate_geno-ind_performance

Calculate $GF_{geno, ind}$ performance at the individual level using genotype models from 01_individual_runs

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[00_calc_ind-averaged_performance](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/03_calculate_geno-ind_performance/00_calc_geno-ind_performance.ipynb)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;see how averaging across individuals affects perceived performance.

### 04_main_questions

Answer main questions outlined in the manuscript.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[01_Q1_effect_of_marker_set_size](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/04_main_questions/01_Q1_effect_of_marker_set_size.ipynb)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Answer Q1 of the manuscript: How does the number of markers used as input affect performance?

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[02_Q2_Q3_effect_of_genetic_source](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/04_main_questions/02_Q2_Q3_effect_of_genetic_source.ipynb)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Answer Q2 and Q3 of the manuscript:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  Q2 How does the format of evaluation data affect performance?

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  Q3 How does the format of the genetic training data affect performance?


&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[03_Q5_computational_requirements](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/04_main_questions/03_Q5_computational_requirements.ipynb)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Answer Q5 from the manuscript: How does the size of the dataset affect computational time and memory requirements?

### 05_supplemental

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[01_pca_loading_corrs.ipynb](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/05_supplement/01_pca_loading_corrs.ipynb)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;why do GF runs using 500 markers do about as well as runs using 10k-20k markers? 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[03_all_compare_workflows](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/05_supplement/03_all_compare_workflows.ipynb)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;see which datasets differed the most between AF and genotype models

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[04_check_overlap_of_loci](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/05_supplement/04_check_overlap_of_loci.ipynb)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;check the overlap of loci used by *GO<sub>geno,ind</sub>* and *GO<sub>AF,pop</sub>* models

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[05_explore_r2_geno_ind_loci](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/05_supplement/05_explore_r2_geno_ind_loci.ipynb)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;explore differences in R2 from loci used by $GF_{geno}$ models but not $GF_{AF}$ models

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[06_determine_levels_of_failed_replicates](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/05_supplement/06_determine_levels_of_failed_replicates.ipynb)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;figure out if there is any commonality regarding the simulation levels for the 269 replicates that died when using 20k loci encoded as individual genotypes

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[07_conceptual_fig](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/05_supplement/07_conceptual_fig.ipynb)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;create a 4x4 conceptual figure that shows a scatter plot between offset and fitness for seed 1231422. The 4x4 is evaluation level (rows - ind or pop) and input data (columns - AF or geno)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[08_lotterhos_env_figs](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/05_supplement/08_lotterhos_env_figs.ipynb)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;recreate the environmental figures from Lotterhos 2023

### 06_calc_af_ind

Calulate $GF_{AF, ind}$ performance

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[00_af_ind_performance.ipynb](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/06_calc_af_ind/00_af_ind_performance.ipynb)

### 07_continuous_space_sims

Set up datasets from the continuous space simulation for GF evaluation: $GF_{geno, ind(ind-env)}$ and $GF_{geno, ind(pop-env)}$

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[00_train_GF_ind-level_envs_and_pop-level_envs.ipynb](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/07_continuous_space_sims/00_train_GF_ind-level_envs_and_pop-level_envs.ipynb)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;set up the files that will be used to train gradientForests where the individuals of the same 'population' are not all assigned the same environmental values
submit GF training jobs

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[01_estimate_fitness.ipynb](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/07_continuous_space_sims/01_estimate_fitness.ipynb)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;estimate fitness of individuals in the common garden on the spatially continuous simulation landscape

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[02_validate_continuous_sims_offset.ipynb](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/07_continuous_space_sims/02_validate_continuous_sims_offset.ipynb)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;validate (calculate kendall's tau) between off set and fitness in the common gardens on the spatially continuous landscape

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[03_env_PCA_map.ipynb](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/07_continuous_space_sims/03_env_PCA_map.ipynb)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;create figure of environmental values for continuous space simulation by color-coding the loadings of the first three principal component axes of PCA-transformed environmental values.

---

`pythonimports` in notebooks can be found here: https://github.com/brandonlind/pythonimports















