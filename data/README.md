# spatially continuous evaluation dataframes

|    column      | python data type | description |
|:---------|:--------|:--------|
| num_loci | int64   | number of loci provided to train the model (['500' '5000' '10000' '20000']) |
| level    | object  | 'ind-env' if $GF_{geno, ind(ind-env)}$ or 'pop-env' if $GF_{geno,ind(ind-env)}$ |
| garden   | int64   | garden ID on the landscape {1..100} |
| score    | float64 | kendall's tau performance |
| program  | object  | always 'GF' |

- spatially_continuous/validation.txt.gz
    - evaluation results for the continuous space simulation. This includes results from both individual-level environmental training data $GF_{geno, ind}(ind-env)$, as well as population-level environmental training data, $GF_{geno, ind}(pop-env)$
    - this file was created and saved in [07_continuous_space_sims/02_validate_continuous_sims_offset.ipynb](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/07_continuous_space_sims/02_validate_continuous_sims_offset.ipynb)

# spatially discrete evaluation dataframes
All of the following columns are found within the .txt.gz files found within spatially_dicrete subdirectories. The exception is that the ['demography', 'final_la_bin', 'seed_garden'] columns are not found within $GF_{geno, ind}$ or $GF_{AF, ind}$ results dataframes.

### column descriptions

| column name | python data type | description |
|:-----------------|:--------|:---|
| garden           | int64   | garden ID on the landscape {1..100} |
| score            | float64 | kendall's tau performance |
| final_LA         | float64 | $LA_\Delta SA$ - the degree of local adaptation in the meta population |
| glevel           | object  | polygenicity level (['highly-polygenic', 'mod-polygenic', 'oligogenic']) |
| plevel           | object  | number of traits under selection (['2-trait']) |
| pleio            | object  | pleiotropy (['no pleiotropy', 'pleiotropy']) |
| slevel           | object  | equality of selection strength on both traits (['equal-S' 'unequal-S']) |
| landscape        | object  | spatially discrete environmental arrangement ['Est-Clines' 'SS-Clines' 'SS-Mtn'] |
| popsize          | object  | pattern of population size clines (['N-cline-center-to-edge' 'N-cline-N-to-S' 'N-equal' 'N-variable']) |
| migration        | object  | pattern of migration (['m-constant' 'm-breaks' 'm-variable']) |
| simulation_level | object  | underscore-joined glevel_plevel_pleio_landscape_migration |
| rep              | object  | name for the 1st, 2nd, and 3rd group of 180 simulation levels ['0-225' '225-450' '450-675'] |
| num_loci         | object  | number of loci provided to train the model (['500' '5000' '10000' '20000']) |
| seed             | int64   | unique simulation ID, N = 540 = 180 * 3 |
| demography       | object  | underscore-joined popsize_migration |
| source           | object  | 'ind' if $GF_{geno, pop}$, 'pooled' if $GF_{AF, pop}$, 'af,ind' if $GF_{AF, ind}$, 'ind-avg' if $GF_{geno, ind}$ |
| offset_level     | object  | underscore-joined simulation_level_rep_garden |
| seed_garden      | object  | underscore-joined seed_garden |
| final_la_bin     | object  | bin of final_LA (['0.42 < LA ≤ 0.58' '0.27 < LA ≤ 0.42'])  |

- spatially_discrete/af_ind/af_ind_results.txt.gz
    - evaluation results for the spatially discrete evaluation using allele frequencies and individual-level fitnesses, $GF_{AF, ind}$
    - this file was created and saved in [06_calc_af_ind/00_af_ind_performance.ipynb](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/06_calc_af_ind/00_af_ind_performance.ipynb)
- spatially_discrete/af_pop/pooled_performance.txt.gz
    - evaluation results for the spatially discrete evaluation using allele frequencies and population-mean fitnesses, $GF_{AF, pop}$
    - this file was created and saved in [02_pooled_runs/03_gather_pooled_scores.ipynb](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/02_pooled_runs/03_gather_pooled_scores.ipynb)
    - within notebooks this file is often loaded, concatenated with `ind-averaged_results.txt`, via the `runtime_API.load_results` function.
- spatially_discrete/geno_ind/ind-averaged_results.txt.gz
    - evaluation results for the spatially discrete evaluation using individual genotypes and individual-level fitnesses, $GF_{geno, ind}$
    - this file was created and saved in [03_calculate_geno-ind_performance/00_calc_geno-ind_performance.ipynb](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/03_calculate_geno-ind_performance/00_calc_geno-ind_performance.ipynb)
- spatially_discrete/geno_pop/ind_performance.txt.gz
    - evaluation results for the spatially discrete evaluation using individual genotypes and individual-level fitnesses, $GF_{geno, pop}$
    - this file was created and saved in [01_individual_runs/03_gather_individual_scores.ipynb](https://nbviewer.org/github/brandonlind/geno_af_gradient_forests/blob/main/01_individual_runs/03_gather_individual_scores.ipynb)
    - within notebooks this file is often loaded, concatenated with `pooled_performance.txt`, via the `runtime_API.load_results` function.
