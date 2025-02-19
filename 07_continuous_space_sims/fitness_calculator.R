#------------------------------------------------------------------------------------------------#
# Calculate fitness for all populations to a garden environment with the six optima from
#    continuous space simulations.
#
# Usage
# -----
# conda activate MVP_env_R4.0.3
# Rscript fitness_calculator.R seed output_file opt0 opt1 opt2 opt3 opt4 opt5
# 
# Parameters
# ----------
# seed - not used, kept for comparing to previous script from Lind & Lotterhos 2024 Mol Ecol Res
# output_file - where to save fitness calculation
# opt1 - temp or MTWetQ optimum for the common garden
# opt0 - (optional) Env2 or MAT optimum for the common garden
# opt2 - MTDQ
# opt3 - PDM
# opt4 - PwarmQ
# opt5 - PWM
#
# Notes
# -----
# `individual_data.txt` is available in the tutorial from Lotterhos 2023 doi: 10.1073/pnas.2220313120
#	the original data, formatted and edited in 07_continuous_space_sims/00_train_GF_ind-level_envs_and_pop-level_envs.ipynb,
#	can be found here:  https://marineomics.github.io/RDAtraitPredictionTutorial.html
# 
#------------------------------------------------------------------------------------------------#
library('mvtnorm')
library(progress)
len = length

print(sessionInfo())

args = commandArgs(trailingOnly=TRUE)

seed = as.character(args[1])
output_file = args[2]
opt0 = as.numeric(args[3])  # mat
opt1 = as.numeric(args[4])  # MTWetQ
opt2 = as.numeric(args[5])  # MTDQ
opt3 = as.numeric(args[6])  # PDM
opt4 = as.numeric(args[7])  # PwarmQ
opt5 = as.numeric(args[8])  # PWM

# find phenos from complex sims (Lotterhos 2023 PNAS)
    # from 06_run_time_project/07_continuous_space_sims/00_train_GF_ind-level_envs_and_pop-level_envs.ipynb
phenodata = read.table('/work/lotterhos/brandon/continuous_space_runtime/individual_data.txt',
                       header=T)
rownames(phenodata) = phenodata[,'samp']

traits = c("phenotype1_mat",
           "phenotype2_MTWetQ",
           "phenotype3_MTDQ",
           "phenotype4_PDM",
           "phenotype5_PwarmQ",
           "phenotype6_PWM")

# create empty dataframe (1 row because input args #2 and #3 specify optima for a single garden ID)
fitness = data.frame(matrix(nrow=1, ncol=1000))
colnames(fitness) = phenodata[, 'samp']

# fill in empty dataframe with fitness of each transplant into garden with optima opt0 [opt1] etc
# for (transplant_ID in 1:100){
    # complex sims, six selective environments
#         phenos = cbind(phenodata[phenodata$subpopID == transplant_ID, traits[1]],
#                        phenodata[phenodata$subpopID == transplant_ID, traits[2]],
#                        phenodata[phenodata$subpopID == transplant_ID, traits[3]],
#                        phenodata[phenodata$subpopID == transplant_ID, traits[4]],
#                        phenodata[phenodata$subpopID == transplant_ID, traits[5]],
#                        phenodata[phenodata$subpopID == transplant_ID, traits[6]])
        
phenos = phenodata[, traits]  # get just the phenotype data
stopifnot(nrow(phenos) == 1000)  # expect 1000 pops

print(c('dim(phenos) = ', dim(phenos)))

print(phenos)

sigma_k = 2
fitness_varcov = matrix(
    c(sigma_k, 0, 0, 0, 0, 0,
      0, sigma_k, 0, 0, 0, 0,
      0, 0, sigma_k, 0, 0, 0,
      0, 0, 0, sigma_k, 0, 0,
      0, 0, 0, 0, sigma_k, 0,
      0, 0, 0, 0, 0, sigma_k
     ),
    nrow=6,
    ncol=6
)

fitness_norm = dmvnorm(c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                       c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                       fitness_varcov)
print('fitness_norm')
print(fitness_norm)
print(c(opt0, opt1, opt2, opt3, opt4, opt5))

fits = dmvnorm(
    phenos,
    c(opt0, opt1, opt2, opt3, opt4, opt5),
    sigma=fitness_varcov
) / fitness_norm

# }

#     fits = round(fits, 2)
    
fitness[1, ] = fits
    
# }

write.table(fitness, output_file, sep='\t', row.names=T, col.names=T)

cat(sprintf('wrote fitness to %s', output_file))

