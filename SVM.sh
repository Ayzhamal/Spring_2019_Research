#!/bin/bash
#SBATCH --job-name=test_job_name
#SBATCH --output=test.o_SVM.txt #output of your pogram prints here
#SBATCH --mail-user=zhamangaraea1@gator.uhd.edu #email
#SBATCH --error=test.e_SVM.txt #file where any error will be written
#SBATCH --mail-type=ALL

#Rscript Logistic_Regression_Angie_cluster.R
#Rscript Random_Forest_Angie_cluster.R
Rscript SVM_Angie_cluster.R
#Rscript Boosted_Tree_Angie_cluster.R
#Rscript 
#python test.py

