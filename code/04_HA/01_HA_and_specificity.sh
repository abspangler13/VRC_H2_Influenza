#!/bin/bash -e

#$ -N HA
#$ -cwd
#$ -l mem_free=20G,h_vmem=20G
#$ -m be
#$ -M abby.spangler@nih.gov
#$ -o ./analysis/code/04_HA/01_HA_and_specificity.txt
#$ -e ./analysis/code/04_HA/01_HA_and_specificity.txt


echo "**** Job starts ****"
date

echo "**** LOCUS info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
pwd

module load R/4.2.1

Rscript ./analysis/code/04_HA/01_HA_and_specificity.R

echo "**** Job ends ****"
date