#!/bin/bash -e

#$ -N DE
#$ -cwd
#$ -l mem_free=50G,h_vmem=50G
#$ -m be
#$ -M abby.spangler@nih.gov
#$ -o ./analysis/code/06_DE/01_DE_all_genes.txt
#$ -e ./analysis/code/06_DE/01_DE_all_genes.txt


echo "**** Job starts ****"
date

echo "**** LOCUS info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
pwd

module load R/4.2.1

Rscript ./analysis/code/06_DE/01_DE_all_genes.R

echo "**** Job ends ****"
date
