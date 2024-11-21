#!/bin/bash -e

#$ -N find_bcells_all
#$ -cwd
#$ -l mem_free=100G,h_vmem=100G
#$ -m be
#$ -M abby.spangler@nih.gov
#$ -o ./analysis/code/02_dsb_normalization/02_find_bcells_all.txt
#$ -e ./analysis/code/02_dsb_normalization/02_find_bcells_all.txt


echo "**** Job starts ****"
date

echo "**** LOCUS info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
pwd

module load R/4.1.3

Rscript ./analysis/code/02_dsb_normalization/02_find_bcells_all.R

echo "**** Job ends ****"
date