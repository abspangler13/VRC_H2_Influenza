#!/bin/bash -e

#$ -N vdj2
#$ -cwd
#$ -l mem_free=100G,h_vmem=100G
#$ -m be
#$ -M abby.spangler@nih.gov
#$ -o ./analysis/code/03_vdj/02_vdj_all.txt
#$ -e ./analysis/code/03_vdj/02_vdj_all.txt


echo "**** Job starts ****"
date

echo "**** LOCUS info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
pwd

module load R/4.2.1

Rscript ./analysis/code/03_vdj/02_vdj_all.R

echo "**** Job ends ****"
date