#!/bin/bash -e

#$ -N GO
#$ -cwd
#$ -l mem_free=50G,h_vmem=50G
#$ -m be
#$ -M abby.spangler@nih.gov
#$ -o ./analysis/code/06_DE/02_GO.txt
#$ -e ./analysis/code/06_DE/02_GO.txt


echo "**** Job starts ****"
date

echo "**** LOCUS info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
pwd

module load R/4.2.1

Rscript ./analysis/code/06_DE/02_GO.R

echo "**** Job ends ****"
date
