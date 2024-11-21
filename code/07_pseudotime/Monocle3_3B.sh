#!/bin/bash -e

#$ -N pseudotime3B
#$ -cwd
#$ -l mem_free=50G,h_vmem=50G
#$ -m be
#$ -M abby.spangler@nih.gov
#$ -o ./analysis/code/07_pseudotime/Monocle3_3B.txt
#$ -e ./analysis/code/07_pseudotime/Monocle3_3B.txt


echo "**** Job starts ****"
date

echo "**** LOCUS info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
pwd

module load R/4.2.1

Rscript ./analysis/code/07_pseudotime/Monocle3_3B.R

echo "**** Job ends ****"
date
