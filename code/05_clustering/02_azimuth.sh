#!/bin/bash -e

#$ -N azimuth
#$ -cwd
#$ -l mem_free=50G,h_vmem=50G
#$ -m be
#$ -M abby.spangler@nih.gov
#$ -o ./analysis/code/05_clustering/02_azimuth.txt
#$ -e ./analysis/code/05_clustering/02_azimuth.txt


echo "**** Job starts ****"
date

echo "**** LOCUS info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
pwd

module unload R/4.2.1
module load R/4.1.3

Rscript ./analysis/code/05_clustering/02_azimuth.R

echo "**** Job ends ****"
date