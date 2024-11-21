#!/bin/bash -e

#$ -N clonal_cluster_overlap
#$ -cwd
#$ -l mem_free=50G,h_vmem=50G
#$ -m be
#$ -M abby.spangler@nih.gov
#$ -o ./analysis/code/08_cluster_clonal_overlap/01_clonal_cluster_overlap.txt
#$ -e ./analysis/code/08_cluster_clonal_overlap/01_clonal_cluster_overlap.txt


echo "**** Job starts ****"
date

echo "**** LOCUS info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
pwd

module load R/4.2.1

Rscript ./analysis/code/08_cluster_clonal_overlap/clonal_cluster_overlap.R

echo "**** Job ends ****"
date
