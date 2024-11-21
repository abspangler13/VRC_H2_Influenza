#!/bin/bash -e

#$ -N build_seurat_run4
#$ -cwd
#$ -l h_vmem=20G
#$ -m be
#$ -M abby.spangler@nih.gov
#$ -o ./analysis/code/01_build_seurat/build_seurat_run4.txt
#$ -e ./analysis/code/01_build_seurat/build_seurat_run4.txt


echo "**** Job starts ****"
date

echo "**** LOCUS info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
pwd

module load R/4.1.3

Rscript ./analysis/code/01_build_seurat/01_build_seurat_run4.R

echo "**** Job ends ****"
date