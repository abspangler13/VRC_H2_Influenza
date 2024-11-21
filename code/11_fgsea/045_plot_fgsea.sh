#!/bin/bash -e

#$ -N plot_fgsea
#$ -cwd
#$ -l mem_free=50G,h_vmem=50G
#$ -m be
#$ -M abby.spangler@nih.gov
#$ -o ./analysis/code/11_fgsea/045_plot_fgsea.txt
#$ -e ./analysis/code/11_fgsea/045_plot_fgsea.txt


echo "**** Job starts ****"
date

echo "**** LOCUS info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
pwd

module load R/4.1.3

Rscript ./analysis/code/11_fgsea/045_plot_fgsea.R

echo "**** Job ends ****"
date
