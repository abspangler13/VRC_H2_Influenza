#!/bin/bash -e

#$ -N fgsea
#$ -cwd
#$ -l mem_free=50G,h_vmem=50G
#$ -m be
#$ -M abby.spangler@nih.gov
#$ -o ./analysis/code/13_pb_fgsea/04_fgsea_loop.txt
#$ -e ./analysis/code/13_pb_fgsea/04_fgsea_loop.txt


echo "**** Job starts ****"
date

echo "**** LOCUS info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
pwd

module load R/4.1.3

Rscript ./analysis/code/13_pb_fgsea/04_fgsea_loop.R

echo "**** Job ends ****"
date
