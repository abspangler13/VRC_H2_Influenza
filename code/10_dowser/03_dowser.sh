#!/bin/bash -e

#$ -N dowser
#$ -cwd
#$ -l mem_free=50G,h_vmem=50G
#$ -m be
#$ -M abby.spangler@nih.gov
#$ -o analysis/code/10_dowser/logs/dowser.$TASK_ID.txt
#$ -e analysis/code/10_dowser/logs/dowser.$TASK_ID.txt
#$ -t 1-5


echo "**** Job starts ****"
date

echo "**** LOCUS info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
pwd

module load immcantation/4.4.0

singularity exec -B /hpcdata,/nethome,/sysapps /sysapps/cluster/software/immcantation/4.4.0/immcantation_suite-4.4.0.sif R -f /hpcdata/vrc_vip/Abby/Experiment_316/analysis/code/10_dowser/03_dowser.R

echo "**** Job ends ****"
date
