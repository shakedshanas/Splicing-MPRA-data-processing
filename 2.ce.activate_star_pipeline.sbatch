#!/bin/bash
 
#SBATCH --job-name=pm455b
#SBATCH --array=0-9608
#SBATCH --ntasks=8
#SBATCH --output=/lustre1/home/martinmikl/mmikl/data/Admera/August2024/CEpromoters/pm455b/%A_%a.output
#SBATCH --partition=hive1d

read1=/lustre1/home/martinmikl/sshanas/demult_Admera_August2024/22007FL-07-01-80_S239_L007_R2_001.fastq.gz
read2=/lustre1/home/martinmikl/sshanas/demult_Admera_August2024/22007FL-07-01-80_S239_L007_R1_001.fastq.gz

python ce.pipeline.py --libseq $SLURM_ARRAY_TASK_ID -r1 $read1 -r2 $read2 --pwrdir /lustre1/home/martinmikl/mmikl/data/Admera/August2024/CEpromoters/pm455b --caslib no

