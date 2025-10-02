#!/bin/bash

read1=/lustre1/home/martinmikl/sshanas/demult_Admera_August2024/22007FL-07-01-80_S239_L007_R2_001.fastq.gz
read2=/lustre1/home/martinmikl/sshanas/demult_Admera_August2024/22007FL-07-01-80_S239_L007_R1_001.fastq.gz

sbatch -J pm455b_bc -o pm455b_bc.out -e pm455b_bc.err --wrap "python 1.create_bcread_dict.py -r1 $read1 -r2 $read2 -o pm455b"
