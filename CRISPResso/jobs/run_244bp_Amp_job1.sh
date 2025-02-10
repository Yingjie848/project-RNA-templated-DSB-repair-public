#!/bin/bash

fastq_r1=(
  fastq_files/30-645893503/30-645893503/00_fastq/cleaned/C185-Hek293T-siNT-C_R1_001.fastq.gz
  fastq_files/30-645893503/30-645893503/00_fastq/cleaned/C185-Hek293T-siNT-pMJ119-2_R1_001.fastq.gz
  fastq_files/30-645893503/30-645893503/00_fastq/cleaned/C195a-Hek293T-siNT-pMJ119-1_R1_001.fastq.gz
  fastq_files/30-645893503/30-645893503/00_fastq/cleaned/C195a-Hek293T-siNT-pMJ119-2_R1_001.fastq.gz
  fastq_files/30-645893503/30-645893503/00_fastq/cleaned/C197a-Hek293T-siNT-pMJ119-1_R1_001.fastq.gz
  fastq_files/30-645893503/30-645893503/00_fastq/cleaned/C197a-Hek293T-siNT-pMJ119-2_R1_001.fastq.gz
)

for f in ${fastq_r1[@]}; do
	echo "sh run_CRISPResso.244bp_Amp.sh $f output_30-645893503_cleaned" | bash
done

