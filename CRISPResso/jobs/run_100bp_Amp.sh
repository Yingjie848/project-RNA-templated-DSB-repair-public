#!/bin/bash

fastq_r1=(
        fastq_files_juber/fastq_files/cleaned/C172_Hek293T_C_IGO_10465_8_S87_R1_001.fastq.gz
        fastq_files_juber/fastq_files/cleaned/C182_Hek293T_C_IGO_10489_8_S74_R1_001.fastq.gz
        fastq_files_juber/fastq_files/cleaned/C184_Hek293T_C_IGO_10489_20_S63_R1_001.fastq.gz
)

for f in ${fastq_r1[@]}; do
	echo "sh run_CRISPResso.100bp_Amp.sh $f output_juber_cleaned_fastq_100bp_window_size_20bp" | bash
done