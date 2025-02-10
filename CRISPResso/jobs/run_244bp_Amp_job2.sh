
#!/bin/bash

fastq_r1=(
    fastq_files_juber/fastq_files/cleaned/C195a-Hek293T-siNT-C_R1_001.fastq.gz
    fastq_files_juber/fastq_files/cleaned/C197a-Hek293T-siNT-C_R1_001.fastq.gz
    fastq_files_juber/fastq_files_cleaned/C279a_Hek293T_siNT_C_IGO_15060_1_R1_001.fastq.gz
    fastq_files_juber/fastq_files_cleaned/C279a_Hek293T_siNT_pMJ119_IGO_15060_3_R1_001.fastq.gz
    fastq_files_juber/fastq_files_cleaned/C279a_Hek293T_siRev3_pMJ119_IGO_15060_6_R1_001.fastq.gz
    fastq_files_juber/fastq_files_cleaned/C285a_Hek293T_siNT_C_IGO_15140_9_R1_001.fastq.gz
    fastq_files_juber/fastq_files_cleaned/C285a_Hek293T_siNT_pMJ119_IGO_15140_10_R1_001.fastq.gz
    fastq_files_juber/fastq_files_cleaned/C285a_Hek293T_siRev3_pMJ119_IGO_15140_12_R1_001.fastq.gz
    fastq_files_juber/fastq_files_cleaned/C288a_Hek293T_siNT_C_IGO_15155_9_R1_001.fastq.gz
    fastq_files_juber/fastq_files_cleaned/C288a_Hek293T_siNT_pMJ119_IGO_15155_10_R1_001.fastq.gz
    fastq_files_juber/fastq_files_cleaned/C288a_Hek293T_siRev3_pMJ119_IGO_15155_12_R1_001.fastq.gz
)

for f in ${fastq_r1[@]}; do
	echo "sh run_CRISPResso.244bp_Amp.sh $f output_juber_fastq_cleaned_244bp" | bash
done