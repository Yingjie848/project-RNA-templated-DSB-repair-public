#!/bin/bash

source activate ~/my_apps/miniconda3
conda activate crispresso2_env

fastq_r1=$1
fastq_r2=$(echo $fastq_r1 | sed 's/_R1_/_R2_/')
name=$(echo $fastq_r1 | sed 's/.*\///' | sed 's/_R1_.*//')

echo $name
echo $fastq_r1
echo $fastq_r2

output=$2

echo $output

mkdir -p $output log

amplicon_seq="GGCCTAAGGATGGGGCTTTTCTGTCACCAATCCTGTCCCTAGTGGCCCCACTGTGGGGTGGAGGGGACAGATAAAAGTACCCAGAACCAGAGCCACATTA"

sub_now.sh 12:00 100 1 log/log_$name.err log/log_$name.out "CRISPResso --fastq_r1 $fastq_r1 --amplicon_seq $amplicon_seq --guide_seq ATCCTGTCCCTAGTGGCCCC --quantification_window_center -10 --quantification_window_size 20 --plot_window_size 20 --base_editor_output --ignore_substitutions --write_detailed_allele_table --min_frequency_alleles_around_cut_to_plot 0.01 --allele_plot_pcts_only_for_assigned_reference --annotate_wildtype_allele \"**\" --max_rows_alleles_around_cut_to_plot 200 -o $output --name $name"
