#!/bin/bash

# make AF files from CRISPRESSO files
workdir=$1
summary=$2

for d in `ls -d $workdir/CRISPResso_on*/`; 
do 
	echo $d;
	name=${d/_IGO_*/}; 
	name=${name/$workdir\/CRISPResso_on_/}; 
	name=${name/\//}; 
	echo $name; 

	# create allele frequency table
	# read in Alleles_frequency_table_around_sgRNA_ATCCTGTCCCTAGTGGCCCC.txt
	# ignore first line
	# if Unedited column is True, assign Read and Reference columns to Non-Indel and Non-Indel as key, else, assing Aligned_Sequence to Read, Reference_Sequence to Reference as key
	# count number of reads for each key
	# count total number of reads for each sample
	awk -v sample=$name 'BEGIN{FS=OFS="\t"; PROCINFO["sorted_in"]="@val_num_desc";}{if(NR==1) next; if($3=="True") key="Non-Indel\tNon-Indel"; else key=$1"\t"$2; reads[key] += $7; total += $7; if(key!="Non-Indel\tNon-Indel") mutated += $7}END{print "Sample\tRead\tReference\tNumber\tTotal\tAF\tTotal-Indels\tFraction of Indel Reads"; for(key in reads) print sample, key, reads[key], total, reads[key]/total, mutated, reads[key]/mutated}' $d/Alleles_frequency_table_around_sgRNA_ATCCTGTCCCTAGTGGCCCC.txt > $workdir/$name-AF-table.txt; 

done


# get the Non-Indel, NHEJ, MMEJ, CTG Insertion, GAT Insertion, ATG Insertion
for f in `ls $workdir/*-AF-table.txt`;
do
        awk 'BEGIN{FS=OFS="\t"}{if(NR==1) next; if($2=="Non-Indel" && $3=="Non-Indel") print $0, "Non-Indel"; else if($2=="TCTGTCACCAATC-TGTCCCTAGTGGCCCCACTGTGGGGT" && $3=="TCTGTCACCAATCCTGTCCCTAGTGGCCCCACTGTGGGGT") print $0, "NHEJ"; else if($2=="TCTGTC------------CCTAGTGGCCCCACTGTGGGGT" && $3=="TCTGTCACCAATCCTGTCCCTAGTGGCCCCACTGTGGGGT") print $0, "MMEJ"; else if($2=="GTCACCAATCGATCTGTCCCTAGTGGCCCCACTGTGGGGT" && $3=="GTCACCAATC---CTGTCCCTAGTGGCCCCACTGTGGGGT") print $0, "GAT Insertion"; else if($2=="GTCACCAATCCTGCTGTCCCTAGTGGCCCCACTGTGGGGT" && $3=="GTCACCAATCCTG---TCCCTAGTGGCCCCACTGTGGGGT") print $0, "CTG Insertion"; else if($2=="GTCACCAATCATGCTGTCCCTAGTGGCCCCACTGTGGGGT" && $3=="GTCACCAATC---CTGTCCCTAGTGGCCCCACTGTGGGGT") print $0, "ATG Insertion";}' $f >> $summary;

done
