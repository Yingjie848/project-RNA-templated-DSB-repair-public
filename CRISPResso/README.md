# Processing steps:
1. To examine sample barcode, first run fetch_sample_barcode.sh to get the top 1 barcode as the sample barcode, then run RemoveBadBarcodes.py to remove reads not matched with sample barcode. 

2. Run CRISPResso

Run run_CRISPResso.PE150.sh for paired-end 150bp sequencing reads, and run_CRISPResso.PE100.sh for paired-end 100bp sequencing reads. For paired-end 100bp sequencing reads, the middle part of amplicon sequence is not covered, therefore, CRISPResso only uses the first read.

3. Create allele fraction table and count the Non-Indel, NHEJ, MMEJ, GAT Insertion, ATG Insertion using create_AF_table.PE150.sh and create_AF_table.PE100.sh 

4. Calculate proportion of 3bp insertions, 2bp insertions, and 4bp insertions around cutting site: analyze_2bp_3bp_4bp_insertions.R
