# Processing steps:
1. To examine sample barcode, first run fetch_sample_barcode.sh to get the top 1 barcode as the sample barcode, then run RemoveBadBarcodes.py to remove reads not matched with sample barcode. 

2. Run CRISPResso

Run run_CRISPResso.244bp_Amp.sh for 244bp amplicon sequence.

Run run_CRISPResso.100bp_Amp.sh for 100bp amplicon sequence, because the read length in these samples is shorter.

3. Create allele fraction table and count the Non-Indel, NHEJ, MMEJ, CTG Insertion, GAT Insertion, ATG Insertion using create_AF_table.100bp.sh and create_AF_table.244bp.sh 

4. Calculate proportion of 3bp insertions, 2bp insertions, and 4bp insertions around cutting site: analyze_2bp_3bp_4bp_insertions.R
