# Description
A Python program (CategorizeDeletions.py) was used to identify whole intron deletions (WID). 

The python package 'intervaltree' is a dependency for this program. The program was run using Python 3.10. The paths in the main() function should be adjusted to reflect local paths.

It takes 2 files as input: 1. gene list bed file defining exon boundaries, 2. maf file with deletions called in the IMPACT cohort

The program builds trees of exon boundaries and then intersects the deletion with those trees to see if there is an overlap with any exon. Whole intron deletion is identified by matching exon boundaries and the start and end of the deletion, with some acceptable margin. The ouput is a maf file with additional columns idnefying WIDs.



