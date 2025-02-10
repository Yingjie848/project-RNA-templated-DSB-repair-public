# pylint: skip-file
'''
This script categorizes deletions based on their location in the genome.

Created on Apr 22, 2022

@author: Juber Patel
'''

#import portion as p

import os
import sys
import intervaltree as it



def main():
    
    d = "change to the directory where the input files are located"    
    os.chdir(d)

    exonsFile = "input/hg19-gene-list.bed"
    deletionsFile = "input/data_mutations_extended.txt"
    outFile = deletionsFile.replace(".txt", "-categorized.txt")
    
    genes = set()
    max = {}
    min = {}
    chr = {}
    
    intervalTrees = {}
    geneTrees = {}
    
    # create intervalTrees
    for i in range(1, 23):
        intervalTrees[str(i)] = it.IntervalTree()
        geneTrees[str(i)] = it.IntervalTree()
        
    intervalTrees["X"] = it.IntervalTree()
    intervalTrees["Y"] = it.IntervalTree()
    
    geneTrees["X"] = it.IntervalTree()
    geneTrees["Y"] = it.IntervalTree()
    
    
    print("Building interval trees...")

    seen = set()
    counter = 0
    intervalPadding = 12

    # read exons file
    with open(exonsFile) as f:
        for line in f:
            tokens = line.rstrip("\n").split("\t")
            
            if tokens[4].startswith("snp"):
                continue
            
            if tokens[4].startswith("MIR"):
                continue
            
            if "exon" not in tokens[4]:
                continue
            
            if tokens[0] not in intervalTrees:
                continue
            
            # add the gene to the gene set
            info = tokens[4].split("_")
            gene = info[0]
            name = info[0] + "_" + info[1]

            
            if name in seen:
                continue

            seen.add(name)
            
            genes.add(gene)
            chr[gene] = tokens[0]
            
            start = int(tokens[1]) - intervalPadding
            end = int(tokens[2]) + intervalPadding
            intervalTrees[tokens[0]][start:end] = name
            
            # update gene range
            if gene not in min or start < min[gene]:
                min[gene] = start
                
            if gene not in max or end > max[gene]:
                max[gene] = end
            
            counter += 1
            if counter % 100000 == 0:
                print(counter)
            
    
    f.close()

    seen.clear()
    
    
    # build gene trees
    for gene in min:
        c = chr[gene]
        geneTrees[c][min[gene]:max[gene]] = gene
        
    
    
    '''
    #13    26975485    26975602
    c = "13"
    position = 26984149
    overlap = sorted(geneTrees[c][position])
    
    print(overlap)
    sys.exit()
    '''
    
    print("Categorizing deletions...")
    
    
    # open output file
    w = open(outFile, "w")
    
    margin = 10
    counter = 0
    # read the deletions file
    with open(deletionsFile) as f:
        for line in f:
            
            if counter == 0:
                w.write(line.rstrip("\n") + "\tDeletion_Category\tLeft_Margin\tRight_Margin\n")
                counter += 1
                continue
            
            
            category = "NA\tNA\tNA"
            tokens = line.rstrip("\n").split("\t")
            delStart = int(tokens[5])
            delEnd = int(tokens[6])
            length = (delEnd - delStart) + 1
            
            if length < 10:
                category = "NA\tNA\tNA"
                w.write(line.rstrip("\n") + "\t" + category + "\n")
                counter += 1
                if counter % 10000 == 0:
                    print(counter)
                continue

            t = intervalTrees[tokens[4]]
            
            
            # intron deletion
            found = False
            for i in range(margin*-1, margin+1):
                if found:
                    break
                
                for j in range(margin*-1, margin+1):
                    start = delStart + i
                    end = delEnd - j + 1
                    overlap = sorted(t[start:end])
                    #print(str(start) + "\t" + str(end) + "\t" + str(overlap))
                    
                    if len(overlap) != 2:
                        continue
                    
                    interval = overlap[0]
                    interval1 = overlap[1]
                    
                    leftMargin = delStart-(interval.end-(intervalPadding-1))
                    rightMargin = delEnd-(interval1.begin+(intervalPadding-1))

                    # found a whole intron deletion!!!
                    if abs(leftMargin) <= margin and abs(rightMargin) <= margin:
                        category = "intron deletion\t" + str(leftMargin) + "\t" + str(rightMargin)
                        w.write(line.rstrip("\n") + "\t" + category + "\n")
                        counter += 1
                        found = True
                        break
            
            if not found:
                category = "NA\tNA\tNA"
                counter += 1
                w.write(line.rstrip("\n") + "\t" + category + "\n")

            if counter % 10000 == 0:
                print(counter)
            


                    
            '''
            # exon deletion
            found = False
            for i in range(-2, 3):
                if found:
                    break
                
                for j in range(-2, 3):
                    start = int(tokens[5]) + i
                    end = int(tokens[6]) - j + 1
                    overlap = sorted(t[start:end])
                    #print(str(start) + "\t" + str(end) + "\t" + str(overlap))
                    
                    if len(overlap) != 1:
                        continue
                    
                    interval = overlap[0]
                    if abs(int(tokens[5])-interval.begin) <= 2 and abs(int(tokens[6])-interval.end) <= 2:
                        category = "exon deletion"
                        w.write(line.rstrip("\n") + "\t" + category + "\n")
                        counter += 1
                        found = True
                        break
            if found:
                continue
            
            
            # genes absent in the annotation: assume exonic variants
            if tokens[0] not in genes:
                category = "exonic"
                w.write(line.rstrip("\n") + "\t" + category + "\n")
                counter += 1
                continue
                
            
            start = int(tokens[5])
            end = int(tokens[6]) + 1
            overlap = sorted(t[start:end])
            
             
            # if the overlap is empty, it's an intronic or intergenic deletion
            if len(overlap) == 0:
                
                # decide if intronic or intergenic
                p = sorted(geneTrees[tokens[4]][start])
                
                if len(p) == 0:
                    category = "intergenic"
                else:
                    category = "intronic"
                
                w.write(line.rstrip("\n") + "\t" + category + "\n")
                counter += 1
                continue
                
            
            interval = overlap[0]
            # if contained within the first interval, it's an exonic deletion
            if interval.begin < start and interval.end > end:
                category = "exonic"
                w.write(line.rstrip("\n") + "\t" + category + "\n")
                counter += 1
                continue
            
            # if bigger than first interval, it is intron-exon spanning
            if start <= interval.begin or end >= interval.end:
                category = "intron-exon spanning"
                w.write(line.rstrip("\n") + "\t" + category + "\n")
                counter += 1
                continue
            
            '''
            
               
                
            #print(tokens[0] + "\t" + tokens[4] + "\t" + tokens[5] + "\t" + tokens[6])
            #print(overlap)
            
            
                
            
    f.close()
    w.close()
    
    
    
    
    
    
    
    



if __name__ == '__main__':
    main()
