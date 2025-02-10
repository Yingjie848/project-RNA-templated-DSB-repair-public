# calculate proportion of 3bp insertions, 2bp insertions, and 4bp insertions around cutting site
rm(list=ls())
library(data.table)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(tidyverse)

# due to alignment issue, gaps are added to alternative positions
# at cutting site
refseq_244bp_3bp          <- "TTTCTGTCACCAATC---CTGTCCCTAGTGGCCCCACTGTGGGGTGGAGG"
refseq_244bp_2bp          <- "TTTTCTGTCACCAATC--CTGTCCCTAGTGGCCCCACTGTGGGGTGGAGG"
refseq_244bp_4bp          <- "TTCTGTCACCAATC----CTGTCCCTAGTGGCCCCACTGTGGGGTGGAGG"

refseq_100bp_3bp          <- "GTCACCAATC---CTGTCCCTAGTGGCCCCACTGTGGGGT"
refseq_100bp_2bp          <- "TGTCACCAATC--CTGTCCCTAGTGGCCCCACTGTGGGGT"
refseq_100bp_4bp          <- "TCACCAATC----CTGTCCCTAGTGGCCCCACTGTGGGGT"

## shift nucleotide on the right to the left

# shift CTG 3bp to the left
refseq_244bp_3bp_shiftCTG <- "TTTCTGTCACCAATCCTG---TCCCTAGTGGCCCCACTGTGGGGTGGAGG"
refseq_100bp_3bp_shiftCTG <- "GTCACCAATCCTG---TCCCTAGTGGCCCCACTGTGGGGT"

# shift C 3bp to the left
refseq_244bp_3bp_shiftC   <- "TTTCTGTCACCAATCC---TGTCCCTAGTGGCCCCACTGTGGGGTGGAGG"
refseq_100bp_3bp_shiftC   <- "GTCACCAATCC---TGTCCCTAGTGGCCCCACTGTGGGGT"

# shift CT 3bp to the left
refseq_244bp_3bp_shiftCT   <- "TTTCTGTCACCAATCCT---GTCCCTAGTGGCCCCACTGTGGGGTGGAGG"
refseq_100bp_3bp_shiftCT   <- "GTCACCAATCCT---GTCCCTAGTGGCCCCACTGTGGGGT"

# shift C 2bp to the left
refseq_244bp_2bp_shiftC          <- "TTTTCTGTCACCAATCC--TGTCCCTAGTGGCCCCACTGTGGGGTGGAGG"
refseq_100bp_2bp_shiftC          <- "TGTCACCAATCC--TGTCCCTAGTGGCCCCACTGTGGGGT"

# shift CT 2bp to the left
refseq_244bp_2bp_shiftCT         <- "TTTTCTGTCACCAATCCT--GTCCCTAGTGGCCCCACTGTGGGGTGGAGG"
refseq_100bp_2bp_shiftCT         <- "TGTCACCAATCCT--GTCCCTAGTGGCCCCACTGTGGGGT"

# shift C, CT, CTG, CTGT 4bp to the left
refseq_244bp_4bp_shiftC          <- "TTCTGTCACCAATCC----TGTCCCTAGTGGCCCCACTGTGGGGTGGAGG"
refseq_244bp_4bp_shiftCT         <- "TTCTGTCACCAATCCT----GTCCCTAGTGGCCCCACTGTGGGGTGGAGG"
refseq_244bp_4bp_shiftCTG        <- "TTCTGTCACCAATCCTG----TCCCTAGTGGCCCCACTGTGGGGTGGAGG"
refseq_244bp_4bp_shiftCTGT       <- "TTCTGTCACCAATCCTGT----CCCTAGTGGCCCCACTGTGGGGTGGAGG"
refseq_100bp_4bp_shiftC          <- "TCACCAATCC----TGTCCCTAGTGGCCCCACTGTGGGGT"
refseq_100bp_4bp_shiftCT         <- "TCACCAATCCT----GTCCCTAGTGGCCCCACTGTGGGGT"
refseq_100bp_4bp_shiftCTG        <- "TCACCAATCCTG----TCCCTAGTGGCCCCACTGTGGGGT"
refseq_100bp_4bp_shiftCTGT       <- "TCACCAATCCTGT----CCCTAGTGGCCCCACTGTGGGGT"


## shift nucleotide on the left to the right

#refseq_244bp_3bp                <- "TTTCTGTCACCAATC---CTGTCCCTAGTGGCCCCACTGTGGGGTGGAGG"

# shift TC 3bp to the right
refseq_244bp_3bp_shiftTCright    <- "TTTCTGTCACCAA---TCCTGTCCCTAGTGGCCCCACTGTGGGGTGGAGG"
refseq_100bp_3bp_shiftTCright    <- "GTCACCAA---TCCTGTCCCTAGTGGCCCCACTGTGGGGT"

# shift C 3bp to the right
refseq_244bp_3bp_shiftCright    <- "TTTCTGTCACCAAT---CCTGTCCCTAGTGGCCCCACTGTGGGGTGGAGG"
refseq_100bp_3bp_shiftCright    <- "GTCACCAAT---CCTGTCCCTAGTGGCCCCACTGTGGGGT"

# shift ATC 3bp to the right
refseq_244bp_3bp_shiftATCright    <- "TTTCTGTCACCA---ATCCTGTCCCTAGTGGCCCCACTGTGGGGTGGAGG"
refseq_100bp_3bp_shiftATCright    <- "GTCACCA---ATCCTGTCCCTAGTGGCCCCACTGTGGGGT"


# shift AATC 4bp to the right
refseq_244bp_4bp_shiftAATCright    <- "TTCTGTCACC----AATCCTGTCCCTAGTGGCCCCACTGTGGGGTGGAGG"
refseq_100bp_4bp_shiftAATCright    <- "TCACC----AATCCTGTCCCTAGTGGCCCCACTGTGGGGT"

# shift ATC 4bp to the right
refseq_244bp_4bp_shiftATCright     <- "TTCTGTCACCA----ATCCTGTCCCTAGTGGCCCCACTGTGGGGTGGAGG"
refseq_100bp_4bp_shiftATCright     <- "TCACCA----ATCCTGTCCCTAGTGGCCCCACTGTGGGGT"

# shift TC 4bp to the right
refseq_244bp_4bp_shiftTCright      <- "TTCTGTCACCAA----TCCTGTCCCTAGTGGCCCCACTGTGGGGTGGAGG"
refseq_100bp_4bp_shiftTCright      <- "TCACCAA----TCCTGTCCCTAGTGGCCCCACTGTGGGGT"

# shift C 4bp to the right
refseq_244bp_4bp_shiftCright       <- "TTCTGTCACCAAT----CCTGTCCCTAGTGGCCCCACTGTGGGGTGGAGG"
refseq_100bp_4bp_shiftCright       <- "TCACCAAT----CCTGTCCCTAGTGGCCCCACTGTGGGGT"


# shift TC 2bp to the right
refseq_244bp_2bp_shiftTCright      <- "TTTTCTGTCACCAA--TCCTGTCCCTAGTGGCCCCACTGTGGGGTGGAGG"
refseq_100bp_2bp_shiftTCright      <- "TGTCACCAA--TCCTGTCCCTAGTGGCCCCACTGTGGGGT"

# shift C 2bp to the right
refseq_244bp_2bp_shiftCright       <- "TTTTCTGTCACCAAT--CCTGTCCCTAGTGGCCCCACTGTGGGGTGGAGG"
refseq_100bp_2bp_shiftCright       <- "TGTCACCAAT--CCTGTCCCTAGTGGCCCCACTGTGGGGT"


# get frequencies of 3bp insertions around cutting site
# cutting site: 
# GGGGCCACTAGGGACAG”GATTGG ( “-cut site)
# CCAATC"CTGTCCCTAGTGGCCCC
# 1) filter AF table by reference sequence with 3bp deletion at the cutting site
# 2) add Event column to the AF table, which indicates the inserted sequences. Upstream and downstream sequences should be the same as the reference sequence.
# 3) save the AF table with Event column to a file, and return an annotated AF table
get_3bp_insertions_244bpAmp <- function(d, sample, reference_sequence){

    #reference_sequence <- "TTTCTGTCACCAATC---CTGTCCCTAGTGGCCCCACTGTGGGGTGGAGG"
    d_inserted_3bp <- d %>% dplyr::filter(Reference==reference_sequence)
    reads_3bp_insertion <- sum(d_inserted_3bp$Number)
    print(reads_3bp_insertion)

    if(reads_3bp_insertion==0){
        return(NULL)
    }

    # add Event column
    assign_event <- function(aligned_sequence, reference_sequence){
        # find gap index in reference sequence
        indexes <- gregexpr("-", reference_sequence)[[1]]

        # find inserted sequences in aligned sequence
        inserted_sequences <- sapply(indexes, function(i){
            substr(aligned_sequence, i, i)
        })
        inserted_sequences <- paste(inserted_sequences, collapse="")

        # get upstream sequences
        ref_upstream_sequences <- substr(reference_sequence, 1, min(indexes)-1)
        ref_upstream_sequences <- paste(ref_upstream_sequences, collapse="")

        upstream_sequences <- substr(aligned_sequence, 1, min(indexes)-1)
        upstream_sequences <- paste(upstream_sequences, collapse="")

        # get downstream sequences
        ref_downstream_sequences <- substr(reference_sequence, max(indexes)+1, nchar(reference_sequence))
        ref_downstream_sequences <- paste(ref_downstream_sequences, collapse="")

        downstream_sequences <- substr(aligned_sequence, max(indexes)+1, nchar(aligned_sequence))
        downstream_sequences <- paste(downstream_sequences, collapse="")

        # compare upstream and downstream sequences
        if(ref_upstream_sequences==upstream_sequences & ref_downstream_sequences==downstream_sequences){
            return(paste(inserted_sequences, " Insertion", sep=""))
        }else{
            return("Other events")
        }

    }

    d_inserted_3bp$Event <- sapply(d_inserted_3bp$Read, function(x){
        Event <- assign_event(x, reference_sequence)
    })

    dir.create(file.path("myoutput",sample,reference_sequence), recursive=TRUE)

    d_inserted_3bp <- d_inserted_3bp %>% dplyr::mutate(Class='3bp Insertion', reference_sequence=reference_sequence)

    d_inserted_3bp %>% fwrite(file.path("myoutput",sample,reference_sequence,"3bp_insertions_AF_table.txt"), sep="\t")

    d_inserted_3bp

}

# get frequencies of 2bp insertions around cutting site
get_2bp_insertions_244bpAmp <- function(d, sample, reference_sequence){

    #reference_sequence <- "TTTTCTGTCACCAATC--CTGTCCCTAGTGGCCCCACTGTGGGGTGGAGG"
    d_inserted_2bp <- d %>% dplyr::filter(Reference==reference_sequence)
    reads_2bp_insertion <- sum(d_inserted_2bp$Number)
    print(reads_2bp_insertion)

    if(reads_2bp_insertion==0){
        return(NULL)
    }

    # add Event column
    assign_event <- function(aligned_sequence, reference_sequence){
        # find gap index in reference sequence
        indexes <- gregexpr("-", reference_sequence)[[1]]

        # find inserted sequences in aligned sequence
        inserted_sequences <- sapply(indexes, function(i){
            substr(aligned_sequence, i, i)
        })
        inserted_sequences <- paste(inserted_sequences, collapse="")

        # get upstream sequences
        ref_upstream_sequences <- substr(reference_sequence, 1, min(indexes)-1)
        ref_upstream_sequences <- paste(ref_upstream_sequences, collapse="")

        upstream_sequences <- substr(aligned_sequence, 1, min(indexes)-1)
        upstream_sequences <- paste(upstream_sequences, collapse="")

        # get downstream sequences
        ref_downstream_sequences <- substr(reference_sequence, max(indexes)+1, nchar(reference_sequence))
        ref_downstream_sequences <- paste(ref_downstream_sequences, collapse="")

        downstream_sequences <- substr(aligned_sequence, max(indexes)+1, nchar(aligned_sequence))
        downstream_sequences <- paste(downstream_sequences, collapse="")

        # compare upstream and downstream sequences
        if(ref_upstream_sequences==upstream_sequences & ref_downstream_sequences==downstream_sequences){
            return(paste(inserted_sequences, " Insertion", sep=""))
        }else{
            return("Other events")
        }

    }

    d_inserted_2bp$Event <- sapply(d_inserted_2bp$Read, function(x){
        Event <- assign_event(x, reference_sequence)
    })

    dir.create(file.path("myoutput",sample,reference_sequence), recursive=TRUE)

    d_inserted_2bp <- d_inserted_2bp %>% dplyr::mutate(Class='2bp Insertion', reference_sequence=reference_sequence)

    d_inserted_2bp %>% fwrite(file.path("myoutput",sample,reference_sequence,"2bp_insertions_AF_table.txt"), sep="\t")

    d_inserted_2bp

}

# get frequencies of 4bp insertions around cutting site
get_4bp_insertions_244bpAmp <- function(d, sample, reference_sequence){

    #reference_sequence <- "TTCTGTCACCAATC----CTGTCCCTAGTGGCCCCACTGTGGGGTGGAGG"
    d_inserted_4bp <- d %>% dplyr::filter(Reference==reference_sequence)
    reads_4bp_insertion <- sum(d_inserted_4bp$Number)
    print(reads_4bp_insertion)

    if(reads_4bp_insertion==0){
        return(NULL)
    }

    # add Event column
    assign_event <- function(aligned_sequence, reference_sequence){
        # find gap index in reference sequence
        indexes <- gregexpr("-", reference_sequence)[[1]]

        # find inserted sequences in aligned sequence
        inserted_sequences <- sapply(indexes, function(i){
            substr(aligned_sequence, i, i)
        })
        inserted_sequences <- paste(inserted_sequences, collapse="")

        # get upstream sequences
        ref_upstream_sequences <- substr(reference_sequence, 1, min(indexes)-1)
        ref_upstream_sequences <- paste(ref_upstream_sequences, collapse="")

        upstream_sequences <- substr(aligned_sequence, 1, min(indexes)-1)
        upstream_sequences <- paste(upstream_sequences, collapse="")

        # get downstream sequences
        ref_downstream_sequences <- substr(reference_sequence, max(indexes)+1, nchar(reference_sequence))
        ref_downstream_sequences <- paste(ref_downstream_sequences, collapse="")

        downstream_sequences <- substr(aligned_sequence, max(indexes)+1, nchar(aligned_sequence))
        downstream_sequences <- paste(downstream_sequences, collapse="")

        # compare upstream and downstream sequences
        if(ref_upstream_sequences==upstream_sequences & ref_downstream_sequences==downstream_sequences){
            return(paste(inserted_sequences, " Insertion", sep=""))
        }else{
            return("Other events")
        }

    }

    d_inserted_4bp$Event <- sapply(d_inserted_4bp$Read, function(x){
        Event <- assign_event(x, reference_sequence)
    })

    dir.create(file.path("myoutput",sample,reference_sequence), recursive=TRUE)

    d_inserted_4bp <- d_inserted_4bp %>% dplyr::mutate(Class='4bp Insertion', reference_sequence=reference_sequence)

    d_inserted_4bp %>% fwrite(file.path("myoutput",sample,reference_sequence,"4bp_insertions_AF_table.txt"), sep="\t")

    d_inserted_4bp

}

# get frequencies of 3bp, 2bp, and 4bp insertions around cutting site for 244bp amplicon samples
get_insertions_244bpAmp <- function(f){
    d <- fread(f)
    sample=d[1,]$Sample
    print(sample)
    d_inserted_3bp <- get_3bp_insertions_244bpAmp(d, sample, refseq_244bp_3bp); if(!is.null(d_inserted_3bp)){ d_inserted_3bp <- d_inserted_3bp %>% mutate(correction="at cutting site") }
    d_inserted_2bp <- get_2bp_insertions_244bpAmp(d, sample, refseq_244bp_2bp); if(!is.null(d_inserted_2bp)){ d_inserted_2bp <- d_inserted_2bp %>% mutate(correction="at cutting site") }
    d_inserted_4bp <- get_4bp_insertions_244bpAmp(d, sample, refseq_244bp_4bp); if(!is.null(d_inserted_4bp)){ d_inserted_4bp <- d_inserted_4bp  %>% mutate(correction="at cutting site") }

    # shift the right nucleotide to the left

    d_inserted_3bp_shiftCTG <- get_3bp_insertions_244bpAmp(d, sample, refseq_244bp_3bp_shiftCTG); if(!is.null(d_inserted_3bp_shiftCTG)){ d_inserted_3bp_shiftCTG <- d_inserted_3bp_shiftCTG %>% mutate(correction="shift CTG 3bp to the left") }

    d_inserted_3bp_shiftC <- get_3bp_insertions_244bpAmp(d, sample, refseq_244bp_3bp_shiftC); if(!is.null(d_inserted_3bp_shiftC)){ d_inserted_3bp_shiftC <- d_inserted_3bp_shiftC %>% mutate(correction="shift C 3bp to the left") }

    d_inserted_3bp_shiftCT <- get_3bp_insertions_244bpAmp(d, sample, refseq_244bp_3bp_shiftCT); if(!is.null(d_inserted_3bp_shiftCT)){ d_inserted_3bp_shiftCT <- d_inserted_3bp_shiftCT %>% mutate(correction="shift CT 3bp to the left") }

    d_inserted_2bp_shiftC <- get_2bp_insertions_244bpAmp(d, sample, refseq_244bp_2bp_shiftC); if(!is.null(d_inserted_2bp_shiftC)){ d_inserted_2bp_shiftC <- d_inserted_2bp_shiftC %>% mutate(correction="shift C 2bp to the left") }

    d_inserted_2bp_shiftCT <- get_2bp_insertions_244bpAmp(d, sample, refseq_244bp_2bp_shiftCT); if(!is.null(d_inserted_2bp_shiftCT)){ d_inserted_2bp_shiftCT <- d_inserted_2bp_shiftCT %>% mutate(correction="shift CT 2bp to the left") }

    d_inserted_4bp_shiftC <- get_4bp_insertions_244bpAmp(d, sample, refseq_244bp_4bp_shiftC); if(!is.null(d_inserted_4bp_shiftC)){ d_inserted_4bp_shiftC <- d_inserted_4bp_shiftC %>% mutate(correction="shift C 4bp to the left") }
    d_inserted_4bp_shiftCT <- get_4bp_insertions_244bpAmp(d, sample, refseq_244bp_4bp_shiftCT); if(!is.null(d_inserted_4bp_shiftCT)){ d_inserted_4bp_shiftCT <- d_inserted_4bp_shiftCT %>% mutate(correction="shift CT 4bp to the left") }
    d_inserted_4bp_shiftCTG <- get_4bp_insertions_244bpAmp(d, sample, refseq_244bp_4bp_shiftCTG); if(!is.null(d_inserted_4bp_shiftCTG)){ d_inserted_4bp_shiftCTG <- d_inserted_4bp_shiftCTG %>% mutate(correction="shift CTG 4bp to the left") }
    d_inserted_4bp_shiftCTGT <- get_4bp_insertions_244bpAmp(d, sample, refseq_244bp_4bp_shiftCTGT); if(!is.null(d_inserted_4bp_shiftCTGT)){ d_inserted_4bp_shiftCTGT <- d_inserted_4bp_shiftCTGT %>% mutate(correction="shift CTGT 4bp to the left") }

    # shift the left nucleotide to the right
    d_inserted_3bp_shiftTCright <- get_3bp_insertions_244bpAmp(d, sample, refseq_244bp_3bp_shiftTCright); if(!is.null(d_inserted_3bp_shiftTCright)){ d_inserted_3bp_shiftTCright <- d_inserted_3bp_shiftTCright %>% mutate(correction="shift TC 3bp to the right") }
    d_inserted_3bp_shiftCright <- get_3bp_insertions_244bpAmp(d, sample, refseq_244bp_3bp_shiftCright); if(!is.null(d_inserted_3bp_shiftCright)){ d_inserted_3bp_shiftCright <- d_inserted_3bp_shiftCright %>% mutate(correction="shift C 3bp to the right") }
    d_inserted_3bp_shiftATCright <- get_3bp_insertions_244bpAmp(d, sample, refseq_244bp_3bp_shiftATCright); if(!is.null(d_inserted_3bp_shiftATCright)){ d_inserted_3bp_shiftATCright <- d_inserted_3bp_shiftATCright %>% mutate(correction="shift ATC 3bp to the right") }

    d_inserted_4bp_shiftAATCright <- get_4bp_insertions_244bpAmp(d, sample, refseq_244bp_4bp_shiftAATCright); if(!is.null(d_inserted_4bp_shiftAATCright)){ d_inserted_4bp_shiftAATCright <- d_inserted_4bp_shiftAATCright %>% mutate(correction="shift AATC 4bp to the right") }
    d_inserted_4bp_shiftATCright <- get_4bp_insertions_244bpAmp(d, sample, refseq_244bp_4bp_shiftATCright); if(!is.null(d_inserted_4bp_shiftATCright)){ d_inserted_4bp_shiftATCright <- d_inserted_4bp_shiftATCright %>% mutate(correction="shift ATC 4bp to the right") }
    d_inserted_4bp_shiftTCright <- get_4bp_insertions_244bpAmp(d, sample, refseq_244bp_4bp_shiftTCright); if(!is.null(d_inserted_4bp_shiftTCright)){ d_inserted_4bp_shiftTCright <- d_inserted_4bp_shiftTCright %>% mutate(correction="shift TC 4bp to the right") }
    d_inserted_4bp_shiftCright <- get_4bp_insertions_244bpAmp(d, sample, refseq_244bp_4bp_shiftCright); if(!is.null(d_inserted_4bp_shiftCright)){ d_inserted_4bp_shiftCright <- d_inserted_4bp_shiftCright %>% mutate(correction="shift C 4bp to the right") }

    d_inserted_2bp_shiftTCright <- get_2bp_insertions_244bpAmp(d, sample, refseq_244bp_2bp_shiftTCright); if(!is.null(d_inserted_2bp_shiftTCright)){ d_inserted_2bp_shiftTCright <- d_inserted_2bp_shiftTCright %>% mutate(correction="shift TC 2bp to the right") }
    d_inserted_2bp_shiftCright <- get_2bp_insertions_244bpAmp(d, sample, refseq_244bp_2bp_shiftCright); if(!is.null(d_inserted_2bp_shiftCright)){ d_inserted_2bp_shiftCright <- d_inserted_2bp_shiftCright %>% mutate(correction="shift C 2bp to the right") }

    inserted <- rbind(d_inserted_3bp, d_inserted_2bp, d_inserted_4bp, 
        d_inserted_3bp_shiftCTG,d_inserted_3bp_shiftC,d_inserted_3bp_shiftCT,
        d_inserted_2bp_shiftC,d_inserted_2bp_shiftCT,
        d_inserted_4bp_shiftC,d_inserted_4bp_shiftCT,d_inserted_4bp_shiftCTG,d_inserted_4bp_shiftCTGT,
        d_inserted_3bp_shiftTCright,d_inserted_3bp_shiftCright,d_inserted_3bp_shiftATCright,
        d_inserted_4bp_shiftAATCright,d_inserted_4bp_shiftATCright,d_inserted_4bp_shiftTCright,d_inserted_4bp_shiftCright,
        d_inserted_2bp_shiftTCright,d_inserted_2bp_shiftCright)

}


# get frequencies of 2bp insertions around cutting site for 100bp amplicon samples
get_3bp_insertions_100bpAmp <- function(d, sample, reference_sequence){

    #reference_sequence <- "GTCACCAATC---CTGTCCCTAGTGGCCCCACTGTGGGGT"
    d_inserted_3bp <- d %>% dplyr::filter(Reference==reference_sequence)
    reads_3bp_insertion <- sum(d_inserted_3bp$Number)
    print(reads_3bp_insertion)

    if(reads_3bp_insertion==0){
        return(NULL)
    }

    # add Event column
    assign_event <- function(aligned_sequence, reference_sequence){
        # find gap index in reference sequence
        indexes <- gregexpr("-", reference_sequence)[[1]]

        # find inserted sequences in aligned sequence
        inserted_sequences <- sapply(indexes, function(i){
            substr(aligned_sequence, i, i)
        })
        inserted_sequences <- paste(inserted_sequences, collapse="")

        # get upstream sequences
        ref_upstream_sequences <- substr(reference_sequence, 1, min(indexes)-1)
        ref_upstream_sequences <- paste(ref_upstream_sequences, collapse="")

        upstream_sequences <- substr(aligned_sequence, 1, min(indexes)-1)
        upstream_sequences <- paste(upstream_sequences, collapse="")

        # get downstream sequences
        ref_downstream_sequences <- substr(reference_sequence, max(indexes)+1, nchar(reference_sequence))
        ref_downstream_sequences <- paste(ref_downstream_sequences, collapse="")

        downstream_sequences <- substr(aligned_sequence, max(indexes)+1, nchar(aligned_sequence))
        downstream_sequences <- paste(downstream_sequences, collapse="")

        # compare upstream and downstream sequences
        if(ref_upstream_sequences==upstream_sequences & ref_downstream_sequences==downstream_sequences){
            return(paste(inserted_sequences, " Insertion", sep=""))
        }else{
            return("Other events")
        }

    }

    d_inserted_3bp$Event <- sapply(d_inserted_3bp$Read, function(x){
        Event <- assign_event(x, reference_sequence)
    })

    dir.create(file.path("myoutput",sample,reference_sequence), recursive=TRUE)

    d_inserted_3bp <- d_inserted_3bp %>% dplyr::mutate(Class='3bp Insertion', reference_sequence=reference_sequence)

    d_inserted_3bp %>% fwrite(file.path("myoutput",sample,reference_sequence,"3bp_insertions_AF_table.txt"), sep="\t")

    d_inserted_3bp

}

# get frequencies of 2bp insertions around cutting site for 100bp amplicon samples
get_2bp_insertions_100bpAmp <- function(d, sample, reference_sequence){

    #reference_sequence <- "TGTCACCAATC--CTGTCCCTAGTGGCCCCACTGTGGGGT"
    d_inserted_2bp <- d %>% dplyr::filter(Reference==reference_sequence)
    reads_2bp_insertion <- sum(d_inserted_2bp$Number)
    print(reads_2bp_insertion)

    if(reads_2bp_insertion==0){
        return(NULL)
    }

    # add Event column
    assign_event <- function(aligned_sequence, reference_sequence){
        # find gap index in reference sequence
        indexes <- gregexpr("-", reference_sequence)[[1]]

        # find inserted sequences in aligned sequence
        inserted_sequences <- sapply(indexes, function(i){
            substr(aligned_sequence, i, i)
        })
        inserted_sequences <- paste(inserted_sequences, collapse="")

        # get upstream sequences
        ref_upstream_sequences <- substr(reference_sequence, 1, min(indexes)-1)
        ref_upstream_sequences <- paste(ref_upstream_sequences, collapse="")

        upstream_sequences <- substr(aligned_sequence, 1, min(indexes)-1)
        upstream_sequences <- paste(upstream_sequences, collapse="")

        # get downstream sequences
        ref_downstream_sequences <- substr(reference_sequence, max(indexes)+1, nchar(reference_sequence))
        ref_downstream_sequences <- paste(ref_downstream_sequences, collapse="")

        downstream_sequences <- substr(aligned_sequence, max(indexes)+1, nchar(aligned_sequence))
        downstream_sequences <- paste(downstream_sequences, collapse="")

        # compare upstream and downstream sequences
        if(ref_upstream_sequences==upstream_sequences & ref_downstream_sequences==downstream_sequences){
            return(paste(inserted_sequences, " Insertion", sep=""))
        }else{
            return("Other events")
        }

    }

    d_inserted_2bp$Event <- sapply(d_inserted_2bp$Read, function(x){
        Event <- assign_event(x, reference_sequence)
    })

    dir.create(file.path("myoutput",sample,reference_sequence), recursive=TRUE)

    d_inserted_2bp <- d_inserted_2bp %>% dplyr::mutate(Class='2bp Insertion', reference_sequence=reference_sequence)

    d_inserted_2bp %>% fwrite(file.path("myoutput",sample,reference_sequence,"2bp_insertions_AF_table.txt"), sep="\t")

    d_inserted_2bp

}

# get frequencies of 4bp insertions around cutting site for 100bp amplicon samples
get_4bp_insertions_100bpAmp <- function(d, sample, reference_sequence){

    #reference_sequence <- "TCACCAATC----CTGTCCCTAGTGGCCCCACTGTGGGGT"
    d_inserted_4bp <- d %>% dplyr::filter(Reference==reference_sequence)
    reads_4bp_insertion <- sum(d_inserted_4bp$Number)
    print(reads_4bp_insertion)

    if(reads_4bp_insertion==0){
        return(NULL)
    }

    # add Event column
    assign_event <- function(aligned_sequence, reference_sequence){
        # find gap index in reference sequence
        indexes <- gregexpr("-", reference_sequence)[[1]]

        # find inserted sequences in aligned sequence
        inserted_sequences <- sapply(indexes, function(i){
            substr(aligned_sequence, i, i)
        })
        inserted_sequences <- paste(inserted_sequences, collapse="")

        # get upstream sequences
        ref_upstream_sequences <- substr(reference_sequence, 1, min(indexes)-1)
        ref_upstream_sequences <- paste(ref_upstream_sequences, collapse="")

        upstream_sequences <- substr(aligned_sequence, 1, min(indexes)-1)
        upstream_sequences <- paste(upstream_sequences, collapse="")

        # get downstream sequences
        ref_downstream_sequences <- substr(reference_sequence, max(indexes)+1, nchar(reference_sequence))
        ref_downstream_sequences <- paste(ref_downstream_sequences, collapse="")

        downstream_sequences <- substr(aligned_sequence, max(indexes)+1, nchar(aligned_sequence))
        downstream_sequences <- paste(downstream_sequences, collapse="")

        # compare upstream and downstream sequences
        if(ref_upstream_sequences==upstream_sequences & ref_downstream_sequences==downstream_sequences){
            return(paste(inserted_sequences, " Insertion", sep=""))
        }else{
            return("Other events")
        }

    }

    d_inserted_4bp$Event <- sapply(d_inserted_4bp$Read, function(x){
        Event <- assign_event(x, reference_sequence)
    })

    dir.create(file.path("myoutput",sample,reference_sequence), recursive=TRUE)

    d_inserted_4bp <- d_inserted_4bp %>% dplyr::mutate(Class='4bp Insertion', reference_sequence=reference_sequence)

    d_inserted_4bp %>% fwrite(file.path("myoutput",sample,reference_sequence,"4bp_insertions_AF_table.txt"), sep="\t")

    d_inserted_4bp

}

# get frequencies of 3bp, 2bp, and 4bp insertions around cutting site for 100bp amplicon samples
get_insertions_100bpAmp <- function(f){
    d <- fread(f)
    sample=d[1,]$Sample
    print(sample)
    d_inserted_3bp <- get_3bp_insertions_100bpAmp(d, sample, refseq_100bp_3bp); if(!is.null(d_inserted_3bp)){ d_inserted_3bp <- d_inserted_3bp %>% mutate(correction="at cutting site") }
    d_inserted_2bp <- get_2bp_insertions_100bpAmp(d, sample, refseq_100bp_2bp); if(!is.null(d_inserted_2bp)){ d_inserted_2bp <- d_inserted_2bp %>% mutate(correction="at cutting site") }
    d_inserted_4bp <- get_4bp_insertions_100bpAmp(d, sample, refseq_100bp_4bp); if(!is.null(d_inserted_4bp)){ d_inserted_4bp <- d_inserted_4bp  %>% mutate(correction="at cutting site") }

    # shift the right nucleotide to the left

    d_inserted_3bp_shiftCTG <- get_3bp_insertions_100bpAmp(d, sample, refseq_100bp_3bp_shiftCTG); if(!is.null(d_inserted_3bp_shiftCTG)){ d_inserted_3bp_shiftCTG <- d_inserted_3bp_shiftCTG %>% mutate(correction="shift CTG 3bp to the left") }

    d_inserted_3bp_shiftC <- get_3bp_insertions_100bpAmp(d, sample, refseq_100bp_3bp_shiftC); if(!is.null(d_inserted_3bp_shiftC)){ d_inserted_3bp_shiftC <- d_inserted_3bp_shiftC %>% mutate(correction="shift C 3bp to the left") }

    d_inserted_3bp_shiftCT <- get_3bp_insertions_100bpAmp(d, sample, refseq_100bp_3bp_shiftCT); if(!is.null(d_inserted_3bp_shiftCT)){ d_inserted_3bp_shiftCT <- d_inserted_3bp_shiftCT %>% mutate(correction="shift CT 3bp to the left") }

    d_inserted_2bp_shiftC <- get_2bp_insertions_100bpAmp(d, sample, refseq_100bp_2bp_shiftC); if(!is.null(d_inserted_2bp_shiftC)){ d_inserted_2bp_shiftC <- d_inserted_2bp_shiftC %>% mutate(correction="shift C 2bp to the left") }

    d_inserted_2bp_shiftCT <- get_2bp_insertions_100bpAmp(d, sample, refseq_100bp_2bp_shiftCT); if(!is.null(d_inserted_2bp_shiftCT)){ d_inserted_2bp_shiftCT <- d_inserted_2bp_shiftCT %>% mutate(correction="shift CT 2bp to the left") }

    d_inserted_4bp_shiftC <- get_4bp_insertions_100bpAmp(d, sample, refseq_100bp_4bp_shiftC); if(!is.null(d_inserted_4bp_shiftC)){ d_inserted_4bp_shiftC <- d_inserted_4bp_shiftC %>% mutate(correction="shift C 4bp to the left") }
    d_inserted_4bp_shiftCT <- get_4bp_insertions_100bpAmp(d, sample, refseq_100bp_4bp_shiftCT); if(!is.null(d_inserted_4bp_shiftCT)){ d_inserted_4bp_shiftCT <- d_inserted_4bp_shiftCT %>% mutate(correction="shift CT 4bp to the left") }
    d_inserted_4bp_shiftCTG <- get_4bp_insertions_100bpAmp(d, sample, refseq_100bp_4bp_shiftCTG); if(!is.null(d_inserted_4bp_shiftCTG)){ d_inserted_4bp_shiftCTG <- d_inserted_4bp_shiftCTG %>% mutate(correction="shift CTG 4bp to the left") }
    d_inserted_4bp_shiftCTGT <- get_4bp_insertions_100bpAmp(d, sample, refseq_100bp_4bp_shiftCTGT); if(!is.null(d_inserted_4bp_shiftCTGT)){ d_inserted_4bp_shiftCTGT <- d_inserted_4bp_shiftCTGT %>% mutate(correction="shift CTGT 4bp to the left") }

    # shift the left nucleotide to the right
    d_inserted_3bp_shiftTCright <- get_3bp_insertions_100bpAmp(d, sample, refseq_100bp_3bp_shiftTCright); if(!is.null(d_inserted_3bp_shiftTCright)){ d_inserted_3bp_shiftTCright <- d_inserted_3bp_shiftTCright %>% mutate(correction="shift TC 3bp to the right") }
    d_inserted_3bp_shiftCright <- get_3bp_insertions_100bpAmp(d, sample, refseq_100bp_3bp_shiftCright); if(!is.null(d_inserted_3bp_shiftCright)){ d_inserted_3bp_shiftCright <- d_inserted_3bp_shiftCright %>% mutate(correction="shift C 3bp to the right") }
    d_inserted_3bp_shiftATCright <- get_3bp_insertions_100bpAmp(d, sample, refseq_100bp_3bp_shiftATCright); if(!is.null(d_inserted_3bp_shiftATCright)){ d_inserted_3bp_shiftATCright <- d_inserted_3bp_shiftATCright %>% mutate(correction="shift ATC 3bp to the right") }

    d_inserted_4bp_shiftAATCright <- get_4bp_insertions_100bpAmp(d, sample, refseq_100bp_4bp_shiftAATCright); if(!is.null(d_inserted_4bp_shiftAATCright)){ d_inserted_4bp_shiftAATCright <- d_inserted_4bp_shiftAATCright %>% mutate(correction="shift AATC 4bp to the right") }
    d_inserted_4bp_shiftATCright <- get_4bp_insertions_100bpAmp(d, sample, refseq_100bp_4bp_shiftATCright); if(!is.null(d_inserted_4bp_shiftATCright)){ d_inserted_4bp_shiftATCright <- d_inserted_4bp_shiftATCright %>% mutate(correction="shift ATC 4bp to the right") }
    d_inserted_4bp_shiftTCright <- get_4bp_insertions_100bpAmp(d, sample, refseq_100bp_4bp_shiftTCright); if(!is.null(d_inserted_4bp_shiftTCright)){ d_inserted_4bp_shiftTCright <- d_inserted_4bp_shiftTCright %>% mutate(correction="shift TC 4bp to the right") }
    d_inserted_4bp_shiftCright <- get_4bp_insertions_100bpAmp(d, sample, refseq_100bp_4bp_shiftCright); if(!is.null(d_inserted_4bp_shiftCright)){ d_inserted_4bp_shiftCright <- d_inserted_4bp_shiftCright %>% mutate(correction="shift C 4bp to the right") }

    d_inserted_2bp_shiftTCright <- get_2bp_insertions_100bpAmp(d, sample, refseq_100bp_2bp_shiftTCright); if(!is.null(d_inserted_2bp_shiftTCright)){ d_inserted_2bp_shiftTCright <- d_inserted_2bp_shiftTCright %>% mutate(correction="shift TC 2bp to the right") }
    d_inserted_2bp_shiftCright <- get_2bp_insertions_100bpAmp(d, sample, refseq_100bp_2bp_shiftCright); if(!is.null(d_inserted_2bp_shiftCright)){ d_inserted_2bp_shiftCright <- d_inserted_2bp_shiftCright %>% mutate(correction="shift C 2bp to the right") }


    inserted <- rbind(d_inserted_3bp, d_inserted_2bp, d_inserted_4bp, 
        d_inserted_3bp_shiftCTG, 
        d_inserted_3bp_shiftC,
        d_inserted_3bp_shiftCT, d_inserted_2bp_shiftC, d_inserted_2bp_shiftCT,
        d_inserted_4bp_shiftC, d_inserted_4bp_shiftCT, d_inserted_4bp_shiftCTG, d_inserted_4bp_shiftCTGT,
        d_inserted_3bp_shiftTCright, d_inserted_3bp_shiftCright, d_inserted_3bp_shiftATCright,
        d_inserted_4bp_shiftAATCright, d_inserted_4bp_shiftATCright, d_inserted_4bp_shiftTCright, d_inserted_4bp_shiftCright,
        d_inserted_2bp_shiftTCright, d_inserted_2bp_shiftCright)
}


main <- function(){

    setwd("/home/zhuy1/my_projects_nrlab/Manisha_RNA_templated_repair/analysis/crispr_3bp_2bp_4bp/")

    # an example
    #f <- "output_juber_fastq_cleaned_244bp/C279a_Hek293T_siNT_C-AF-table.txt"
    #inserted <- get_insertions_244bpAmp(f)

    # there is an issue on Event variable, some are list, for example, this one
    # Found the issue, some are empty, so Event variable goes to a list
    #f="/lila/data/riazlab/projects/zhuy1/my_projects/Manisha_RNA_templated_repair/analysis/crispr_3bp_2bp_4bp/output_juber_fastq_cleaned_244bp/C285a_Hek293T_siNT_pMJ119-AF-table.txt"
    #inserted <- get_insertions_244bpAmp(f)

    # get insertions for all samples with 244bp amplicon
    # used 3 references: refseq_244bp, refseq_244bp_shiftCTG, refseq_244bp_shiftC
    filelist <- fread("AF_table_244bpAmp_path.txt", header=FALSE)
    inserted_1 <- do.call(rbind,lapply(filelist$V1, get_insertions_244bpAmp))
 
    # get insertions for all samples with 100bp amplicon
    # used 3 references: refseq_100bp, refseq_100bp_shiftCTG, refseq_100bp_shiftC
    filelist <- fread("AF_table_100bpAmp_path.txt", header=FALSE)
    inserted_2 <- do.call(rbind,lapply(filelist$V1, get_insertions_100bpAmp))

    results <- rbind(inserted_1, inserted_2) %>% arrange(reference_sequence, Sample, Class, Event) %>%
        dplyr::mutate(design=ifelse(grepl('_C$',Sample) | grepl('-C$',Sample),'C',
            ifelse(grepl('siNT',Sample) & grepl('pMJ119',Sample),'siNT pMJ119',
            ifelse(grepl('siRev3',Sample),'siRev3 pMJ119','other')))) %>%
        dplyr::mutate(Event.uncorrected=Event)

    # since the inserted sequence is shifted in alignment, so we need to correct the inserted sequence
    correct_insertion_events <- function(results){

        # correct 3bp insertions shifted to the left
        results <- results %>% dplyr::mutate(Event=ifelse((reference_sequence==refseq_244bp_3bp_shiftC | reference_sequence==refseq_100bp_3bp_shiftC) & substr(Event,3,3)=='C',
            paste0('C',substr(Event,1,2)," Insertion"),
            ifelse((reference_sequence==refseq_244bp_3bp_shiftC | reference_sequence==refseq_100bp_3bp_shiftC) & substr(Event,3,3)!='C',"Other events",Event)))

        results <- results %>% dplyr::mutate(Event=ifelse((reference_sequence==refseq_244bp_3bp_shiftCT | reference_sequence==refseq_100bp_3bp_shiftCT) & substr(Event,2,3)=='CT',
            paste0('CT',substr(Event,1,1)," Insertion"),
            ifelse((reference_sequence==refseq_244bp_3bp_shiftCT | reference_sequence==refseq_100bp_3bp_shiftCT) & substr(Event,2,3)!='CT',"Other events",Event)))

        results <- results %>% dplyr::mutate(Event=ifelse((reference_sequence==refseq_244bp_3bp_shiftCTG | reference_sequence==refseq_100bp_3bp_shiftCTG) & Event!="CTG Insertion",'Other events',Event))  
        
        # correct 2bp insertions shifted to the left
        results <- results %>% dplyr::mutate(Event=ifelse((reference_sequence==refseq_244bp_2bp_shiftC | reference_sequence==refseq_100bp_2bp_shiftC) & substr(Event,2,2)=='C',
            paste0('C',substr(Event,1,1)," Insertion"),
            ifelse((reference_sequence==refseq_244bp_2bp_shiftC | reference_sequence==refseq_100bp_2bp_shiftC) & substr(Event,2,2)!='C',"Other events",Event)))

        results <- results %>% dplyr::mutate(Event=ifelse((reference_sequence==refseq_244bp_2bp_shiftCT | reference_sequence==refseq_100bp_2bp_shiftCT) & Event!="CT Insertion",'Other events',Event))  
        
        # correct 4bp insertions shifted to the left
        results <- results %>% dplyr::mutate(Event=ifelse((reference_sequence==refseq_244bp_4bp_shiftC | reference_sequence==refseq_100bp_4bp_shiftC) & substr(Event,4,4)=='C',
            paste0('C',substr(Event,1,3)," Insertion"),
            ifelse((reference_sequence==refseq_244bp_4bp_shiftC | reference_sequence==refseq_100bp_4bp_shiftC) & substr(Event,4,4)!='C',"Other events",Event))) 

        results <- results %>% dplyr::mutate(Event=ifelse((reference_sequence==refseq_244bp_4bp_shiftCT | reference_sequence==refseq_100bp_4bp_shiftCT) & substr(Event,3,4)=='CT',
            paste0('CT',substr(Event,1,2)," Insertion"),
            ifelse((reference_sequence==refseq_244bp_4bp_shiftCT | reference_sequence==refseq_100bp_4bp_shiftCT) & substr(Event,3,4)!='CT',"Other events",Event))) 

        results <- results %>% dplyr::mutate(Event=ifelse((reference_sequence==refseq_244bp_4bp_shiftCTG | reference_sequence==refseq_100bp_4bp_shiftCTG) & substr(Event,2,4)=='CTG',
            paste0('CTG',substr(Event,1,1)," Insertion"),
            ifelse((reference_sequence==refseq_244bp_4bp_shiftCTG | reference_sequence==refseq_100bp_4bp_shiftCTG) & substr(Event,2,4)!='CTG',"Other events",Event))) 

        results <- results %>% dplyr::mutate(Event=ifelse((reference_sequence==refseq_244bp_4bp_shiftCTGT | reference_sequence==refseq_100bp_4bp_shiftCTGT) & Event!="CTGT Insertion",'Other events',Event))  


        # correct 3bp insertions shifted to the right
        results <- results %>% dplyr::mutate(Event=ifelse((reference_sequence==refseq_244bp_3bp_shiftTCright | reference_sequence==refseq_100bp_3bp_shiftTCright) & substr(Event,1,2)=='TC',
            paste0("TC",substr(Event,3,3)," Insertion"),
            ifelse((reference_sequence==refseq_244bp_3bp_shiftTCright | reference_sequence==refseq_100bp_3bp_shiftTCright) & substr(Event,1,2)!='TC',"Other events",Event)))

        results <- results %>% dplyr::mutate(Event=ifelse((reference_sequence==refseq_244bp_3bp_shiftCright | reference_sequence==refseq_100bp_3bp_shiftCright) & substr(Event,1,1)=='C',
            paste0("C",substr(Event,2,3)," Insertion"),
            ifelse((reference_sequence==refseq_244bp_3bp_shiftCright | reference_sequence==refseq_100bp_3bp_shiftCright) & substr(Event,1,1)!='C',"Other events",Event)))

        results <- results %>% dplyr::mutate(Event=ifelse((reference_sequence==refseq_244bp_3bp_shiftATCright | reference_sequence==refseq_100bp_3bp_shiftATCright) & Event!='ATC',"Other events",Event))

        # correct 4bp insertions shifted to the right
        results <- results %>% dplyr::mutate(Event=ifelse((reference_sequence==refseq_244bp_4bp_shiftAATCright | reference_sequence==refseq_100bp_4bp_shiftAATCright) & Event!='AATC',"Other events",Event))

        results <- results %>% dplyr::mutate(Event=ifelse((reference_sequence==refseq_244bp_4bp_shiftATCright | reference_sequence==refseq_100bp_4bp_shiftATCright) & substr(Event,1,3)=='ATC',
            paste0("ATC",substr(Event,4,4)," Insertion"),
            ifelse((reference_sequence==refseq_244bp_4bp_shiftATCright | reference_sequence==refseq_100bp_4bp_shiftATCright) & substr(Event,1,3)!='ATC',"Other events",Event)))

        results <- results %>% dplyr::mutate(Event=ifelse((reference_sequence==refseq_244bp_4bp_shiftTCright | reference_sequence==refseq_100bp_4bp_shiftTCright) & substr(Event,1,2)=='TC',
            paste0("TC",substr(Event,3,4)," Insertion"),
            ifelse((reference_sequence==refseq_244bp_4bp_shiftTCright | reference_sequence==refseq_100bp_4bp_shiftTCright) & substr(Event,1,2)!='TC',"Other events",Event)))

        results <- results %>% dplyr::mutate(Event=ifelse((reference_sequence==refseq_244bp_4bp_shiftCright | reference_sequence==refseq_100bp_4bp_shiftCright) & substr(Event,1,1)=='C',
            paste0("C",substr(Event,2,4)," Insertion"),
            ifelse((reference_sequence==refseq_244bp_4bp_shiftCright | reference_sequence==refseq_100bp_4bp_shiftCright) & substr(Event,1,1)!='C',"Other events",Event)))

        # correct 2bp insertions shifted to the right
        results <- results %>% dplyr::mutate(Event=ifelse((reference_sequence==refseq_244bp_2bp_shiftTCright | reference_sequence==refseq_100bp_2bp_shiftTCright) & Event!='TC',"Other events",Event))

        results <- results %>% dplyr::mutate(Event=ifelse((reference_sequence==refseq_244bp_2bp_shiftCright | reference_sequence==refseq_100bp_2bp_shiftCright) & substr(Event,1,1)=='C',
            paste0("C",substr(Event,2,2)," Insertion"),
            ifelse((reference_sequence==refseq_244bp_2bp_shiftCright | reference_sequence==refseq_100bp_2bp_shiftCright) & substr(Event,1,1)!='C',"Other events",Event)))

    }

    results <- correct_insertion_events(results)

    # check which events are corrected
    results %>% dplyr::filter(Event!=Event.uncorrected) %>% dplyr::select(reference_sequence, Class, Event.uncorrected, Event, correction) %>% distinct %>% arrange(Class,Event) %>% fwrite("myoutput/insertions_corrected.txt", sep="\t")

    # save
    results %>% fwrite("myoutput/insertions_AF_table.txt", sep="\t")
    results %>% dplyr::filter(Event!="Other events",design != 'siRev3 pMJ119') %>% fwrite("myoutput/insertions_AF_table_filtered.txt", sep="\t")

    # generate Insertions x Samples matrix
    results_3bp <- results %>% dplyr::filter(Event!="Other events",design != 'siRev3 pMJ119',Class=="3bp Insertion") %>% dplyr::select(design, Sample, Event, `Fraction of Indel Reads`) %>% arrange(design,Sample,Event)
    results_3bp %>% fwrite("myoutput/insertions_AF_table_filtered_3bp.txt", sep="\t")
    mat_3bp <- results_3bp %>% reshape2::acast(Event ~ Sample, value.var="Fraction of Indel Reads") %>% as.data.frame %>% rownames_to_column(var="Event") %>% dplyr::mutate(Event=gsub(" Insertion","",Event))


    # count how many Events for each reference sequence
    results %>% dplyr::group_by(reference_sequence, design, Sample, Class, Event) %>% dplyr::summarize(reads=sum(Number)) %>% fwrite("myoutput/insertions_count.txt", sep="\t")


    ## -------------------------------------------------------------------------------------------------------------------
    # make plots for each type of insertions

    # for insertions at the cutting site
    make_plots_at_cutting_site <- function(res, outDir){

        dir.create(outDir, recursive=TRUE)

        # filter out siRev3 pMJ119
        res <- res %>% dplyr::filter(design != 'siRev3 pMJ119') %>% arrange(Class,design,Sample) # let's focus on siNT pMJ119 and C

        # filter out "Other events" and "N" in Event
        results_3bp <- res %>% dplyr::filter(Class=="3bp Insertion", Event!="Other events", !grepl('N',Event))
        results_2bp <- res %>% dplyr::filter(Class=="2bp Insertion", Event!="Other events", !grepl('N',Event))
        results_4bp <- res %>% dplyr::filter(Class=="4bp Insertion", Event!="Other events", !grepl('N',Event))

        # convert Event, Sample, and "Fraction of Indel Reads" columns to a matrix
        results_3bp_matrix <- results_3bp %>% reshape2::acast(Event ~ Sample, value.var="Fraction of Indel Reads")
        results_2bp_matrix <- results_2bp %>% reshape2::acast(Event ~ Sample, value.var="Fraction of Indel Reads")
        results_4bp_matrix <- results_4bp %>% reshape2::acast(Event ~ Sample, value.var="Fraction of Indel Reads")

        # filing NA with 0
        results_3bp_matrix[is.na(results_3bp_matrix)] <- 0
        results_2bp_matrix[is.na(results_2bp_matrix)] <- 0
        results_4bp_matrix[is.na(results_4bp_matrix)] <- 0

        # order samples by design
        samples_3bp_C    <- results_3bp %>% dplyr::filter(design=="C") %>% dplyr::select(Sample) %>% dplyr::pull() %>% unique()
        samples_3bp_siNT <- results_3bp %>% dplyr::filter(design=="siNT pMJ119") %>% dplyr::select(Sample) %>% dplyr::pull() %>% unique()
        samples_2bp_C    <- results_2bp %>% dplyr::filter(design=="C") %>% dplyr::select(Sample) %>% dplyr::pull() %>% unique()
        samples_2bp_siNT <- results_2bp %>% dplyr::filter(design=="siNT pMJ119") %>% dplyr::select(Sample) %>% dplyr::pull() %>% unique()
        samples_4bp_C    <- results_4bp %>% dplyr::filter(design=="C") %>% dplyr::select(Sample) %>% dplyr::pull() %>% unique()
        samples_4bp_siNT <- results_4bp %>% dplyr::filter(design=="siNT pMJ119") %>% dplyr::select(Sample) %>% dplyr::pull() %>% unique()

        results_3bp_matrix <- cbind(results_3bp_matrix[,samples_3bp_C], results_3bp_matrix[,samples_3bp_siNT])
        results_2bp_matrix <- cbind(results_2bp_matrix[,samples_2bp_C], results_2bp_matrix[,samples_2bp_siNT])
        results_4bp_matrix <- cbind(results_4bp_matrix[,samples_4bp_C], results_4bp_matrix[,samples_4bp_siNT])

        # make heatmap for 3bp insertions
        ht <- Heatmap(results_3bp_matrix, name="Fraction of Indel Reads", col = colorRamp2(c(0,0.001, 0.01, 0.05), c("white", "orange","red","blue")), 
            cluster_rows = FALSE, cluster_columns = FALSE, border=TRUE)

        pdf(paste0(outDir,"/3bp_insertions_heatmap.pdf"), width=10, height=10)
        draw(ht)
        dev.off()

        # make heatmap for 2bp insertions
        ht <- Heatmap(results_2bp_matrix, name="Fraction of Indel Reads", col = colorRamp2(c(0,0.001, 0.01, 0.01), c("white", "orange","red","blue")), 
            cluster_rows = FALSE, cluster_columns = FALSE, border=TRUE)

        pdf(paste0(outDir,"/2bp_insertions_heatmap.pdf"), width=10, height=10)
        draw(ht)
        dev.off()

        # make heatmap for 4bp insertions
        ht <- Heatmap(results_4bp_matrix, name="Fraction of Indel Reads", col = colorRamp2(c(0,0.0001, 0.001), c("white", "orange","red")), 
            cluster_rows = FALSE, cluster_columns = FALSE, border=TRUE)

        pdf(paste0(outDir,"/4bp_insertions_heatmap.pdf"), width=10, height=10)
        draw(ht)
        dev.off()


        ## -------------------------------------------------------------------------------------------------------------------
        # for each type of insertions, make barplot for mean of "Fraction of Indel Reads" 
        results_3bp_mean <- results_3bp %>% dplyr::group_by(design,Event) %>% dplyr::summarize(mean=mean(`Fraction of Indel Reads`))
        results_2bp_mean <- results_2bp %>% dplyr::group_by(design,Event) %>% dplyr::summarize(mean=mean(`Fraction of Indel Reads`))
        results_4bp_mean <- results_4bp %>% dplyr::group_by(design,Event) %>% dplyr::summarize(mean=mean(`Fraction of Indel Reads`))

        p <- ggplot(results_3bp_mean, aes(x=Event, y=mean)) + geom_bar(stat="identity") + theme_bw() + ylab("Fraction of Indel Reads") + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + facet_wrap(~design, ncol=1)
        ggsave(paste0(outDir,"/3bp_insertions_barplot.mean.pdf"),width=7,height=8)

        p <- ggplot(results_2bp_mean, aes(x=Event, y=mean)) + geom_bar(stat="identity") + theme_bw() + ylab("Fraction of Indel Reads") + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + facet_wrap(~design, ncol=1)
        ggsave(paste0(outDir,"/2bp_insertions_barplot.mean.pdf"),width=7,height=8)

        p <- ggplot(results_4bp_mean, aes(x=Event, y=mean)) + geom_bar(stat="identity") + theme_bw() + ylab("Fraction of Indel Reads") + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + facet_wrap(~design,ncol=1)
        ggsave(paste0(outDir,"/4bp_insertions_barplot.mean.pdf"),width=7,height=8)

        ## -------------------------------------------------------------------------------------------------------------------
        # for each type of insertions, make barplot for sum of "Fraction of Indel Reads" 
        results_3bp_sum <- results_3bp %>% dplyr::group_by(design,Event) %>% dplyr::summarize(sum=sum(`Fraction of Indel Reads`))
        results_2bp_sum <- results_2bp %>% dplyr::group_by(design,Event) %>% dplyr::summarize(sum=sum(`Fraction of Indel Reads`))
        results_4bp_sum <- results_4bp %>% dplyr::group_by(design,Event) %>% dplyr::summarize(sum=sum(`Fraction of Indel Reads`))

        p <- ggplot(results_3bp_sum, aes(x=Event, y=sum)) + geom_bar(stat="identity") + theme_bw() + ylab("Fraction of Indel Reads") + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + facet_wrap(~design, ncol=1)
        ggsave(paste0(outDir,"/3bp_insertions_barplot.sum.pdf"),width=7,height=8)

        p <- ggplot(results_2bp_sum, aes(x=Event, y=sum)) + geom_bar(stat="identity") + theme_bw() + ylab("Fraction of Indel Reads") + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + facet_wrap(~design, ncol=1)
        ggsave(paste0(outDir,"/2bp_insertions_barplot.sum.pdf"),width=7,height=8)

        p <- ggplot(results_4bp_sum, aes(x=Event, y=sum)) + geom_bar(stat="identity") + theme_bw() + ylab("Fraction of Indel Reads") + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + facet_wrap(~design,ncol=1)
        ggsave(paste0(outDir,"/4bp_insertions_barplot.sum.pdf"),width=7,height=8)

    }

    res <- results %>% dplyr::filter(reference_sequence %in% c(refseq_244bp_3bp, refseq_244bp_2bp, refseq_244bp_4bp, refseq_100bp_3bp, refseq_100bp_2bp, refseq_100bp_4bp))
    #make_plots_at_cutting_site(res, "myoutput/insertions_at_cutting_site")

    # for 3bp insertions shifted CTG 3bp to the left
    make_plots_shifted_CTG_3bp_to_left <- function(res, outDir){

        dir.create(outDir, recursive=TRUE)

        # filter out siRev3 pMJ119
        res <- res %>% dplyr::filter(design != 'siRev3 pMJ119') %>% arrange(Class,design,Sample) # let's focus on siNT pMJ119 and C

        # filter out "Other events" and "N" in Event
        results_3bp <- res %>% dplyr::filter(Class=="3bp Insertion", Event!="Other events", !grepl('N',Event))

        # convert Event, Sample, and "Fraction of Indel Reads" columns to a matrix
        results_3bp_matrix <- results_3bp %>% reshape2::acast(Event ~ Sample, value.var="Fraction of Indel Reads")

        # filing NA with 0
        results_3bp_matrix[is.na(results_3bp_matrix)] <- 0

        # order samples by design
        samples_3bp_C    <- results_3bp %>% dplyr::filter(design=="C") %>% dplyr::select(Sample) %>% dplyr::pull() %>% unique()
        samples_3bp_siNT <- results_3bp %>% dplyr::filter(design=="siNT pMJ119") %>% dplyr::select(Sample) %>% dplyr::pull() %>% unique()

        results_3bp_matrix <- cbind(results_3bp_matrix[,samples_3bp_C], results_3bp_matrix[,samples_3bp_siNT])

        # make heatmap for 3bp insertions
        ht <- Heatmap(results_3bp_matrix, name="Fraction of Indel Reads", col = colorRamp2(c(0,0.001, 0.01,0.05), c("white", "orange","red","blue")), 
            cluster_rows = FALSE, cluster_columns = FALSE, border=TRUE)

        pdf(paste0(outDir,"/3bp_insertions_heatmap.pdf"), width=10, height=10)
        draw(ht)
        dev.off()


        ## -------------------------------------------------------------------------------------------------------------------
        # for each type of insertions, make barplot for mean of "Fraction of Indel Reads" 
        results_3bp_mean <- results_3bp %>% dplyr::group_by(design,Event) %>% dplyr::summarize(mean=mean(`Fraction of Indel Reads`))

        p <- ggplot(results_3bp_mean, aes(x=Event, y=mean)) + geom_bar(stat="identity") + theme_bw() + ylab("Fraction of Indel Reads") + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + facet_wrap(~design, ncol=1)
        ggsave(paste0(outDir,"/3bp_insertions_barplot.mean.pdf"),width=7,height=8)

        ## -------------------------------------------------------------------------------------------------------------------
        # for each type of insertions, make barplot for sum of "Fraction of Indel Reads" 
        results_3bp_sum <- results_3bp %>% dplyr::group_by(design,Event) %>% dplyr::summarize(sum=sum(`Fraction of Indel Reads`))

        p <- ggplot(results_3bp_sum, aes(x=Event, y=sum)) + geom_bar(stat="identity") + theme_bw() + ylab("Fraction of Indel Reads") + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + facet_wrap(~design, ncol=1)
        ggsave(paste0(outDir,"/3bp_insertions_barplot.sum.pdf"),width=7,height=8)

    }
    res <- results %>% dplyr::filter(reference_sequence %in% c(refseq_244bp_3bp_shiftCTG, refseq_100bp_3bp_shiftCTG))
    #make_plots_shifted_CTG_3bp_to_left(res, "myoutput/insertions_shifted_CTG_3bp_to_left")

    # for 3bp insertions shifted C 3bp to the left
    make_plots_shifted_C_3bp_to_left <- function(res, outDir){

        dir.create(outDir, recursive=TRUE)

        # filter out siRev3 pMJ119
        res <- res %>% dplyr::filter(design != 'siRev3 pMJ119') %>% arrange(Class,design,Sample) # let's focus on siNT pMJ119 and C

        # filter out "Other events" and "N" in Event
        results_3bp <- res %>% dplyr::filter(Class=="3bp Insertion", Event!="Other events", !grepl('N',Event))

        # convert Event, Sample, and "Fraction of Indel Reads" columns to a matrix
        results_3bp_matrix <- results_3bp %>% reshape2::acast(Event ~ Sample, value.var="Fraction of Indel Reads")

        # filing NA with 0
        results_3bp_matrix[is.na(results_3bp_matrix)] <- 0

        # order samples by design
        samples_3bp_C    <- results_3bp %>% dplyr::filter(design=="C") %>% dplyr::select(Sample) %>% dplyr::pull() %>% unique()
        samples_3bp_siNT <- results_3bp %>% dplyr::filter(design=="siNT pMJ119") %>% dplyr::select(Sample) %>% dplyr::pull() %>% unique()

        results_3bp_matrix <- cbind(results_3bp_matrix[,samples_3bp_C], results_3bp_matrix[,samples_3bp_siNT])

        # make heatmap for 3bp insertions
        ht <- Heatmap(results_3bp_matrix, name="Fraction of Indel Reads", col = colorRamp2(c(0,0.001, 0.01), c("white", "orange","red")), 
            cluster_rows = FALSE, cluster_columns = FALSE, border=TRUE)

        pdf(paste0(outDir,"/3bp_insertions_heatmap.pdf"), width=10, height=10)
        draw(ht)
        dev.off()


        ## -------------------------------------------------------------------------------------------------------------------
        # for each type of insertions, make barplot for mean of "Fraction of Indel Reads" 
        results_3bp_mean <- results_3bp %>% dplyr::group_by(design,Event) %>% dplyr::summarize(mean=mean(`Fraction of Indel Reads`))

        p <- ggplot(results_3bp_mean, aes(x=Event, y=mean)) + geom_bar(stat="identity") + theme_bw() + ylab("Fraction of Indel Reads") + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + facet_wrap(~design, ncol=1)
        ggsave(paste0(outDir,"/3bp_insertions_barplot.mean.pdf"),width=7,height=8)

        ## -------------------------------------------------------------------------------------------------------------------
        # for each type of insertions, make barplot for sum of "Fraction of Indel Reads" 
        results_3bp_sum <- results_3bp %>% dplyr::group_by(design,Event) %>% dplyr::summarize(sum=sum(`Fraction of Indel Reads`))

        p <- ggplot(results_3bp_sum, aes(x=Event, y=sum)) + geom_bar(stat="identity") + theme_bw() + ylab("Fraction of Indel Reads") + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + facet_wrap(~design, ncol=1)
        ggsave(paste0(outDir,"/3bp_insertions_barplot.sum.pdf"),width=7,height=8)

    }
    res <- results %>% dplyr::filter(reference_sequence %in% c(refseq_244bp_3bp_shiftC, refseq_100bp_3bp_shiftC))
    #make_plots_shifted_C_3bp_to_left(res, "myoutput/insertions_shifted_C_3bp_to_left")

    # for 3bp insertions shifted CT 3bp to the left
    make_plots_shifted_CT_3bp_to_left <- function(res, outDir){

        dir.create(outDir, recursive=TRUE)

        # filter out siRev3 pMJ119
        res <- res %>% dplyr::filter(design != 'siRev3 pMJ119') %>% arrange(Class,design,Sample) # let's focus on siNT pMJ119 and C

        # filter out "Other events" and "N" in Event
        results_3bp <- res %>% dplyr::filter(Class=="3bp Insertion", Event!="Other events", !grepl('N',Event))

        # convert Event, Sample, and "Fraction of Indel Reads" columns to a matrix
        results_3bp_matrix <- results_3bp %>% reshape2::acast(Event ~ Sample, value.var="Fraction of Indel Reads")

        # filing NA with 0
        results_3bp_matrix[is.na(results_3bp_matrix)] <- 0

        # order samples by design
        samples_3bp_C    <- results_3bp %>% dplyr::filter(design=="C") %>% dplyr::select(Sample) %>% dplyr::pull() %>% unique()
        samples_3bp_siNT <- results_3bp %>% dplyr::filter(design=="siNT pMJ119") %>% dplyr::select(Sample) %>% dplyr::pull() %>% unique()

        results_3bp_matrix <- cbind(results_3bp_matrix[,samples_3bp_C], results_3bp_matrix[,samples_3bp_siNT])

        # make heatmap for 3bp insertions
        ht <- Heatmap(results_3bp_matrix, name="Fraction of Indel Reads", col = colorRamp2(c(0,0.001, 0.01), c("white", "orange","red")), 
            cluster_rows = FALSE, cluster_columns = FALSE, border=TRUE)

        pdf(paste0(outDir,"/3bp_insertions_heatmap.pdf"), width=10, height=10)
        draw(ht)
        dev.off()


        ## -------------------------------------------------------------------------------------------------------------------
        # for each type of insertions, make barplot for mean of "Fraction of Indel Reads" 
        results_3bp_mean <- results_3bp %>% dplyr::group_by(design,Event) %>% dplyr::summarize(mean=mean(`Fraction of Indel Reads`))

        p <- ggplot(results_3bp_mean, aes(x=Event, y=mean)) + geom_bar(stat="identity") + theme_bw() + ylab("Fraction of Indel Reads") + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + facet_wrap(~design, ncol=1)
        ggsave(paste0(outDir,"/3bp_insertions_barplot.mean.pdf"),width=7,height=8)

        ## -------------------------------------------------------------------------------------------------------------------
        # for each type of insertions, make barplot for sum of "Fraction of Indel Reads" 
        results_3bp_sum <- results_3bp %>% dplyr::group_by(design,Event) %>% dplyr::summarize(sum=sum(`Fraction of Indel Reads`))

        p <- ggplot(results_3bp_sum, aes(x=Event, y=sum)) + geom_bar(stat="identity") + theme_bw() + ylab("Fraction of Indel Reads") + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + facet_wrap(~design, ncol=1)
        ggsave(paste0(outDir,"/3bp_insertions_barplot.sum.pdf"),width=7,height=8)

    }
    res <- results %>% dplyr::filter(reference_sequence %in% c(refseq_244bp_3bp_shiftCT, refseq_100bp_3bp_shiftCT))
    #make_plots_shifted_CT_3bp_to_left(res, "myoutput/insertions_shifted_CT_3bp_to_left")


    # for 3bp insertions all together
    make_plots_3bp <- function(res, outDir){

        dir.create(outDir, recursive=TRUE)

        # filter out siRev3 pMJ119
        res <- res %>% dplyr::filter(design != 'siRev3 pMJ119') %>% arrange(Class,design,Sample) # let's focus on siNT pMJ119 and C

        # filter out "Other events" and "N" in Event
        results_3bp <- res %>% dplyr::filter(Class=="3bp Insertion", Event!="Other events", !grepl('N',Event))
        results_3bp %>% group_by(Event,Sample) %>% summarise(n=n()) %>% arrange(-n) %>% fwrite(paste0(outDir,"/3bp_insertions_count.tsv"), sep="\t")
        results_3bp %>% fwrite(paste0(outDir,"/3bp_insertions.tsv"), sep="\t")

        # convert Event, Sample, and "Fraction of Indel Reads" columns to a matrix
        results_3bp_matrix <- results_3bp %>% reshape2::acast(Event ~ Sample, value.var="Fraction of Indel Reads")

        # filing NA with 0
        results_3bp_matrix[is.na(results_3bp_matrix)] <- 0

        # order samples by design
        samples_3bp_C    <- results_3bp %>% dplyr::filter(design=="C") %>% dplyr::select(Sample) %>% dplyr::pull() %>% unique()
        samples_3bp_siNT <- results_3bp %>% dplyr::filter(design=="siNT pMJ119") %>% dplyr::select(Sample) %>% dplyr::pull() %>% unique()

        results_3bp_matrix <- cbind(results_3bp_matrix[,samples_3bp_C], results_3bp_matrix[,samples_3bp_siNT])

        results_3bp_matrix %>% as.data.frame %>% rownames_to_column(var="Event") %>% dplyr::mutate(Insertion=gsub(" .*","",Event)) %>% dplyr::select(-Event) %>%
            dplyr::select(Insertion,everything()) %>% 
            fwrite(paste0(outDir,"/3bp_insertions_matrix.tsv"), sep="\t")


        # make heatmap for 3bp insertions
        ht <- Heatmap(results_3bp_matrix, name="Fraction of Indel Reads", col = colorRamp2(c(0,0.001, 0.01,0.1), c("white", "orange","red","blue")), 
            cluster_rows = FALSE, cluster_columns = FALSE, border=TRUE)

        pdf(paste0(outDir,"/3bp_insertions_heatmap.pdf"), width=10, height=10)
        draw(ht)
        dev.off()


        ## -------------------------------------------------------------------------------------------------------------------
        # for each type of insertions, make barplot for mean of "Fraction of Indel Reads" 
        # Method1: take mean of replicates for each insertion, then sum up all insertion events
        results_3bp_mean <- results_3bp %>% dplyr::group_by(design,Event) %>% dplyr::summarize(mean=mean(`Fraction of Indel Reads`),sd=sd(`Fraction of Indel Reads`),n=n())
        results_3bp_mean %>% fwrite(paste0(outDir,"/3bp_insertions_barplot.mean.tsv"), sep="\t")

        p <- ggplot(results_3bp_mean, aes(x=Event, y=mean)) + geom_bar(stat="identity") + theme_bw() + ylab("Fraction of Indel Reads") + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
            geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2, position=position_dodge(.9)) + 
            scale_y_continuous(trans=scales::pseudo_log_trans(sigma=1e-6,base=10),breaks=c(-1e4,-1e3,-1e2,-1e1,-1e0,-1e-1,-1e-2,-1e-3,-1e-4,-1e-5,0,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7)) +
            facet_wrap(~design, ncol=1)
        ggsave(paste0(outDir,"/3bp_insertions_barplot.mean.log10.pdf"),width=10,height=8)

        p <- ggplot(results_3bp_mean, aes(x=Event, y=mean)) + geom_bar(stat="identity") + theme_bw() + ylab("Fraction of Indel Reads") + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + 
            geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2, position=position_dodge(.9)) +
            facet_wrap(~design, ncol=1)
        ggsave(paste0(outDir,"/3bp_insertions_barplot.mean.pdf"),width=10,height=8)

        # excluding donor "GAT Insertion", sum up all combinations for C and siNT pMJ119
        results_3bp_mean_summary <- results_3bp_mean %>% dplyr::filter(Event != 'GAT Insertion') %>% dplyr::group_by(design) %>% dplyr::summarize(sum=sum(mean),mean=mean(mean),sd=sd(mean))

        p <- ggplot(results_3bp_mean_summary, aes(x=design, y=sum)) + geom_bar(stat="identity") + theme_bw() + ggtitle("3bp") + ylab("Fraction of Indel Reads") + xlab("") + geom_text(aes(label=round(sum,4)), vjust=-0.5, size=3)
        ggsave(paste0(outDir,"/3bp_insertions_barplot.m1_sum_for_all_combinations.pdf"),width=5,height=5)

        p <- ggplot(results_3bp_mean_summary, aes(x=design, y=mean)) + geom_bar(stat="identity") + theme_bw() + ggtitle("3bp") + ylab("Fraction of Indel Reads") + xlab("") + geom_text(aes(label=round(mean,4)), vjust=-0.5, size=3)
        ggsave(paste0(outDir,"/3bp_insertions_barplot.m1_mean_for_all_combinations.pdf"),width=5,height=5)


        ## -------------------------------------------------------------------------------------------------------------------
        # for each type of insertions, make barplot for mean of "Fraction of Indel Reads" 
        # Method2: take mean of insertion events for each replicate, then sum up all replicates
        results_3bp_mean_m2 <- results_3bp %>% dplyr::filter(Event != 'GAT Insertion') %>% dplyr::group_by(design,Sample) %>% dplyr::summarize(mean=mean(`Fraction of Indel Reads`),sd=sd(`Fraction of Indel Reads`),n=n())
        results_3bp_mean_m2 %>% fwrite(paste0(outDir,"/3bp_insertions_barplot.mean.m2.tsv"), sep="\t")

        # excluding donor "GAT Insertion", sum up all combinations for C and siNT pMJ119
        results_3bp_mean_m2_summary <- results_3bp_mean_m2 %>% dplyr::group_by(design) %>% dplyr::summarize(sum=sum(mean),mean=mean(mean),sd=sd(mean))

        p <- ggplot(results_3bp_mean_m2_summary, aes(x=design, y=sum)) + geom_bar(stat="identity") + theme_bw() + ggtitle("3bp") + ylab("Fraction of Indel Reads") + xlab("") + geom_text(aes(label=round(sum,4)), vjust=-0.5, size=3)
        ggsave(paste0(outDir,"/3bp_insertions_barplot.m2_sum_for_all_combinations.pdf"),width=5,height=5)

        p <- ggplot(results_3bp_mean_m2_summary, aes(x=design, y=mean)) + geom_bar(stat="identity") + theme_bw() + ggtitle("3bp") + ylab("Fraction of Indel Reads") + xlab("") + geom_text(aes(label=round(mean,4)), vjust=-0.5, size=3)
        ggsave(paste0(outDir,"/3bp_insertions_barplot.m2_mean_for_all_combinations.pdf"),width=5,height=5)

    }
    res <- results
    make_plots_3bp(res, "myoutput/insertions_3bp")


    # for 2bp insertions all together
    make_plots_2bp <- function(res, outDir){

        dir.create(outDir, recursive=TRUE)

        # filter out siRev3 pMJ119
        res <- res %>% dplyr::filter(design != 'siRev3 pMJ119') %>% arrange(Class,design,Sample) # let's focus on siNT pMJ119 and C

        # filter out "Other events" and "N" in Event
        results_2bp <- res %>% dplyr::filter(Class=="2bp Insertion", Event!="Other events", !grepl('N',Event))

        # convert Event, Sample, and "Fraction of Indel Reads" columns to a matrix
        results_2bp_matrix <- results_2bp %>% reshape2::acast(Event ~ Sample, value.var="Fraction of Indel Reads")

        # filing NA with 0
        results_2bp_matrix[is.na(results_2bp_matrix)] <- 0

        # order samples by design
        samples_2bp_C    <- results_2bp %>% dplyr::filter(design=="C") %>% dplyr::select(Sample) %>% dplyr::pull() %>% unique()
        samples_2bp_siNT <- results_2bp %>% dplyr::filter(design=="siNT pMJ119") %>% dplyr::select(Sample) %>% dplyr::pull() %>% unique()

        results_2bp_matrix <- cbind(results_2bp_matrix[,samples_2bp_C], results_2bp_matrix[,samples_2bp_siNT])

        results_2bp_matrix %>% as.data.frame %>% rownames_to_column(var="Event") %>% dplyr::mutate(Insertion=gsub(" .*","",Event)) %>% dplyr::select(-Event) %>%
            dplyr::select(Insertion,everything()) %>% 
            fwrite(paste0(outDir,"/2bp_insertions_matrix.tsv"), sep="\t")


        # make heatmap for 2bp insertions
        print(summary(as.vector(results_2bp_matrix)))
        ht <- Heatmap(results_2bp_matrix, name="Fraction of Indel Reads", col = colorRamp2(c(0,0.01, 0.1,1), c("white", "orange","red",'blue')), 
            cluster_rows = FALSE, cluster_columns = FALSE, border=TRUE)

        pdf(paste0(outDir,"/2bp_insertions_heatmap.pdf"), width=10, height=10)
        draw(ht)
        dev.off()


        ## -------------------------------------------------------------------------------------------------------------------
        # for each type of insertions, make barplot for mean of "Fraction of Indel Reads" 
        # Method1: take mean of replicates for each insertion, then sum up all insertion events
        results_2bp_mean <- results_2bp %>% dplyr::group_by(design,Event) %>% dplyr::summarize(mean=mean(`Fraction of Indel Reads`), sd=sd(`Fraction of Indel Reads`), n=n())
        results_2bp_mean %>% fwrite(paste0(outDir,"/2bp_insertions_barplot.mean.tsv"), sep="\t")

        p <- ggplot(results_2bp_mean, aes(x=Event, y=mean)) + geom_bar(stat="identity") + theme_bw() + ylab("Fraction of Indel Reads") + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + 
            geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2, position=position_dodge(.9)) + 
            scale_y_continuous(trans=scales::pseudo_log_trans(sigma=1e-6,base=10),breaks=c(-1e4,-1e3,-1e2,-1e1,-1e0,-1e-1,-1e-2,-1e-3,-1e-4,-1e-5,0,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7)) +
            facet_wrap(~design, ncol=1)
        ggsave(paste0(outDir,"/2bp_insertions_barplot.mean.log10.pdf"),width=10,height=8)

        p <- ggplot(results_2bp_mean, aes(x=Event, y=mean)) + geom_bar(stat="identity") + theme_bw() + ylab("Fraction of Indel Reads") + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + 
            geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2, position=position_dodge(.9)) + 
            facet_wrap(~design, ncol=1)
        ggsave(paste0(outDir,"/2bp_insertions_barplot.mean.pdf"),width=10,height=8)

        # excluding donor "GAT Insertion", sum up all combinations for C and siNT pMJ119
        results_2bp_mean_summary <- results_2bp_mean %>% dplyr::group_by(design) %>% dplyr::summarize(sum=sum(mean),mean=mean(mean),sd=sd(mean))

        p <- ggplot(results_2bp_mean_summary, aes(x=design, y=sum)) + geom_bar(stat="identity") + theme_bw() + ggtitle("2bp") + ylab("Fraction of Indel Reads") + xlab("") + geom_text(aes(label=round(sum,4)), vjust=-0.5, size=3)
        ggsave(paste0(outDir,"/2bp_insertions_barplot.m1_sum_for_all_combinations.pdf"),width=5,height=5)

        p <- ggplot(results_2bp_mean_summary, aes(x=design, y=mean)) + geom_bar(stat="identity") + theme_bw() + ggtitle("2bp") + ylab("Fraction of Indel Reads") + xlab("") + geom_text(aes(label=round(mean,4)), vjust=-0.5, size=3)
        ggsave(paste0(outDir,"/2bp_insertions_barplot.m1_mean_for_all_combinations.pdf"),width=5,height=5)


        ## -------------------------------------------------------------------------------------------------------------------
        # for each type of insertions, make barplot for mean of "Fraction of Indel Reads" 
        # Method2: take mean of insertion events for each replicate, then sum up all replicates
        results_2bp_mean_m2 <- results_2bp %>% dplyr::group_by(design,Sample) %>% dplyr::summarize(mean=mean(`Fraction of Indel Reads`),sd=sd(`Fraction of Indel Reads`),n=n())
        results_2bp_mean_m2 %>% fwrite(paste0(outDir,"/2bp_insertions_barplot.mean.m2.tsv"), sep="\t")

        # excluding donor "GAT Insertion", sum up all combinations for C and siNT pMJ119
        results_2bp_mean_m2_summary <- results_2bp_mean_m2 %>% dplyr::group_by(design) %>% dplyr::summarize(sum=sum(mean),mean=mean(mean),sd=sd(mean))

        p <- ggplot(results_2bp_mean_m2_summary, aes(x=design, y=sum)) + geom_bar(stat="identity") + theme_bw() + ggtitle("2bp") + ylab("Fraction of Indel Reads") + xlab("") + geom_text(aes(label=round(sum,4)), vjust=-0.5, size=3)
        ggsave(paste0(outDir,"/2bp_insertions_barplot.m2_sum_for_all_combinations.pdf"),width=5,height=5)

        p <- ggplot(results_2bp_mean_m2_summary, aes(x=design, y=mean)) + geom_bar(stat="identity") + theme_bw() + ggtitle("2bp") + ylab("Fraction of Indel Reads") + xlab("") + geom_text(aes(label=round(mean,4)), vjust=-0.5, size=3)
        ggsave(paste0(outDir,"/2bp_insertions_barplot.m2_mean_for_all_combinations.pdf"),width=5,height=5)


    }
    res <- results
    make_plots_2bp(res, "myoutput/insertions_2bp")


    # for 4bp insertions all together
    make_plots_4bp <- function(res, outDir){

        dir.create(outDir, recursive=TRUE)

        # filter out siRev3 pMJ119
        res <- res %>% dplyr::filter(design != 'siRev3 pMJ119') %>% arrange(Class,design,Sample) # let's focus on siNT pMJ119 and C

        # filter out "Other events" and "N" in Event
        results_4bp <- res %>% dplyr::filter(Class=="4bp Insertion", Event!="Other events", !grepl('N',Event))

        # convert Event, Sample, and "Fraction of Indel Reads" columns to a matrix
        results_4bp_matrix <- results_4bp %>% reshape2::acast(Event ~ Sample, value.var="Fraction of Indel Reads")

        # filing NA with 0
        results_4bp_matrix[is.na(results_4bp_matrix)] <- 0

        # order samples by design
        samples_4bp_C    <- results_4bp %>% dplyr::filter(design=="C") %>% dplyr::select(Sample) %>% dplyr::pull() %>% unique()
        samples_4bp_siNT <- results_4bp %>% dplyr::filter(design=="siNT pMJ119") %>% dplyr::select(Sample) %>% dplyr::pull() %>% unique()

        results_4bp_matrix <- cbind(results_4bp_matrix[,samples_4bp_C], results_4bp_matrix[,samples_4bp_siNT])

        results_4bp_matrix %>% as.data.frame %>% rownames_to_column(var="Event") %>% dplyr::mutate(Insertion=gsub(" .*","",Event)) %>% dplyr::select(-Event) %>%
            dplyr::select(Insertion,everything()) %>% 
            fwrite(paste0(outDir,"/4bp_insertions_matrix.tsv"), sep="\t")


        # make heatmap for 4bp insertions
        print(summary(as.vector(results_4bp_matrix)))
        ht <- Heatmap(results_4bp_matrix, name="Fraction of Indel Reads", col = colorRamp2(c(0,0.00001,0.0001,0.001), c("white", "orange","red","blue")), 
            cluster_rows = FALSE, cluster_columns = FALSE, border=TRUE)

        pdf(paste0(outDir,"/4bp_insertions_heatmap.pdf"), width=10, height=10)
        draw(ht)
        dev.off()


        ## -------------------------------------------------------------------------------------------------------------------
        # for each type of insertions, make barplot for mean of "Fraction of Indel Reads" 
        # Method1: take mean of replicates for each insertion, then sum up all insertion events
        results_4bp_mean <- results_4bp %>% dplyr::group_by(design,Event) %>% dplyr::summarize(mean=mean(`Fraction of Indel Reads`), sd=sd(`Fraction of Indel Reads`), n=n())
        results_4bp_mean %>% fwrite(paste0(outDir,"/4bp_insertions_barplot.mean.tsv"), sep="\t")

        p <- ggplot(results_4bp_mean, aes(x=Event, y=mean)) + geom_bar(stat="identity") + theme_bw() + ylab("Fraction of Indel Reads") + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + 
            geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2, position=position_dodge(.9)) +
            scale_y_continuous(trans=scales::pseudo_log_trans(sigma=1e-6,base=10),breaks=c(-1e4,-1e3,-1e2,-1e1,-1e0,-1e-1,-1e-2,-1e-3,-1e-4,-1e-5,0,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7)) +
            facet_wrap(~design, ncol=1)
        ggsave(paste0(outDir,"/4bp_insertions_barplot.mean.log10.pdf"),width=10,height=8)

        p <- ggplot(results_4bp_mean, aes(x=Event, y=mean)) + geom_bar(stat="identity") + theme_bw() + ylab("Fraction of Indel Reads") + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + 
            geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2, position=position_dodge(.9)) +
            facet_wrap(~design, ncol=1)
        ggsave(paste0(outDir,"/4bp_insertions_barplot.mean.pdf"),width=10,height=8)

        # excluding donor "GAT Insertion", sum up all combinations for C and siNT pMJ119
        results_4bp_mean_summary <- results_4bp_mean %>% dplyr::group_by(design) %>% dplyr::summarize(sum=sum(mean),mean=mean(mean),sd=sd(mean))

        p <- ggplot(results_4bp_mean_summary, aes(x=design, y=sum)) + geom_bar(stat="identity") + theme_bw() + ggtitle("4bp") + ylab("Fraction of Indel Reads") + xlab("") + geom_text(aes(label=round(sum,4)), vjust=-0.5, size=3)
        ggsave(paste0(outDir,"/4bp_insertions_barplot.sum_for_all_combinations.pdf"),width=5,height=5)

        p <- ggplot(results_4bp_mean_summary, aes(x=design, y=mean)) + geom_bar(stat="identity") + theme_bw() + ggtitle("4bp") + ylab("Fraction of Indel Reads") + xlab("") + geom_text(aes(label=round(mean,5)), vjust=-0.5, size=3)
        ggsave(paste0(outDir,"/4bp_insertions_barplot.mean_for_all_combinations.pdf"),width=5,height=5)


        ## -------------------------------------------------------------------------------------------------------------------
        # for each type of insertions, make barplot for mean of "Fraction of Indel Reads" 
        # Method2: take mean of insertion events for each replicate, then sum up all replicates
        results_4bp_mean_m2 <- results_4bp %>% dplyr::group_by(design,Sample) %>% dplyr::summarize(mean=mean(`Fraction of Indel Reads`),sd=sd(`Fraction of Indel Reads`),n=n())
        results_4bp_mean_m2 %>% fwrite(paste0(outDir,"/4bp_insertions_barplot.mean.m2.tsv"), sep="\t")

        # excluding donor "GAT Insertion", sum up all combinations for C and siNT pMJ119
        results_4bp_mean_m2_summary <- results_4bp_mean_m2 %>% dplyr::group_by(design) %>% dplyr::summarize(sum=sum(mean),mean=mean(mean),sd=sd(mean))

        p <- ggplot(results_4bp_mean_m2_summary, aes(x=design, y=sum)) + geom_bar(stat="identity") + theme_bw() + ggtitle("4bp") + ylab("Fraction of Indel Reads") + xlab("") + geom_text(aes(label=round(sum,4)), vjust=-0.5, size=3)
        ggsave(paste0(outDir,"/4bp_insertions_barplot.m2_sum_for_all_combinations.pdf"),width=5,height=5)

        p <- ggplot(results_4bp_mean_m2_summary, aes(x=design, y=mean)) + geom_bar(stat="identity") + theme_bw() + ggtitle("4bp") + ylab("Fraction of Indel Reads") + xlab("") + geom_text(aes(label=round(mean,5)), vjust=-0.5, size=3)
        ggsave(paste0(outDir,"/4bp_insertions_barplot.m2_mean_for_all_combinations.pdf"),width=5,height=5)


    }
    res <- results
    make_plots_4bp(res, "myoutput/insertions_4bp")


}

main()