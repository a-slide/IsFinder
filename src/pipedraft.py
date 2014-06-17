#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio import pairwise2, SeqIO
from Utilities import file_basename
import gzip


def quality_control ():

    R1_path = "../test/fastq/med_R1.fastq.gz"
    R2_path = "../test/fastq/med_R2.fastq.gz"
    adapter_list = ["gcgcgcg", "atattaat", "tttaggcg", "aataccgattc"]
    min_size = 50
    min_qual = 31
    qual_scale = "fastq-sanger"

    try:
        print("Uncompressing and extracting data")
        # Input fastq files
        input_R1 = gzip.open(R1_path, "r")
        input_R2 = gzip.open(R2_path, "r")
        # Output fastq files
        output_R1 = gzip.open(file_basename(R1_path)+"_trim.fastq.gz", 'w')
        output_R2 = gzip.open(file_basename(R2_path)+"_trim.fastq.gz", 'w')
        # Initialize a generator to iterate over fastq files
        genR1 = SeqIO.parse(input_R1, qual_scale)
        genR2 = SeqIO.parse(input_R2, qual_scale)
        # Initialize variables
        ListR1 = []
        ListR2 = []

        print("Parsing files and filtering sequences")
        # Parsing files and filtering sequences matching quality requirements
        while True:
            seqR1 = next(genR1, None)
            seqR2 = next(genR2, None)

            # End of file
            if not seqR1 or not seqR2:
                #print_stat ()
                break

            #TOTAL += 1

            # Quality filter
            if _rec_qual(seqR1) < min_qual or _rec_qual(seqR1) < min_qual:
                continue
            #QUAL_PASSED += 1

            # Adapter trimming and size selection
            seqR1 = trim_adapter(seqR1, adapter_list, min_size)
            seqR2 = trim_adapter(seqR2, adapter_list, min_size)
            if not seqR1 or not seqR2:
                continue
            #LEN_PASSED += 1

            # If all test passed = append the sequence to the list
            output_R1.write (seqR1.format (qual_scale))
            output_R2.write (seqR2.format (qual_scale))

        input_R1.close()
        input_R2.close()
        output_R1.close()
        output_R2.close()

    except IOError as E:
        print (E)
        exit (0)
    return 1

def _rec_qual(record):
    """
    Compute mean quality score
    """
    return sum(record.letter_annotations['phred_quality'])/len(record)

def trim_adapter (record, adapter_list, min_size):

    match_list = []
    for adapter in adapter_list:
        match_list.extend(align(adapter, str(record.seq).lower()))

    if not match_list:
        return record

    #print ("{} adapter hit(s) found in {}".format(len(match_list), record.id))
    start, end = longer_interval (match_list, len(record))
    #TRIMMED_SEQ +=1
    #TRIMMED_BASE += (len(record)-(end-start))

    return record[start:end] if end-start >= min_size else None

def align (adapter, sequence):
    """
    Pairwise alignment using the buildin Biopython function pairwise2
    """
    pairwise2.MAX_ALIGNMENTS = 5
    match = 2
    mismatch = -2
    open_gap = -2
    extend_gap = -2
    cutoff = 1.5

    match_list = pairwise2.align.localms(adapter, sequence, match, mismatch, open_gap, extend_gap)

    #for seqA, seqB, score, begin, end in match_list:
        #if score/len(adapter) > cutoff:
            #print("{}\n{}\nScore:{}\tBegin:{}\tEnd:{}".format(seqA, seqB, score, begin, end))

    return [[begin, end] for seqA, seqB, score, begin, end in match_list if score/len(adapter) > cutoff]


def longer_interval(match_list, maxsize):
    """
    Find the larger interval that do not overlapp an adapter location
    """

    coverage = [0 for i in range(maxsize)]

    for begin, end in match_list:
        for i in range (begin, end):
            coverage[i] = 1

    start_max = end_max = inter_max = start = inter = 0

    for i in range(maxsize):
        if coverage[i]:
            start = i+1
            inter = 0

        else:
            inter += 1
            if inter > inter_max:
                inter_max = inter
                start_max = start
                end_max = i
    #print ("Longer interval = {} [{}:{}]".format(inter_max, start_max+1, end_max-1))

    return start_max, end_max

#def print_stat ():
    #print ("Total sequence = {}".format(TOTAL))
    #print ("Passed quality filter = {}".format(QUAL_PASSED))
    #print ("Untrimmed sequences = {}".format(QUAL_PASSED - TRIMMED_SEQ))
    #print ("Trimmed sequences = {}".format(TRIMMED_SEQ))
    #print ("Trimmed bases = {}".format(TRIMMED_BASE))
    #print ("Passed length filter after adapter trimming = {}".format(LEN_PASSED))
    #return 1

if __name__ == '__main__':
    TOTAL = 0
    QUAL_PASSED = 0
    TRIMMED_SEQ = 0
    TRIMMED_BASE = 0
    LEN_PASSED = 0

    quality_control()


############# FUCKING STRANGE BEHAVIOUR WITH GLOBAL VARIABLES ?????
