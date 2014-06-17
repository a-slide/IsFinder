#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio import pairwise2, SeqIO
from multiprocessing import Pool, cpu_count
from functools import partial
from pprint import pprint as pp
import gzip

MATCH_LIST = []

def main ():
    adapter_list = ["gcgcg", "gcgcgcg", "atatat", "taggcg", "tgtaa"]
    min_size = 50
    min_qual = 32

    print("Uncompressing and extracting data")
    R1 = gzip.open("../test/fastq/small.R1.fastq.gz", "r")
    R2 = gzip.open("../test/fastq/small.R2.fastq.gz", "r")

    genR1 = SeqIO.parse(R1, "fastq")
    genR2 = SeqIO.parse(R2, "fastq")

    ListR1 = []
    ListR2 = []

    Total = 0
    Qual_passed = 0
    Trim_passed = 0

    # Parsing file and filtering sequences matching quality requirements
    while True:
        seqR1 = next(genR1, None)
        seqR2 = next(genR2, None)

        # End of file
        if not seqR1 or not seqR2:
            print ("Total sequence = {}".format(Total))
            print ("Passed quality filter = {}".format(Qual_passed))
            print ("Passed adapter trimming = {}".format(Trim_passed))
            break
        Total +=1

        # Quality filter
        if _rec_qual(seqR1) < min_qual or _rec_qual(seqR1) < min_qual:
            continue
        Qual_passed +=1

        # Adapter trimming and size selection
        seqR1 = trim_adapter(seqR1, adapter_list, min_size)
        seqR2 = trim_adapter(seqR2, adapter_list, min_size)
        if not seqR1 or not seqR2:
            continue
        Trim_passed +=1

        # If all test passed = append the sequence to the list
        ListR1.append(seqR1)
        ListR1.append(seqR2)

    R1.close()
    R2.close()

    for fastq in ListR1:
        print ("{} : {}\n".format(fastq.id, fastq.seq))

    for fastq in ListR1:
        print ("{} : {}\n".format(fastq.id, fastq.seq))

def _rec_qual(record):
    """
    Compute mean quality score
    """
    return sum(record.letter_annotations['phred_quality'])/len(record)

def trim_adapter (record, adapter_list, min_size):

    MATCH_LIST = []
    # Define the number of parallel thread to use
    p = Pool(cpu_count()+1)
    # Define an helper function with partial to allow dealing with multi arg
    partial_align = partial(align, str(record.seq).lower())
    # Parallel processing of parwise alignments using Pool map
    match_list = p.map(partial_align, adapter_list)
    pp(match_list)

    coverage = [0 for i in range(maxsize)]

    for match in match_list:
        for seqA, seqB, score, begin, end in match:
            for i in range (begin, end):
                coverage[i] = 1

    if not match_list:
        return record

    start, end = longer_interval (match_list, len(record))
    return record[start:end] if end-start >= min_size else None

def align (adapter, sequence):
    """
    Pairwise alignment using the buildin Biopython function pairwise2
    """
    match_list = pairwise2.align.localms(sequence, adapter, 2, -1, -1, -.1)

    if match_list:
        for seqA, seqB, score, begin, end in match_list
            MATCH_LIST.append([begin, end]):


def longer_interval(match_list, maxsize):
    """
    Find the larger interval that do not overlapp an adapter location
    """
    coverage = [0 for i in range(maxsize)]

    for seqA, seqB, score, begin, end in match_list:
        print begin, end
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
    print ("Longer interval = {} [{}:{}]".format(inter_max, start_max+1, end_max-1))

    return start_max, end_max

if __name__ == '__main__':
    main()
