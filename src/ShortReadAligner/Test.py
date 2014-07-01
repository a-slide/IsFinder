# Starndard library imports
from os import remove

# Local imports

#from AlignerSplitter import AlignerSplitter
#from SamSpliter import SamSpliter
from FastaReader import FastaReader
from BwaIndexer import GenerateIndex, ExistingIndex
from BwaAligner import BwaAligner


def test ():
    
    bwa_path = "bwa"
    ref1_path = "../../test/references/rAAV_genome.fa"
    ref2_path = "../../test/references/chr19_HRDEL.fa.gz"
    index_path = "./bwa_index/tmp1OQD4s.gz"
    R1_path = "../../test/fastq/100seq_R1.fastq.gz"
    R2_path = "../../test/fastq/100seq_R2.fastq.gz"
    index_opt = {}
    # Option -M = mark secondary hits as not primary aligned in sam flag
    align_opt = {"M":None}
    spliter_opt = {}
    main_opt = {}
    
    
    try:
        if not index_path:
            raise Exception("No index provided. An index will be generated")
        print("Existing index provided")
        FastaRef = FastaReader(ref1_path, ref2_path, False)
        Index = ExistingIndex(bwa_path, index_path)
        
    except Exception as E:
        print (E)
        print("Generating index...")
        FastaRef = FastaReader(ref1_path, ref2_path, True)
        Index = GenerateIndex(bwa_path, FastaRef.merge_ref, index_opt)
        remove (FastaRef.merge_ref)
    
    
    Aligner = BwaAligner(bwa_path, Index, align_opt)
    
    return (Aligner.align_fifo(R1_path, R2_path))
    
    #Spliter = SamSpliter(FastaRef.ref1_list, FastaRef.ref2_list, spliter_opt)
    #Main = (Aligner, Spliter, main_opt)

    #Main.align_and_split()
