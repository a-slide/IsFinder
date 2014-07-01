#from BwaAligner import BwaAligner
#from AlignerSplitter import AlignerSplitter
#from BwaIndexer import GenerateIndex, ExistingIndex
#from SamSpliter import SamSpliter
from FastaReader import FastaReader

def test ():
    
    bwa_path = "bwa"
    ref1_path = "../../test/references/rAAV_genome.fa"
    ref2_path = "../../test/references/"
    index_path = ""
    R1_path = ""
    R2_path = ""
    index_opt = {}
    align_opt = {}
    spliter_opt = {}
    main_opt = {}

    if index_path:
        FastaRef = FastaReader(ref1, ref2, False)
        print repr((FastaRef))
        #Index = ExistingIndex(index_path)

    else:
        FastaRef = FastaReader(ref1_path, ref2_path, True)
        print repr((FastaRef))
        #Index = GenerateIndex(bwa_path, FastaRef.merge_ref, index_opt)

    #Aligner = BwaAligner(bwa_path, Index)
    #Spliter = SamSpliter(FastaRef.ref1_list, FastaRef.ref2_list, spliter_opt)
    #Main = (Aligner, Spliter, main_opt)

    #Main.align_and_split()
