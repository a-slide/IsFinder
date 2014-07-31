"""
@package AlignerSpliter
@brief ...
@copyright [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

#~~~~~~~GLOBAL IMPORTS~~~~~~~#

# Standard library packages import
from os import remove, path

# Local library packages import
#from SamSpliter import SamSpliter
from FastaReader import FastaReader
from BWA.BwaIndexer import GenerateIndex, ExistingIndex
from BWA.BwaAligner import BwaAligner


def align  (source_list,
            R1,
            R2=None,
            index = None,
            ref_outdir = "./reference/"
            bwa_mem = "bwa mem",
            align_opt=None,
            align_outdir= "./bwa_align/",
            bwa_index = "bwa index",
            index_opt=None,
            index_outdir = "./bwa_index/")
    """
    Main function of RefMasker that integrate database creation, blast and homology masking
    * Instantiate Blast database and blastn object
    * Perform iterative blasts of query sequences against the subject database and create a list of
    hits.
    """

    # Try to validate a index from an existing one
    try:
        if not index_path:
            raise Exception("No index provided. An index will be generated")

        print("Existing index provided")
        FastaRef = FastaReader(ref1_path, ref2_path, write_merge=False)
        Index = ExistingIndex(bwa_path, index_path)

    # If no index or if an error occured during validation of the existing index = create a new one
    except Exception as E:
        print (E)

        print("Merge References...")
        mkdir(ref_outdir)

        FastaRef = FastaReader([ref1_path,ref2_path], write_merge=True, output="merged.fa")

        print("Generating index...")
        mkdir(db_outdir)
        Index = GenerateIndex(bwa_path, FastaRef.merge_ref, index_opt)
        remove (FastaRef.merge_ref)

    Aligner = BwaAligner(bwa_path, Index, align_opt)

    return (Aligner.align(R1_path, R2_path))

    #Spliter = SamSpliter(FastaRef.ref1_list, FastaRef.ref2_list, spliter_opt)
    #Main = (Aligner, Spliter, main_opt)

    #Main.align_and_split()
