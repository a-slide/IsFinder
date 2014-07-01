
#~~~~~~~GLOBAL IMPORTS~~~~~~~#
# Standard library packages import
import gzip
from tempfile import mkstemp
from os.path import getsize, isfile
from sys import exit as sys_exit
from os import remove

# Third party package import
from Bio import SeqIO

# Local Package import
from Utilities import file_extension

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class FastaReader(object):
    """
    @class  FastaReader
    @brief
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __repr__(self):
        msg = self.__str__()
        msg += "Sequences in Reference 1:\n\t{}\n".format("\n\t".join([ref for ref in self.ref1_list]))
        msg += "Sequences in Reference 2:\n\t{}\n".format("\n\t".join([ref for ref in self.ref2_list]))
        if self.merge_ref:
            msg += "Path of the merged reference:\n\t{}\n".format(self.merge_ref)
        return msg

    def __str__(self):
        return "\n<Instance of {} from {}>\n".format(self.__class__.__name__, self.__module__)

    def __init__ (self, ref1, ref2, write_merge=False):
        """
        @param ref1 Path of the fasta file corresponding to the first reference
        @param ref2 Path of the fasta file corresponding to the second reference
        @param write_merge If True a merged fasta file containing both references will be
        generated in a temporary file
        """
        self.ref1 = ref1
        self.ref2 = ref2
        self.merge_ref = ""
        self.ref1_list = []
        self.ref2_list = []

        if write_merge:
            self._read_and_merge()
        else:
            self._simple_read()

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _read_and_merge(self):
        """
        Read fasta file ref1 and ref2, store names of sequences in lists and write a merged fasta
        reference in a temporary file
        """
        try:
            # Create a temp file for storing merged fasta
            OSlevel, output = mkstemp(suffix=".gz")
            outfile = open(output, "w")

            for ref, ref_list in [[self.ref1, self.ref1_list], [self.ref2, self.ref2_list]]:
                print ("Start parsing {}...".format(ref))
                # Open gzipped fasta or uncompressed fasta
                if file_extension(ref) == "gz" or file_extension(ref) == "GZ":
                    infile = gzip.open(ref, "r")
                else:
                    infile = open(ref, "r")
                    
                # Parse the file
                for sequence in SeqIO.parse(infile, "fasta"):
                    print("\t"+sequence.id)
                    # Fill the list and write the sequence in the merged reference
                    self._ckeck_duplicate (sequence.id)
                    ref_list.append(sequence.id)
                    outfile.write(sequence.format("fasta"))
                               
            # Close and verify the output file
            outfile.close()
            # Verify if the list was filled if the file was generated
            assert self.ref1_list, "No sequences were imported from {}".format(ref1)
            assert self.ref2_list, "No sequences were imported from {}".format(ref2)
            assert isfile(output), "No output file found"
            assert getsize(output) > 0, "The output file is empty"
            
        except Exception as E:
            print(E)
            print("Removing temporary file")
            try:
                if isfile (output):
                    remove (output)
            except Exception:
                pass
            sys_exit()
        
        # Store the path of the tempory file in the object variable : merge_ref
        self.merge_ref = output
        return True
        

    def _simple_read(self):
        """
        Read fasta file ref1 and ref2 and store names of sequences in lists
        """
        try:

            # Parse fastq file from ref1 and ref2
            for ref, ref_list in [[self.ref1, self.ref1_list], [self.ref2, self.ref2_list]]:
                print ("Start parsing {}...".format(ref))
                # Open gzipped fasta or uncompressed fasta
                if file_extension(ref) == "gz" or file_extension(ref) == "GZ":
                    infile = gzip.open(ref, "r")
                else:
                    infile = open(ref, "r")
                    
                # Parse the file
                for sequence in SeqIO.parse(infile, "fasta"):
                    print("\t"+sequence.id)
                    # Fill the list
                    self._ckeck_duplicate (sequence.id)
                    ref_list.append(sequence.id)
                
            # Verify if the list was filled
            assert self.ref1_list, "No sequences were imported from {}".format(ref1)
            assert self.ref2_list, "No sequences were imported from {}".format(ref2)
            
        except Exception as E:
            print(E)
            sys_exit()
        
        return True
    
    
    def _ckeck_duplicate (self, seqid):
        """
        Raise an error if a duplicated sequence identifier is found
        """
        if seqid in self.ref1_list:
            raise ValueError("{} is was already found in {}\n".format(seqid, self.ref1))
        if seqid in self.ref2_list:
            raise ValueError("{} is was already found in {}\n".format(seqid, self.ref2))
        
        return True
