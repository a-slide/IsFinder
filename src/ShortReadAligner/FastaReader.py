
#~~~~~~~GLOBAL IMPORTS~~~~~~~#
# Standard library packages import
import gzip
from tempfile import mkstemp


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
        return "{}\n".format(self.__str__())

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def __init__ (self, ref1, ref2, write_merge=False):
        """
        @param ref1
        @param ref2
        @param write_merge
        """
        self.ref1 = ref1
        self.ref2 = ref2
        self.merge_ref = None
        self.ref1_list = []
        self.ref2_list = []

        if write_merge:
            self._read_merge_fasta()
        else:
            self._read_fasta()

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _read_and_merge():

        try:
            # Create a temp file for storing merged fasta
            OSlevel, output = mkstemp(suffix=".gz")
            outfile = gzip.open(output, "wb")

            # Parse fastq file from ref1 and ref2
            for ref, ref_list in [[self.ref1, self.ref1_list], [self.ref2, self.ref2_list]]

                # Open gzipped fasta or uncompressed fasta
                try:
                    if file_extension(ref) == "gz" or file_extension(ref) == "GZ":
                        infile = gzip.open(ref, "r")
                    else:
                        infile = open(ref, "r")
                    # Parse the file
                    for sequence in SeqIO.parse(infile, "fasta"):
                        ref_list.append(sequence.id)
                        outfile.write(sequence.format("fasta"))

                except:










    def _simple_read():
