"""
@package    FastqQualityControl
@brief
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

#~~~~~~~GLOBAL IMPORTS~~~~~~~#
# Standard library packages import
import gzip
from time import time

# Third party package import
from Bio import SeqIO
from Bio import pairwise2

# Local Package import
from Utilities import import_seq
from Utilities import fill_between_graph
from Utilities import file_basename

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class FastqFilter(object):
    """
    @class  FastqFilter
    @brief Require the third party package Biopython
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __repr__(self):
        return self.__str__() + self.get_report()

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def __init__ (self, quality_filter, adapter_trimmer, input_qual):
        """
        """
        # Declare and initialize object variables
        self.total = 0

        # Store parameters in object variables
        self.qual = quality_filter
        self.adapt = adapter_trimmer
        self.input_qual = input_qual

    #~~~~~~~CLASS METHODS~~~~~~~#

    def filter(self, R1_path, R2_path):
        """
        """
        # Start a time counter
        start_time = time()

        try:
            print("Uncompressing and extracting data")
            # Open fastq and initialize a generator to iterate over files
            in_R1 = gzip.open(R1_path, "r")
            in_R2 = gzip.open(R2_path, "r")
            genR1 = SeqIO.parse(in_R1, self.input_qual)
            genR2 = SeqIO.parse(in_R2, self.input_qual)

            # Open output stream to write in fastq.gz files
            out_name_R1 = file_basename(R1_path) + "_trimmed.fastq.gz"
            out_name_R2 = file_basename(R2_path) + "_trimmed.fastq.gz"
            out_R1 = gzip.open(out_name_R1, "w")
            out_R2 = gzip.open(out_name_R2, "w")

            print("Parsing files and filtering sequences\n")
            print("Number of sequence analyzed :\n")
            # Parsing files and filtering sequences matching quality requirements
            while True:

                # Import a new pair of fastq until end of file
                seqR1 = next(genR1, None)
                seqR2 = next(genR2, None)
                if not seqR1 or not seqR2:
                    break
                self.total +=1
                if self.total%1000 == 0:
                    print (self.total)

                # Quality filtering
                if self.qual:
                    seqR1 = self.qual.filter(seqR1)
                    seqR2 = self.qual.filter(seqR2)
                    if not seqR1 or not seqR2:
                        continue

                # Adapter trimming and size filtering
                if self.adapt:
                    seqR1 = self.adapt.trimmer(seqR1)
                    seqR2 = self.adapt.trimmer(seqR2)
                    if not seqR1 or not seqR2:
                        continue

                # If all test passed = write fastq to the files
                out_R1.write(seqR1.format("fastq-sanger"))
                out_R2.write(seqR2.format("fastq-sanger"))

            in_R1.close()
            in_R2.close()
            out_R1.close()
            out_R2.close()

        except IOError as E:
            print (E)
            exit (0)

        return (self.get_report(time()-start_time, R1_path, R2_path))


    def get_report (self, exec_time=None, R1_path=None, R2_path=None):
        """
        Return a report
        """
        report = "====== FASTQFILTER QUALITY CONTROL ======\n\n"
        if exec_time:
            report += "  Execution time : {} s\n".format(exec_time)
            report += "  R1 file : {}\n".format (R1_path)
            report += "  R2 file : {}\n".format (R2_path)
        report += "  Total sequences processed : {}\n".format(self.total)
        report += "  Input quality score : {}\n\n".format (self.input_qual)
        if self.qual:
            report += self.qual.get_report()
        if self.adapt:
            report += self.adapt.get_report()

        return report

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class QualityFilter(object):
    """
    @class  QualityFilter
    @brief
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __repr__(self):
        return self.__str__() + self.get_report()

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def __init__(self, min_qual):
        """
        """
        # Init object variables
        self.min_qual = min_qual
        self.qual_pass = 0
        self.qual_fail = 0

        # Count each of the mean quality in a dictionary
        self.qdict = { i : 0 for i in range (41)}

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def filter(self, record):
        """
        Compute mean quality score and compare to the minimal quality required
        """
        # Compute the mean quality
        mean = sum(record.letter_annotations['phred_quality'])/len(record)
        # Fill the quality dictionary
        self.qdict[mean] += 1

        # Return the record if its quality is high enought
        if mean >= self.min_qual:
            self.qual_pass += 1
            return record
        else:
            self.qual_fail += 1
            return None

    def get_report (self):
        """
        Return a report
        """
        report = "====== QUALITY FILTER ======\n\n"
        report += "  Quality Threshold : {}\n".format(self.min_qual)
        report += "  Fail quality filter : {}\n".format(self.qual_fail)
        report += "  Pass quality filter : {}\n".format(self.qual_pass)

        if self.qual_fail+self.qual_pass != 0:
            # Simplify the dict by removing empty quality means
            reduce_qdict = {key: value for key, value in self.qdict.items() if value !=0}
            # Calculate simple probability
            mean_qual = sum([key*value for key, value in reduce_qdict.items()])/sum(reduce_qdict.values())
            min_qual = reduce_qdict.keys()[0]
            max_qual = reduce_qdict.keys()[-1]

            report += "  Mean quality : {}\n".format(mean_qual)
            report += "  Minimal quality : {}\n".format(min_qual)
            report += "  Maximal quality : {}\n".format(max_qual)

        report += "\n"
        return report

    def trace_graph (self):
        """

        """
        print ("\tCreating a graphical output of mean quality distribution")

        X = self.qdict.keys()
        Y = self.qdict.values()
        basename = "Mean_read_quality"
        img_type = "svg"
        title = "Distribution of mean read quality"
        xlabel = "Mean quality"
        ylabel = "Number of reads"

        fill_between_graph(X, Y, basename, img_type, title, xlabel, ylabel)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class AdapterTrimmer(object):
    """
    @class  AdapterTrimmer
    @brief
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __repr__(self):
        return self.__str__() + self.get_report()

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def __init__ (self, aligner, adapter_path, min_size, read_len):
        """
        """
        # Init object variables
        self.min_size = min_size
        self.read_len = read_len
        self.aligner = aligner

        # Import a list of adapters and add the reverse complements of adapters to the list
        self.adapter_list = import_seq(adapter_path, "list", "fasta")
        assert self.adapter_list, "The adater list is empty"

        adapter_rc = []
        for i in self.adapter_list:
            rc = i.reverse_complement()
            rc.id = i.id + "_RC"
            adapter_rc.append(rc)
        self.adapter_list.extend(adapter_rc)

        # Initialize description field to 0 to serve as a counter of match found
        for i in self.adapter_list:
            i.description = 0

        # Initialize generic counters
        self.seq_untrimmed = 0
        self.seq_trimmed = 0
        self.base_trimmed = 0
        self.len_pass = 0
        self.len_fail = 0

        # Dict to accumulate remaining coverage after trimming
        self.coverage = { i : 0 for i in range (self.read_len)}

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def trimmer (self, record):
        """
        """
        #stime = time()
        match_list = []

        for adapter in self.adapter_list:
            # Find matches of the adapter along the current read
            adapter_match = self.aligner.find_match(str(adapter.seq).lower(), str(record.seq).lower())
            # Update the counter of the adapter with the number of matches found
            adapter.description += (len(adapter_match))
            # Add to the list of
            match_list.extend(adapter_match)

        #print ("FIND MATCH = {} ms".format((time()-stime)/1000))
        #stime = time()

        # In case no match were found, the sequence doesn't need to be modify
        if not match_list:
            self.seq_untrimmed += 1
            return record

        # Else find the longer interval without adaptor matches
        start, end = self._longer_interval (match_list, len(record))
        assert 0 <= start < end <= len(record), "Invalid interval returned"

        #print ("LONGER INTERVAL = {} ms".format((time()-stime)/1000))
        #stime = time()

        # Update counters
        self.seq_trimmed += 1
        self.base_trimmed += (len(record)-(end-start))
        self._update_coverage (start, end)

        #print ("COUNTUP = {} ms".format((time()-stime)/1000))

        # Return a slice of the reccord corresponding to the longer interval
        if end-start >= self.min_size:
            self.len_pass +=1
            return record[start:end]
        # Or None if smaller than min_size
        else:
            self.len_fail +=1
            return None

    def get_report (self):
        """
        """
        report = "====== ADAPTER TRIMMER ======\n\n"
        report += "  List of adapters imported for trimming\n"
        for i in self.adapter_list:
            report += "    Name : {}\tSequence : {}\tTimes trimmed : {}\n".format(
                i.id, i.seq, i.description)
        report += "\n  Sequences untrimmed : {}\n".format(self.seq_untrimmed)
        report += "  Sequences trimmed : {}\n".format(self.seq_trimmed)
        report += "  DNA base trimmed : {}\n".format(self.base_trimmed)
        report += "  Fail len filtering: {}\n".format(self.len_fail)
        report += "  Pass len filtering : {}\n".format(self.len_pass)
        report += "  Total pass : {}\n\n".format(self.len_pass+self.seq_untrimmed)
        report += self.aligner.get_report()
        return report

    def trace_graph (self):
        """

        """
        print ("\tCreating a graphical output of read coverage after trimming")

        X = self.coverage.keys()
        Y = self.coverage.values()
        basename = "Read_coverage_trimming"
        img_type = "svg"
        title = "Coverage of reads after trimming"
        xlabel = "Position"
        ylabel = "Number of overlapping reads"

        fill_between_graph(X, Y, basename, img_type, title, xlabel, ylabel)

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _longer_interval(self, match_list, maxsize):
        """
        Find the first larger interval that do not overlapp any match in match list
        This strategy allow to use an unsorted list of match
        """
        # Initialize a list of boolean to False of the same size as the read
        coverage = [False for i in range(maxsize)]

        # Flag positions overlapped by a read by changing the boolean to True
        for begin, end in match_list:
            for i in range (begin, end):
                coverage[i] = True

        # Read through the list to find the longer inteval between True flags
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

    def _update_coverage (self, start, end):
        """
        """
        for i in range(start, end):
            self.coverage[i] += 1


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class PairwiseAligner(object):
    """
    @class  PairwiseAligner
    @brief  Require the third party package Biopython
    """
    # TODO write a wrapper for an external faster SW aligner
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __repr__(self):
        return self.__str__() + self.get_report()

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def __init__ (self, max_alignment=5, match=2, mismatch=-2, open_gap=-2, extend_gap=-2, cutoff=1.5):
        """
        @param max_alignment Maximal number of alignement per read to output
        @param match Gain value in case match
        @param mismatch Penality value in case mismatch
        @param open_gap Penality value in case opening a gap
        @param extend_gap Penality value in case extending a gap
        @param cutoff Raw SW score divided by the lenght of the adapter
        """

        # Store parameters in object variables
        pairwise2.MAX_ALIGNMENTS = max_alignment
        self.match = match
        self.mismatch = mismatch
        self.open_gap = open_gap
        self.extend_gap = extend_gap
        self.cutoff = cutoff

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def find_match (self, adapter, sequence):
        """
        Pairwise alignment using the buildin Biopython function pairwise2
        """

        # Perform a local SW alignment with specific penalties
        match_list = pairwise2.align.localms(
            adapter,
            sequence,
            self.match,
            self.mismatch,
            self.open_gap,
            self.extend_gap)

        #for seqA, seqB, score, begin, end in match_list:
            #if score/len(adapter) > self.cutoff:
                #print("{}\n{}\nScore:{}\tBegin:{}\tEnd:{}".format(seqA, seqB, score, begin, end))

        # Return begin and end position if a match has a score higher than the cutoff value
        return [[begin, end] for seqA, seqB, score, begin, end in match_list if score/len(adapter) > self.cutoff]

    def get_report (self):
        """
        """
        report = "====== PAIRWISE ALIGNER ======\n\n"
        report += "  Parameters of pairwise2\n"
        report += "  Maximal number of alignement : {}\n".format(pairwise2.MAX_ALIGNMENTS)
        report += "  Match bonus : {}\n".format(self.match)
        report += "  Mismatch Penality : {}\n".format(self.mismatch)
        report += "  Gap open Penality : {}\n".format(self.open_gap)
        report += "  Gap extend Penality : {}\n".format(self.extend_gap)
        report += "  Score cutoff : {}\n".format(self.cutoff)

        return report