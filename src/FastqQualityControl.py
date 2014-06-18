"""
@package    FastqFilter
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

# Standard library packages
import gzip
from time import time

# Third party package
from Bio import pairwise2, SeqIO

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class FastqFilter(object):
    """
    @class  FastqFilter
    @brief  
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    
    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#
    
    def __repr__(self):
        pass
    
    def __str__(self):
        return "<Instance of {} from {} >".format(self.__class__.__name__, self.__module__)
    
    def __init_ (self, quality_filter, adapter_trimmer, input_qual, output_qual):
        """
        
        """
        # Declare and initialize object variables
        self.total = 0
        
        # Store instances of QualityFilter and AdapterTrimmer
        self.qual = quality_filter
        self.adapt = adapter_trimmer
        
        # Store parameters in object variables
        self.input_qual = input_qual
        self.output_qual = output_qual
    
    #~~~~~~~CLASS METHODS~~~~~~~#
    
    def filter(self)
        """
        bla bla bla
        """
        # Start a time counter
        start_time = time()
        # Declare list to store valid fastq sequences
        ListR1 = []
        ListR2 = []
        
        try:
            print("Uncompressing and extracting data")
            # Input fastq files
            input_R1 = gzip.open(R1_path, "r")
            input_R2 = gzip.open(R2_path, "r")
            
            # Initialize a generator to iterate over fastq files
            genR1 = SeqIO.parse(input_R1, qual_scale)
            genR2 = SeqIO.parse(input_R2, qual_scale)

            print("Parsing files and filtering sequences")
            # Parsing files and filtering sequences matching quality requirements
            while True:
                # Import a new pair of fastq until end of file
                seqR1 = next(genR1, None)
                seqR2 = next(genR2, None)
                if not seqR1 or not seqR2:
                    break
                self.total +=1
                
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
                
                # If all test passed = append Sequence to the list
                
                ListR1.append(seqR1)
                ListR2.append(seqR2)
            
            input_R1.close()
            input_R2.close()
        
        except IOError as E:
            print (E)
            exit (0)
        
        print (self.get_report(R1_path, R2_path, time()-start_time))
        assert ListR1 =! [] and ListR2 =! [], , "All reads were filtered out by parameters"
        return ListR1, ListR2
    
    def get_report (self, R1_path, R2_path, exec_time):
        """
        Return a report
        """
        report = "====== FASTQFILTER QUALITY CONTROL REPORT ======\n\n" 
        report =+ "Execution time : {} s".format(exec_time)
        report += "  Total sequences processed : {}\n".format(self.total)
        report += "  R1 file : {}\n".format (R1_path)
        report += "  R2 file : {}\n".format (R2_path)
        report += "  Input quality score : {}\n".format (self.input_qual)
        report += "  Output quality score : {}\n\n".format (self.input_qual)
        report = "====== QUALITY FILTER ======\n\n" 
        if self.qual:
            report += self.qual.get_report()
        else:
            "  No quality fitering\n\n"
        report = "====== ADAPTER TRIMMING ======\n\n" 
        if self.adapt:
            report += self.adapt.get_report()
        else:
            "  No adapter trimming\n\n"
            
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
        pass
    
    def __str__(self):
        return "<Instance of {} from {} >".format(self.__class__.__name__, self.__module__)

    def __init__(self, min_qual):
        """
        """
        # Init object variables
        self.min_qual = min_qual
        self.passed = 0
        self.fitered = 0
        
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
            self.passed += 1
            return record
        else:
            self.filtered += 1
            return None
        
    def get_report (self):
        """
        Return a report
        """
        # Simplify the dict by removing empy quality means
        reduce_qdict = {key: value for key, value in self.qdict.items() if value !=0}
        # Calculate simple probability
        mean_qual = sum([key*value for key, value in reduce_qdict()])/sum(reduce_qdict.values())
        min_qual = reduce_qdict.keys()[0]
        max_qual = reduce_qdict.keys()[-1]
        
        report = "  Count filtered : {}\n".format(self.fitered)
        report += "  Cout passed : {}\n".format(self.passed)
        report += "  Mean quality : {}\n".format(mean_qual)
        report += "  Minimal quality : {}\n".format(min_qual)
        report += "  Maximal quality : {}\n\n".format(max_qual)
        
        return report
        
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class AdapterTrimmer(object):
    """
    @class  AdapterTrimmer
    @brief  
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    
    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#
    
    def __repr__(self):
        pass
    
    def __str__(self):
        return "<Instance of {} from {} >".format(self.__class__.__name__, self.__module__)

    def __init__ (self, adapter_list, min_size, read_len):
        """
        """
        # Init object variables
        self.min_size = min_size
        self.adapter_list = adapter_list
        self.seq_untrimmed = 0
        self.seq_trimmed = 0
        self.base_trimmed = 0
        
        # Dict to accumulate remaining coverage after trimming
        self.coverage = { i : 0 for i in range (read_len)}
        
    #~~~~~~~PUBLIC METHODS~~~~~~~#
    
    def trim_adapter (record, ):

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
        
    ###############################################################################################
    
    def get_report (self):
        """
        Return a report
        """
        pass
        
    #~~~~~~~PRIVATE METHODS~~~~~~~#
    
    def _align (adapter, sequence):
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


    def _longer_interval(match_list, maxsize):
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


