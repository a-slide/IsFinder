# Standard library packages
from subprocess import Popen, PIPE
from multiprocessing import cpu_count

####################################################################################################

class Blast(object):
    """Singleton class to manage blast db creation and blastn
    """
####################################################################################################
    
    ####    PUBLIC METHODS     ####
    
    @ classmethod
    def blast (self, query, subject, makedb=True, evalue=10):
        """ 
        """
        # If needed a blast database will be generated
        if makedb:
            cmd = "makeblastdb -in {0} -dbtype nucl -input_type fasta".format(subject)
            print ("Creating subject database") 
            print(self._cmd_return(cmd))
    
        cmd = "blastn -task blastn -outfmt 6 -dust no -num_threads {0} -evalue {1} -query {2} -db {3}".format(
            cpu_count(), evalue, query, subject)
    
        print ("Blast query against subject database")
        for h in self._parse_blast(cmd):
            hit = BlastHit(h[0], h[1] , h[2], h[3], h[4], h[5], h[6], h[7], h[8], h[9], h[10], h[11])
            print (hit)
        
        BlastHit.count_total()
        BlastHit.count_per_ref()
        BlastHit.get()
        BlastHit.get_ref()
        
    ####    PRIVATE METHODS     ####
    
    @ classmethod
    def _parse_blast (self, cmd):
        return [i.split() for i in self._split_lines (cmd)]
    
    @ classmethod
    def _split_lines (self, cmd):
        return [i for i in self._cmd_return(cmd).splitlines()]
    
    @ classmethod
    def _cmd_return(self, command):
        try:
            return (Popen(command, shell=True, stdout=PIPE).stdout.read())
        except (OSError, ValueError) as E:
            print (E)
            exit (0)
            
####################################################################################################

class BlastHit(object):
    """Class description
    """
####################################################################################################
    
    ####    CLASS FIELD    ####
    
    Instances = [] # Class field used for instance tracking 
    
    ####    CLASS METHODS    ####
    
    @ classmethod
    def count_total (self):
        return (len(self.Instances))
    
    @ classmethod
    def count_per_ref (self):
        d = {}
        for hit in Instances:
            if hit.s_id in d:
                d[ref] +=1
            else:
                d[ref] = 0
        return d
        
    @ classmethod
    def get (self):
        return Instances
        
    @ classmethod
    def get_ref (self, ref):
        return [hit for hit in Instances if hit.s_id == "ref"]
        
    ####    FONDAMENTAL METHODS     ####
   
    def __init__(self, q_id, s_id, identity, length, mis, gap, q_start, q_end, s_start, s_end, evalue, bscore):
        """Object constructor"""
        
        # Name of the reference seq in the query and in the subject were a hit was found
        self.q_id = q_id
        self.s_id = s_id
        # Percentage of identity
        self.identity = float(identity)
        # Length of the alignement
        self.length = int(length)
        # Number of mismatches and gaps
        self.mis = int(mis)
        self.gap = int(gap)
        # Start and end position along the query sequence
        self.q_orient = int(q_start) < int(q_end) # True if orientation is positive
        self.q_start = int(q_start) if self.q_orient else int(q_end)
        self.q_end = int(q_end) if self.q_orient else int(q_start)
        # Start and end position along the subject sequence
        self.s_orient = int(s_start) < int(s_end) # True if orientation is positive
        self.s_start = int(s_start) if self.s_orient else int(s_end)
        self.s_end = int(s_end) if self.s_orient else int(s_start)
        # Evalue and bit score of the match
        self.evalue = float(evalue)
        self.bscore = bscore
        # Add the instance to the class instance tracking list
        self.Instances.append(self)
    
    def __repr__(self):
        """Long representation"""
        msg = "{}\n".format(self.__str__())
        msg += "\tQuery\t{}:{}-{}({})\n".format(self.q_id, self.q_start, self.s_end, self.q_orient)
        msg += "\tSubject\t{}:{}-{}({})\n".format(self.q_id, self.q_start, self.s_end, self.q_orient)
        msg += "\tLenght : {}\tIdentity : {}%\tEvalue : {}\tBit score : {}\t".format(self.length, self.identity, self.evalue, self.bscore)
        
    def __str__(self):
        """Short representation"""
        return "<Instance of " + self.__module__ + ">"

####################################################################################################

class BlastException(Exception):
    """Custom Exception Class to handle error during blast
    """
####################################################################################################

    def __init__(self, msg): # Object constructor initialized with a custom user message
        self.msg = msg

    def __str__(self): # String returned by print
        err_msg = "BlastHitsException : An error occured while trying to blast ref1 against ref2\n"
        return err_msg
