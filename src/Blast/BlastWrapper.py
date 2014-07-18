#~~~~~~~GLOBAL IMPORTS~~~~~~~#

# Standard library packages
from multiprocessing import cpu_count

# Local library packages
from Utilities import run_command, file_basename, make_cmd_str

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Blastn(object):
    """
    @class  Blast
    @brief  Singleton class allowing to create a blast database and to perform de blast of a query
    against a blast database. Blast+ 2.8+ needs to be install and correctly added to the path
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __repr__(self):
        msg = self.__str__()
        msg += "Blastn path : {}\n".format(self.blastn)
        msg += "Options :\n"
        for i, j in self.blastn_opt.items():
            msg += "\tFlag : {}\tValue : {}\n".format(i,j)
        msg += "BlastDB : {}\n".format(repr(self.Blastdb))
        return msg

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def __init__ (self, Blastdb, blastn_opt=None, blastn="blastn"):
        """
        Initialize the object and index the reference genome if necessary
        """
        # Creating object variables
        self.blastn = blastn
        self.Blastdb = Blastdb

        # init an option dict and attribute defaut options
        self.blastn_opt = blastn_opt if blastn_opt else {}
        self.blastn_opt["db"] = self.Blastdb.db_path
        if "num_threads" not in self.blastn_opt:
            self.blastn_opt["num_threads"] = cpu_count()
        if "task" not in self.blastn_opt:
            self.blastn_opt["task"] = "blastn"
        if "outfmt" not in self.blastn_opt:
            self.blastn_opt["outfmt"] = 6
        if "dust" not in self.blastn_opt:
            self.blastn_opt["dust"] = "no"


    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def align (self, query, evalue=0.1):
        """
        Blast query against a subject database and return a list of BlastHit object
        @param  query Path to a fasta file containing the query sequences
        @param  evalue  Cutoff used in blast to select valid hits
        @return A list of BlastHit objects if at least one hit was found
        @exception (SystemError,OSerror) May be returned by run_command in case of invalid command line
        """
        # Build the command line string
        self.blastn_opt["evalue"] = evalue
        self.blastn_opt["query"] = query
        cmd = make_cmd_str(self.blastn, self.blastn_opt)

        # Execute blastn (Can raise a SystemError) and create BlastHit objects
        print ("\nBlast {} against {} database".format(file_basename(query), file_basename (self.Blastdb.db_path)))
        try:
            # Run the command line without stdin and asking only stdout
            blast_lines = run_command(cmd, None, False, True).splitlines()
            for line in blast_lines:
                # Parse each result lines and create a BlastHit object
                h = line.split()
                BlastHit(h[0], h[1] , h[2], h[3], h[4], h[5], h[6], h[7], h[8], h[9], h[10], h[11])

            # Sumarize the hit count in the different references
            print ("\t{} hits found in the subject database".format(BlastHit.count_total()))
            for ref, val in BlastHit.stat_per_ref().items():
                print ("\t* {} hit(s) in ref {}\t({} pb)".format(val[0], ref, val[1]))

            hits_list = BlastHit.get()
            BlastHit.reset_list()
            return hits_list

        except Exception as E:
            print(E)
            raise

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class BlastHit(object):
    """
    @class  BlastHit
    @brief  Object oriented class containing information of one blast hit
    The following instance field are accessible :
    * q_id : Query sequence name
    * s_id : Subject sequence name
    * identity : % of identity in the hit
    * length : length of the hit
    * mis : Number of mismatch in the hit
    * gap : Number of gap in the hit
    * q_orient : Orientation of the query along the hit
    * q_start : Hit start position of the query
    * q_end : Hit end position of the query
    * s_orient : Orientation of the subject along the hit
    * s_start : Hit start position of the subject
    * s_end : Hit end position of the subject
    * evalue : E value of the alignement
    * bscore : Bit score of the alignement
    A class list is used to track all instances generated.
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~CLASS FIELDS~~~~~~~#

    Instances = [] # Class field used for instance tracking

    #~~~~~~~CLASS METHODS~~~~~~~#

    @ classmethod
    def count_total (self):
        """
        @return Overall number of BlastHit object in Instance list
        """
        return (len(self.Instances))

    @ classmethod
    def stat_per_ref (self):
        """
        @return Number of BlastHit object in Instance list sorted by reference subject sequence
        """
        d = {}
        for hit in self.Instances:
            if hit.s_id in d:
                d[hit.s_id][0] += 1
                d[hit.s_id][1] += hit.length
            else:
                d[hit.s_id] = [1, hit.length]
        return d

    @ classmethod
    def get (self):
        """
        @return The list of all BlastHit object generated
        """
        return self.Instances

    @ classmethod
    def get_ref (self, ref):
        """
        @param ref Name of a reference sequence in the subject database
        @return The list of all BlastHit object generated for this reference
        """
        return [hit for hit in self.Instances if hit.s_id == "ref"]

    @ classmethod
    def reset_list (self):
        """
        Reset the instance tracking list (Usefull after
        """
        self.Instances = []

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__(self, q_id, s_id, identity, length, mis, gap, q_start, q_end, s_start, s_end, evalue, bscore):
        """
        Create a BlastHit object which is automatically added to the class tracking instance list
        The object with the following parameters are required for object initialisation
        @param  q_id    Query sequence name
        @param  s_id    Subject sequence name
        @param  identity    % of identity in the hit
        @param  length  length of the hit
        @param  mis Number of mismatch in the hit
        @param  gap Number of gap in the hit
        @param  q_start Hit start position of the query
        @param  q_end   Hit end position of the query
        @param  s_start Hit start position of the subject
        @param  s_end   Hit end position of the subject
        @param  evalue  E value of the alignement
        @param  bscore Bit score of the alignement
        """
        self.q_id = q_id
        self.s_id = s_id
        self.identity = float(identity)
        self.length = int(length)
        self.mis = int(mis)
        self.gap = int(gap)
        self.evalue = float(evalue)
        self.bscore = float(bscore)

        # Autoadapt start and end so that start is always smaller than end
        self.q_start = int(q_start) if int(q_start) < int(q_end) else int(q_end)
        self.q_end = int(q_end) if int(q_start) < int(q_end) else int(q_start)
        self.s_start = int(s_start) if int(s_start) < int(s_end) else int(s_end)
        self.s_end = int(s_end) if int(s_start) < int(s_end) else int(s_start)

        # Orientation of the query and subject along the hit. True if positive
        self.q_orient = int(q_start) < int(q_end)
        self.s_orient = int(s_start) < int(s_end)

        # Add the instance to the class instance tracking list
        self.Instances.append(self)

    def __repr__(self):
        msg = "{}\n".format(self.__str__())
        msg += "\tQuery\t{}:{}-{}({})\n".format(self.q_id, self.q_start, self.q_end, "+" if self.q_orient else "-")
        msg += "\tSubject\t{}:{}-{}({})\n".format(self.s_id, self.s_start, self.s_end, "+" if self.q_orient else "-")
        msg += "\tLenght : {}\tIdentity : {}%\tEvalue : {}\tBit score : {}\n".format(self.length, self.identity, self.evalue, self.bscore)
        return (msg)

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)
