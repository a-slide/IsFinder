"""
@package    HomologyMasker
@brief      Compared a list of query DNA sequence with a subject and mask the eventual homologies in the subject sequence with Ns
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger <adrien.leger@gmail.com>
"""

# Standard library packages
from multiprocessing import cpu_count

# Local library packages
from Utilities import run_shell, mkdir, import_fasta, file_basename, file_name

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

class RefMasker(object):
    """
    @class  RefMasker
    @brief  Singleton class allowing to find DNA sequences homologies between a list of query
    and a subject and to write a modified version of the subject sequence.
    First a blast database is created from a reference fasta file if it is not provided by
    the user. Then, query sequences fasta files are blasted against the newly created subject database.
    Finally, if blast hits are found (homologies) the subject genome is imported, hit locations are
    masked with "Ns" and a new masked reference is written
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~PUBLIC CLASS METHODS~~~~~~~#

    @ classmethod
    def masker (self, query_list, subject, evalue, subjectdb=''):
        """
        Main function of RefMasker that integrate database creation, blast and homology masking
        @param  query_list List of paths indicating fasta files containing query sequences. Fasta can contains multiple sequences
        @param  subject Path to a fasta file containing the reference subject sequences
        @param  evalue  Cutoff used in blast to select valid hits
        @param  subjectdb   Facultative paramater. Path of the blastn database of the subject's basename created by "makeblastdb"
        @return If the sequence was edited the path of the edited reference is indicated, else False
        """
        if not subjectdb:
            # Create a directory to store the database
            mkdir("blastdb")
            subjectdb = "blastdb/" + file_basename(subject)
            # Create the database
            print (Blast.makedb (subject, subjectdb))

        # Generate a list of hit containing hits of all sequence in query list in subject
        hit_list = self._list_homologies (query_list, subjectdb, evalue)

        # If needed the subject sequence will be imported masked and rewritten on disk
        if not hit_list:
            print ("No hits found. The original subject sequences will not be edited")
            return 0
        else:
            mkdir("references")
            masked_subject = "references/masked_" + file_name(subject)
            self._mask_homologies(hit_list, subject, masked_subject)
            return masked_subject

     #~~~~~~~PRIVATE CLASS METHODS~~~~~~~#

    @ classmethod
    def _list_homologies (self, query_list, subjectdb, evalue):
        """
        Perform iterative blasts of query sequences against the subject database and create a list of hits
        @param  query_list List of paths indicating fasta files containing query sequences. Fasta can contains multiple sequences
        @param  subjectdb Path of the blastn database of the subject's basename created by "makeblastdb"
        @param  evalue  Cutoff used in blast to select valid hits
        @return A list of BlastHit objects
        """
        # Empty list to store hits
        hit_list = []

        # Append the list of hits for each query in a bigger list.
        for query in query_list:
            hit_list += Blast.do_blast(query, subjectdb, evalue)

        return hit_list

    @ classmethod
    def _mask_homologies (self, hit_list, subject, output):
        """
        Import the reference subject genome, edit it and rewrite an edited version
        @param  hit_list List BlastHit objects descibing the position of hits in the subject sequences
        @param  subject Path to a fasta file containing the reference subject sequences
        @param  output Path name of the output subject fasta file
        @return A list of BlastHit objects
        """
        # Importing sequences from a fasta file in a dictionnary of SeqRecord
        print ("\nImporting subject sequences for hard masking of homologies")
        ref_dict = import_fasta(subject, "dict")

        # Casting Seq type to MutableSeq Type to allow string editing
        for record in ref_dict.values():
            record.seq = record.seq.tomutable()

        print ("\nEditing subject sequences sequence")
        # For each hit in the list of blast found in the reference sequences
        for hit in hit_list:

            # For all position between start and end coordinates modify the base by N
            for position in range (hit.s_start, hit.s_end):
                ref_dict[hit.s_id].seq[position]= 'N'

        # Write the new reference in fasta format
        print ("\nWriting new version of {} in which homologies with the query list are masked".format(file_basename (subject)))
        with open(output, 'w') as f:
            for ref in ref_dict.values():
                f.write(ref.format("fasta"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

class Blast(object):
    """
    @class  Blast
    @brief  Singleton class allowing to create a blast database and to perform de blast of a query
    against a blast database. Blast+ 2.8+ needs to be install and correctly added to the path
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~PUBLIC CLASS METHODS~~~~~~~#

    @ classmethod
    def makedb (self, subject, output):
        """
        Create a blastn database in subjectdb using makeblastdb
        @param  subject Path to a fasta file containing the reference subject sequences
        @param  output Path to the output database basename
        @return Raw output of makeblastdb
        @exception (SystemError,OSerror) May be returned by run_shell in case of invalid command line
        """
        print ("\nCreating {} database".format(file_basename (subject)))

        # Build the command line string
        cmd = "makeblastdb -in {0} -out {1} -dbtype nucl -input_type fasta".format(subject, output)
        # Create the database
        return(run_shell(cmd))

    @ classmethod
    def do_blast (self, query, subjectdb, evalue):
        """
        Blast query against a subject database and return a list of BlastHit object
        @param  query Path to a fasta file containing the query sequences
        @param  subjectdb Path to the subject blast database basename
        @param  evalue  Cutoff used in blast to select valid hits
        @return A list of BlastHit objects if at least one hit was found
        @exception (SystemError,OSerror) May be returned by run_shell in case of invalid command line
        """
        # Build the command line string
        cmd = "blastn -task blastn -outfmt 6 -dust no -num_threads {0} -evalue {1} -query {2} -db {3}".format(
            cpu_count(), evalue, query, subjectdb)

        # Execute blastn (Can raise a SystemError) and create BlastHit objects
        print ("\nBlast {} against {} database".format(file_basename(query), file_basename (subjectdb)))
        try:
            for h in self._parse_blast(cmd):
                BlastHit(h[0], h[1] , h[2], h[3], h[4], h[5], h[6], h[7], h[8], h[9], h[10], h[11])

            # Sumarize the hit count in the different references
            print ("\t{} hits found in the subject database".format(BlastHit.count_total()))
            for ref, count in BlastHit.count_per_ref().items():
                print ("\t* {} {} in ref {}".format(count, "hit" if count == 1 else "hits", ref))

            hits_list = BlastHit.get()
            BlastHit.reset_list()
            return hits_list

        except AttributeError:
            print ("\t no hit found in the subject database")
            return []

    #~~~~~~~PRIVATE CLASS METHODS~~~~~~~#

    @ classmethod
    def _parse_blast (self, cmd):
        return [i.split() for i in self._split_lines (cmd)]

    @ classmethod
    def _split_lines (self, cmd):
        return [i for i in run_shell(cmd).splitlines()]

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
    * b_score : Bit score of the alignement
    A class list is used to track all instances generated.
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~CLASS FIELD~~~~~~~#

    Instances = [] # Class field used for instance tracking

    #~~~~~~~CLASS FIELD~~~~~~~#

    @ classmethod
    def count_total (self):
        """
        @return Overall number of BlastHit object in Instance list
        """
        return (len(self.Instances))

    @ classmethod
    def count_per_ref (self):
        """
        @return Number of BlastHit object in Instance list sorted by reference subject sequence
        """
        d = {}
        for hit in self.Instances:
            if hit.s_id in d:
                d[hit.s_id] +=1
            else:
                d[hit.s_id] = 1
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
        return "<Instance of {} from package {} >".format(self.__class__.__name__, self.__module__)
