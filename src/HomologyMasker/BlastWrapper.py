#~~~~~~~GLOBAL IMPORTS~~~~~~~#

# Standard library packages
from multiprocessing import cpu_count
# Local library packages
from Utilities import run_command, file_basename

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class BlastWrapper(object):
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
        @exception (SystemError,OSerror) May be returned by run_command in case of invalid command line
        """
        print ("\nCreating {} database".format(file_basename (subject)))

        # Build the command line string
        cmd = "makeblastdb -in {0} -out {1} -dbtype nucl -input_type fasta".format(subject, output)
        # Create the database
        return(run_command(cmd))

    @ classmethod
    def do_blast (self, query, subjectdb, evalue):
        """
        Blast query against a subject database and return a list of BlastHit object
        @param  query Path to a fasta file containing the query sequences
        @param  subjectdb Path to the subject blast database basename
        @param  evalue  Cutoff used in blast to select valid hits
        @return A list of BlastHit objects if at least one hit was found
        @exception (SystemError,OSerror) May be returned by run_command in case of invalid command line
        """

        print query
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
        return [i for i in run_command(cmd).splitlines()]
