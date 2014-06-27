#~~~~~~~GLOBAL IMPORTS~~~~~~~#

# Local library packages
from Utilities import mkdir, import_seq, file_basename, file_name

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
    def masker (self, query_list, subject, evalue, subjectdb=None):
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
        ref_dict = import_seq(subject, "dict", "fasta")

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