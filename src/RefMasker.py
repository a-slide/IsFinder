"""
@package RefMasker
@brief Find DNA sequences homologies between a list of query and a subject and to write a modified
version of the subject sequence.
* First a blast database is created from a reference fasta file if it is not provided by the user.
* Then, query sequences fasta files are blasted against the newly created subject database.
* Finally, if blast hits are found (homologies) the subject genome is imported, hit locations are
masked with "Ns" and a new masked reference is written
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
from tempfile import mkstemp
import gzip

# Local library packages
from Utilities import mkdir, import_seq, file_basename, file_name, file_extension, fgunzip
from BlastnWrapper import Blastn, NewDB, ExistingDB

#~~~~~~~MAIN METHODS~~~~~~~#

def mask   (query_list,
            subject,
            evalue=0.1,
            db_path=None,
            blastn_opt=None,
            makeblastdb_opt=None,
            db_outdir="./blastdb/",
            maskref_outdir="./references/",
            blastn="blastn",
            makeblastdb="makeblastdb"):
    """
    Main function of RefMasker that integrate database creation, blast and homology masking
    * Instantiate Blast database and blastn object
    * Perform iterative blasts of query sequences against the subject database and create a list of
    hits.
    * If hits ar found the refercne sequence will be imported masked and rewritten on disk
    @param  query_list List of paths indicating fasta files containing query sequences. Fasta can
    contains multiple sequences
    @param  subject Path to a fasta file containing the reference subject sequences (can be gziped)
    @param  evalue  Cutoff used in blast to select valid hits
    @param  db_path   Facultative parameter. Path of the blastn database of the subject's basename
    created by "makeblastdb"
    @param blastn_opt Dict of options for blastn (see make_cmd_str from Utilities library)
    @param makeblastdb_opt Dict of options for makeblastdb (see make_cmd_str from Utilities library)
    @param db_outdir Output directory for the database is requested
    @param maskref_outdir Output directory for the masked reference is requested
    @param blastn Path ot the blastn executable.
    @param makeblastdb Path ot the makeblastdb executable.
    @return If the sequence was edited the path of the edited reference is indicated, else False
    """

    # Try to validate a database from an existing one
    try:
        if not db_path:
            raise Exception("No Blast database was provided")
        else:
            BlastDB = ExistingDB(db_path)

    # If no DB or if an error occured during validation of the existing DB = create a new db
    except Exception as E:
        print (E)
        print ("Try to create a database from the subject reference sequence")

        # TO DO remove the folder in case of error
        mkdir(db_outdir)
        db_outname = path.join (db_outdir, file_basename(subject))

        # If the fastq file is compressed = extract the file in a tempory file
        if file_extension(subject) in ["gz","GZ"]:
            tmphandle, tmpname = mkstemp(suffix='.fa')
            fgunzip (subject, tmpname)
            BlastDB = NewDB(tmpname, db_outname, makeblastdb_opt, makeblastdb)
            remove(tmpname)
        else:
            BlastDB = NewDB(subject, db_outname, makeblastdb_opt, makeblastdb)

    # Initialise a Blastn object
    Blaster = Blastn(BlastDB, blastn_opt, blastn)

    # Generate a list of hit containing hits of all sequence in query list in subject
    hit_list = []
    # Extend the list of hits for each query in a bigger list.
    for query in query_list:
        hit_list.extend(Blaster.align(query, evalue))

    # If needed the subject sequence will be imported as biopython seqRecord,
    # masked and rewritten on disk
    if not hit_list:
        print ("No hits found. The original subject sequences will not be edited")
        return subject
    else:
        mkdir(maskref_outdir)
        maskref_outname =  path.join (maskref_outdir + "masked_{}.fa.gz".format(file_basename(subject)))
        mask_homologies(hit_list, subject, maskref_outname)
        return maskref_outname

def mask_homologies (hit_list, subject, outname):
    """
    Import the reference subject genome, edit it and rewrite an edited version
    @param  hit_list List BlastHit objects descibing the position of hits in the subject sequences
    @param  subject Path to a fasta file containing the reference subject sequences
    @param  outname Path name of the output subject fasta file
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
    with gzip.open(outname, 'w') as f:
        for ref in ref_dict.values():
            f.write(ref.format("fasta"))
