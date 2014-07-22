"""
@package    Utilities
@brief      Contains several usefull functions to interact with OS environement and to parse files
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

#~~~~~~~COMMAND LINE UTILITIES~~~~~~~#

def run_command(cmd, stdin=None, ret_stderr=False, ret_stdout=True):
    """
    Run a command line in the default shell and return the standard output
    @param  cmd A command line string formated as a string
    @param  stdinput    Facultative parameters to redirect an object to the standard input
    @param  ret_stderr  If True the standard error output will be returned
    @param  ret_stdout  If True the standard output will be returned
    @note If ret_stderr and ret_stdout are True a tuple will be returned and if both are False
    None will be returned
    @return If no standard error return the standard output as a string
    @exception  OSError Raise if a message is return on the standard error output
    @exception  (ValueError,OSError) May be raise by Popen
    """
    # Function specific imports
    from subprocess import Popen, PIPE

    # Execute the command line in the default shell
    if stdin:
        proc = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        stdout, stderr = proc.communicate(input=stdin)
    else:
        proc = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        stdout, stderr = proc.communicate()

    if proc.returncode == 1:
        msg = "An error occured during execution of following command :\n"
        msg += "COMMAND : {}\n".format(cmd)
        msg += "STDERR : {}\n".format(stderr)
        raise Exception (msg)

    # Else return data according to user choices is returned
    if ret_stdout and ret_stderr:
        return stdout, stderr
    elif ret_stdout:
        return stdout
    elif ret_stderr:
        return stderr
    else:
        return None

def make_cmd_str(prog_name, opt_dict={}, opt_list=[]):
    """
    Create a Unix like command line string from a
    @param prog_name Name (if added to the system path) or path of the programm
    @param opt_dict Dictionnary of option arguments such as "-t 5". The option flag have to
    be the key (without "-") and the the option value in the dictionnary value. If no value is
    requested after the option flag "None" had to be asigned to the value field.
    @param opt_list List of simple command line arguments
    @exemple make_cmd_str("bwa", {"b":None, t":6, "i":"../idx/seq.fa"}, ["../read1", "../read2"])
    """

    # Start the string by the name of the program
    cmd = "{} ".format(prog_name)

    # Add options arguments from opt_dict
    if opt_dict:
        for key, value in opt_dict.items():
            if value:
                cmd += "-{} {} ".format(key, value)
            else:
                cmd += "-{} ".format(key)

    # Add arguments from opt_list
    if opt_list:
        for value in opt_list:
            cmd += "{} ".format(value)

    return cmd

#~~~~~~~FILE MANIPULATION~~~~~~~#

def fgzip(in_name, out_name=None):
    """
    @param in_name Path of the input uncompressed file
    @param out_name Path of the output compressed file (facultative)
    @exception  OSError Can be raise by open
    """
    # Function specific imports
    import gzip
    from os import remove, path

    # Generate a automatic name if none is given
    if not out_name:
        out_name = in_name +".gz"

    # Try to initialize handle for
    try:
        in_handle = open(in_name, "r")
        out_handle = gzip.open(out_name, "w")
        # Write input file in output file
        print ("Compressing {}".format(in_name))
        out_handle.write (in_handle.read())
        # Close both files
        in_handle.close()
        out_handle.close()
        return path.abspath(out_name)

    except IOError as E:
        print(E)
        if path.isfile (out_name):
            try:
                remove (out_name)
            except OSError:
                print "Can't remove {}".format(out_name)

def fgunzip(in_name, out_name=None):
    """
    @param in_name Path of the input compressed file
    @param out_name Path of the output uncompressed file (facultative)
    @exception  OSError Can be raise by open
    """
    # Function specific imports
    import gzip
    from os import remove, path

    # Generate a automatic name without .gz extension if none is given
    if not out_name:
        out_name = in_name[0:-3]

    try:
        # Try to initialize handle for
        in_handle = gzip.open(in_name, "r")
        out_handle = open(out_name, "w")
        # Write input file in output file
        print ("Uncompressing {}".format(in_name))
        out_handle.write (in_handle.read())
        # Close both files
        out_handle.close()
        in_handle.close()
        return path.abspath(out_name)

    except IOError as E:
        print(E)
        if path.isfile (out_name):
            try:
                remove (out_name)
            except OSError:
                print "Can't remove {}".format(out_name)


def mkdir(fp):
    """
    Create a directory at the indicated path\n
    Reproduce the ability of UNIX "mkdir -p" command
    (ie if the path already exits no exception will be raised).
    @param  fp path name where the folder should be created
    @exception  OSError Can be raise by os.mkdir
    """
    # Function specific imports
    from os import mkdir, path

    if path.exists(fp) and path.isdir(fp):
        print ("'{}' already exist in the current directory".format(fp))
    else:
        print ("Creating '{}' in the current directory".format(fp))
        mkdir(fp)

def merge_files (input_list, output=None, compress_output=True):
    """
    @param input_list List of files to merge
    @param output Destination file
    """
    # Standard library import
    import gzip
    from os.path import getsize, isfile, abspath
    from tempfile import mkstemp

    # Creating and storing a file for writting output
    if not output:
        if compress_output:
            OSlevel, output = mkstemp(suffix=".gz")
        else:
            OSlevel, output = mkstemp()
    else:
        output = abspath(output)

    try:
        # Opening the outpout file for writting
        if compress_output:
            outfile = gzip.open(output, "wb")
        else:
            outfile = open(output, "wb")

        # For all files in the input list
        for input in input_list:

            # In case the file is gzipped
            if file_extension(input) in ["gz","GZ"]:
                infile = gzip.open(input, "rb")
            # In case the file is not compressed
            else:
                infile = open(input, "rb")

            # writting infile in outfile
            outfile.write(infile.read())
            infile.close()

        # Close and verify the output file
        outfile.close()
        assert isfile(output), "No output file found"
        assert getsize(output) > 0, "The output file is empty"

    except (IOError,AssertionError)  as E:
        print(E)
        exit()

    return output


def file_basename (path):
    """
    Return the basename of a file without folder location and extension
    @param path Filepath as a string
    """
    return path.rpartition('/')[2].partition('.')[0]

def file_extension (path):
    """
    Return the extension of a file.
    @param path Filepath as a string
    """
    return path.rpartition(".")[2]

def file_name (path):
    """
    Return the complete name of a file with the extension but without folder location
    @param path Filepath as a string
    """
    return path.rpartition("/")[2]

def dir_name (path):
    """
    Return the complete path where is located the file without the file name
    @param path Filepath as a string
    """
    return path.rpartition("/")[0].rpartition("/")[2]

#~~~~~~~FASTA UTILITIES~~~~~~~#

def import_seq(filename, col_type="dict", seq_type="fasta"):
    """
    Import sequences from a fasta files in a list of biopython SeqRecord
    @param filename Valid path to a fasta file. Can contains several sequences and can be gzipped
    @param col_type Type of the collection where SeqReccord entries will be added ("list" or "dict").
    @param seq_type Type of the sequence file to parse (see Biopython seqIO for supported format)
    @return A list or a dictionnary containing all seqReccord objects from the fastq file
    @exception IOError  Raise if the path in invalid or unreadeable
    """
    # Require the Third party package Biopython
    from Bio import SeqIO
    import gzip

    # Try to open the file first gz compressed and uncompressed
    try:

        # Verify if the type of the input sequence is valid
        seq_type = seq_type.lower()
        allowed_seq = ["fasta", "genbank", "gb", "fastq-illumina", "fastq-solexa" , "fastq",
        "fastq-sanger", "embl ", "abi ", "seqxml", "sff", "uniprot-xml"]
        assert seq_type in allowed_seq , "The input file format have to be in the following list : "+ ", ".join(allowed_seq)

        # Verify if the type of the output collection is valid
        col_type = col_type.lower()
        allowed_types = ["dict", "list"]
        assert col_type  in allowed_types, "The output collection type have to be in the following list : "+ ", ".join(allowed_types)

        # Open gzipped or uncompressed file
        if file_extension(filename) in ["gz","GZ"]:
            #print("\tUncompressing and extracting data")
            handle = gzip.open(filename, "r")
        else:
            #print("\tExtracting data")
            handle = open(filename, "r")

        # Create the collection
        if col_type == "list":
            seq_col = list(SeqIO.parse(handle, seq_type))
        else:
            seq_col = SeqIO.to_dict(SeqIO.parse(handle, seq_type))

        # Close file, verify if the collection is filled and returned it
        handle.close()
        assert seq_col, 'The collection contains no SeqRecord after file parsing. Exit'
        return seq_col

    except IOError as E:
        print(E)
        exit()

    except AssertionError as E:
        print (E)
        exit()

#~~~~~~~GRAPHICAL UTILIES~~~~~~~#

def fill_between_graph (X, Y, basename="out", img_type="png", title=None, xlabel=None, ylabel=None, baseline=0):
    """
    Trace a generic fill between graph with matplotlib pyplot
    @param X List of values for x axis
    @param Y List of values for y axis
    @param title Title of graph (facultative)
    @param xlabel Label for x axis (facultative)
    @param ylabel Label for y axis (facultative)
    @param basename Output basename of the image file (Default "out")
    @param img_type Type of the image file (Default "png")
    """

    # Require the Third party package Biopython
    from matplotlib import pyplot as plt

    # Create a figure object and adding details
    fig = plt.figure(figsize=(15, 10), dpi=100)
    if title:
        plt.title(title)
    if xlabel:
        plt.xlabel(xlabel)
    if ylabel:
        plt.ylabel(ylabel)

    # Plot an area representing the coverage depth
    plt.fill_between(X, Y, baseline, facecolor='green', alpha=0.5)

    # Tweak spacing to prevent clipping of ylabel
    plt.subplots_adjust(left=0.15)

    # Export figure to file
    try:
        fig.savefig(basename+"."+img_type, format = img_type)
    except ValueError as E:
        print (E)
        print ("Saving file as png")
        fig.savefig(basename+".png", format = "png")
