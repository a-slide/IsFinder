#~~~~~~~GLOBAL IMPORTS~~~~~~~#
# Standard library packages import
from multiprocessing import cpu_count
from sys import exit as sys_exit
from subprocess import Popen, PIPE
from os import symlink, remove, path

# Third party package import
import pysam

# Local Package import
from Utilities import mkdir, file_name, make_cmd_str

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class BwaWrapper(object):
    """
    @class BwaWrapper
    @brief Wrapper for bwa index and bwa mem aligner. Initialise an object with a reference index
    and can be call through align to
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#


    def __repr__(self):
        return "{}\nSTANDARD OUTPUT\n{}\nSTANDARD ERROR\n{}".format(
            self.__str__(),
            self.index_report,
            self.align_report)

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def __init__ (self, bwa_path, ref_path, make_index=True, align_opt={}, index_opt={}):
        """
        Initialize the object and index the reference genome if necessary
        """
        # Store path to aligner and indexer programs
        self.indexer = bwa_path+" index"
        self.aligner = bwa_path+" mem"
        self.index_opt = index_opt
        self.align_opt = align_opt
        self.index_report = ""
        self.align_report = ""

        # Create an index for the reference genome
        if make_index:
            # Complete the index option dict with default values
            if "a" not in self.index_opt:
                self.index_opt["a"] = "bwtsw"
            if "p" not in self.index_opt:
                try:
                    out_dir = "./bwa_index/"
                    mkdir(out_dir)
                except Exception:
                    out_dir = "./"
                finally:
                    self.index_opt["p"] = out_dir + file_name(ref_path)

            # Make index with bwa index and store the path in an object variable
            self._make_index(ref_path)
            self.index_path = self.index_opt["p"]

        # Else just attribute the user index to the object variable
        else:
            self.index_path = ref_path

        # Complete the aligner option dict with default values
        if "t" not in self.align_opt:
            self.align_opt["t"] = cpu_count()

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def align(self, R1_path, R2_path, out_prefix="out"):
        """
        """
        # Build the command line thanks to make_cmd_str
        #out_name = "> {}.sam".format(out_prefix)
        opt_list = [self.index_path, R1_path, R2_path]#, out_name]
        cmd = make_cmd_str(self.aligner, self.align_opt, opt_list)

        print "Aligning reads against {} with {}".format(file_name(self.index_path), self.aligner)

        # Try to execute the command and verify if stdout is not None
        try:
            proc = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
            pysam.Samfile("-", "rb")
            # Launch the process
            stdout, stderr = proc.communicate()
            # Verify the return code
            if proc.returncode == 1:
                raise Exception ("An error occured during alignment with bwa mem" + stderr)

            ###### TO TRANSFORM INTO A GENERATOR THAT CAN YIELD READS..... proc.stdout.next()

            print ("[bwa mem] Finished alignment\n")
            print (stdout)
            self.align_report = stderr

        except Exception as E:
            print (E)
            sys_exit()

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _make_index(self, ref_path):
        """
        """
        # Build the command line thanks to make_cmd_str
        opt_list = [ref_path]
        cmd = make_cmd_str(self.indexer, self.index_opt, opt_list)
        print cmd
        print "Indexing {}.".format(file_name(ref_path))

        # Try to execute the command and verify if stdout is not None
        try:
            proc = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
            # Launch the process
            stdout, stderr = proc.communicate()
            # Verify the return code
            if proc.returncode == 1:
                raise Exception ("An error occured during indexing with bwa index" + stderr)

            print (stdout)
            self.index_report = stderr

        except Exception  as E:
            for ext in ["amb", "ann", "bwt", "pac", "sa"]:
                f = "{}.{}".format(self.index_opt["p"], ext)
                print f
                try:
                    if path.isfile (f):
                        remove (f)
                except Exception:
                    print "Unable to remove file : {}".format(f)

            print (E)
            sys_exit()
