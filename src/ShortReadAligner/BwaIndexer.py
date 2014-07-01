#~~~~~~~GLOBAL IMPORTS~~~~~~~#
# Standard library packages import
from multiprocessing import cpu_count
from sys import exit as sys_exit
from subprocess import Popen, PIPE
from os import remove, path, rmdir
from tempfile import NamedTemporaryFile as tempfile

# Third party package import
import pysam

# Local Package import
from Utilities import mkdir, file_name, make_cmd_str

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class GenerateIndex(object):
    """
    @class GenerateIndex
    @brief Wrapper for bwa index. Create a reference index from a fasta file
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __repr__(self):
        msg = self.__str__()
        msg += "Bwa Index Report:\n{}\n".format(self.index_report)
        msg += "Bwa Index Path:\n{}\n".format(self.index_path)
        return msg

    def __str__(self):
        return "\n<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def __init__ (self, bwa_path, ref_path, index_opt=None):
        """
        Initialize the object and index the reference genome if necessary
        """
        # Creating object variables
        self.indexer = bwa_path+" index"
        self.ref_path = ref_path
        self.index_opt = index_opt if index_opt else {}
        self.index_report = ""

        # By default bwtsw is selected if no algorithm is indicated in index_opt
        if "a" not in self.index_opt:
            self.index_opt["a"] = "bwtsw"
        
        # By default a local directory named bwa_index is created to store index files if 
        # no output directory is already indicated in index_opt
        if "p" not in self.index_opt:
            try:
                self.out_dir = "./bwa_index/"
                mkdir(self.out_dir)
            except Exception:
                self.out_dir = "./"
            finally:
                self.index_opt["p"] = self.out_dir + file_name(self.ref_path)

        # Make index with bwa index and store the path in an object variable
        self._make_index()
        self.index_path = self.index_opt["p"]

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _make_index(self):
        """
        """
        # Build the command line thanks to make_cmd_str
        opt_list = [self.ref_path]
        cmd = make_cmd_str(self.indexer, self.index_opt, opt_list)
        print "Indexing {}".format(file_name(self.ref_path))

        # Try to execute the command
        try:
            proc = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
            # Launch the process
            stdout, stderr = proc.communicate()
            # Verify the return code
            if proc.returncode == 1:
                raise Exception ("An error occured during indexing with bwa index" + stderr)

            print (stderr)
            self.index_report = stderr
        
        # In case of error the files created during indexing are deleted
        except Exception  as E:
            print (E)
            print("Remove index files and folder")
            
            # Removing Index file
            for ext in ["amb", "ann", "bwt", "pac", "sa"]:
                f = "{}.{}".format(self.index_opt["p"], ext)
                try:
                    if path.isfile (f):
                        remove (f)
                except Exception:
                    pass
                    
            # Removing Index folder if empty
            try:
                rmdir(self.out_dir)
            except Exception:
                    pass
            print ("All clear")
            sys_exit()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class ExistingIndex(object):
    """
    @class GenerateIndex
    @brief Wrapper for bwa index. Create a reference index from a fasta file
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    def __repr__(self):
        msg = self.__str__()
        msg += "Bwa Index Path:\n{}\n".format(self.index_path)
        return msg

    def __str__(self):
        return "\n<Instance of {} from {}>\n".format(self.__class__.__name__, self.__module__)
        
    def __init__ (self, bwa_path, index_path):
        """
        Initialize the object and index the reference genome if necessary
        """
        # Creating object variables
        self.index_path = index_path
        
        print ("Checking index files") 
        # Checking if all index files needed by bwa are 
        for ext in ["amb", "ann", "bwt", "pac", "sa"]:
            index_file = self.index_path+"."+ext
            
            if not path.isfile (index_file):
                raise Exception ("{} does not exist".format(index_file))
            
            if path.getsize(index_file) == 0:
                raise Exception ("{} is empty".format(index_file))
        
        print ("All index files are valid")
