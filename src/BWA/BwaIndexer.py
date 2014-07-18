#~~~~~~~GLOBAL IMPORTS~~~~~~~#

# Standard library packages import
from os import path, remove

# Local Package import
from Utilities import file_basename, make_cmd_str, run_command

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class CreateIndex(object):
    """
    @class CreateIndex
    @brief Wrapper for bwa index. Create a reference index from a fasta file
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __repr__(self):
        msg = self.__str__()
        msg += "bwa index path : {}\n".format(self.indexer)
        msg += "Blastn database path : {}\n".format(self.index_path)
        msg += "Options :\n"
        for i, j in self.index_opt.items():
            msg += "\tFlag : {}\tValue : {}\n".format(i,j)
        return msg

    def __str__(self):
        return "\n<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def __init__ (self, ref_path, out_dir="./", index_opt=None, bwa_path = "bwa"):
        """
        Initialize the object and index the reference genome if necessary
        """
        # Creating object variables
        self.indexer = bwa_path+" index"
        self.ref_path = ref_path
        self.index_path = path.join(out_dir + file_basename(self.ref_path))

        # init an option dict and attribute defaut options
        self.index_opt = index_opt if index_opt else {}

        if "a" not in self.index_opt:
            self.index_opt["a"] = "bwtsw"
        if "p" not in self.index_opt:
            self.index_opt["p"] = self.out_dir + file_name(self.ref_path)

        # Make index with bwa index and store the path in an object variable
        self._make_index()
        print ("\tIndex ready")

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _make_index(self):
        """
        Create a bwa index from ref_path using bwa index
        @exception (SystemError,OSerror) May be returned by run_command in case of invalid command line
        """

        # Build the command line thanks to make_cmd_str
        opt_list = [self.ref_path]
        cmd = make_cmd_str(self.indexer, self.index_opt, opt_list)
        print "Creating a BWA index for {}".format(file_basename(self.ref_path))

        # Try to execute the command
        try:

            # Run the command line without stdin and asking both stdout and stderr
            stdout, stderr = run_command(cmd, None, True, True)

            # Verify the output
            if not stdout:
                raise Exception ("Error, no data received from standard output\n"+stderr)

            print (stdout)
            return True

        # In case of error the files created during indexing are deleted
        except Exception  as E:
            print (E)
            print("Remove index files and folder")

            # Removing Index file
            for ext in ["amb", "ann", "bwt", "pac", "sa"]:
                f = "{}.{}".format(self.index_path, ext)
                try:
                    if path.isfile (f):
                        remove (f)
                except Exception:
                    print "Can't remove {}".format(f)
                    pass
            exit()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class ExistingIndex(object):
    """
    @class ExistingIndex
    @brief
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    def __repr__(self):
        msg = self.__str__()
        msg += "Bwa Index Path : {}\n".format(self.index_path)
        return msg

    def __str__(self):
        return "\n<Instance of {} from {}>\n".format(self.__class__.__name__, self.__module__)

    def __init__ (self, index_path):
        """
        """
        # Creating object variables
        self.index_path = index_path

        print ("Checking index files")
        # Checking if all index files needed by bwa are
        for ext in ["amb", "ann", "bwt", "pac", "sa"]:
            f = "{}.{}".format(self.index_path, ext)

            if not path.isfile (f):
                raise Exception ("{} does not exist".format(f))
            if path.getsize(f) == 0:
                raise Exception ("{} is empty".format(f)

        print ("All index files are valid")
