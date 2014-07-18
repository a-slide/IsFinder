#~~~~~~~GLOBAL IMPORTS~~~~~~~#

# Standard library packages import
from os import path, remove

# Local Package import
from Utilities import file_basename, make_cmd_str, run_command

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class CreateDB(object):
    """
    @class CreateDB
    @brief Wrapper for makeblastdb index. Create a subject database from a fasta file
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __repr__(self):
        msg = self.__str__()
        msg += "Makeblastdb path : {}\n".format(self.makeblastdb)
        msg += "Blastn database path : {}\n".format(self.db_path)
        msg += "Options :\n"
        for i, j in self.makeblastdb_opt.items():
            msg += "\tFlag : {}\tValue : {}\n".format(i,j)
        return msg

    def __str__(self):
        return "\n<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def __init__ (self, ref_path, out_dir="./", makeblastdb_opt=None, makeblastdb="makeblastdb"):
        """
        Initialize the object and index the reference genome if necessary
        """
        # Creating object variables
        self.makeblastdb = makeblastdb
        self.ref_path = ref_path
        self.db_path = path.join(out_dir + file_basename(self.ref_path))

        # init an option dict and attribute defaut options
        self.makeblastdb_opt = index_opt if makeblastdb_opt else {}

        if "in" not in self.makeblastdb_opt:
            self.makeblastdb_opt["in"] = ref_path
        if "out" not in self.makeblastdb_opt:
            self.makeblastdb_opt["out"] = self.db_path
        if "dbtype" not in self.makeblastdb_opt:
            self.makeblastdb_opt["dbtype"] = "nucl"
        if "input_type" not in self.makeblastdb_opt:
            self.makeblastdb_opt["input_type"] = "fasta"

        # Make index with bwa index and store the path in an object variable
        self._make_db()
        print ("\tDatabase ready")

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _make_db(self):
        """
        Create a blastn database from ref_path using makeblastdb
        @exception (SystemError,OSerror) May be returned by run_command in case of invalid command line
        """

        # Build the command line thanks to make_cmd_str
        cmd = make_cmd_str(self.makeblastdb, self.makeblastdb_opt)
        print "Creating a blastn database for {}".format(file_basename(self.ref_path))

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
            print("Remove database files")

            # Removing Index file
            for ext in ["00.nhr", "nhr", "00.nin", "nin", "00.nsq", "nsq"]:
                f = "{}.{}".format(self.db_path, ext)
                try:
                    if path.isfile (f):
                        remove (f)
                except Exception:
                    print "Can't remove {}".format(f)
                    pass
            exit()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class ExistingDB(object):
    """
    @class ExistingDB
    @brief
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    def __repr__(self):
        msg = self.__str__()
        msg += "Blastn database path : {}\n".format(self.db_path)
        return msg

    def __str__(self):
        return "\n<Instance of {} from {}>\n".format(self.__class__.__name__, self.__module__)

    def __init__ (self, db_path):
        """
        """
        # Creating object variables
        self.db_path = db_path

        print ("Checking db files")
        # Checking if all index files needed by bwa are
        for ext in ["nhr", "nin", "nsq"]:
            f = "{}.{}".format(self.db_path, ext)

            if not path.isfile (f):
                raise Exception ("Invalid database : {} does not exist".format(f))
            if path.getsize(f) == 0:
                raise Exception ("Invalid database : {} is empty".format(f))

        print ("All index files are valid")
