#~~~~~~~GLOBAL IMPORTS~~~~~~~#
# Standard library packages import
from multiprocessing import cpu_count
from sys import exit as sys_exit

# Third party package import

# Local Package import
from Utilities import mkdir, file_name, make_cmd_str, run_command

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class BwaWrapper(object):
    """
    @class BwaWrapper
    @brief Wrapper for bwa index and bwa mem aligner
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#


    def __repr__(self):
        return "{}\nSTANDARD OUTPUT\n{}\nSTANDARD ERROR\n{}".format(
            self.__str__(),
            self.index_stdout,
            self.index_stderr)

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def __init__ (self, bwa_path, ref_path, make_index=True, align_opt={}, index_opt={}):
        """
        Initialize the object by indexing the reference genome if necessary
        """
        # Store path to aligner and indexer programs
        self.indexer = bwa_path+" index"
        self.aligner = bwa_path+" mem"

        # Create an index for the reference genome
        if make_index:
            # Complete the index option dict with default values
            if "a" not in index_opt:
                index_opt["a"] = "bwtsw"
            if "p" not in index_opt:
                mkdir("./bwa_index/")
                index_opt["p"] = "./bwa_index/"+file_name(ref_path)

            # Make index with bwa index and store the path in an object variable
            self.index_stdout, self.index_stderr = self._make_index(ref_path, index_opt)
            self.index_path = index_opt["p"]

        # Else just attribute the user index to the object variable
        else:
            self.index_path = ref_path

        # Complete the aligner option dict with default values
        if "t" not in align_opt:
            align_opt["t"] = cpu_count()

    #~~~~~~~CLASS METHODS~~~~~~~#

    def align(self, R1_path, R2_path, opt_dict):
        """
        """
       # Build the command line thanks to make_cmd_str
        opt_list = [self.index_path, R1_path, R2_path, "> outfile.sam"]
        cmd = make_cmd_str(self.aligner, opt_dict, opt_list)

        print "Aligning reads against {} with {}".format(file_name(self.index_path, self.aligner))
        print "This may last several hours\n"

        try:
            stdout, stderr = run_command(cmd, None, False, True, True)
            assert stdout
            return stdout, stderr

        except AssertionError as E:
            print "An error occured during reference indexing with bwa index"
            print (stderr)
            sys_exit()


    def _make_index(self, ref_path, opt_dict):
        """
        """
        # Build the command line thanks to make_cmd_str
        opt_list = [ref_path]
        cmd = make_cmd_str(self.indexer, opt_dict, opt_list)

        print "Indexing {}.".format(file_name(ref_path))
        print "This may last several hours depending of the size of the reference genome\n"

        try:
            stdout, stderr = run_command(cmd, None, False, True, True)
            assert stdout
            return stdout, stderr

        except AssertionError as E:
            print "An error occured during reference indexing with bwa index"
            print (stderr)
            sys_exit()
