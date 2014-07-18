#~~~~~~~GLOBAL IMPORTS~~~~~~~#
# Standard library packages import
from multiprocessing import cpu_count
from os import remove, path

# Local Package import
from Utilities import mkdir, file_name, make_cmd_str

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class BwaWrapper(object):
    """
    @class BwaAligner
    @brief
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __repr__(self):
        msg = self.__str__()
        msg += "Aligner path : {}\n".format(self.aligner)
        msg += "Index : {}\n".format(self.index)
        msg += "Options :\n"
        for i, j in self.align_opt.items():
            msg += "\tFlag : {}\tValue : {}\n".format(i,j)
        return msg

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def __init__ (self, Index, align_opt=None, bwa_path = "bwa"):
        """
        Initialize the object and index the reference genome if necessary
        """
        # Creating object variables
        self.aligner = bwa_path+" mem"
        self.Index = Index
        self.align_opt = align_opt if align_opt else {}

        # By default the option t (number of thread to use) is setted to the max number of
        # available threads
        if "t" not in self.align_opt:
            self.align_opt["t"] = cpu_count()

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def align(self, R1_path, R2_path, out_path="out.sam"):
        """
        """
        # Build the command line with make_cmd_str
        opt_list = [self.Index.index_path, R1_path, R2_path, "> "+out_path]
        cmd = make_cmd_str(self.aligner, self.align_opt, opt_list)
        print cmd
        print "Aligning reads with {}".format(self.aligner)
        # Try to execute the command and verify if stdout is not None
        try:
            # Run the command line without stdin and asking both stdout and stderr
            stdout, stderr = run_command(cmd, None, True, True)

            # In bwa stderr return a report of alignment
            print (stderr)

        except Exception as E:
            print (E)
            raise
