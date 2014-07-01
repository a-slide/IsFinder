#~~~~~~~GLOBAL IMPORTS~~~~~~~#
# Standard library packages import
from multiprocessing import cpu_count
from sys import exit as sys_exit
from subprocess import Popen, PIPE
from os import remove, path, rmdir, mkfifo
from tempfile import mkstemp
from tempfile import mkdtemp
from shutil import rmtree

# Third party package import
import pysam

# Local Package import
from Utilities import mkdir, file_name, make_cmd_str

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class BwaAligner(object):
    """
    @class BwaAligner
    @brief 
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#


    def __repr__(self):
        msg = self.__str__()
        return msg

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def __init__ (self, bwa_path, Index, align_opt=None):
        """
        Initialize the object and index the reference genome if necessary
        """
        # Creating object variables
        self.aligner = bwa_path+" mem"
        self.Index = Index
        self.align_opt = align_opt if align_opt else {}
        self.align_report = ""
        
        # By default the option t (number of thread to use) is setted to the max number of
        # available threads
        if "t" not in self.align_opt:
            self.align_opt["t"] = cpu_count()

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def align_tempfile(self, R1_path, R2_path, out_prefix="out"):
        """
        """
        # Build the command line with make_cmd_str
        tmp_out = mkstemp()
        opt_list = [self.Index.index_path, R1_path, R2_path, "> "+tmp_out]
        cmd = make_cmd_str(self.aligner, self.align_opt, opt_list)
         
        print "Aligning reads against {} with {}".format(file_name(self.index_path), self.aligner)
        # Try to execute the command and verify if stdout is not None
        try:
            proc = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
            # Launch the process
            stdout, stderr = proc.communicate()
            
            # Verify the return code
            if proc.returncode == 1:
                raise Exception ("An error occured during alignment with bwa mem" + stderr)
            
            print (stderr)
            return (tmp_out.name)
            
        except Exception as E:
            print (E)
            sys_exit()
    
    
    def align_fifo(self, R1_path, R2_path, out_prefix="out"):
        """
        """
        # Create a fifo
        tmpdir = mkdtemp()
        out_name = "{}/{}.sam".format(tmpdir, out_prefix)
        fifo = out_name 
        mkfifo(fifo) # It's not a file anymore it's a fifo
        print fifo
        
        # Build the command line with make_cmd_str
        opt_list = [self.Index.index_path, R1_path, R2_path, "> "+out_name]
        cmd = make_cmd_str(self.aligner, self.align_opt, opt_list)
        print cmd
        print "Aligning reads with {}".format(self.aligner)
        # Try to execute the command and verify if stdout is not None
        try:
            proc = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
            #stdout, stderr = proc.communicate()
            
            #if proc.returncode == 1:
            #    raise Exception ("An error occured during alignment with bwa mem" + stderr)
            
            # Read from the named pipe
            samfile = pysam.Samfile(fifo, "r")

            # Print out the names of each record
            for read in samfile:
                print read
            
            # Clean up the named pipe and associated temp directory
            rmtree(tmpdir)
            
        except Exception as E:
            print (E)
            sys_exit()
