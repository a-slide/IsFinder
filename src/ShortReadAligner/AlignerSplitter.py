
#~~~~~~~GLOBAL IMPORTS~~~~~~~#
# Standard library packages import
from time import time

# Third party package import

# Local Package import

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class AlignerSplitter(object):
    """
    @class  ShortReadAligner
    @brief
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __repr__(self):
        return self.__str__() + self.get_report()

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def __init__ (self, FastaReader, aligner):
        """
        """
        # Store parameters in object variables
        self.aligner = aligner
        self.spliter = spliter

    #~~~~~~~CLASS METHODS~~~~~~~#

    def align(self, R1_path, R2_path, align_opt=None, split_opt=None):
        """
        """
        pass



    def get_report (self, exec_time=None, R1_path=None, R2_path=None):
        """
        Return a report
        """
        pass
