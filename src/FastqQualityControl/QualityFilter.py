#~~~~~~~GLOBAL IMPORTS~~~~~~~#
# Local Package import
from Utilities import fill_between_graph

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class QualityFilter(object):
    """
    @class  QualityFilter
    @brief
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __repr__(self):
        return self.__str__() + self.get_report()

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def __init__(self, min_qual):
        """
        """
        # Init object variables
        self.min_qual = min_qual
        self.qual_pass = 0
        self.qual_fail = 0

        # Count each of the mean quality in a dictionary
        self.qdict = { i : 0 for i in range (41)}

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def filter(self, record):
        """
        Compute mean quality score and compare to the minimal quality required
        """
        # Compute the mean quality
        mean = sum(record.letter_annotations['phred_quality'])/len(record)
        # Fill the quality dictionary
        self.qdict[mean] += 1

        # Return the record if its quality is high enought
        if mean >= self.min_qual:
            self.qual_pass += 1
            return record
        else:
            self.qual_fail += 1
            return None

    def get_report (self):
        """
        Return a report
        """
        report = "====== QUALITY FILTER ======\n\n"
        report += "  Quality Threshold : {}\n".format(self.min_qual)
        report += "  Fail quality filter : {}\n".format(self.qual_fail)
        report += "  Pass quality filter : {}\n".format(self.qual_pass)

        if self.qual_fail+self.qual_pass != 0:
            # Simplify the dict by removing empty quality means
            reduce_qdict = {key: value for key, value in self.qdict.items() if value !=0}
            # Calculate simple probability
            mean_qual = sum([key*value for key, value in reduce_qdict.items()])/sum(reduce_qdict.values())
            min_qual = reduce_qdict.keys()[0]
            max_qual = reduce_qdict.keys()[-1]

            report += "  Mean quality : {}\n".format(mean_qual)
            report += "  Minimal quality : {}\n".format(min_qual)
            report += "  Maximal quality : {}\n".format(max_qual)

        report += "\n"
        return report

    def trace_graph (self):
        """

        """
        print ("\tCreating a graphical output of mean quality distribution")

        X = self.qdict.keys()
        Y = self.qdict.values()
        basename = "Mean_read_quality"
        img_type = "svg"
        title = "Distribution of mean read quality"
        xlabel = "Mean quality"
        ylabel = "Number of reads"

        fill_between_graph(X, Y, basename, img_type, title, xlabel, ylabel)
