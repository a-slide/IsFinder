#~~~~~~~GLOBAL IMPORTS~~~~~~~#

# Local Package import
from Utilities import import_seq, fill_between_graph

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class AdapterTrimmer(object):
    """
    @class  AdapterTrimmer
    @brief
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __repr__(self):
        return self.__str__() + self.get_report()

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def __init__ (self, aligner, adapter_path, min_size, read_len):
        """
        """
        # Init object variables
        self.min_size = min_size
        self.read_len = read_len
        self.aligner = aligner

        # Import a list of adapters and add the reverse complements of adapters to the list
        self.adapter_list = import_seq(adapter_path, "list", "fasta")
        assert self.adapter_list, "The adater list is empty"

        adapter_rc = []
        for i in self.adapter_list:
            rc = i.reverse_complement()
            rc.id = i.id + "_RC"
            adapter_rc.append(rc)
        self.adapter_list.extend(adapter_rc)

        # Initialize description field to 0 to serve as a counter of match found
        for i in self.adapter_list:
            i.description = 0

        # Initialize generic counters
        self.seq_untrimmed = 0
        self.seq_trimmed = 0
        self.base_trimmed = 0
        self.len_pass = 0
        self.len_fail = 0

        # Dict to accumulate remaining coverage after trimming
        self.coverage = { i : 0 for i in range (self.read_len)}

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def trimmer (self, record):
        """
        """
        #stime = time()
        match_list = []

        for adapter in self.adapter_list:
            # Find matches of the adapter along the current read
            adapter_match = self.aligner.find_match(str(adapter.seq).lower(), str(record.seq).lower())
            # Update the counter of the adapter with the number of matches found
            adapter.description += (len(adapter_match))
            # Add to the list of
            match_list.extend(adapter_match)

        #print ("FIND MATCH = {} ms".format((time()-stime)/1000))
        #stime = time()

        # In case no match were found, the sequence doesn't need to be modify
        if not match_list:
            self.seq_untrimmed += 1
            return record

        # Else find the longer interval without adaptor matches
        start, end = self._longer_interval (match_list, len(record))
        assert 0 <= start < end <= len(record), "Invalid interval returned"

        #print ("LONGER INTERVAL = {} ms".format((time()-stime)/1000))
        #stime = time()

        # Update counters
        self.seq_trimmed += 1
        self.base_trimmed += (len(record)-(end-start))
        self._update_coverage (start, end)

        #print ("COUNTUP = {} ms".format((time()-stime)/1000))

        # Return a slice of the reccord corresponding to the longer interval
        if end-start >= self.min_size:
            self.len_pass +=1
            return record[start:end]
        # Or None if smaller than min_size
        else:
            self.len_fail +=1
            return None

    def get_report (self):
        """
        """
        report = "====== ADAPTER TRIMMER ======\n\n"
        report += "  List of adapters imported for trimming\n"
        for i in self.adapter_list:
            report += "    Name : {}\tSequence : {}\tTimes trimmed : {}\n".format(
                i.id, i.seq, i.description)
        report += "\n  Sequences untrimmed : {}\n".format(self.seq_untrimmed)
        report += "  Sequences trimmed : {}\n".format(self.seq_trimmed)
        report += "  DNA base trimmed : {}\n".format(self.base_trimmed)
        report += "  Fail len filtering: {}\n".format(self.len_fail)
        report += "  Pass len filtering : {}\n".format(self.len_pass)
        report += "  Total pass : {}\n\n".format(self.len_pass+self.seq_untrimmed)
        report += self.aligner.get_report()
        return report

    def trace_graph (self):
        """

        """
        print ("\tCreating a graphical output of read coverage after trimming")

        X = self.coverage.keys()
        Y = self.coverage.values()
        basename = "Read_coverage_trimming"
        img_type = "svg"
        title = "Coverage of reads after trimming"
        xlabel = "Position"
        ylabel = "Number of overlapping reads"

        fill_between_graph(X, Y, basename, img_type, title, xlabel, ylabel)

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _longer_interval(self, match_list, maxsize):
        """
        Find the first larger interval that do not overlapp any match in match list
        This strategy allow to use an unsorted list of match
        """
        # Initialize a list of boolean to False of the same size as the read
        coverage = [False for i in range(maxsize)]

        # Flag positions overlapped by a read by changing the boolean to True
        for begin, end in match_list:
            for i in range (begin, end):
                coverage[i] = True

        # Read through the list to find the longer inteval between True flags
        start_max = end_max = inter_max = start = inter = 0
        for i in range(maxsize):
            if coverage[i]:
                start = i+1
                inter = 0
            else:
                inter += 1
                if inter > inter_max:
                    inter_max = inter
                    start_max = start
                    end_max = i

        #print ("Longer interval = {} [{}:{}]".format(inter_max, start_max+1, end_max-1))
        return start_max, end_max

    def _update_coverage (self, start, end):
        """
        """
        for i in range(start, end):
            self.coverage[i] += 1
