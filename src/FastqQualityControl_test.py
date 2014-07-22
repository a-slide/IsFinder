from FastqQualityControl.QualityFilter import QualityFilter
from FastqQualityControl.AdapterTrimmer import AdapterTrimmer
from FastqQualityControl.FastqFilterPP import FastqFilterPP
from FastqQualityControl.FastqFilter import FastqFilter
from ssw_wrap import Aligner

def test ():

    filter_quality = True
    min_qual = 25

    trim_adapter = True
    adapter_path = "/home/adrien/Programming/Python/IsFinder/test/adapter.fa"

    R1 = "/home/adrien/Programming/Python/IsFinder/test/fastq/all_AAV_R1.fastq.gz"
    R2 = "/home/adrien/Programming/Python/IsFinder/test/fastq/all_AAV_R2.fastq.gz"

    # Define a quality filter object
    if filter_quality:
        q_filter = QualityFilter (min_qual)
    else:
        q_filter = None

    # Define a adapter trimmer object
    if trim_adapter:
        a = Aligner()
        trimmer = AdapterTrimmer(a, adapter_path)
    else:
        trimmer = None

    # Define the global fastq filter
    F = FastqFilterPP(R1, R2, quality_filter=q_filter, adapter_trimmer=trimmer, input_qual="fastq-sanger")

    print (repr(F))
    print (repr(q_filter))
    print (repr(a))
    print (repr(trimmer))


def test2 ():

    filter_quality = True
    min_qual = 25

    trim_adapter = True
    adapter_path = "/home/adrien/Programming/Python/IsFinder/test/adapter.fa"

    R1_path = "/home/adrien/Programming/Python/IsFinder/test/fastq/all_AAV_R1.fastq.gz"
    R2_path = "/home/adrien/Programming/Python/IsFinder/test/fastq/all_AAV_R2.fastq.gz"

    # Define a quality filter object
    if filter_quality:
        q_filter = QualityFilter (min_qual)
    else:
        q_filter = None

    # Define a adapter trimmer object
    if trim_adapter:
        a = Aligner(report_cigar=True)
        trimmer = AdapterTrimmer(a, adapter_path)
    else:
        trimmer = None

    # Define the global fastq filter
    f_filter = FastqFilter(q_filter, trimmer)

    print(f_filter.filter(R1_path, R2_path))
    print (repr(q_filter))
    print (repr(a))
    print (repr(trimmer))
