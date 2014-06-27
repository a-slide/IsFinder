min_qual = 31
min_size = 100
input_qual = "fastq-sanger"
R1_path = "../test/fastq/small_R1.fastq.gz"
R2_path = "../test/fastq/small_R2.fastq.gz"
adapter_path = "../test/adapter.fa"
read_len = 150


from FastqQualityControl import FastqFilter, QualityFilter, AdapterTrimmer, PairwiseAligner
a = PairwiseAligner() # VALID
b = AdapterTrimmer(a, adapter_path, min_size, read_len)  # VALID
c = QualityFilter(min_qual)# VALID
d = FastqFilter(R1_path, R2_path, c, b, input_qual)# VALID

d.filter()
