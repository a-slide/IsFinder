min_qual = 31
adapter_list = ["gcgcgcg", "atattaat", "tttaggcg", "aataccgattc"]
min_size = 50
input_qual = "fastq-sanger"
output_qual = "fastq-sanger"
R1_path = "../test/fastq/med_R1.fastq.gz"
R2_path = "../test/fastq/med_R2.fastq.gz"


from FastqQualityControl import FastqFilter, QualityFilter, AdapterTrimmer, PairwiseAligner

a = FastqQualityControl.PairwiseAligner() # VALID
b = FastqQualityControl.AdapterTrimmer(a, "adapter.fa", 100, 150)  # VALID
c = QualityFilter(min_qual)
d = FastqFilter(b,c...
