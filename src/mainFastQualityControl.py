min_qual = 31
adapter_list = ["gcgcgcg", "atattaat", "tttaggcg", "aataccgattc"]
min_size = 50
input_qual = "fastq-sanger"
output_qual = "fastq-sanger"
R1_path = "../test/fastq/med_R1.fastq.gz"
R2_path = "../test/fastq/med_R2.fastq.gz"


from FastqQualityControl import FastqFilter, QualityFilter, AdapterTrimmer

qual = QualityFilter(min_qual)
adapt = AdapterTrimmer(adapter_list, min_size)
fast = FastqFilter(qual, adapt, input_qual, output_qual)

fast.filter(R1_path, R2_path)
