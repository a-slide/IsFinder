from ShortReadAligner.BwaWrapper import BwaWrapper
from pprint import pprint
import pysam

bwa = BwaWrapper("bwa", "../../../../bioinformatics/NGS/1_INDEX/mm10_AAV_backbone/Bwa_mm10_AAV_pAAV/mm10_AAV_pAAV.fa", False )
tmp = bwa.align("fastq/100seq_R1.fastq.gz", "fastq/100seq_R2.fastq.gz")

with open (tmp, "r") as f:
    print(f.read())
 


sam = pysam.Samfile(tmp)
#pprint(sam.header)

prop_read = []
prop_pair = []
chim_pair = []
chim_read = []
unmapped = []

while True:
    try:
        read1 = sam.next()
        read2 = sam.next()
       
        if read1.mapq < 30 or read2.mapq < 30 or read1.is_unmapped or read2.is_unmapped:
            unmapped.extend([read1,read2])
            continue
                            
        if read1.tid != read2.tid:
            chim_pair.extend([read1,read2])
            continue

        SA = [value for tag, value in read1.tags if tag == "SA"]
        if SA:
            SA_name =  SA[0].split(',')[0]
            if sam.getrname(read1.tid) != SA_name:
                chim_read.extend([read1,read2])
                continue

        SA = [value for tag, value in read2.tags if tag == "SA"]
        if SA:
            SA_name =  SA[0].split(',')[0]
            if sam.getrname(read2.tid) != SA_name:
                chim_read.extend([read1,read2])
                continue
        
        if read1.tid == read2.tid:
            prop_pair.extend([read1,read2])
            continue
        
    except Exception as E:
        print "EOF"
        print E
        break
