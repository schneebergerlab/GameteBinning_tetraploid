samtools view test_bam.bam | T10xbam2fq - test
# you would get an error in this test as read number is not paired in R1 and R2; which should not happen for paired-end readsets.
