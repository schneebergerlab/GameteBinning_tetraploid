# A test on separating reads/cells with 10x barcode

minimumRP=5000 
sample=A
R1=/path/to/4279_A_run615_SI-GA-D4_S1_L003_R1_001.fastq.gz
R2=/path/to/4279_A_run615_SI-GA-D4_S1_L003_R2_001.fastq.gz

asCellseparator 16 ${R1} ${R2} ${minimumRP} 500000 asCellseparator_sep_cells

# here 500000 is minimum number of read pairs collected for all tmp barcodes to output as a part -- you can change it.

