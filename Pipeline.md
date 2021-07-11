Pipeline
=
This is the pipeline explaining how gamete binning in tetraploid (potato) works.

##### Step 0. Prepare data

All data are from a tetraploid (potato cultivar) of interest, including PacBio HiFi reads from somatic tissue and short reads from
single-cell sequencing of sufficient (e.g., hundreds of) gamete genomes. In this example, we use the following (available on NCBI Bioproject: PRJNA726019):

* PacBio HiFi: long_reads_raw.fa
* 10x Genomics+Illumina (sc): e.g., A_seq4414plus4431_R1.fastq.gz and A_seq4414plus4431_R2.fastq.gz (from 4414_A_run633_SI-GA-A1,4461_A_run636_SI-GA-A1), B_seq4414plus4431_R1.fastq.gz and B_seq4414plus4431_R2.fastq.gz (4414_B_run633_SI-GA-B1,4461_B_run636_SI-GA-B1)
* 10x Genomics+Illumina (sm): e.g., C_seq2806_R1.fastq.gz, C_seq2806_R2.fastq.gz (for depth analysis)

Suppose all these raw data are collected in the path below, and for convenience, create softlinks for fastq files (Note, 10x Genomics tools need the full name, so we use both namings),

    wd=/path/to/reads/
    cd ${wd}
    
    ln -s  A_seq4414plus4431_R1.fastq.gz gamete_libA_R1.fastq.gz
    ln -s  A_seq4414plus4431_R2.fastq.gz gamete_libA_R2.fastq.gz
    
    ln -s  B_seq4414plus4431_R1.fastq.gz gamete_libB_R1.fastq.gz
    ln -s  B_seq4414plus4431_R2.fastq.gz gamete_libB_R2.fastq.gz

##### Step 1. Trim reads

Trim 16 bp barcodes off R1's (10x Genomics library setting, including hexamer so 22 bp trimmed off),

    wd=/path/to/reads/
    cd ${wd}
    
    T10X_barcode_trimmer gamete_libA_R1.fastq.gz gamete_libA_R2.fastq.gz
    T10X_barcode_trimmer gamete_libB_R1.fastq.gz gamete_libB_R2.fastq.gz
    T10X_barcode_trimmer C_seq2806_R1.fastq.gz C_seq2806_R1.fastq.gz

This leads to

* gamete_libA_R1_clean.fastq.gz, gamete_libA_R2_clean.fastq.gz 
* gamete_libB_R1_clean.fastq.gz, gamete_libB_R2_clean.fastq.gz
* C_seq2806_R1_clean.fastq.gz, C_seq2806_R2_clean.fastq.gz

##### Step 2. Preliminary assembly

    wd=/path/to/pre_assembly/
    cd ${wd}
    
    otavahifi=/path/to/reads/long_reads_raw.fa
    hifiasm -t 10 -o otava ${otavahifi} >hifiasm.log

This leads to preliminary assembly (we select the utg-level) by 

    cat otava.p_utg.gfa | grep '^S' | cut -f2,3 | awk '{print ">"$1"\n"$2}' > otava.p_utg.fasta

* otava.p_utg.fasta

(Polish the preliminary assembly - not necessary, so details not given here)

##### Step 3. Curation of assembly using read alignment depth (or, purge redundant contigs representing the same genomic regions)    

This leads to a version of manually curated assembly (please refer to manuscript supplementary information, section "Initial tetraploid genome assembly, polishing and purging" for details)

    wd=/path/to/curated_asm/
    cd ${wd}
    
* HiFiasm_ref_6366long_ctgs_selected.fasta
* HiFiasm_ref_6366long_ctgs_selected.chrsizes (this is the contig size file with two tab-separated columns: contig_id	contig_size)

Note, it is a mixture of four haplotypes (with potentially collapsed homozygous regions between any more than one haplotype).

Index the sequence as reference for later steps,

    refgenome=HiFiasm_ref_6366long_ctgs_selected.fasta
    bowtie2-build -f ${refgenome} ${refgenome} --threads 4

##### Step 4. Read alignment of various sequencings for haplo/diplo/triplo/tetraplotig marker identification (note we used sc, sm and HiFi depth to collect more power to differentiate different types of contigs)

    wd=/path/to/marker_creation/
    cd ${wd}

Align pooled gamete reads to the reference

    refgenome=/path/to/curated_asm/HiFiasm_ref_6366long_ctgs_selected.fasta
    bowtie2 -x ${refgenome} -1 /path/to/reads/gamete_libA_R1_clean.fastq.gz,/path/to/reads/gamete_libB_R1_clean.fastq.gz -2 /path/to/reads/gamete_libA_R2_clean.fastq.gz,/path/to/reads/gamete_libB_R2_clean.fastq.gz -p 20 | samtools view -@ 20 -bS - | samtools sort -@ 20 -o gamete_ManualCurated.bam -
    
Remove duplicates and get position-wise depth

    java -jar picard.jar MarkDuplicates I=gamete_ManualCurated.bam O=gamete_ManualCurated_markeduplicates.bam M=gamete_ManualCurated_marked_dup_metrics.txt
    samtools index gamete_ManualCurated_markeduplicates.bam
    
Similarly, align sm related reads and get bam file 

    refgenome=/path/to/curated_asm/HiFiasm_ref_6366long_ctgs_selected.fasta
    bowtie2 -x ${refgenome} -1 /path/to/reads/C_seq2806_R1_clean.fastq.gz -2 /path/to/reads/C_seq2806_R2_clean.fastq.gz -p 20 | samtools view -@ 20 -bS - | samtools sort -@ 20 -o sm_ManualCurated.bam -
    java -jar picard.jar MarkDuplicates I=sm_ManualCurated.bam O=sm_ManualCurated_markeduplicates.bam M=sm_ManualCurated_marked_dup_metrics.txt
    samtools index sm_ManualCurated_markeduplicates.bam
    
Align HiFi reads, remove non-primary alignments and get position-wise depth 
    
    otavahifi=/path/to/reads/long_reads_raw.fa
    refgenome=/path/to/curated_asm/HiFiasm_ref_6366long_ctgs_selected.fasta
    minimap2 -ax map-pb -t 20 -N 1 --secondary=no ${refgenome} ${otavahifi} | samtools view -@ 20 -bS - | samtools sort -@ 20 -o HiFi_ManualCurated.bam -
    samtools view -h -F 3840 -bS HiFi_ManualCurated.bam | samtools sort -o HiFi_ManualCurated_clean.bam - 
    samtools depth -Q 1 -a HiFi_ManualCurated_clean.bam > HiFi_ManualCurated_depth_clean.txt

Get position-wise sequencing depth 

    samtools depth -Q 1 -a gamete_ManualCurated.bam sm_ManualCurated.bam HiFi_ManualCurated_clean.bam > otava_Hifiasm_ref_pb_illu_depth.txt
    awk '{printf ("%s\t%s\t%s\n", $1, $2, $3+$4+$5)}' otava_Hifiasm_ref_pb_illu_depth.txt > otava_Hifiasm_ref_pb_illu_depth_sum.txt
    rm otava_Hifiasm_ref_pb_illu_depth.txt
    
##### step 5: find distribution of average depth at non-overlapping windows: winstep = winsize

    chrsizes=/path/to/HiFiasm_ref_6366long_ctgs_selected.ctgsizes
    winsize=10000
    
    # => cnv_winsize10000_step10000_hq.txt
    samplecol=1
    avgdepth=113 # 137*0.85
    CNV_HQ_v3 ${chrsizes} otava_Hifiasm_ref_pb_illu_depth_sum.txt ${winsize} ${winsize} ${samplecol} ${avgdepth}
    
    # => cnv_winsize10000_step10000_hq_HiFi_only.txt
    samplecol=1
    avgdepth=30 # HiFi only
    CNV_HQ_v3 ${chrsizes} HiFi_ManualCurated_depth_clean.txt ${winsize} ${winsize} ${samplecol} ${avgdepth}
    
    
##### step 6: new window marker generation => cnv_winsize10000_step10000_hq_markers_20210714_wsize50kb_final.txt

    paste cnv_winsize10000_step10000_hq.txt cnv_winsize10000_step10000_hq_HiFi_only.txt | cut -f1-7,12,13 > cnv_winsize10000_step10000_hq_merged_vs_hifi.txt
    tig_marker_finder cnv_winsize10000_step10000_hq_merged_vs_hifi.txt 113 30 50000 20210714 > tig_marker_finder.log

##### Step 7. 10x Genomics barcode correction and nuclei separation

    wd=/path/to/individual_nuclei_extraction/
    cd ${wd}

We use DM assembly at chr-level - potato_dm_v404_all_pm_un.fasta - as reference (download: http://solanaceae.plantbiology.msu.edu/pgsc_download.shtml), and prepare a chrsizes file - potato_dm_v404_all_pm_un_modified.chrsizes (two columns, tab-separated with chr_id chr_size):

    sed -i 's/chr00/chrX/' potato_dm_v404_all_pm_un.fasta
    sed -i 's/ChrUn/chrY/' potato_dm_v404_all_pm_un.fasta

Correspondingly, we prepare "a JSON file - /file_aux/contig_defs.json - describing primary contigs" (with below), 

    {
            "species_prefixes": [""],
            "primary_contigs": [
            "chr01", "chr02", "chr03","chr04", "chr05", "chr06", "chr07", "chr08", "chr09", "chr10", "chr11", "chr12", "chrX", "chrY"
            ],
            "sex_chromosomes": {
                    "_male": {
                            "chrX": 1,
                            "chrY": 1
                     },
                     "_female": {
                            "chrX": 2,
                            "chrY": 0
                     }
            },
            "non_nuclear_contigs": ["chrM"]
    }

and then index the genome with cellranger-dna,

    refgenome=potato_dm_v404_all_pm_un.fa
    cellranger-dna mkref ${refgenome} /path/to/file_aux/contig_defs.json

This will create a new reference folder of "/refdata-potato_dm_v404_all_pm_un/".

Correct 10x Genomics barcodes (note, if there are multiple libraries, this step needs to be done library by library, as same barcodes might be shared across libraries. However, different runs of the same library can be run together by setting option --sample=libx-run-1,libx-run-2\[,...\]),

    cellranger-dna cnv --id=A_all3runs --reference=/path/to/refdata-potato_dm_v404_all_pm_un/ --fastq=/path/to/reads/ --sample=4414_A_run633_SI-GA-A1,4461_A_run636_SI-GA-A1 --localcores=40 --localmem=30
    cellranger-dna cnv --id=B_all3runs --reference=/path/to/refdata-potato_dm_v404_all_pm_un/ --fastq=/path/to/reads/ --sample=4414_B_run633_SI-GA-B1,4461_B_run636_SI-GA-B1 --localcores=40 --localmem=30

Sort the bam (from the above, which is with corrected barcode information) with read name

    bam=/path/to/A_all3runs/outs/possorted_bam.bam
    samtools sort -n ${bam} -o RNsorted_bam.bam

Create reads ligated with corrected barcode,

    bam=/path/to/A_all3runs/outs/RNsorted_bam.bam
    samtools view ${bam} | T10xbam2fq - A

This leads to

    A_fqfrom10xBam_bxCorrected_R1.fq.gz
    A_fqfrom10xBam_bxCorrected_R2.fq.gz

Extract individual nuclei

    R1=A_fqfrom10xBam_bxCorrected_R1.fq.gz
    R2=A_fqfrom10xBam_bxCorrected_R2.fq.gz
    barcode_len=16
    minimumRP=40000
    min_readpair=20000000
    asCellseparator ${barcode_len} ${R1} ${R2} ${minimumRP} ${min_readpair} cells_sep

Get list of good barcodes,

    cd /path/to/cells_sep
    awk '$2>=40000' asCellseparator_intermediate_raw_barcode_stat.txt > barcode_over_40000rpairs.list
   
Similarly, extract cell-wise sequencings from library B.

##### step 8. read alignment for each pollen genome sequencing

    wd=/path/to/individual_nuclei_read_align/
    cd ${wd}
    
    cellpath=/path/to/individual_nuclei_extraction/
    cd ${cellpath}/
    refgenome=/path/to/HiFiasm_ref_6366long_ctgs_selected.fasta
    for sample in A B; do
        cd ${cellpath}/sample_${sample}_asCellseparator_40krp/
        # get individual barcode list
        echo */ | sed 's/ /\n/g' | sed -e 's/\///g' > ${sample}_this_barcode_list
        # mapping within barcode list
        while read bc; do 
            cd ${cellpath}/sample_${sample}_asCellseparator_40krp/${bc}
            R1=`ls part*_${bc}_R1.fq.gz`
            R1=${R1//[[:space:]]/,}
            R2=`ls part*_${bc}_R2.fq.gz`
            R2=${R2//[[:space:]]/,}
            bowtie2 -p 1 -x ${refgenome} -1 \${R1} -2 \${R2} 2> longctg_bowtie2.err | samtools view -@ 1 -bS - | samtools sort -@ 1 -o longctg_${bc}.bam -
            samtools index longctg_${bc}.bam
            java -jar picard.jar MarkDuplicates I=longctg_${bc}.bam O=longctg_${bc}_markeduplicates.bam M=longctg_${bc}_marked_dup_metrics.txt
            samtools index longctg_${bc}_markeduplicates.bam
        done < ${sample}_this_barcode_list
        cd ${cellpath}
    done
    
##### step 9. genotype each pollen at each contig marker - note 1: "-F 3840" == "-F 256 -F 512 -F 1024 -F 2048", excluding all kinds of non-primary alignment! - note 2: MQ=5 is too stringent that it was observed some pollen lost coverage at haplotig markers (where there were primary-reads with MQ=1); need to use MQ1: found with IGV

    wd=/path/to/individual_nuclei_read_align/
    cd ${wd}
    
    cut -f1-3 /path/to/cnv_winsize10000_step10000_hq_markers_20210714_wsize50kb_final.txt > cnv_winsize10000_step10000_hq_markers_20210714_wsize50kb_final.bed
    MQ=1
    cellpath=/path/to/individual_nuclei_extraction/
    for sample in A B; do
        while read r; do 
            cd ${cellpath}/sample_${sample}_asCellseparator_40krp/${r}
            samtools view -h -F 3840 -q ${MQ} -bS longctg_${r}_markeduplicates.bam | samtools sort -o longctg_${r}_markeduplicates_MQ${MQ}_clean.bam -
            samtools index longctg_${r}_markeduplicates_MQ${MQ}_clean.bam
            bedtools coverage -counts -a /path/to/cnv_winsize10000_step10000_hq_markers_20210714_wsize50kb_final.bed -b longctg_${r}_markeduplicates_MQ${MQ}_clean.bam -bed > longctg_${r}_win_marker_read_count_MQ${MQ}.bed            
        done < ${cellpath}/sample_${sample}_asCellseparator_40krp/${sample}_this_barcode_list
    done
    
##### step 10 add more info on markers to read cnt file 

    cellpath=/path/to/individual_nuclei_extraction/
    cd ${cellpath}
    marker=/path/to/cnv_winsize10000_step10000_hq_markers_20210714_wsize50kb_final.txt
    MQ=1
    for sample in A B; do
        while read r; do 
            echo -e ${sample}"\t"${r}
            cd ${cellpath}/sample_${sample}_asCellseparator_40krp/${r}
            paste longctg_${r}_win_marker_read_count_MQ${MQ}.bed ${marker} > longctg_${r}_win_marker_read_count_MQ${MQ}_updated.bed
        done < ${cellpath}/sample_${sample}_asCellseparator_40krp/${sample}_this_barcode_list
    done
    
Note, we finally selected 717 nuclei to perform gamete binning, given under: "/aux_data/", i.e., longctg_s4p3_selected_717_good_nucei_lib[A|B].txt

##### step 11 prepare nuclei depth data

    wd=/path/to/s5_selected_long_contigs_sc_read_coverage_genotype_v2/
    ls /path/to/individual_nuclei_extraction/sample_*_asCellseparator_40krp/*/longctg_*_win_marker_read_count_MQ1_updated.bed > longctg_list_bed_files.txt
    #
    >longctg_list_bed_files_selected717.txt
    for sample in A B; do 
        while read bc; do 
             grep ${bc} longctg_list_bed_files.txt >> longctg_list_bed_files_selected717.txt
        done < ../s4_selected_long_contigs_sc_read_coverage_calculation/longctg_s4p3_selected_717_good_nucei_lib${sample}.txt
    done
    
##### step 12 CORE FUNCTION: find linkage groups and match "homologous" linkage groups
 
Input 1. bed files on read counts for each pollen
    selectedbed=longctg_list_bed_files_selected717.txt
Input 2. minimum correlation to build contact graph of haplotig markers    <= best 1
    cor=0.55
Input 3. contig sizes
    ctgsize=/path/to/HiFiasm_ref_6366long_ctgs_selected.ctgsizes
Input 4. expected number of linkage groups
    nLG=48
Input 5. maximum correlation to find "homologous haplotype linkage groups" <= best -1
    ncor=-0.25
Input 6. minimum contig size to select initial haplotig markers to build the backbone of linkage groups. 15000 bp worked, but might cause consufison when integrating other contigs! so here we are 100 kb to build the initial linkage group of haplotigs
    min_hapctg_size=100000
Check purpose: re-calculate intra-LG inter-CTG correlation for all related markers!
    recalc="_recalc"

    gamete_binning_tetra ${selectedbed} ${cor} ${ctgsize} ${nLG} ${ncor} ${min_hapctg_size} gamete_binning_selected717_cor${cor}_ncor${ncor}_minHap${min_hapctg_size}_ncorminus${recalc} >gamete_binning_selected717_cor${cor}_ncor${ncor}_ncorminus_minHap${min_hapctg_size}bp${recalc}.log
    
    
    
