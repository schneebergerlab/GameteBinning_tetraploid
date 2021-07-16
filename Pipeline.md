Pipeline
=
This is the pipeline explaining how gamete binning in tetraploid (potato) works. Intermediate data (for running gamete_binning_tetra) are provided [here](https://mega.nz/folder/OhkkyDDQ#cb5j-u8BxyCIfxQ-rBAkfQ), if you want to have a test run.

##### Step 0. Prepare data

All data are from a tetraploid (potato cultivar) of interest, including PacBio HiFi reads from somatic tissue and short reads from
single-cell sequencing of sufficient gamete genomes; note that here we also included additional Illumina sequencing (of 10x single-molecule libraries) for contig depth analysis (to define contig markers more reliably). Single-molecule libraries can be replaced with any other normal Illumina short read sequencings in future applications. 

In this example pipeline, suppose all raw data (available on NCBI Bioproject PRJNA726019) are collected in the path below (where fastq file names might be indicatively different from downloaded ones),

    wd=/path/to/s0_reads/
    cd ${wd}

PacBio HiFi: 

    4396_A_CCS.fastq
    
10x Genomics+Illumina (sc: single-cell): 

Raw:

    4414_A_run633_SI-GA-A1_S1_L001_R1_001.fastq.gz
    4414_A_run633_SI-GA-A1_S1_L001_R2_001.fastq.gz
    4461_A_run636_SI-GA-A1_S3_L005_R1_001.fastq.gz
    4461_A_run636_SI-GA-A1_S3_L005_R2_001.fastq.gz
    4461_A_run636_SI-GA-A1_S3_L006_R1_001.fastq.gz
    4461_A_run636_SI-GA-A1_S3_L006_R2_001.fastq.gz
    4414_B_run633_SI-GA-B1_S2_L001_R1_001.fastq.gz
    4414_B_run633_SI-GA-B1_S2_L001_R2_001.fastq.gz
    4461_B_run636_SI-GA-B1_S2_L003_R1_001.fastq.gz
    4461_B_run636_SI-GA-B1_S2_L003_R2_001.fastq.gz
    4461_B_run636_SI-GA-B1_S2_L007_R1_001.fastq.gz
    4461_B_run636_SI-GA-B1_S2_L007_R2_001.fastq.gz


Combined:

    A_seq4414plus4461_R1.fastq.gz
    A_seq4414plus4461_R2.fastq.gz (from 4414_A_run633_SI-GA-A1,4461_A_run636_SI-GA-A1)
    B_seq4414plus4461_R1.fastq.gz
    B_seq4414plus4461_R2.fastq.gz (from 4414_B_run633_SI-GA-B1,4461_B_run636_SI-GA-B1)
    
10x Genomics+Illumina (sm: single-molecule): 

    C_seq2806_R1.fastq.gz
    C_seq2806_R2.fastq.gz (combined from all runs)

##### Step 1. Trim barcodes off short reads from 10x libraries - this is for pooled read alignment (for depth analysis and pollen genotyping)

Trim 16 bp barcodes off R1's (considering an additional hexamer, 22 bp would be trimmed off),

    wd=/path/to/s0_reads/
    cd ${wd}
    
    T10X_barcode_trimmer A_seq4414plus4461_R1.fastq.gz A_seq4414plus4461_R2.fastq.gz
    T10X_barcode_trimmer B_seq4414plus4461_R1.fastq.gz B_seq4414plus4461_R2.fastq.gz
    T10X_barcode_trimmer C_seq2806_R1.fastq.gz C_seq2806_R1.fastq.gz

This leads to

    trimmed_A_seq4414plus4461_R1.fastq.gz
    trimmed_A_seq4414plus4461_R2.fastq.gz 
    #
    trimmed_B_seq4414plus4461_R1.fastq.gz
    trimmed_B_seq4414plus4461_R2.fastq.gz
    #
    trimmed_C_seq2806_R1.fastq.gz
    trimmed_C_seq2806_R2.fastq.gz

##### Step 2. Preliminary assembly - markers will be defined from these contigs (after depth analysis)

    wd=/path/to/s2_pre_assembly/
    cd ${wd}
    
    otavahifi=/path/to/s0_reads/4396_A_CCS.fastq
    hifiasm -t 10 -o otava ${otavahifi} >hifiasm.log

This leads to preliminary assembly (we select the utg-level) by 

    cat otava.p_utg.gfa | grep '^S' | cut -f2,3 | awk '{print ">"$1"\n"$2}' > otava.p_utg.fasta

(Polish the preliminary assembly - not necessary, so details not given here)

##### Step 3. Curation of assembly using read alignment depth (or, purge redundant contigs representing the same genomic regions)    

This leads to a version of manually curated assembly (please refer to manuscript supplementary information: section "Initial tetraploid genome assembly, polishing and purging" for details)

    wd=/path/to/s3_curated_asm/
    cd ${wd}

We rename the purged assembly (with corresponding contig size information - two tab-separated columns: contig_id	contig_size) as below ([available here](https://mega.nz/folder/GktXEYCR#F3I8uTKvKO0Fu8VY2yc2WA)),
    
    HiFiasm_ref_6366long_ctgs_selected.fasta
    HiFiasm_ref_6366long_ctgs_selected.chrsizes

Note, it is a mixture of four haplotypes (with potentially collapsed homozygous regions between any more than one haplotype).

Index the sequence as reference for later steps,

    refgenome=HiFiasm_ref_6366long_ctgs_selected.fasta
    bowtie2-build -f ${refgenome} ${refgenome} --threads 4

##### Step 4. Read alignment of various sequencings for haplo/diplo/triplo/tetraplotig marker identification (note we used sc, sm and HiFi depth to collect more power to differentiate different types of contigs)

    wd=/path/to/s4_marker_creation/
    cd ${wd}

Align pooled gamete reads to the reference

    refgenome=/path/to/s3_curated_asm/HiFiasm_ref_6366long_ctgs_selected.fasta
    bowtie2 -x ${refgenome} -1 /path/to/s0_reads/trimmed_A_seq4414plus4461_R1.fastq.gz,/path/to/s0_reads/trimmed_B_seq4414plus4461_R1.fastq.gz -2 /path/to/s0_reads/trimmed_A_seq4414plus4461_R2.fastq.gz,/path/to/s0_reads/trimmed_B_seq4414plus4461_R2.fastq.gz -p 20 | samtools view -@ 20 -bS - | samtools sort -@ 20 -o gamete_ManualCurated.bam -
    
Remove duplicates and get position-wise depth

    java -jar picard.jar MarkDuplicates I=gamete_ManualCurated.bam O=gamete_ManualCurated_markeduplicates.bam M=gamete_ManualCurated_marked_dup_metrics.txt
    samtools index gamete_ManualCurated_markeduplicates.bam
    
Similarly, align sm related reads and get bam file 

    refgenome=/path/to/s3_curated_asm/HiFiasm_ref_6366long_ctgs_selected.fasta
    bowtie2 -x ${refgenome} -1 /path/to/s0_reads/C_seq2806_R1_clean.fastq.gz -2 /path/to/s0_reads/C_seq2806_R2_clean.fastq.gz -p 20 | samtools view -@ 20 -bS - | samtools sort -@ 20 -o sm_ManualCurated.bam -
    java -jar picard.jar MarkDuplicates I=sm_ManualCurated.bam O=sm_ManualCurated_markeduplicates.bam M=sm_ManualCurated_marked_dup_metrics.txt
    samtools index sm_ManualCurated_markeduplicates.bam
    
Align HiFi reads, remove non-primary alignments and get position-wise depth 
    
    otavahifi=/path/to/s0_reads/4396_A_CCS.fastq
    refgenome=/path/to/s3_curated_asm/HiFiasm_ref_6366long_ctgs_selected.fasta
    minimap2 -ax map-pb -t 20 -N 1 --secondary=no ${refgenome} ${otavahifi} | samtools view -@ 20 -bS - | samtools sort -@ 20 -o HiFi_ManualCurated.bam -
    samtools view -h -F 3840 -bS HiFi_ManualCurated.bam | samtools sort -o HiFi_ManualCurated_clean.bam - 
    samtools depth -Q 1 -a HiFi_ManualCurated_clean.bam > HiFi_ManualCurated_depth_clean.txt

Get position-wise sequencing depth 

    samtools depth -Q 1 -a gamete_ManualCurated.bam sm_ManualCurated.bam HiFi_ManualCurated_clean.bam > otava_Hifiasm_ref_pb_illu_depth.txt
    awk '{printf ("%s\t%s\t%s\n", $1, $2, $3+$4+$5)}' otava_Hifiasm_ref_pb_illu_depth.txt > otava_Hifiasm_ref_pb_illu_depth_sum.txt
    rm otava_Hifiasm_ref_pb_illu_depth.txt
    
##### step 5: find distribution of average depth at non-overlapping windows: winstep = winsize

    wd=/path/to/s4_marker_creation/
    cd ${wd}
    chrsizes=/path/to/s3_curated_asm/HiFiasm_ref_6366long_ctgs_selected.ctgsizes
    winsize=10000
    
    # => cnv_winsize10000_step10000_hq.txt
    samplecol=1
    avgdepth=113 # 137*0.85
    CNV_HQ_v3 ${chrsizes} otava_Hifiasm_ref_pb_illu_depth_sum.txt ${winsize} ${winsize} ${samplecol} ${avgdepth}
    
    # => cnv_winsize10000_step10000_hq_HiFi_only.txt
    samplecol=1
    avgdepth=30 # HiFi only
    CNV_HQ_v3 ${chrsizes} HiFi_ManualCurated_depth_clean.txt ${winsize} ${winsize} ${samplecol} ${avgdepth}
    
    
##### step 6: new window marker generation => cnv_winsize10000_step10000_hq_markers_20210712_wsize50kb_final.txt

    wd=/path/to/s4_marker_creation/
    cd ${wd}
    #
    paste cnv_winsize10000_step10000_hq.txt cnv_winsize10000_step10000_hq_HiFi_only.txt | cut -f1-7,12,13 > cnv_winsize10000_step10000_hq_merged_vs_hifi.txt
    tig_marker_finder cnv_winsize10000_step10000_hq_merged_vs_hifi.txt 113 30 50000 20210712 > tig_marker_finder.log

##### Step 7. 10x Genomics barcode correction and nuclei separation

    wd=/path/to/s7_individual_nuclei_extraction/
    cd ${wd}

We use DM assembly at chr-level - potato_dm_v404_all_pm_un.fasta - as reference (download: http://solanaceae.plantbiology.msu.edu/pgsc_download.shtml), and prepare a chrsizes file - potato_dm_v404_all_pm_un_modified.chrsizes (two columns, tab-separated with chr_id chr_size).

Rename some chromsome IDs (naming required by cellranger),

    sed -i 's/chr00/chrX/' potato_dm_v404_all_pm_un.fasta
    sed -i 's/ChrUn/chrY/' potato_dm_v404_all_pm_un.fasta

Correspondingly, we prepare "a JSON file - available under /file_aux/: contig_defs.json - describing primary contigs" (with below), 

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

This will create a new reference folder under current folder.

    refdata-potato_dm_v404_all_pm_un

Correct 10x Genomics barcodes (note, if there are multiple libraries, this step needs to be done library by library, as same barcodes might be shared across libraries. However, different runs of the same library can be run together by setting option --sample=libx-run-1,libx-run-2\[,...\]),

    cellranger-dna cnv --id=A_all3runs --reference=/path/to/refdata-potato_dm_v404_all_pm_un/ --fastq=/path/to/s0_reads/ --sample=4414_A_run633_SI-GA-A1,4461_A_run636_SI-GA-A1 --localcores=40 --localmem=30

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
    mkdir sample_${sample}_asCellseparator_40krp
    cd sample_${sample}_asCellseparator_40krp
    asCellseparator ${barcode_len} ${R1} ${R2} ${minimumRP} ${min_readpair} ./
    
Barcode-specific fastq file(s) will be collected under 

    sample_${sample}_asCellseparator_40krp/XXXbarcodeXXX/

Get list of good barcodes,

    awk '$2>=40000' asCellseparator_intermediate_raw_barcode_stat.txt > barcode_over_40000rpairs.list
   
Similarly, extract cell-wise sequencings from library B.

##### step 8. read alignment for each pollen genome sequencing

    wd=/path/to/s7_individual_nuclei_extraction/
    cd ${wd}
    
    cellpath=/path/to/s7_individual_nuclei_extraction/
    cd ${cellpath}/
    refgenome=/path/to/s3_curated_asm/HiFiasm_ref_6366long_ctgs_selected.fasta
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
            bowtie2 -p 1 -x ${refgenome} -1 ${R1} -2 ${R2} 2> longctg_bowtie2.err | samtools view -@ 1 -bS - | samtools sort -@ 1 -o longctg_${bc}.bam -
            samtools index longctg_${bc}.bam
            java -jar picard.jar MarkDuplicates I=longctg_${bc}.bam O=longctg_${bc}_markeduplicates.bam M=longctg_${bc}_marked_dup_metrics.txt
            samtools index longctg_${bc}_markeduplicates.bam
        done < ${sample}_this_barcode_list
        cd ${cellpath}
    done
    
##### step 9. genotype each pollen at each contig marker

Note 1: "-F 3840" == "-F 256 -F 512 -F 1024 -F 2048", excluding all kinds of non-primary alignment!

Note 2: MQ=5 is too stringent that it was observed some pollen lost coverage at haplotig markers (where there were primary-reads with MQ=1); need to use MQ1: found with IGV.

    wd=/path/to/s7_individual_nuclei_extraction/
    cd ${wd}
    
    cut -f1-3 /path/to/s4_marker_creation/cnv_winsize10000_step10000_hq_markers_20210712_wsize50kb_final.txt > cnv_winsize10000_step10000_hq_markers_20210712_wsize50kb_final.bed
    MQ=1
    cellpath=/path/to/s7_individual_nuclei_extraction/
    for sample in A B; do
        while read r; do 
            cd ${cellpath}/sample_${sample}_asCellseparator_40krp/${r}
            samtools view -h -F 3840 -q ${MQ} -bS longctg_${r}_markeduplicates.bam | samtools sort -o longctg_${r}_markeduplicates_MQ${MQ}_clean.bam -
            samtools index longctg_${r}_markeduplicates_MQ${MQ}_clean.bam
            bedtools coverage -counts -a /path/to/s7_individual_nuclei_extraction/cnv_winsize10000_step10000_hq_markers_20210712_wsize50kb_final.bed -b longctg_${r}_markeduplicates_MQ${MQ}_clean.bam -bed > longctg_${r}_win_marker_read_count_MQ${MQ}.bed            
        done < ${cellpath}/sample_${sample}_asCellseparator_40krp/${sample}_this_barcode_list
    done
    
##### step 10 add more info on markers to read cnt file of each cell 

    cellpath=/path/to/s7_individual_nuclei_extraction/
    cd ${cellpath}
    marker=/path/to/s4_marker_creation/cnv_winsize10000_step10000_hq_markers_20210712_wsize50kb_final.txt
    MQ=1
    for sample in A B; do
        while read r; do 
            echo -e ${sample}"\t"${r}
            cd ${cellpath}/sample_${sample}_asCellseparator_40krp/${r}
            paste longctg_${r}_win_marker_read_count_MQ${MQ}.bed ${marker} > longctg_${r}_win_marker_read_count_MQ${MQ}_updated.bed
        done < ${cellpath}/sample_${sample}_asCellseparator_40krp/${sample}_this_barcode_list
    done
    
##### step 11 prepare nuclei depth data

Note, we finally selected 717 nuclei to perform gamete binning (for details, please check Supplementary information: Section "Selection of single-cell sequencings"), for which barcodes are available under: "/aux_data/", i.e., longctg_s4p3_selected_717_good_nucei_lib[A|B].txt

    wd=/path/to/s11_selected_long_contigs_sc_read_coverage_genotype_v2/
    cd ${wd}
    ls /path/to/s7_individual_nuclei_extraction/sample_*_asCellseparator_40krp/*/longctg_*_win_marker_read_count_MQ1_updated.bed > longctg_list_bed_files.txt
    #
    >longctg_list_bed_files_selected717.txt
    for sample in A B; do 
        while read bc; do 
             grep ${bc} longctg_list_bed_files.txt >> longctg_list_bed_files_selected717.txt
        done < ../s4_selected_long_contigs_sc_read_coverage_calculation/longctg_s4p3_selected_717_good_nucei_lib${sample}.txt
    done
    
##### step 12 CORE FUNCTION of Gamete Binning: find linkage groups and match "homologous" linkage groups
 
    wd=/path/to/s11_selected_long_contigs_sc_read_coverage_genotype_v2/
    cd ${wd} 
 
Input 1. bed files on read counts for each pollen

    selectedbed=longctg_list_bed_files_selected717.txt

Input 2. minimum correlation to build contact graph of haplotig markers    <= best 1

    cor=0.55

Input 3. contig sizes

    ctgsize=/path/to/s3_curated_asm/HiFiasm_ref_6366long_ctgs_selected.ctgsizes

Input 4. expected number of linkage groups

    nLG=48

Input 5. maximum correlation to find "homologous haplotype linkage groups" <= best -1

    ncor=-0.25

Input 6. minimum contig size to select initial haplotig markers to build the backbone of linkage groups. 15000 bp worked, but might cause consufison when integrating other contigs! so here we are 100 kb to build the initial linkage group of haplotigs

    min_hapctg_size=100000

Check purpose: re-calculate intra-LG inter-CTG correlation for all related markers

    recalc="_recalc"

Run:

    gamete_binning_tetra ${selectedbed} ${cor} ${ctgsize} ${nLG} ${ncor} ${min_hapctg_size} gamete_binning_selected717_cor${cor}_ncor${ncor}_minHap${min_hapctg_size}_ncorminus${recalc} >gamete_binning_selected717_cor${cor}_ncor${ncor}_ncorminus_minHap${min_hapctg_size}bp${recalc}.log
    
##### step 13 some checkings

    wd=/path/to/s11_selected_long_contigs_sc_read_coverage_genotype_v2/
    cd ${wd}
    min_hapctg_size=100000
    cd s4_gamete_binning_selected717_cor0.55_ncor-0.25_minHap100000_ncorminus_recalc_tmp_integrating_all_ctg_markers_to_LGs/
    sort -k1,1 -k2,2n s4p6_refine_grouping_final_window_markers.txt > s4p6_refine_grouping_final_window_markers_sorted.txt
    # check LG-wise marker sizes 
    >LG_total_sizes_final.txt
    for i in {1..48}; do 
        echo -n -e "LG"$i"\t" >>LG_total_sizes_final.txt
        awk -v var="$i" '$5==var' s4p6_refine_grouping_final_window_markers_sorted.txt | awk '{s+=$3-$2+1} END {print s}' >> LG_total_sizes_final.txt
    done 
    # check inferred total assembly size 
    awk '{s+=$2} END {print s}' LG_total_sizes_final.txt
    # 
    
##### step 14 group HiFi reads to 1-48 linkage groups according to marker phasing/grouping at step 12-13.

    wd=/path/to/s14_HiFi_separation/
    cd ${wd}
    bam=/path/to/s4_marker_creation/HiFi_ManualCurated_clean.bam
    marker=/path/to/s11_selected_long_contigs_sc_read_coverage_genotype_v2/s4_gamete_binning_selected717_cor0.55_ncor-0.25_minHap100000_ncorminus_recalc_tmp_integrating_all_ctg_markers_to_LGs/s4p6_refine_grouping_final_window_markers_sorted.txt
    samtools view ${bam} | long_read_separator - ${marker} hifi_separation_20210712 > hifi_separation.log

check how reads are grouped

    tail -n 62 hifi_separation.log

you would see something like below:
    
    Warning: there are a1=5968 alignments, totaling v1=0.0697008 Gb  without explicit CIAGR info -- collected in unmapped file.
    Warning: there are a2=109882 alignments being secondary/supplementary alignment, skipped. 
    Info: in total a4=7037721 reads from all=7153571 aligment lines collected (<-header line not counted; raw alignment including secondary etc - skipped in read sep). 
    Info: number of pb alignment seqs WITHOUT linkage Info: 681266, totaling v2=9.90746 Gb
          (u0=681266 unique reads in this cluster of no lg or not grouped. )
 
    Info: distribution of pacbio reads in linkage groups: 
 
 	homLG_10_LG_2_reads.fa	132631	132631
	homLG_10_LG_35_reads.fa	113291	113291
	homLG_10_LG_3_reads.fa	134721	134721
	homLG_10_LG_47_reads.fa	111901	111901
	homLG_11_LG_10_reads.fa	106071	106071
	homLG_11_LG_30_reads.fa	95197	95197
	homLG_11_LG_32_reads.fa	101950	101950
	homLG_11_LG_46_reads.fa	104772	104772
	homLG_12_LG_1_reads.fa	163015	163015
	homLG_12_LG_33_reads.fa	159507	159507
	homLG_12_LG_34_reads.fa	133039	133039
	homLG_12_LG_44_reads.fa	133466	133466
	homLG_1_LG_11_reads.fa	90819	90819
	homLG_1_LG_42_reads.fa	89757	89757
	homLG_1_LG_45_reads.fa	83221	83221
	homLG_1_LG_5_reads.fa	75891	75891
	homLG_2_LG_12_reads.fa	113501	113501
	homLG_2_LG_27_reads.fa	101865	101865
	homLG_2_LG_31_reads.fa	114803	114803
	homLG_2_LG_8_reads.fa	114322	114322
	homLG_3_LG_13_reads.fa	149412	149412
	homLG_3_LG_14_reads.fa	137305	137305
	homLG_3_LG_18_reads.fa	142087	142087
	homLG_3_LG_21_reads.fa	126934	126934
	homLG_4_LG_15_reads.fa	121480	121480
	homLG_4_LG_17_reads.fa	110504	110504
	homLG_4_LG_23_reads.fa	118533	118533
	homLG_4_LG_36_reads.fa	120575	120575
	homLG_5_LG_19_reads.fa	162159	162159
	homLG_5_LG_29_reads.fa	153538	153538
	homLG_5_LG_43_reads.fa	132764	132764
	homLG_5_LG_9_reads.fa	153981	153981
	homLG_6_LG_20_reads.fa	214290	214290
	homLG_6_LG_26_reads.fa	182949	182949
	homLG_6_LG_48_reads.fa	200294	200294
	homLG_6_LG_6_reads.fa	208353	208353
	homLG_7_LG_22_reads.fa	141740	141740
	homLG_7_LG_24_reads.fa	147419	147419
	homLG_7_LG_39_reads.fa	121814	121814
	homLG_7_LG_7_reads.fa	135765	135765
	homLG_8_LG_25_reads.fa	107689	107689
	homLG_8_LG_37_reads.fa	127179	127179
	homLG_8_LG_38_reads.fa	107652	107652
	homLG_8_LG_4_reads.fa	126167	126167
	homLG_9_LG_16_reads.fa	168754	168754
	homLG_9_LG_28_reads.fa	168063	168063
	homLG_9_LG_40_reads.fa	148658	148658
	homLG_9_LG_41_reads.fa	146657	146657
	mapped_no_lg_pb_reads.fa	681266	681266
	unmapped_starcigar_pb_reads.fa	5968	5968
    Info: extract reads from bam/sam into linkage done. 
    Time consumed: 13030.3 seconds

##### step 15 LG-wise assembly 

    wd=/path/to/s15_asm_version_20210712/
    cd ${wd} 
    tail -n 62 /path/to/s14_HiFi_separation/hifi_separation.log | grep 'homLG' | sed 's/\t//g' | sed 's/_reads/\t/g' | cut -f1 > > fa_to_run.list
    
The file "fa_to_run.list" should contain haplotype chromsome ids at each line: homLG_10_LG_2, homLG_10_LG_35..., etc.
    
    while read fa; do
        cd ${wd}
        mkdir hifiasm_${fa}
        cd hifiasm_${fa}        
        otavahifi=/path/to/s14_HiFi_separation/hifi_separation_20210712_window_marker_separated_reads/${fa}_reads.fa
        ll ${otavahifi}
        # assemble haplotype: hifiasm v0.7
        hifiasm -t 10 -o ${fa} ${otavahifi} > otava_hifiasm.log
        # extract sequence
        cat ${fa}.p_ctg.gfa | grep '^S' | cut -f2,3 | awk '{print ">"$1"\n"$2}' > ${fa}.p_ctg.fasta
        # calculate N50 etc, perl: v5.28.1
        perl /path/to/src_calc_N50/calc_CN50.pl ${fa}.p_ctg.fasta 100000000 1 > ${fa}.p_ctg.N50.calc.result.txt
        cd ..
    done < ../fa_to_run.list
    
    
    
    
    
    
    
    
    
    
    
    
    
