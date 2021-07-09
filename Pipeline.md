Pipeline
=
This is the pipeline explaining how gamete binning in tetraploid (potato) works.

##### Step 0. Prepare data

All data are from a tetraploid (potato cultivar) of interest, including PacBio HiFi reads from somatic tissue and short reads from
single-cell sequencing of sufficient (e.g., hundreds of) gamete genomes. In this example, we use the following:

* PacBio: long_reads_raw.fa
* 10x Genomics+Illumina: e.g., A_seq4414plus4431_R1.fastq.gz, A_seq4414plus4431_R2.fastq.gz (available on NCBI Bioproject: PRJNA726019)

Suppose all these raw data are collected the path below, and for convenience, make some softlinks (Note, 10x Genomics tools need the full name, so we use both namings),

    wd=/path/to/reads/
    cd ${wd}
    
    ln -s  A_seq4414plus4431_R1.fastq.gz gamete_libx_R1.fastq.gz
    ln -s  A_seq4414plus4431_R2.fastq.gz gamete_libx_R2.fastq.gz

##### Step 1. Trim reads

Trim 16 bp barcodes off R1's (10x Genomics library setting, including hexamer so 22 bp trimmed off),

    wd=/path/to/reads/
    cd ${wd}
    
    T10X_barcode_trimmer gamete_libx_R1.fastq.gz gamete_libx_R2.fastq.gz

This leads to

* gamete_libx_R1_clean.fastq.gz, gamete_libx_R2_clean.fastq.gz

##### Step 2. Genome size estimation with cleaned reads

Count k-mers (k=21),

    wd=/path/to/kmer_analysis/
    cd ${wd}
    zcat /path/to/reads/gamete_libx_R1_clean.fastq.gz /path/to/reads/gamete_libx_R2_clean.fastq.gz | jellyfish count /dev/fd/0  -C -o gamete_21mer_trimmed -m 21 -t 20 -s 5G
    jellyfish histo -h 200000 -o gamete_21mer_trimmed.histo gamete_21mer_trimmed

Estimate genome size,

    R
    library("findGSE")
    findGSE(histo="gamete_21mer_trimmed.histo", sizek=21, outdir=".", exp_hom=200)

In apricot, we got ~ 242.5 Mb.

##### Step 3. Preliminary assembly

    wd=/path/to/pre_assembly/
    cd ${wd}
    
    canu -p preasm -d canu_preasm useGrid=false genomeSize=242500000 corMhapSensitivity=high corMinCoverage=0 corOutCoverage=100 correctedErrorRate=0.105 -pacbio-raw long_reads_raw.fa executiveThreads=20 >canu_preasm.log

or,

    flye --pacbio-raw long_reads_raw.fa --genome-size 243m --out-dir flye_preasm --threads 4

This leads to preliminary assembly

* pre_asm.fasta

Polish the preliminary assembly using PacBio reads (using arrow/pilon), leading to the polished assembly

* pre_asm_pilon.fasta

Index it with bowtie2 for later read alignment (4 threads used)

    bowtie2-build -f pre_asm_pilon.fasta pre_asm_pilon.fasta --threads 4

Note, reference contig id might include "|", which is not accepted by many existing tools. Replace it with sed 's/|arrow|/_/g'.


##### Step 4. Curation of assembly (with purge_haplotigs pipeline)

Align Illumina reads to the preliminary assembly


    wd=/path/to/curated_asm/
    cd ${wd}

    bowtie2 -x /path/to/pre_assembly/pre_asm_pilon.fasta -1 /path/to/reads/gamete_libx_R1_clean.fastq.gz -2 /path/to/reads/gamete_libx_R2_clean.fastq.gz -p 20 | samtools view -@ 20 -bS - | samtools sort -@ 20 -o RP_PE_sorted.bam -

Generate a coverage histogram

    purge_haplotigs readhist -b RP_PE_sorted.bam -g /path/to/pre_assembly/pre_asm_pilon.fasta -t 4

Select coverage cutoffs (./tmp_purge_haplotigs/MISC/aligned_pe.bam.histogram.csv) - here we have het-peak at 84x, hom-peak at 171x:

    lowcutoff=40   -- this is the value at the valley on the left of het-peak
    midcutoff=121  -- this is the value at the valley between het- and hom-peaks
    highcutoff=342 -- this is ~2*hom-peak

Analyse the coverage on a contig by contig basis. This script produces a contig coverage stats csv file with suspect contigs flagged for further analysis or removal.

    purge_haplotigs contigcov -i RP_PE_sorted.bam.gencov -l ${lowcutoff} -m ${midcutoff} -h ${highcutoff} -o coverage_stats.csv -j 95 -s 80

Run a BEDTools windowed coverage analysis (note, this would lead to curated.fasta and curated.haplotigs.fasta)

    ABAM=RP_PE_sorted.bam
    genome=/path/to/pre_assembly/pre_asm_pilon.fasta
    purge_haplotigs purge -g ${genome} -c coverage_stats.csv -t 16 -o curated -d -b ${ABAM} -wind_min 1000 -wind_nmax 250 -v

Re-check haplotigs (blast them with the curated genome, i.e., the one built up with selected primary contigs). If a defined haplotig is not covered by more than 50% (by primary contigs), correct it as a primary contig, and merge it with curated.fasta. (Note, visualization of the blast result with R_scripts_aux/visualize_blast_haplotig_against_purged.R -- necessary settings on paths needed).

    db=curated.fasta
    makeblastdb -in ${db} -dbtype nucl > formatdb.log
    contigpath=/path/to/curated_asm/tmp_purge_haplotigs/CONTIGS/
    grep '>' curated.haplotigs.fasta | sed 's/ /\t/g' | cut -f1 | sed 's/>//' > haplotigs_as_predicted_by_purgeHaplotig.list
    while read ctg; do q=${ctg}.fasta; blastn -query ${contigpath}/${q} -db ${db} -out ${q}.oblast -outfmt 7 > blastall_blastn.log;done < ../haplotigs_as_predicted_by_purgeHaplotig.list

This leads to a version of manually curated assembly

* manually_curated.fasta
* manually_curated.chrsizes (this is contig size file of the assembly with format: contig_id	contig_size, tab-separated)

Note, although this is supposed to be a haploid assembly, it is built up with a high mixture of both haplotypes.

Index the sequence as reference for later steps,

    refgenome=manually_curated.fasta
    bowtie2-build -f ${refgenome} ${refgenome} --threads 4

##### Step 5. Read alignment of pooled gamete nuclei for SNP marker definition

    wd=/path/to/marker_creation/
    cd ${wd}

Align pooled gamete reads to the reference

    refgenome=/path/to/curated_asm/manually_curated.fasta
    bowtie2 -x ${refgenome} -1 /path/to/reads/gamete_libx_R1_clean.fastq.gz -2 /path/to/reads/gamete_libx_R2_clean.fastq.gz -p 20 | samtools view -@ 20 -bS - | samtools sort -@ 20 -o gamete_ManualCurated.bam -; samtools depth gamete_ManualCurated.bam > gamete_ManualCurated.depth.txt; samtools index gamete_ManualCurated.bam

Get vcf and variants (a file of ploidy need, columns are '\t'-separated)

    echo "*	*	*	*	2" > ploidy
    bcftools mpileup -d 800 -f ${refgenome} gamete_ManualCurated.bam |bcftools call -m -v -Ov --ploidy-file ${ploidy} > gamete_ManualCurated.vcf

Convert the variant format to plain text

    SHOREmap convert --marker gamete_ManualCurated.vcf --folder shoremap_converted --indel-size 5 --min-AF 0.1 -runid 20200426

This leads to 

    20200426_converted_variant.txt

Filtering for allelic snps. This needs to tune quality and coverage according to the specific data. For allele frequency, as we are looking for allelic snps, it should be around 0.5; while the coverage on the alt allele should be around half the genome wide average. In the meanwhile, the total coverage should not be too small or too large. In my case, the avg is around 171x (and half:84x), and we can try the cutoffs below:

Note, $6 is mapping quality; $7 is coverage of alt allele, we can try with

    awk '$6>=100 && $7>=60 && $7<=140 && $7/$8>=120 && $7/$8<=280 && $8>=0.38 && $8<=0.62' /path/to/20200426_converted_variant.txt > final_snp_markers.txt

##### Step 6. 10x Genomics barcode correction and nuclei separation

    wd=/path/to/individual_nuclei_extraction/
    cd ${wd}

We use an artifical reference from almond at chr-level (in other cases, one can select a closely-related species with chr-level assembly):

    wget http://getentry.ddbj.nig.ac.jp/getentry?database=na&accession_number=AP019297-AP019304&filetype=gz&limit=1000&format=fasta
    gunzip fasta_na.AP019297-AP019304.txt.gz
    sed 's/AP019297|AP019297\.1/chr1/g' fasta_na.AP019297-AP019304.txt| sed 's/AP019298|AP019298\.1/chr2/g' | sed 's/AP019299|AP019299\.1/chr3/g' | sed 's/AP019300|AP019300\.1/chr4/g' | sed 's/AP019301|AP019301\.1/chr5/g' | sed 's/AP019302|AP019302\.1/chr6/g' | sed 's/AP019303|AP019303\.1/chrX/g' | sed 's/AP019304|AP019304\.1/chrY/g' | sed 's/\./ /g' | sed 's/,//g'  | sed 's/ Prunus dulcis DNA pseudomolecule /_/g' > almond_genome.fa

Correspondingly, we prepare "a JSON file - /file_aux/contig_defs.json - describing primary contigs", and then index the genome with cellranger-dna,

    refgenome=almond_genome.fa
    cellranger-dna mkref ${refgenome} /path/to/file_aux/contig_defs.json

This will create a new reference folder of "/refdata-almond_genome/".

Correct 10x Genomics barcodes (note, if there are multiple libraries, this step needs to be done library by library, as same barcodes might be shared across libraries. However, different runs of the same library can be run together by setting option --sample=libx-run-1,libx-run-2\[,...\]),

    cellranger-dna cnv --id=4279_A_run615_cellranger --reference=/path/to/refdata-almond_genome/ --fastq=/path/to/reads/ --sample=4279_A_run615_SI-GA-D4 --localcores=20 --localmem=30

Sort the bam (from the above, which is with corrected barcode information) with read name

    bam=/path/to/4279_A_run615_cellranger/outs/possorted_bam.bam
    samtools sort -n ${bam} -o RNsorted_bam.bam

Create reads ligated with corrected barcode,

    bam=/path/to/4279_A_run615_cellranger/outs/RNsorted_bam.bam
    samtools view ${bam} | T10xbam2fq - 4279_A

This leads to

    4279_A_fqfrom10xBam_bxCorrected_R1.fq.gz
    4279_A_fqfrom10xBam_bxCorrected_R2.fq.gz

Extract individual nuclei

    R1=4279_A_fqfrom10xBam_bxCorrected_R1.fq.gz
    R2=4279_A_fqfrom10xBam_bxCorrected_R2.fq.gz
    barcode_len=16
    minimumRP=5000
    min_readpair=10000000
    asCellseparator ${barcode_len} ${R1} ${R2} ${minimumRP} ${min_readpair} cells_sep

This would report something like below,

    Info: R1 file: 4279_A_fqfrom10xBam_bxCorrected_R1.fq.gz
    Info: total number of barcodes observed in file1: 8656307
    Info: total number of reads    observed in file1: 124447727
    	among the above, good with read pairs >= 5000:
    	number of barcodes observed in file1: 401
    	number of reads	   observed in file1: 103317664
    Info: number of effective barcodes: 401

Get list of good barcodes,

    cd /path/to/cells_sep
    awk '$2>=5000' asCellseparator_intermediate_raw_barcode_stat.txt > barcode_over_5000rpairs.list

Merge parts of read extraction

    while read bc; do
    	cd /path/to/cells_sep/${bc}
    	cat *R1* > ${bc}_R1.fastq.gz
    	cat *R2* > ${bc}_R2.fastq.gz
    	cd ..
    done < barcode_over_5000rpairs.list

##### step 7. Read alignment and variant calling for each gamete nuclei

    wd=/path/to/individual_nuclei_read_align_var_calling/
    cd ${wd}

Prepare ploidy file,

    echo "*	*	*	*	1" > ploidy1

Align individual nuclei to manually curated genome (see Step 4.),

    refgenome=/path/to/curated_asm/manually_curated.fasta
    while read bc; do
        cd /path/to/individual_nuclei_extraction/cells_sep/${bc}
        bowtie2 -p 1 -x ${refgenome} -1 ${bc}_R1.fastq.gz -2 ${bc}_R2.fastq.gz 2> bowtie2.err | samtools view -@ 1 -bS - | samtools sort -@ 1 -o part1_${bc}.bam -
        bcftools mpileup -Oz -o bt2_${bc}_PE_mpileup.gz -f ${refgenome} part1_${bc}.bam
        bcftools call -A -m -Ov --ploidy-file /path/to/ploidy1 bt2_${bc}_PE_mpileup.gz > bt2_${bc}_PE.vcf
        SHOREmap convert --marker bt2_${bc}_PE.vcf --folder shoremap_converted -runid 20200426 -no-r >convert.log
        bcftools view -Ob bt2_${bc}_PE.vcf -o bt2_${bc}_PE.bcf
        rm bt2_${bc}_PE_mpileup.gz;
    done < barcode_over_5000rpairs.list

Note, at this step, there is possiblity to figure out which nuclei could be related to contamination according to alignment rate.

##### Step 8. Extract allele count (individual gamete nuclei) at SNP markers (from Step 5.)

    wd=/path/to/individual_nuclei_read_align_var_calling/
    cd ${wd}

    while read bc; do
        cd /path/to/individual_nuclei_extraction/cells_sep/${bc}
        SHOREmap extract --marker /path/to/marker_creation/final_snp_markers.txt --chrsizes /path/to/curated_asm/manually_curated.chrsizes --folder . --consen 20200426_converted_consen.txt
    done < barcode_over_5000rpairs.list

##### Step 9. Phasing SNPs within gamete genomes.

    wd=/path/to/snp_phasing/
    cd ${wd}

Prepare a meta-file of subset662_consen_cells.txt, where each line points to a consensus file of a specific barcode: /path/to/cells_sep/barcodeX/shoremap_converted/extracted_consensus_0.txt

    date=20200426
    marker=/path/to/marker_creation/final_snp_markers.txt
    cells=/path/to/subset662_consen_cells.txt
    sizes=/path/to/curated_asm/manually_curated.chrsizes
    asPollinator_1.0 --marker ${marker} --pollen ${cells} -o z${date}_phasing_with_correction_XXX_samples_full_markerSet_scorep81 --corr --ims 0.81 --size ${sizes} > z${date}_phasing_with_correction_XXX_samples_full_markerSet_scorep81.log

Collect the PM pattern of each nuclei at each contig for next step of haploid evaluation,

    ls z${date}_phasing_with_correction_XXX_samples_full_markerSet_scorep81_tmp_pollen_genotypes/s1_genotype_pollen_seq_ctgwise/s1_genotype_pollen_seq*.txt > pattern_nuclei_full_markerSet_list.txt

##### Step 10. haploidy level evaluation

    wd=/path/to/haploid_eval/
    cd ${wd}

Find potential transitions between two genotypes with the PM pattern of each nuclei at each contig from asPollinator_1.0,

    contigPM=/path/to/snp_phasing/pattern_nuclei_full_markerSet_list.txt
    noise_checker ${contigPM} > checkling_full_markerSert.log
    sort -k1,1n haplotype_swaps_observed_in_cells.txt > haplotype_swaps_observed_in_cells_sorted.txt
    rm haplotype_swaps_observed_in_cells.txt

The script "/path/to/R_scripts_aux/visualize_noise_stat_haploidy_evalu.R" can be used to visualize the haploidy evaluation after setting up paths on respective variables. The nuclei with <=5% genotype transitions are considered as haploid cells. With this, suppose we have an updated (and nuclei-reduced) barcode/nuclei list

    final_barcode_over_5000rpairs.list

Correspondingly, we have an updated list of consensus info of haploid nuclei to analyze,

    subset445_consen_cells.txt

##### Step 11. Contig grouping and genetic mapping using JoinMap4.0

    wd=/path/to/ctg_grp_gmap/
    cd ${wd}

This is done with JoinMap (interactive work with the software in windows), leading to

    genetic_map_limited_markers.txt --genetic map with 216 ordered markers, see here /path/to/GameteBinning/file_aux/GMaps/genetic_map_limited_markers.txt
    zphase_contigs_linksage.txt     --linkage groups with 653 non-ordered markers, see here /path/to/GameteBinning/file_aux/LGs/zphase_contigs_linksage.txt

Note, for different species, this would be different. We will later use the 216-ctg genetic map as a backbone (together with the phasing-of-contigs given by the 653-grouping) to insert un-mapped contigs, i.e., to get a more complete genetic map.

##### Step 12. Deletion marker definition (= select large regions without SNP markers)

    wd=/path/to/del_marker_definition/
    cd ${wd}

    marker=/path/to/marker_creation/final_snp_markers.txt
    sizes=/path/to/curated_asm/manually_curated.chrsizes
    mkdir define_del_like_regions
    cd define_del_like_regions
    del_marker_finder --marker ${marker} --min-del-size 2000 --chrsizes ${sizes} -o del_like > define_del_isize2000.log
    mv del_like_isize2000_final_snp_markers.txt del_like_isize2000_final_snp_markers.bed

##### Step 13. Genotype deletion markers with read counts. Note, -a of samtools depth must be on, to use coordiate as key to match regions by del_depth_finder

    wd=/path/to/del_marker_definition/
    cd ${wd}

Get depth for regions in bed,

    bam=/path/to/marker_creation/gamete_ManualCurated.bam
    samtools depth -a -b /path/to/del_like_isize2000_final_snp_markers.bed ${bam} > del_like_isize2000_final_snp_markers.bed.depth.illumina

Define hap- and hom-regions, according to (where 120x hap-hom cutoff needs to be selected based on data)

    region=/path/to/del_marker_definition/del_like_isize2000_final_snp_markers.bed
    depth=/path/to/del_marker_definition/del_like_isize2000_final_snp_markers.bed.depth.illumina
    del_depth_finder --region ${region} --depth ${depth} --hap-hom-cutoff 120 -o del_like_isize2000_final_snp_markers > del_depth_finder_illumina.log

This leads to

    del_like_isize2000_final_snp_markers_del_like_interval_avg_depth_sorted.bed

For each nuclei, count reads in the above hap- and hom- regions (if barcodes have been filtered during some of the above steps, pls replace it with the update one),

    bam=/path/to/marker_creation/gamete_ManualCurated.bam
    bed=/path/to/del_marker_definition/del_like_isize2000_final_snp_markers.bed
    while read bc; do
        cd /path/to/cells_sep/${bc}
        bedtools coverage -counts -a ${bed} -b part1_${bc}.bam -bed > ${bc}_del_like_read_count.bed
    done < /path/to/final_barcode_over_5000rpairs.list

Note, the read counts will be normalized (to RPKM) during del-marker phasing and positioning in genetic map.

Genotype deletion markers,

    date=2020426
    depth=/path/to/del_marker_definition/del_like_isize2000_final_snp_markers_del_like_interval_avg_depth_sorted.bed
    cells=/path/to/snp_phasing/subset445_consen_cells.txt
    del_marker_genotyper --pollen ${cells} --leaf-depth ${leaf_depth} --sample 8 --barcode 9 -o ${date} > ${date}_del_marker_genotyper.log

This leads to

    ${date}_tmp_pollen_del_like_genotypes_addi/s2_genotype_contig_seq_del_like.txt

##### Step 14. Genetic map completing

Complete the genetic map using del-like markers

    wd=/path/to/genetic_map_completing/
    cd ${wd}
    #
    map=/path/to/ctg_grp_gmap/genetic_map_limited_markers.txt - genetic map with 216 ordered markers    are generated by JoinMap, see here /path/to/GameteBinning/file_aux/GMaps/genetic_map_limited_markers.txt
    LG=/path/to/ctg_grp_gmap/zphase_contigs_linksage.txt      - linkage groups with 653 non-ordered markers are generated by JoinMap, see here /path/to/GameteBinning/file_aux/LGs/zphase_contigs_linksage.txt
    PM=/path/to/s2_genotype_contig_seq.txt                    - PM patterns are generated    by asPollinator_v6
    DEL=/path/to/s2_genotype_contig_seq_del_like.txt          - del markers are generated    by del_marker_genotyper
    asCaffolder_v2 --map ${map} --phase ${LG} --marker ${PM}  - marker-del-like ${DEL} -o phased > del_phased.log

This leads to,

    phased_s2_genotype_contig_seq_del_like.txt
    z_genetic_maps_updated_with_PMsimilarity_of_snp_plus_del_like_contigs/upd_map_group[1-8].txt

This folder includes more complete genetic maps by inserting del-like markers which will be used to anchor contigs of final assembly into chr-level.

Find out contigs upd_map_group[1-8].txt with left/right marker NOT next to each other; put them together at one position with minimum genetic distance to nearby markers:

    cat upd_map_group1.txt | cut -f1 | uniq -c | awk '$1==1'

Those showing 1 need moving.

This finally leads to,

    final_manual_upd_map_group[1-8].txt --see example result here: GameteBinning/file_aux/final_GMaps/

##### Step 15. Long read separation

    wd=/path/to/long_read_separation/
    cd ${wd}

Align PacBio reads to the manully curated assembly

    refgenome=/path/to/curated_asm/manually_curated.fasta
    minimap2 -ax map-pb -t 4 ${refgenome} /path/to/long_reads_raw.fa > pacbio_manual_purged_ref.sam

Separate longs with phased snps and dels in each linkage group,

    snps=/path/to/s4_phased_markers.txt ------------------------ phased snps by asPollinator_v6
    dels=/path/to/phased_s2_genotype_contig_seq_del_like.txt --- phased dels by asCaffolder_v2
    LG=/path/to/zphase_contigs_linksage.txt -------------------- grouping/phasing by JoinMap, see here /path/to/GameteBinning/file_aux/LGs/zphase_contigs_linksage.txt
    sam=/path/to/pacbio_manual_purged_ref.sam
    pacbio_genotyper --sam ${sam} --marker ${snps} --marker2 ${dels} --phase ${LG} --ims 0.9 -o pb_separation_with_dels > pb_intermediate_for_checkings.txt

In the end of pb_intermediate_for_checkings.txt (caution this is a large file), there would be a report on how the reads were seprated according to different cases (apricot-case as below),

    tail -n 43 rojopasionpb_intermediate_for_checkings.txt

We would see something like below:

    Warning: there are a1=390453 alignments, totaling v1=0.80742 Gb  without explicit CIAGR info -- collected in unmapped file.
    Warning: there are a2=949982 alignments being secondary/supplementary alignment, skipped.
    Warning: there are a3=1189596 alignments without seq field - secondary/supplementary alignment or not passing filters, skipped.
    Info: in total a4=2196665 reads from all=4726696 aligment lines collected (<-header line not counted).
    Info: number of pb alignment seqs WITHOUT linkage info: 57955, totaling v2=0.493197 Gb
          (u0=57955 unique reads in this cluster of no lg or not grouped. )
    Info: number of pb alignment seqs WITH    linkage info a5=2138710, totaling v3=18.6285 Gb
          (u1=2138710 unique reads in this cluster of lg_mkr+lg_nomkr: it is a sum of reads covering or not covering phased markers); among v3 (with a5),
          a6=801036 covered no phased markers thus cannot determine P/M cluster - alignment seq has been put in both P and M cluster,
          taking a portion of v4=6.24754 Gb (u2=801036 unique reads in this cluster of lg_nomkr)

    Note: u1 included all u2.
    Note: v1+v2+v3    = total raw pacbio data of (full: 19.929) Gb.
    Note: a1+a4       = total raw pacbio read number (full: 2,587,118)
    Note: a1+a2+a3+a4 = all raw alignment num;
          a4 only gives unique readname numbers; 1 readname may have >=2 alignments.)
    Note: u0+u1       = a4.

And, information on distribution of pacbio reads in linkage groups:

    1.txt_MMM_pbreads.fa    150560
    1.txt_PPP_pbreads.fa    158813
    2.txt_MMM_pbreads.fa    84496
    2.txt_PPP_pbreads.fa    86248
    3.txt_MMM_pbreads.fa    125911
    3.txt_PPP_pbreads.fa    122432
    4.txt_MMM_pbreads.fa    105348
    4.txt_PPP_pbreads.fa    114706
    5.txt_MMM_pbreads.fa    99987
    5.txt_PPP_pbreads.fa    102427
    6.txt_MMM_pbreads.fa    132686
    6.txt_PPP_pbreads.fa    138104
    7.txt_MMM_pbreads.fa    140458
    7.txt_PPP_pbreads.fa    148249
    8.txt_MMM_pbreads.fa    204317
    8.txt_PPP_pbreads.fa    223968

##### Step 16. Independent haplotype assemblies within each linkage group

    wd=/path/to/indep_haplotype_assembly/
    cd ${wd}

Here we assembled each haplotye for each linkage group, using flye

    for i in {1..8}; do
       for sample in ${chr}.txt_PPP_pbreads.fa ${chr}.txt_MMM_pbreads.fa; do
           flye --pacbio-raw /path/to/${sample} --genome-size 40m --out-dir flye_${sample} --threads 4
       done
    done

This leads to 16 haplotype-specific assemblies, each representing for one haploid genome in each of the eight linkage groups.

##### Step 17. Scaffolding haplotype-specific assemblies to chromosome-level

    wd=/path/to/indep_chr_scaffolding/
    cd ${wd}

Create a pseudo-reference chromosome for each linkage group with manually curated assembly (from Step 4) and the final complete genetic map (from Step 14: asCaffolder_v2),

    fa=/path/to/curated_asm/manually_curated.fasta
    updatedgm=/path/to/updated_genetic_map_folder/
    >scaffolder.log
    for i in {1..8}; do 
        scaffolder ${updatedgm}/final_manual_upd_map_group${i}.txt ${fa} >> scaffolder.log
    done

This leads to,

    scaffolded_final_manual_upd_map_group[1-8].fa

Scaffod the contigs into chromosome-level assemblies for each haplotype-assembly of each linkage group,

    scaf_dir=/path/to/do/scaffolding/
    cd ${scaf}
    #
    # Caution: only when linking fasta files (from Step 16) at this relative location can Rogoo find them!!!
    #
    for i in {1..8}; do ln -s /path/to/flye_${i}.txt_MMM_pbreads.fa/assembly.fasta ln_asm${i}.txt_MMM_pbreads.fa.contigs.fasta; done
    for i in {1..8}; do ln -s /path/to/flye_${i}.txt_PPP_pbreads.fa/assembly.fasta ln_asm${i}.txt_PPP_pbreads.fa.contigs.fasta; done
    #
    # CAUTION: genetic-map scaffolds - path setting below required by ragoo
    #
    # linkage group 1
    refmapfa=../scaffolded_final_manual_upd_map_group2.fa
    querylab=asm1
    # linkage group 2
    refmapfa=../scaffolded_final_manual_upd_map_group7.fa
    querylab=asm2
    # linkage group 3
    refmapfa=../scaffolded_final_manual_upd_map_group3.fa
    querylab=asm3
    # linkage group 4
    refmapfa=../scaffolded_final_manual_upd_map_group4.fa
    querylab=asm4
    # linkage group 5
    refmapfa=../scaffolded_final_manual_upd_map_group6.fa
    querylab=asm5
    # linkage group 6
    refmapfa=../scaffolded_final_manual_upd_map_group1.fa
    querylab=asm6
    # linkage group 7
    refmapfa=../scaffolded_final_manual_upd_map_group5.fa
    querylab=asm7
    # linkage group 8
    refmapfa=../scaffolded_final_manual_upd_map_group8.fa
    querylab=asm8
    #
    # for each pair of the above eight configurations, do ragoo as below:
    #
    for PM in PPP MMM; do
        mkdir ${querylab}_scaffolding_${PM}
        cd ${querylab}_scaffolding_${PM}
        queryfa=ln_${querylab}.txt_${PM}_pbreads.fa.contigs.fasta # path required by ragoo
        ragoo.py ../${queryfa} ${refmapfa} -C
        cd ..
    done

This leads to final chromosome-level and haplotype-specific assemblies (for 2 haplotypes x eight linkage groups)

    /path/to/do/scaffolding/asm[1-8]_scaffolding_[MMM/PPP]/ragoo_output/ragoo.fasta




