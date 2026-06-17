# Running metilene to get differentially methylated regions

# Why am I doing this? The epidiverse/dmrs pipeline failed to detect parameters X and Y, so I had to run manually the metilene step

# I am using the bedtools output from the epidiverse/dmrs pipeline 

# to get the metilene binary for Linux
wget http://www.bioinf.uni-leipzig.de/Software/metilene/metilene_v02-9.tar.gz
# I had to download it on aa Windows and upload to EBD:Genomics-b via FileZilla. I had to run `make` in EBD:Genomics-b to compile the binary

# Generate an alias for the binary with path (Windows path)
metilene="/mnt/d/Ruben_backup/Postdoc/2022_AlonsoLabPostdoc/Analyses/EcicEcaz_2022/Ecic_WGBS_herb/ecic_herbivory_test/dmr/dmrs/metilene_v02-9/metilene_v0.2-9/metilene"

# To test that metilene runs, just type
$metilene

# Now, for a basic DMR analysis:
#$metilene -a C -b H input.bed > metilene_output.txt
# where -a is group 1 ID (C as in UN) and -b is the group 2 ID (H as in MH)
# output columns are: chr start stop q-value mean_difference_mean_g1-mean_g2 #CpGs p_(MWU) p_(2D_KS) mean_g1 mean_g2

# Interesting options
# option -M: sliding window length (default 300)
# option -m: minimum DMC numbers within sliding window to be considered a DMR (defoult 10)
# option -d: minimum mean methylation differences to consider a DMR (default 0.1)
# option -X/-Y: number of replicates per group a/b that have data to take into account for DMRs (default -1)
# option -t: multithreading
# option -v: valley filter, can split a DMR if there are empty regions between their DMCs, the higher, the more likely it splits the DMR (default 0.7)
# option -B/-f 2: add annotation in BED format (generate promoter, and gene entries)

# A more customized command:
#$metilene -a C -b H -d 0.2 -X 4 -Y 4 -v 0.5 input.bed > metilene_output.txt
# which: compares C vs H, minimum difference per DMR is |0.2|, min number of reps with data per group is 4, stringency for valley filter is 0.5

# In a loop:
projectDir="/mnt/d/Ruben_backup/Postdoc/2022_AlonsoLabPostdoc/Analyses/EcicEcaz_2022/Ecic_WGBS_herb/ecic_herbivory_test/dmr/dmrs"
cd $projectDir/
# detect folders from all methylation contexts
for dir in $(ls -d C[pH][GH]); do
  cd $dir
  for i in $(ls); do
    # Get first character of the folder (grouup 1)
    a=`echo $i | cut -c1-1`
    # Get last character of the folder (group 2)
    b=`echo "${i: -1}"`
    # Get into the folder of interest (group 1 vs group 2)
    cd $a\_vs_$b
    # run metilene with customized groups and generate output with customized name
    $metilene -a $a -b $b -d 0.2 -X 4 -Y 4 -v 0.5 input.bed > metilene_$dir.$a\_vs_$b.txt
    # options from epidiverse/dmr: -M 146 -c 2 -m 10 -d 0.1 -t 8
    cd ..
  done
  cd ..
done

# To annotate the detected DMRs:
# Generate bed from E. cicutarium genome GFF file (metilene generates BED files, so it is better to avoid comparing different formats)
genomeGFF="/mnt/d/Ruben_backup/Postdoc/2022_AlonsoLabPostdoc/Data/genomes/Ecic_genome_pacBio/Gene_prediction_files/AUGUSTUS/Cicuseq2_augustus_v2.gff"
# Generate a BED file only with gene coordinates
grep 'ID=gene_id' $genomeGFF | gff2bed < - > ecic_genome_gene.bed
# To generate 5' promoter regions linked to the genes
# Prepare a tab delimited genome file with "chromosome/contig"\t"length"
bioawk -c fastx '{ print $name, length($seq) }' < /mnt/d/Ruben_backup/Postdoc/2022_AlonsoLabPostdoc/Data/genomes/Ecic_genome_pacBio/MaSuRCA_assembly.fasta.PolcaCorrected.linear.fa > genomeLength.txt
# run bedtools flank with option -l 2000 or -r 2000 to generate 2 Kbp flanking regions before TSS, strand sensitive (if + strand, 2k left; if - strand, 2k right)
grep '+' ecic_genome_gene.bed | bedtools flank -i - -g genomeLength.txt -l 2000 -r 0 > ecic_genome_gene_2kbp_plus.bed
grep -v '+' ecic_genome_gene.bed | bedtools flank -i - -g genomeLength.txt -l 0 -r 2000 > ecic_genome_gene_2kbp_minus.bed
cat ecic_genome_gene_2kbp_*.bed > ecic_genome_gene_2kbp_upstreamTSS.bed
# Generate BED file from RepeatMasker output (filter with grep the categories 'Unknown repeat', 'Simple repeat', and 'Low complexity'):
tail -n +4 /mnt/d/Ruben_backup/Postdoc/2022_AlonsoLabPostdoc/Data/genomes/Ecic_genome_pacBio/Repeats/CicuSeq2_assembly.fasta.noUnknownSimpleRepeatLowComplexity.out | \
  awk 'BEGIN{OFS="\t"} {print $5, $6-1, $7, $10, $11, ".", ($9=="C" ? "-" : "+")}' > TEs.bed
# Asign the produced BED files to variables
genomeBED="/mnt/d/Ruben_backup/Postdoc/2022_AlonsoLabPostdoc/Analyses/EcicEcaz_2022/Ecic_WGBS_herb/ecic_herbivory_test/dmrs/dmrs/ecic_genome_gene.bed"
genomeFLANK="/mnt/d/Ruben_backup/Postdoc/2022_AlonsoLabPostdoc/Analyses/EcicEcaz_2022/Ecic_WGBS_herb/ecic_herbivory_test/dmrs/dmrs/ecic_genome_gene_2kbp_upstreamTSS.bed"
genomeTEs="/mnt/d/Ruben_backup/Postdoc/2022_AlonsoLabPostdoc/Data/genomes/Ecic_genome_pacBio/Repeats/TEs.bed"
# Filter metilene output files by FDR and extract their overlapping gene IDs:
for dir in $(ls -d C[pH][GH]); do
  cd $dir
  for i in $(ls); do
    a=`echo $i | cut -c1-1`
    b=`echo "${i: -1}"`
    cd $a\_vs_$b
    # If the 4th column (FDR) is less than 0.05, print the line
    awk '{ if ($4 <= 0.05) print $0 }' < metilene_$dir.$a\_vs_$b.txt > metilene_$dir.$a\_vs_$b\_FDRfilter.txt
    # Obtain overlaps within genes 
    bedtools intersect -a $genomeBED -b metilene_$dir.$a\_vs_$b\_FDRfilter.txt > metilene_$dir.$a\_vs_$b\_genes.bed
    # Obtain overlaps in 2Kbp flanking region
    bedtools intersect -a $genomeFLANK -b metilene_$dir.$a\_vs_$b\_FDRfilter.txt > metilene_$dir.$a\_vs_$b\_genes_2kbp.bed
    # Obtain overlaps in TE regions
    bedtools intersect -a $genomeTEs -b metilene_$dir.$a\_vs_$b\_FDRfilter.txt > metilene_$dir.$a\_vs_$b\_TEs.bed
    # Obtain intersections between TEs and genes and their flanking regions, 
    cat metilene_$dir.$a\_vs_$b\_genes.bed metilene_$dir.$a\_vs_$b\_genes_2kbp.bed | \
    bedtools intersect -a - -b metilene_$dir.$a\_vs_$b\_TEs.bed > metilene_$dir.$a\_vs_$b\_TEs_overlap.bed
    cd ..
  done
  cd ..
done

# To generate a few reports:
wc -l C*/*/*FDR* > number_DMRs_comparison_report.txt # number of detected DMRs per context and comparison 
wc -l C*/*/*genes.bed > overlap_DMR_genes_report.txt # number of genes overlapping with a DMR
wc -l C*/*/*2kb.bed > overlap_DMR_genes_2kb_report.txt # number of flanking regions ("promoters") overlapping with a DMR
wc -l C*/*/*TEs.bed > overlap_DMR_TEs_report.txt # number of TEs overlapping with a DMR
wc -l C*/*/*TEs_overlap.bed > overlap_DMR_TEs_promoters_report.txt # number of TEs overlapping flanking regions ("promoters") associated with a gene with a DMR  

# Third: Get GO terms and do GO term enrichment analysis
# obtain FDR from ovelapping genes (made in excel as temporary solution)
# run topGO with script TFM_topGO_analysis.R
