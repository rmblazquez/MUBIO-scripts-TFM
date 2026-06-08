#!/usr/bin/env bash

#### SETTING THE PROJECT ####
# Define the path of the project
projectDir="/home/usuario/Windows/Analyses/EcicEcaz_2022/Ecic_rnaseq_herbivory"
cd $projectDir
# Generate all subfolders
mkdir -p $projectDir/src $projectDir/data/genome $projectDir/data/raw-seq \
$projectDir/results/fastqc-rawseq $projectDir/results/fastqc-trimseq  \
$projectDir/results/htseq-count  $projectDir/results/star \
$projectDir/results/trimmomatic $projectDir/results/featureCounts \
$projectDir/results/salmon

# Generate a list of file names from the project
ls -d $projectDir/data/raw-seq/*.fastq.gz | sed -e 's/\_[12].fastq.gz//g' -e "s#$projectDir/data/raw-seq/##g" | sort | uniq > $projectDir/src/ecic_basenames.txt

#### RUN QUALITY CHECK IN RAW SEQUENCES WITH FASTQC ####
# Run FastQC for all the project files
cat $projectDir/src/ecic_basenames.txt | while read line; do
	mkdir $projectDir/results/fastqc-rawseq/$line\_1.fastqc
	fastqc -o $projectDir/results/fastqc-rawseq/$line\_1.fastqc \
	$projectDir/data/raw-seq/$line\_1.fastq.gz;
	mkdir $projectDir/results/fastqc-rawseq/$line\_2.fastqc
	fastqc -o $projectDir/results/fastqc-rawseq/$line\_2.fastqc \
	$projectDir/data/raw-seq/$line\_2.fastq.gz;
done

# Generate a HTML quality chek report for raw reads with MultiQC
cd $projectDir/results/fastqc-rawseq
multiqc $projectDir/results/fastqc-rawseq/*/*.zip
cd $projectDir

#### FILTER READS AND NUCLEOTIDES BY QUALITY WITH TRIMMOMATIC ####
# Define a path for fasta file with Illumina adapters 
# Different adapters are used depending on type of sequencing machine
trimm="/home/usuario/Windows/Analyses/EcicEcaz_2022/Ecic_rnaseq_herbivory/src/adapters.fa" # Illumina TruSeq adapters
# Run Trimmomatic for all the fastq files
cat $projectDir/src/ecic_basenames.txt | while read line; do
  TrimmomaticPE \
  -phred33 \
  $projectDir/data/raw-seq/$line\_1.fastq.gz \
  $projectDir/data/raw-seq/$line\_2.fastq.gz \
  $projectDir/results/trimmomatic/$line\_1.trimmed.paired.fastq.gz \
  $projectDir/results/trimmomatic/$line\_1.trimmed.unpaired.fastq.gz \
  $projectDir/results/trimmomatic/$line\_2.trimmed.paired.fastq.gz \
  $projectDir/results/trimmomatic/$line\_2.trimmed.unpaired.fastq.gz \
  ILLUMINACLIP:$trimm:2:30:10 \
  LEADING:28 TRAILING:28 \
  SLIDINGWINDOW:4:15 MINLEN:30;
  mkdir $projectDir/results/fastqc-trimseq/$line\_1.trimmed.paired.fastqc;
  fastqc -o $projectDir/results/fastqc-trimseq/$line\_1.trimmed.paired.fastqc \
  $projectDir/results/trimmomatic/$line\_1.trimmed.paired.fastq.gz;
  mkdir $projectDir/results/fastqc-trimseq/$line\_2.trimmed.paired.fastqc;
  fastqc -o $projectDir/results/fastqc-trimseq/$line\_2.trimmed.paired.fastqc \
  $projectDir/results/trimmomatic/$line\_2.trimmed.paired.fastq.gz;
done

# Generate a HTML quality chek report for trimmed reads with MultiQC
cd $projectDir/results/fastqc-trimseq
multiqc $projectDir/results/fastqc-trimseq/*/*.zip
cd $projectDir

#### EDIT GENOME FILES TO FIT THE PIPELINE ####
ln -s /home/usuario/Windows/Data/genomes/Ecic_genome_pacBio/MaSuRCA_assembly.fasta.PolcaCorrected.fa $projectDir/data/genome/ecic_genome.fa
ln -s /home/usuario/Windows/Data/genomes/Ecic_genome_pacBio/Gene_prediction_files/AUGUSTUS/CicuSeq2_augustus.gff3 $projectDir/data/genome/ecic_annotation.gff3
grep -v '#' ecic_annotation.gff3 > ecic_annotation.tmp.gff3

#### CONVERT GFF FILES INTO GTF FILES ####
# Use the Rscript gff2gtf_rtracklayer.R to convert the GFF into GTF
# NEED TO USE R version >= 3.5.0
/usr/local/bin/./R --slave --no-restore --file=$projectDir/src/gff2gtf_rtracklayer.R $projectDir/data/genome/ecic_annotation.tmp.gff3 $projectDir/data/genome/ecic_annotation.gtf putEGin_gene_id
# Command above does not run.
# Used the modified Rscript gff2gtf_rtracklayer_windows.R in a Windows OS with Rstudio, copying and pasting the commands.


#### INDEX BIMP GENOME WITH STAR ####
# Create a folder for the index files
mkdir -p $projectDir/data/genome/STAR-2.4.2a-R64_Ecic_Index
# run STAR in index generation mode
STAR --runThreadN 4 \
 --runMode genomeGenerate \
 --genomeDir $projectDir/data/genome/STAR-2.4.2a-R64_Ecic_Index \
 --genomeFastaFiles $projectDir/data/genome/GCA_046563565.1_EBD_Ecic_2.1_genomic.fna.gz \
 --sjdbGTFfile $projectDir/data/genome/GCA_046563565.1_EBD_Ecic_2.1_genomic.gbff.gz \
 --sjdbOverhang 99
 
#### ALIGN THE READS TO BIMP REFERENCE GENOME WITH STAR ####
cat $projectDir/src/ecic_basenames.txt | while read line; do
mkdir -p $projectDir/results/star/$line.STAR;
STAR --runThreadN 4 \
  --genomeDir $projectDir/data/genome/STAR-2.4.2a-R64_Ecic_Index \
  --readFilesIn $projectDir/results/trimmomatic/$line\_1.trimmed.paired.fastq.gz $projectDir/results/trimmomatic/$line\_2.trimmed.paired.fastq.gz \
  --sjdbGTFfile $projectDir/data/genome/GCA_046563565.1_EBD_Ecic_2.1_genomic.gtf \
  --outFileNamePrefix $projectDir/results/star/$line.STAR/$line.STAR \
  --readFilesCommand gunzip -c \
  --quantMode TranscriptomeSAM GeneCounts \
  --twopassMode Basic \
  --outSAMunmapped Within  \
  --outSAMtype BAM SortedByCoordinate;
 # Output your alignment file in BAM format, sorted by coordinate
 # Give your output files a UNIQUE name/prefix (${line} also works)
 # Gene counting for generating gene count tables
 # Two-pass mode for refining splice junctions
done

#### GENERATE MULTIQC REPORT ####
cd $projectDir/results/star/
multiqc $projectDir/results/star/


#### GET MAPPED READ COUNTS WITH FEATURECOUNTS ####
cat $projectDir/src/ecic_basenames.txt | while read line; do
mkdir $projectDir/results/featureCounts/$line.Counts;
featureCounts -a $projectDir/data/genome/GCA_046563565.1_EBD_Ecic_2.1_genomic.gtf \
  -p --countReadPairs -s 0 -t exon -g gene_id \
  -o $projectDir/results/featureCounts/$line.Counts/$line.featurecounts \
  $line.out.bam$projectDir/results/star/$line.STAR/$line.STARAligned.sortedByCoord.out.bam
done

#### GENERATE MULTIQC REPORT ####
cd $projectDir/results/featureCounts/
multiqc $projectDir/results/featureCounts/

