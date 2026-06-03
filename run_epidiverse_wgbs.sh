#######################################
###   RUNNING EPIDIVERSE PIPELINE   ###
#######################################

## set conda to recognize your SSL certificates
conda config --set ssl_verify /etc/ssl/certs/ca-bundle.trust.crt

## make sure anaconda is instaled
curl -O https://repo.anaconda.com/archive/Anaconda3-2020.07-Linux-x86_64.sh
bash ./Anaconda3-2020.07-Linux-x86_64.sh
## browse through the document and accept terms and conditions typing yes
## type yes when asked for conda init to write PATH in .bashrc, and source it
source .bashrc
rm Anaconda3-2020.07-Linux-x86_64.sh

## update conda
#conda update conda
conda update -n base -c defaults conda

## create a new conda environment if needed
containerEnv='epidiv' # switch epidiv with your own conda env
conda create --name $containerEnv

## log in conda environment (you can log in yours/
## set your own from scratch) and update environment
conda activate $containerEnv 
conda update --all

#####################################
#### PREVIOUS STEPS AND TEST RUN ####
#####################################

## these steps are to take if running the analysis 
## in a new conda environment/new computer

## install java v11.0.13 (v11.0.1 gives an error)
conda install -c conda-forge openjdk=11.0.13
java -version # check that v11.0.13 is default java version

## install nextflow
conda install -c bioconda nextflow
nextflow # check that v.24.04.2 is default version

## install r-base v3.4.3 
conda install -c conda-forge r-base=3.4.3
#conda install -c free r-base=3.4.3
#conda install -c pkgs/r r-base=3.4.3 # package inconsistency error!
R --version

## run epidiverse/wgbs test (use nextflow v20.07.1, otherwise, errors)
mkdir test_wgbs && cd test_wgbs
NXF_VER=20.07.1 nextflow run epidiverse/wgbs -profile test,conda # test successful!
## run epidiverse/dmr test
cd .. && mkdir test_dmr && cd test_dmr
NXF_VER=20.07.1 nextflow run epidiverse/dmr -profile test,conda # error: metilene resample rate parameters $X and $Y are empty, thus metilene fails to run when specifying them as variables
## check -nextflow.log for errors. If all is OK, delete test folder 
#cd .. && rm -r test_*


#############################
#### PREPARE INPUT FILES ####
####  FOR AN EXPERIMENT  ####
#############################

## create directory for the reference genome
mkdir -p ~/epidiv/ecic_genome && cd ~/epidiv/

## run the WGBS pipeline with the --INDEX option to generate genome indexes
genome='/home/user/epidiv/ecic_genome'
scp usuario@10.222.7.83:/home/usuario/Windows/Data/genomes/Ecic_genome_pacBio/MaSuRCA_assembly.fasta.PolcaCorrected.fa $genome/ecic_genome.fa

## set folder for reads
readsDir='/home/user/epidiv/reads/herbivory' # change reads path to switch projects
mkdir -p $readsDir

## download reads
# E. cicutarium herbivory in leaves WGBS experiment
scp usuario@10.222.7.83:/media/usuario/One\ Drive/HN00218941/RawFASTQ/*.fastq.gz $readsDir/

## ¡WARNING! Check the integrity of the reads
scp usuario@10.222.7.83:/media/usuario/One\ Drive/HN00218941/md5sum.txt $readsDir/
sed -i 's#RawData/##'g md5sum.txt
md5sum -c md5sum.txt 

## If the md5sum does is not the same, the file might be corrupt.
## Unzip the file, check for lines with corruption in the fastq
## and remove the whole read entry from the file, and its mate pair from the other file
## This has to be done manually
 

#############################
#### RUN PIPELINES: WGBS ####
#############################

## run epidiverse/wgbs pipeline
projectDir='/home/user/epidiv/ecic_herbivory_test' # change project path to switch projects
mkdir -p $projectDir/wgbs && cd $projectDir/wgbs
NXF_VER=20.07.1 nextflow run epidiverse/wgbs \
  -profile conda \
  --input $readsDir \
  --reference $genome/ecic_genome.fa \
  # once BAMs are generated, uncomment this option for downstream reanalysis involving only methylated Cs calling
  #--CALL \
  # perform trimming. Disable when trimmed reads are provided
  --trim \
  # Filter any positions with a phred value below Q30
  --minQual 30 \
  # generate fastQC from the trimmed reads. Disable when disabling --trim option.
  --fastqc \
  # wgbs will calculate bisulfite non conversion rate for each bam file, placing a BisNonConvRate.txt file within the bam directory
  # uncomment this option when calculating bisulfite conversion rate is not needed
  #--noLambda #\
  # activate in case genome indexes are not already generated (erne-bs5 indexing gives an error, prepared ebm index running erne-bs5 index command)
  --INDEX #\
  # in case analysis finishes with errors, uncomment this option after fixing the errors (i.e. incorrect path, deficient dependency ...) to resume the analysis
  #-resume
## oneliner
#NXF_VER=20.07.1 nextflow run epidiverse/wgbs -profile conda --input $readsDir --reference $genome/ecic_genome.fa --trim --minQual 30 --fastqc --INDEX

## obtain report  with all the Non-conversion Rates from each library:
cat $projectDir/wgbs/wgbs/bam/ecicBS*/stats/BisNonConvRate.txt > $projectDir/wgbs/wgbs/BisNonConvRate_summary.txt
# Caution! this report shows the Non-Conversion! Data has to be transformed into conversion percentage with the formula 1 - NonConversionRate * 100


#############################
#### RUN PIPELINES: DMR  ####
#############################

cd $projectDir && mkdir $projectDir/dmrs && cd $projectDir/dmrs
## generate a tab delimited TSV table with samples, experimental groups, and within group replicates:
##sampleA_1    groupA    rep1
##sampleA_2    groupA    rep2
##sampleB_1    groupB    rep1
ls $projectDir/wgbs/wgbs/bam | grep -v bam - > column1.txt
ls $projectDir/wgbs/wgbs/bam | grep -v bam -  | cut -d$'_' -f3 > column2.txt
ls $projectDir/wgbs/wgbs/bam | grep -v bam -  | cut -d$'_' -f1,2 | sed 's/_//g' > column3.txt
paste column* > samples.tsv
rm column*

## Run the DMR script
NXF_VER=20.07.1 nextflow run epidiverse/dmr \
  -profile conda \
  # specify path to BAM files from WGBS analysis
  --input $projectDir/wgbs/wgbs/bedGraph \
  # generate a tab delimited table with sample and group names
  --samples samples.tsv \
  # metilene gives an error about resample variables X and Y not bein integer, so run this option to keep the resample as an integer
  --resample 1
  # designate one of the groups from the former table as control group (comment to allow a H vs S comparison)
  #--control C #\
  # analyze DMPs instead of DMRs
  #--dmp
  ## parameters to consider a DMR/DMP:
  # minimum coverage (default: 5) 
  #--coverage 6 \
  # minimum number of citosines (default: 10, not computed for DMP mode)
  #--CpN 10 \
  # minimum percentage of methylation difference between groups (default: 10)
  #--diff 20

## metilene gives an error:

#Error executing process > 'DMRS:metilene (CpG - C_vs_H)'
#
#Caused by:
#  Process `DMRS:metilene (CpG - C_vs_H)` terminated with an error exit status (1)
#
#Command executed:
#
#  mkdir tmp CpG CpG/C_vs_H
#  bed=inputs/C_vs_H/input.bed
#
#  # define resample rate parameters
#  X=$(printf "%.0f\n" $(echo $(head -1 $bed | grep -o C | wc -l)*0.8 | bc -l)) 
#  # running:$ $(echo $(head -1 $bed | grep -o C | wc -l)*0.8 | bc -l)
#  # gives the error:$ bash: 12*0.8: no se encontró la orden
#  Y=$(printf "%.0f\n" $(echo $(head -1 $bed | grep -o H | wc -l)*0.8 | bc -l))
#
#  # run metilene
#  metilene -X $X -Y $Y -a C -b H \
#  -M 146 -c 2 -m 10 -d 0.1 -t 8 \
#  $bed 1> CpG/C_vs_H/C_vs_H.bed 2> CpG/C_vs_H/C_vs_H.log || exit $?
#
#  # filter metilene output
#  awk 'BEGIN {OFS="\t"} $4 <= 0.05 {len=$3-$2; print $1,$2,$3,$6,$5,$4,len}' CpG/C_vs_H/C_vs_H.bed |
#  sort -k1,1 -k2,2n -T tmp > CpG/C_vs_H/C_vs_H.0.05.bed
#  ln -s C_vs_H/C_vs_H.0.05.bed CpG/C_vs_H.bed
#
#Command exit status:
#  1
#
#Command output:
#  (empty)
#
#Command error:
#  .command.sh: línea 6: printf: 4.8: número inválido
#
#Work dir:
#  /home/mmedrano/epidiv/ecic_herbivory_test/dmr/work/6b/3bf8b427184f55d3ba69c300f7138b
#
#Tip: you can try to figure out what's wrong by changing to the process work dir and showing the script file named `.command.sh`
 
# To complete the epidiverse/dmr step I had to run metilene separately from the pipeline,
# using the bedgraph files originated from epidiverse/dmr, with the script run_metilene.sh
