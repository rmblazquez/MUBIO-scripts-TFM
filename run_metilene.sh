# Running metilene to get differentially methylated regions

# Why I am doing this? The epidiverse/dmrs pipeline failed to detect parameters X and Y, so I had to run manually the metilene step

# I am using the bedtools output from the epidiverse/dmrs pipeline 

# to get the metilene binary for Linux
wget http://www.bioinf.uni-leipzig.de/Software/metilene/metilene_v02-9.tar.gz
# I had to download it on aa Windows and upload to EBD:Genomics-b via FileZilla. I had to run `make` in EBD:Genomics-b to compile the binary

# Generate an alias for the binary with path (Windows path)
metilene="/mnt/d/Ruben_backup/Postdoc/2022_AlonsoLabPostdoc/Analyses/EcicEcaz_2022/Ecic_WGBS_herb/ecic_herbivory_test/dmr/dmrs/metilene_v02-9/metilene_v0.2-9/metilene"

# To test that metilene runs, just type
$metilene

# Now, for a basic DMR analysis:
$metilene -a C -b H input.bed > metilene_output.txt
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
$metilene -a C -b H -d 0.2 -X 4 -Y 4 -v 0.5 input.bed > metilene_output.txt
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

# INTEGRATE ANNOTATION BED FILE
# Investigate "bedtools intersect" and "bedtools annotation" 
