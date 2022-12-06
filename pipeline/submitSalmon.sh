#!/bin/bash
#Script to submit salmon to the job queue

#To be safe
set -eux
#to make the file executable:  chmod +x squeueSalmon.sh 
#inside ~/Git/Master-thesis-project/pipeline call the file since its executable: ./submitSalmon.sh
#we have to specify singularity bindpath -> set the environment variable 
export SINGULARITY_BINDPATH=/mnt:/mnt

proj=u2019016
mail=emanuela.damieto@studenti.unitn.it

#real path to make the code more reproducible
# INDIR=$(realpath ../reference/indices/salmon/Picea-abies-mRNA-temp-merge_without-decoy_salmon-version-1-dot-9-dot0)     #index directory
INDIR=$(realpath ../reference/indices/salmon/Picea-abies-protein-coding_without-decoy_salmon-version-1-dot-9-dot0)
#fq files directory
FQDIR=$(realpath ../data/trimmomatic)
#output directory 
OUTDIR=$(realpath ../results/Salmon)
#take singularity container file from kogia because there are more recent versions of files
SING_SALMON=$(realpath ../singularity/kogia/salmon_1.9.0.sif)


#[[condition]] substitute if condition fi
#! = doesn't exist
#-d directory 
#mkdir -p to create recursively (it works even if result it doesn't exist and creates both)
[[ ! -d $OUTDIR ]] && mkdir -p $OUTDIR

#-A project, -J job name, -o output file, -e error file, --mail.user to specify the mail
for f in $(find $FQDIR -name "*trimmomatic_1.fq.gz"); do 
    #define a variable with the name of the file
    fnam=$(basename ${f/_sortmerna_trimmomatic_1.fq.gz})
    mkdir -p $OUTDIR/$fnam
    sbatch -A $proj -J salmon-${fnam} -o ${OUTDIR}/$fnam.out -e ${OUTDIR}/$fnam.err \
    --mail-user $mail runSalmon.sh $SING_SALMON $INDIR $f ${f/trimmomatic_1.fq.gz/trimmomatic_2.fq.gz} $OUTDIR/$fnam
done 


#Or sbatch from the terminal: sbatch -A u2019016 -J salmon -o ../data/SalmonResults.out -e ../data/SalmonResults.err --mail-user emanuela.damieto@studenti.unitn.it squeueSalmon.sh
#runSalmon.sh $SING_SALMON $INDIR $f ${f/trimmomatic_1.fq.gz/trimmomatic_2.fq.gz} $OUTDIR/$(basename ${f/trimmomatic_1.fq.gz})
