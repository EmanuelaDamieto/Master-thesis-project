#!/bin/bash -l
#Script to submit multiqc to the job queue

## stop on error but be verbose
set -eux

#inside ~/Git/Master-thesis-project/pipeline call the file since its executable: ./submitMultiQC.sh
#we have to specify singularity bindpath -> set the environment variable 
export SINGULARITY_BINDPATH=/mnt:/mnt

proj=u2019016
mail=emanuela.damieto@studenti.unitn.it

#real path to make the code more reproducible
INDIR=$(realpath ../data)
OUTDIR=$(realpath ../data/MultiQC)
#take singularity container file from kogia because there are more recent versions of files
SING_MULTIQC=$(realpath ../singularity/kogia/multiqc_1.10.sif)

#[[condition]] substitute if condition fi
#! = doesn't exist
#-d directory 
#mkdir -p to create recursively (it works even if result it doesn't exist and creates both)
[[ ! -d $OUTDIR ]] && mkdir -p $OUTDIR

#-A project, -J job name, -o output file, -e error file, --mail.user to specify the mail
sbatch -A $proj -J multiqc -o $OUTDIR.out -e $OUTDIR.err \
--mail-user $mail runMultiQC.sh $SING_MULTIQC $INDIR $OUTDIR
