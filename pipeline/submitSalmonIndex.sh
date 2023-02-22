#!/bin/bash
#Script to submit salmon index to the job queue

#To be safe
set -eux

#to make the file executable:  chmod +x squeueSalmonIndex.sh 
#inside ~/Git/Master-thesis-project/pipeline call the file since its executable: ./submitSalmonIndex.sh
#we have to specify singularity bindpath -> set the environment variable 
export SINGULARITY_BINDPATH=/mnt:/mnt

proj=u2019016
mail=emanuela.damieto@studenti.unitn.it

#The first argument is a singularity file
SING_SALMON=$(realpath ../singularity/kogia/salmon_1.9.0.sif)
#The second argument is the input fasta file 
INPUTFILE=$(realpath ../reference/fasta/Picab02_codingAll_mRNA.fa.gz)
#The third argument is an output directory for the indices 
OUTDIR=$(realpath ../reference/indices/salmon/Picea-abies-protein-coding_without-decoy_salmon-version-1-dot-9-dot0)

#[[condition]] substitute if condition fi
#! = doesn't exist
#-d directory 
#mkdir -p to create recursively (it works even if result it doesn't exist and creates both)
[[ ! -d $OUTDIR ]] && mkdir -p $OUTDIR

#-A project, -J job name, -o output file, -e error file, --mail.user to specify the mail
sbatch -A $proj -J salmonIndex -o $OUTDIR.out -e $OUTDIR.err \
--mail-user $mail runSalmonIndex.sh $SING_SALMON $INPUTFILE $OUTDIR


