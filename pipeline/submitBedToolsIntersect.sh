#!/bin/bash
#Script to submit BedToolsIntersect to the job queue
#Perform intersection between annotation and repeats file to see 

#To be safe
set -eux
#to make the file executable:  chmod +x squeueBedToolsIntersect.sh 
#inside ~/Git/Master-thesis-project/pipeline call the file since its executable: ./submitBedToolsIntersect.sh

proj=u2019016
mail=emanuela.damieto@studenti.unitn.it

module load bioinfo-tools
module load BEDTools

#3 arguments: 2 input files and 1 output directory
#Repeats file
INREP=$(realpath ../reference/gff3/pabies-2.0_chromosomes_and_unplaced.fa.repeats.gff)

#Annotation file 
#INANN=$(realpath ../data/gff3/introns.gff3)
INANN=$(realpath ../data/gff3/cds.gff3)

#Specify output directory
OUTDIR=$(realpath ../results/BedToolsIntersect)

#mkdir -p to create recursively (it works even if result it doesn't exist and creates both)
[[ ! -d $OUTDIR ]] && mkdir -p $OUTDIR

#-A project, -J job name, -o output file, -e error file, --mail.user to specify the mail
sbatch -A $proj -J BedToolsIntersect -o $OUTDIR/BedToolsIntersect_$(basename ${INANN/.gff3}).out \
-e $OUTDIR/BedToolsIntersect_$(basename ${INANN/.gff3}).err \
--mail-user $mail runBedToolsIntersect.sh $INREP $INANN $OUTDIR
