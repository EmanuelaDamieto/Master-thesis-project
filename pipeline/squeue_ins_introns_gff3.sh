#!/bin/bash -l

set -eux

export SINGULARITY_BINDPATH=/mnt:/mnt

proj=u2019016
mail=emanuela.damieto@studenti.unitn.it

#real path to make the code more reproducible
infl=$(realpath ../reference/gff3/spruce_pine-extended_master.gff3)
path=$(realpath ../data/gff3)

[[ ! -d $path ]] && mkdir -p $path

#take singularity container file from kogia because there are more recent versions of files
sing_genometools=$(realpath ../singularity/kogia/genometools_1.6.2.sif)

#-A project, -J job name, -o output file, -e error file, --mail.user to specify the mail
sbatch -A $proj -J genometools -o $path/ins_introns.out -e $path/ins_introns.err \
--mail-user $mail ins_introns_gff3.sh $sing_genometools $infl $path/spruce_pine-extended_master_complete.gff3
