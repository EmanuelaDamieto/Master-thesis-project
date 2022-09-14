#!/bin/bash 
#SBATCH -p core
#SBATCH -n 1
#SBATCH --mem=16GB
#SBATCH --mail-type=END,FAIL
#SBATCH -t 12:00:00

# stop on error 
set -eux

#ins_introns_gff3.sh $sing_genometools $infl $path/spruce_pine-extended_master_complete.gff3

#path : cd /mnt/picea/projects/singularity/genometools_1.6.2.sif

## arguments
[[ $# -ne 3 ]] && abort "This script takes three arguments"

[[ ! -f $1 ]] && abort "The first argument needs to be an existing genometools singularity container file."

## enforce singularity
[[ -z ${SINGULARITY_BINDPATH:-} ]] && abort "This function relies on singularity, set the SINGULARITY_BINDPATH environment variable"

[[ ! -f $2 ]] && abort "The second argument needs to be an existing input file."

[[ ! -d $(dirname $3) ]] && abort "The third argument parent directory needs to exist."

# start
singularity exec $1 gt gff3 -force -tidy yes -addids yes \
-fixregionboundaries yes -retainids yes -sortlines yes -checkids yes -addintrons yes -o $3 $2
