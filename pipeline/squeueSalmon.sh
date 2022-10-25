#Script to submit salmon to the job queue

#to make the file executable:  chmod +x squeueSalmon.sh 
#inside ~/Git/Master-thesis-project/pipeline call the file since its executable: ./squeueSalmon.sh
#we have to specify singularity bindpath -> set the environment variable 
export SINGULARITY_BINDPATH=/mnt:/mnt

proj=u2019016
mail=emanuela.damieto@studenti.unitn.it

#real path to make the code more reproducible
INDIR=$(realpath /mnt/picea/storage/reference/Picea-abies/v2.0/indices/salmon)     #index directory
#fq files directory
FQDIR=$(realpath /mnt/picea/projects/spruce/vhurry/drought-stress-roots/preprosessed/trimmomatic)
#output directory 
OUTDIR=$(realpath ../data/SalmonResults)
#take singularity container file from kogia because there are more recent versions of files
SING_SALMON=$(realpath ../singularity/kogia/salmon_1.9.0.sif)

#[[condition]] substitute if condition fi
#! = doesn't exist
#-d directory 
#mkdir -p to create recursively (it works even if result it doesn't exist and creates both)
[[ ! -d $OUTDIR ]] && mkdir -p $OUTDIR

#-A project, -J job name, -o output file, -e error file, --mail.user to specify the mail
i=1
for f in $(find $FQDIR -name "*trimmomatic_1.fq.gz"); do 
    sbatch -A $proj -J salmon_i -o $OUTDIR_$(basename ${f/trimmomatic_1.fq.gz}).out -e $OUTDIR_$(basename ${f/trimmomatic_1.fq.gz}).err \
    --mail-user $mail runSalmon.sh $SING_SALMON $INDIR $f ${f/trimmomatic_1.fq.gz/trimmomatic_2.fq.gz} $OUTDIR/$(basename ${f/trimmomatic_1.fq.gz})
    i++
done 


#Or sbatch from the terminal: sbatch -A u2019016 -J salmon -o ../data/SalmonResults.out -e ../data/SalmonResults.err --mail-user emanuela.damieto@studenti.unitn.it squeueSalmon.sh
#runSalmon.sh $SING_SALMON $INDIR $f ${f/trimmomatic_1.fq.gz/trimmomatic_2.fq.gz} $OUTDIR/$(basename ${f/trimmomatic_1.fq.gz})
