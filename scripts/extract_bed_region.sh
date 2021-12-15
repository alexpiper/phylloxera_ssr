#!/bin/bash
#SBATCH --job-name=ext_region         
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20GB 
#SBATCH --time=640:00:00
#SBATCH --mail-user=alexander.piper@agriculture.vic.gov.au
#SBATCH --mail-type=ALL
#SBATCH --account=pathogens
#SBATCH --export=none
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.out

#--------------------------------------------------------------------------------
#-                                  HEADER		                                -
#--------------------------------------------------------------------------------

# Welcome to the insect SkimSeq Pipeline
# Developed by Alexander Piper 
# Contact: alexander.piper@agriculture.vic.gov.au

# This script searches the genome for the coordinates of all sequences in a reference fast
# then extracts those regions from BAM files and calls variants

if [ -z "$SLURM_ARRAY_TASK_COUNT" ]; then 
  echo SLURM_ARRAY_TASK_COUNT unset; 
  echo You must launch this job as an array
  echo see https://slurm.schedmd.com/job_array.html
  echo for info on how to run arrays
  exit 1
fi

Index=extraction_job_index.txt

# Check sequence index file exists
if [[ ! -f "${Index}" ]]; then
  echo "Error sequence index file ${Index} does not exist"
  exit 1
fi

#--------------------------------------------------------------------------------
#-                                    Preparation                               -
#--------------------------------------------------------------------------------

# Input information 
ReferenceGenome=/group/pathogens/Alexp/genomes/phylloxera/Dv_genome_V4.0.fasta
FullSampleName=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${Index})
Sample=$(basename $FullSampleName | cut -d'.' -f1 )

# Input reference fasta from stdin 1
ref_bed=$1
refname=$(basename $ref_bed | cut -d'.' -f1 )

# set minimum depth, length of non-N bases, and how many extra bases on each side to extract
mindepth=5
minlength=20

# Create output name
outname=$(echo $Sample $refname "dp"$mindepth | sed 's/ /_/g' )

echo ReferenceGenome=${ReferenceGenome}
echo FullSampleName=${FullSampleName}
echo Sample=${Sample}
echo ref_bed=${ref_bed}
echo SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR}
echo SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}
echo mindepth=${mindepth}
echo minlength=${minlength}

# Make directories for outputs
mkdir ${SLURM_SUBMIT_DIR}/extract

# Goto tmp
cd $TMPDIR
tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)
cd ${tmp_dir}
pwd

# Copy data files to temp and decompress
cp ${FullSampleName} .
cp ${FullSampleName}.bai .
cp $(dirname ${ReferenceGenome})/* .
cp ${ReferenceGenome}.fai .
pigz -p ${SLURM_CPUS_PER_TASK} -d ./*.gz

ls

#--------------------------------------------------------------------------------
#-                            	  Extract regions                         		-
#--------------------------------------------------------------------------------
# Load modules
module purge
module load bedops/2.4.35
module load SAMtools/1.12-GCC-9.3.0
module load BCFtools/1.12-GCC-9.3.0

cp ${ref_bed} ${Sample}_alignments.bed

# make a regions file from the bed file for Samtools
awk 'BEGIN { FS="\t" } { print $1":"$2"-"$3}' ${Sample}_alignments.bed > ${Sample}_alignments.regions

# print regions
echo extracting regions:
cat ${Sample}_alignments.regions

# Subset bam to only coordinates
samtools view -@ ${SLURM_CPUS_PER_TASK} -b ${Sample}.bam -L ${Sample}_alignments.bed > ${Sample}_subset.bam

# Index resulting bam file
samtools index -@ ${SLURM_CPUS_PER_TASK} ${Sample}_subset.bam

# Check bam had reads
#samtools tview -p VWMZ01005123.1:17154145 bams/1059c.bam

# Call variants - create a gVCF so all sites are included 
bcftools mpileup -Ou -f ${ReferenceGenome} ${Sample}_subset.bam -C 50 -E -q30 -Q20 --gvcf ${mindepth} | \
bcftools call --gvcf ${mindepth} -Ou -m | \
bcftools norm -f ${ReferenceGenome} -Oz -o ${Sample}.vcf.gz

# Index VCF
tabix ${Sample}.vcf.gz 

mkdir ${outname}
# Call consensus sequence and output as fasta file - Use -a to set absent bases to N
samtools faidx ${ReferenceGenome} -r ${Sample}_alignments.regions  | bcftools consensus ${Sample}.vcf.gz -a N  > ${outname}.fa

# Split multifasta into separate fastas
awk -F "|" '/^>/ {close(F); ID=$1; gsub("^>", "", ID); F=ID".fasta"} {print >> F}' ${outname}.fa

cp -r $(ls | grep .fasta | grep -v $(basename ${ReferenceGenome})) ${outname}

mv ${Sample}_subset.bam ${outname}
mv ${Sample}_subset.bam.bai ${outname}

cp -r ${outname} ${SLURM_SUBMIT_DIR}/extract/.

# Output useful job stats
/usr/local/bin/showJobStats.scr 