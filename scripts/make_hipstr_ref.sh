
cd ~
module purge
module load Python/3.8.2-GCCcore-9.3.0 
source ~/hipstr/bin/activate
git clone https://github.com/HipSTR-Tool/HipSTR-references.git
cd HipSTR-references

# Setup 
mkdir /group/pathogens/Alexp/microsat/trf_reference
cd /group/pathogens/Alexp/microsat/trf_reference
cp -r /home/ap0y/HipSTR-references/scripts .
cp /home/ap0y/HipSTR-references/trf409.legacylinux64 .
chmod 755 trf409.legacylinux64
chmod 755 scripts/run_TRF.sh

# make a regions file from the microsat bed file for Samtools
#awk 'BEGIN { FS="\t" } { print $1":"$2"-"$3}' /group/pathogens/Alexp/microsat/microsats.bed > microsats.regions
#cat /group/pathogens/Alexp/microsat/microsats.bed |  awk 'BEGIN { FS=" " } { print $1}' | sort | uniq > contigs.txt

# Subset the assembly to the regions in bed
#module load SAMtools/1.14-GCC-11.2.0
#samtools faidx /group/pathogens/Alexp/genomes/phylloxera/Dv_genome_V4.0.fasta -r microsats.regions > Dv_genome_V4_subset.fa

# Subset assembly to just the contigs of interest
cat /group/pathogens/Alexp/microsat/microsats.bed |  awk 'BEGIN { FS=" " } { print $1}' | sort | uniq > contigs.txt
module load seqtk/1.3-GCC-8.2.0-2.31.1
seqtk subseq  /group/pathogens/Alexp/genomes/phylloxera/Dv_genome_V4.0.fasta  contigs.txt  > Dv_genome_V4_subset.fa

# Split subset assembly into separate fastas
mkdir split
while read line
do
    if [[ ${line:0:1} == '>' ]]
    then
        outfile=$(echo ${line#>} | sed 's/\:.*$//' ).fa
        header=$(echo $line | sed 's/\:.*$//' )
        echo $header > split/$outfile
    else
        echo $line >> split/$outfile
    fi
done <  Dv_genome_V4_subset.fa

# Run TRF
mkdir trf_results
for chrom in $(cat contigs.txt)
  do
  echo split/$chrom.fa trf_results 5
done | xargs -L 1 -P 30 scripts/run_TRF.sh

# Fix lengths
# NOTE This script gave error: TypeError: 'float' object cannot be interpreted as an integer
# range(len(nmer)/k) changed to range(len(nmer)//k) (added extra /) and now working
mkdir fixed_trf_results
for chrom in $(cat contigs.txt)
do
  echo scripts/fix_trf_output.py trf_results/$chrom.fa fixed_trf_results/$chrom.fa
done | xargs -L 1 -P 40 python

#Filtering the STRs
files=""
for chrom in $(cat contigs.txt)
do
files="$files,fixed_trf_results/$chrom.fa"
done
files=`echo $files | sed "s/,//"`


# Note xrange had to be changed to range in the below script to make it work with python 3
python scripts/trf_parser.py $files >  filtered_repeats.bed

module load BEDTools/2.30.0-GCC-10.2.0
bedtools sort -i filtered_repeats.bed > filtered_repeats.sorted.bed

# Add a dummy line to end of file as the analyze_overlaps.py script drops it
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "contig_dummy" "83" "96" "2" "AG" "7.0" "28" "AGAGAGAGAGAGAG" >> filtered_repeats.sorted.bed

# Merge overlapping STRs
python scripts/analyze_overlaps.py filtered_repeats.sorted.bed pass.dv fail.dv


#Remove any entries within 10bp of a failed merge region
bedtools window -w 10 -a pass.dv -b fail.dv -v > pass.dv.r2

#Extract entries that aren't within 10bp of another entry or are within 10bp of one or more entries that all share the same period
bedtools merge -i pass.dv.r2 -c 4,6 -o collapse -d 10 | grep -v "," > pass.dv.r3
bedtools merge -i pass.dv.r2 -c 4,4,4,6 -o collapse,count_distinct,distinct,collapse -d 10 | grep "," | awk '$5 == 1' | awk -v OFS="\t" '{print $1, $2, $3, $6, $7}' | sed "s/,/\//g"  >> pass.dv.r3

# Construct final references
cat pass.dv.r3 | awk -v OFS="\t" '{print $1, $2, $3, $4, ($3-$2+1)/$4, "STR_"NR, $5}' > pass.dv.bed 
#cat pass.dv.r3 | awk -v OFS="\t" '{print $1, $2, $3, $4, "STR_"NR, $5}' > dv.hipstr_reference.bed

# Sort bed file
bedtools sort -i pass.dv.bed  > pass.dv.bed_sorted

# Subset bed to just those regiosn contained in the original bed file
bedtools intersect -wa -a pass.dv.bed -b /group/pathogens/Alexp/microsat/microsats.bed > dv.hipstr_reference.bed
 
# Join with startign reference list to get the loci names
bedtools sort -i /group/pathogens/Alexp/microsat/microsats.bed > tmp.bed
awk -v OFS="\t" 'FNR==NR{a[NR]=$4;next}{$6=a[FNR]}1' tmp.bed dv.hipstr_reference.bed > final_trf_strs.bed 


# Try intersecting the unfiltered one
bedtools intersect -wa -a trf_results/ -b /group/pathogens/Alexp/microsat/microsats.bed 



