#!/bin/bash

#A simple set of commands to perform the standard 3C mapping operations and 
#generate GRAAL-compatible matrices.
#
#This script must be run in the same directory as fraglist.py as that's what 
#it relies on to perform some of the GRAAL-specific stuff like matrix/fragment data
#file generation.

current_dir="$( cd "$( dirname "$0" )" && pwd )"
arguments=()

trigger_help=0

#Defaults
enzyme=5000
output_dir=$current_dir
quality_min=30
size=0
duplicate=0
clean_up=0
threads=1
minimap=0
circular=0

#Argument parsing
while [[ $# -gt 0 ]]; do

    key="$1"

    case $key in
        -1|--forward)
            reads_for="$2"
            shift 
            shift 
        ;;
        -2|--reverse)
            reads_rev="$2"
            shift 
            shift 
        ;;
        -f|--fasta)
            fasta="$2"
            shift 
            shift 
        ;;
        -e|--enzyme)
            enzyme="$2"
            shift 
            shift 
        ;;
        -o|--output)
            output_dir="$2"
            shift 
            shift 
        ;;
        -q|--quality-min)
            quality_min="$2"
            shift 
            shift 
        ;;
        -s|--size)
            size="$2"
            shift 
            shift 
        ;;
        -t|--threads)
            threads="$2"
            shift 
            shift 
        ;;
        -T|--tmp)
            tmp_dir="$2"
            shift
            shift
        ;;
        -m|--minimap|--minimap2)
            minimap=1
            shift 
        ;;
        -d|--duplicates)
            duplicate=1
            shift
        ;;
        -c|--clean-up)
            clean_up=1
            shift
        ;;
        -C|--circular)
            circular=1
            shift
        ;;
        -h|--help)
            trigger_help=1
            shift
        ;;
        *)    
            arguments+=("$1")
            shift 
        ;;
    esac
done

tmp_dir=$output_dir

set -e
set -o pipefail

set -- "${arguments[@]}"

#Hopefully useful help message if a parameter is missing
if [ -z $reads_for ] || [ -z $reads_rev ] || [ -z $fasta ] || [ $trigger_help -eq 1 ]; then
    echo ""
    echo "    yahcp (Yet Another Hi-C Pipeline) - A simple and relatively painless Hi-C data processing pipeline"
    echo ""
    echo "    Usage: ./yachp.sh -1 reads_forward.fastq -2 reads_reverse.fastq -f genome.fa [-s size] [-o output_directory] [-e enzyme] [-q quality_min] [--duplicates] [--clean-up]"
    echo ""
    echo "    Generate a sparse, GRAAL-compatible contact map from paired-end reads and a reference genome (for more information about GRAAL, see https://github.com/koszullab/GRAAL)."
    echo "    The map can also be easily visualized with HiC-Box (see https://github.com/koszullab/HiC-Box)."
    echo "    Information about fragments and contigs/chromosomes are stored in separate files."
    echo "    The genome can either be partitioned by restriction fragments (specifying the enzyme) or fixed size chunks (specifying a number)."
    echo ""
    echo "    Requires bowtie2, samtools, bedtools and python (with Biopython installed) to run. Optionally, minimap2 can be used instead of bowtie2 if specified."
    echo ""
    echo "    Parameters:"
    echo ""
    echo "        -1 or --forward: Forward FASTQ reads"
    echo "        -2 or --reverse: Reverse FASTQ reads"
    echo "        -f or --fasta: Reference genome to map against in FASTA format"
    echo "        -o or --output: Output directory. Defaults to the current directory."
    echo "        -e or --enzyme: Restriction enzyme if a string, or chunk size (i.e. resolution) if a number. Defaults to 5000 bp chunks."
    echo "        -q or --quality-min: Minimum mapping quality for selecting contacts. Defaults to 30."
    echo "        -d or --duplicates: If enabled, removes adapters and PCR duplicates prior to mapping. Not enabled by default."
    echo "        -s or --size: Minimum size threshold to consider contigs. Defaults to 0 (keep all contigs)."
    echo "        -c or --clean-up: If enabled, removes intermediary BED files after generating the contact map. Not enabled by default."
    echo "        -t or --threads: Number of threads to use for the aligner and samtools. Defaults to 1."
    echo "        -m or --minimap: Use the minimap2 aligner instead of bowtie2. Not enabled by default."
    echo "        -h or --help: Display this help message"
    echo ""
    echo "    The expected files in the output directory will take the form:"
    echo "        -abs_fragments_contacts_weighted.txt: the sparse contact map"
    echo "        -fragments_list.txt: information about restriction fragments (or chunks)"
    echo "        -info_contigs.txt: information about contigs or chromosomes"
    echo ""
    exit 1
fi

#Check everything's in place
if [ ! -f $current_dir/fraglist.py ]; then
    echo "Something went wrong: couldn't detect fraglist.py in the same directory as toolbox.sh"
    exit 1
fi

python -c "import Bio" >/dev/null 2>&1 || { echo "Error! Biopython is missing from your python libraries. Please install it (using either your package manager or pip)"; exit 1; }

if [ $minimap -eq 0 ]; then
    aligner=bowtie2
else
    aligner=minimap2
fi

for tool in $aligner samtools bedtools; do
    command -v $tool >/dev/null 2>&1 || { echo "Error! $tool is needed and could not be found on your machine."; exit 1; }
done

index=${fasta%.fa}
t=$(( $threads/2 < 1 ? 1 : $threads/2 ))

mkdir -p $output_dir
mkdir -p $tmp_dir

#Write fragments_list.txt info_contigs.txt
echo "Writing fragment information..."
python $current_dir/fraglist.py --fasta $fasta --enzyme $enzyme --output-dir $output_dir --size $size

#Build fasta index files

if [ $minimap -eq 0 ]; then
    if [ ! -f ${index}.1.bt2 ]; then
        echo "Building fasta index files..."
        bowtie2-build --quiet $fasta $index
    fi
fi

#Remove adapters and PCR duplicates
if [ $duplicate -eq 1 ]; then
    echo "Removing adapters and PCR duplicates..."

    if [ ! -f $current_dir/pcr_duplicate_Hiseq20.pl ]; then
        echo "Something went wrong: couldn't detect pcr_duplicate_Hiseq20.pl in the same directory as toolbox.sh"
        exit 1
    fi

    perl $current_dir/pcr_duplicate_Hiseq20.pl $reads_for $reads_rev "$reads_for".trimmed "$reads_rev".trimmed
    reads_for=${reads_for}.trimmed
    reads_rev=${reads_rev}.trimmed
fi

#Do the job: 
# 1/ bowtie2/minimap2 alignment
# 2/ filter reads according to flag and mapping quality
# 3/ convert to bed
# 4/ keep relevant fields (chromosome/contig, starting position, end position, read name, orientation) 
# 5/ sort by chromosome/contig (#1 in dict order) then by position (#2 in numerical order)
echo "Performing alignment and generating bed files..."
if [ $minimap -eq 0 ]; then
    alignment_instruction_for="bowtie2 --very-sensitive-local -p $t -x $index -U $reads_for"
    alignment_instruction_rev="bowtie2 --very-sensitive-local -p $t -x $index -U $reads_rev"
else
    alignment_instruction_for="minimap2 -2 -t $t -ax sr $fasta $reads_for"
    alignment_instruction_rev="minimap2 -2 -t $t -ax sr $fasta $reads_rev"
fi

$alignment_instruction_for \
    | samtools view -bS -F 260 -@ $t -q $quality_min - \
    | bedtools bamtobed -i - \
    | awk 'OFS="\t" { print $1,$2,$3,$4,$6 }' \
    > $tmp_dir/unsorted_contacts_for.bed &

$alignment_instruction_rev \
    | samtools view -bS -F 260 -@ $t -q $quality_min - \
    | bedtools bamtobed -i - \
    | awk 'OFS="\t" { print $1,$2,$3,$4,$6 }' \
    > $tmp_dir/unsorted_contacts_rev.bed & 

wait

sort -T $tmp_dir -k1,1d -k2,2n $tmp_dir/unsorted_contacts_for.bed $tmp_dir/unsorted_contacts_rev.bed \
    > $tmp_dir/total_contacts.bed

#Make a bed out of fragments_list.txt
awk 'NR>1 { print $2"\t"$3"\t"$4 }' $tmp_dir/fragments_list.txt > $tmp_dir/fragments_list.bed 

#Intersect fragment list with mapping data
echo "Intersecting bed files..."
bedtools intersect -a $tmp_dir/total_contacts.bed -b $tmp_dir/fragments_list.bed -wa -wb \
    | awk 'OFS="\t" { print $1,$2,$3,$4,$5,$(NF-2),$(NF-1),$NF }' \
    | sort -k4d -T $tmp_dir \
    > $tmp_dir/contact_intersect_sorted.bed

#Write GRAAL matrix out of intersecting bed file
echo "Generating contact map..."
python $current_dir/fraglist.py --intersection $tmp_dir/contact_intersect_sorted.bed --frags $tmp_dir/fragments_list.txt --output-dir $output_dir

if [ $clean_up -eq 1 ]; then
    rm $tmp_dir/fragments_list.bed
    rm $tmp_dir/unsorted_contacts_for.bed
    rm $tmp_dir/unsorted_contacts_for.bed
    rm $tmp_dir/total_contacts.bed
    rm $tmp_dir/contact_intersect_sorted.bed
fi

echo "Finito"
