#!/bin/bash

input_bams_dir=$1 # directory with names of bam files
annotation_file=$2 # file with annotations in GTF format
output_dir=$3 # output directory for count matrix

# activate our conda env
source activate subread

# safety check for input file
if [ -z "$input_bams_dir" ] || [ ! -f "$annotation_file" ] || [ -z "$output_dir" ]; then
    echo "Usage: <input_bam_txt_directory <gtf_file>"
    exit 1
fi

# makes the output dir just in case
mkdir -p $output_dir

# params for featurecounts run
featureCounts -s 0 \
    -T 8 \
    -t exon \
    -g gene_id \
    -a "$annotation_file" \
    -o "$output_dir"/feature_counts.txt \
    "$input_bams_dir"/*.bam

echo "Finished calculating counts"

# params for feature counts
# -t is the feautre, exon counts only coding sequences
# -s tells you strandess of rna library, 0 is for unstranded
# options are, 0 unstranded, 1 stranded, 2 reverse stranded
# -a is for annotation files, GTF
# -g is extremely important, it collapses isoforms to prevent watering down of analysis
# -o is for output
# -T is for adding threads
# -O, allows the usage of multi-mapping reads
# the input bam file is postional