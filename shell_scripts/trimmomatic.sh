#!/bin/bash
################################################################################
# Run the following command to download the adapter file for Single End Illumina
# Modify as needed to obtain the Paired End file as needed
# wget https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq3-SE.fa
#
# Whatever directory you ran the command in will be the input for the adapter dir if no argument is
# provided
#
# The following parameters will be applied
# TrueSeq3-SE adapters removed or Paired End if chosen
# 2 mismatches allowed
# 30 palindrome clip threshold
# 10 simple clip threshold
# leading:3 for the removal of low-qual base pairs from the beginning (below 3)
# trailing:3 for removing low-qual bases from the end (below 3)
# slidingwindow:4:15, scans with a 4 base window and cuts after drop of 15 in qual
# minlen:36 discards reads shorter than 36 bases

# Example usage:
# ./trimmomatic.sh /input/to/files/ /output/dir/ /adapter/dir/ single

command -v trimmomatic >/dev/null 2>&1 || { echo "trimmomatic not found"; exit 1; }

function trim_single_end() {
    # Loops through input directory
    for fastq_file in $input_directory/*;do
        # Comparison test to make sure files are fasta files
        if [[ $fastq_file == *fastq.gz || $fastq_file == *fastq ]]; then
            # Creates base names for output files
            file_name=$(basename $fastq_file)
            output_file=$output_directory/trimmed_$file_name

            echo "now trimming $fastq_file -> $output_file";

            trimmomatic SE -phred33 $fastq_file $output_file \
                ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 \
                LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
                
            # Warning check for errors that occur
            if  [ $? -ne 0 ]; then
                echo "Error Trimming $fastq_file"
            fi
        fi
    done
}

function 

source activate genomics

date 
input_directory=$1 
output_directory=$2
adapters=$3
mode=$4

# Added safety test to ensure all parameters are provided
if [ ! -d "$input_directory" ] || [ ! -d "$output_directory" ] || [ -z "$adapters" ] || [ -z "$mode" ]; then
    echo "Usage: $0 <input_directory> <output_directory> </path/to/adapters>"
    exit 1
fi

echo $input_directory will now be processed

# Creates output directory just in case it does not exist
mkdir -p $output_directory

# Argument check
if [[ "$mode" == "single" ]]; then
    trim_single_end
elif [[ "$mode" == "paired" ]]; then
    trim_paired_end
else
    echo "Unexpected input: $mode\n single or paired"
    exit 1
fi
date