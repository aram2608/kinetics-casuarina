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


source activate genomics

date 
input_directory=$1 
output_directory=$2
adapters=$3
mode=$4

# Added safety test to ensure all parameters are provided
if [ -z $input_directory ] || [ -z $output_directory ] || [ -z $adapters ] || [ -z $mode ]; then
    echo "Usage: $0 <input_directory> <output_directory> </path/to/adapters>"
    exit 1
fi

echo $input_directory will now be processed

# Creates output directory just in case it does not exist
mkdir -p $output_directory

# Argument check
if [ "$mode" == "single" ]; then
    trim_single_end
elif [ $mode == "paired" ]
    trim_paired_end
else
    echo "Unexpected input: $mode\n single or paired"
    exit 1
fi

# Single end function
function trim_single_end() {
    for fastq_file in $input_directory/*;do #loop through files of input directory
        if [[ $fastq_file == *fastq.gz ]]; then #tests for fastq file
            file_name=$(basename $fastq_file) #creates basename for files
            output_file=$output_directory/trimmed_$file_name

            echo now trimming $fastq_file -> $output_file; #prints file to be trimmed

            trimmomatic SE -phred33 $fastq_file $output_file \
                ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 \
                LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
                
            #an error handling step in case a file is not processed for some reason
            if  [ $? -ne 0 ]; then # $? is a special variable that stores the status of the previous command, -ne 0 means a non zero error was thrown
                echo "Error Trimming $fastq_file"
            fi
        fi
    done
}
date