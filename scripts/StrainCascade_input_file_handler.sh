#!/bin/bash

# StrainCascade_input_file_handler.sh - Version 1.2.0
# Author: Sebastian Bruno Ulrich Jordi

if [ "$#" -ne 8 ]; then
    echo "Usage: $0 <script_dir> <logs_dir> <apptainer_images_dir> <input_file> <output_dir> <sequencing_reads_main_abs> <genome_assembly_main_abs> <input_type>"
    exit 1
fi

script_dir=$1
logs_dir=$2
apptainer_images_dir=$3
input_file=$4
output_dir=$5
sequencing_reads_main_abs=$6
genome_assembly_main_abs=$7
input_type=$8

# Load utils from the script directory
utils_file="${script_dir}/utils.sh"
if [ -f "$utils_file" ]; then
  source "$utils_file"
else
  echo "Error: utils.sh not found in $script_dir"
  exit 1
fi

# Find the appropriate Apptainer image
readarray -t matching_files < <(find "$apptainer_images_dir" -name 'straincascade_genome_assembly*.sif' -print)

if [ ${#matching_files[@]} -eq 0 ]; then
    log "$logs_dir" "Input_file_handler.log"  "No matching .sif files found in $apptainer_images_dir."
    exit 1
elif [ ${#matching_files[@]} -gt 1 ]; then
    log "$logs_dir" "Input_file_handler.log" "Warning: Multiple matching .sif files found. Using the first match: ${matching_files[0]}"
fi

straincascade_genome_assembly=${matching_files[0]}

original_input="$input_file"
bam_file="not_available"

# Process input file
if [ -n "$input_file" ]; then
    dir=$(dirname "$input_file")
    base=$(basename "$input_file")
    
    # Handle both single and double extensions (like .fastq.gz)
    filename=$(basename "$base" .gz)  # Remove .gz if present
    extension="${filename##*.}"
    filename="${filename%.*}"
    
    if [[ "$base" == *.fastq.gz ]]; then
        extension="fastq.gz"
    fi

    if [ "$input_type" == "assembly" ]; then
        case "$extension" in
        fasta|fa|fna)
            # Copy the assembly file to the genome_assembly_main_abs directory with .fasta extension
            cp "$input_file" "$genome_assembly_main_abs/${filename}.fasta"
            input_file="$genome_assembly_main_abs/${filename}.fasta"
            ;;
        *)
            log "$logs_dir" "Input_file_handler.log" "Error: Unsupported assembly file extension. Expected .fasta, .fa, or .fna"
            exit 1
            ;;
        esac
    else
        case "$extension" in
        fasta|fa|fna)
            # Copy the file to the sequencing_reads_main_abs directory with .fasta extension
            cp "$input_file" "$sequencing_reads_main_abs/${filename}.fasta"
            input_file="$sequencing_reads_main_abs/${filename}.fasta"
            ;;
        fastq|fastq.gz|bam)
            if [ "$extension" = "bam" ]; then
                bam_file="$input_file"
                # Copy BAM file to sequencing_reads_main_abs
                cp "$input_file" "$sequencing_reads_main_abs/$base"
            else
                # Copy FASTQ file to sequencing_reads_main_abs
                cp "$input_file" "$sequencing_reads_main_abs/$base"
            fi
            
            # Convert to FASTA
            apptainer exec \
                --bind "$sequencing_reads_main_abs":/mnt/input \
                "$straincascade_genome_assembly" \
                /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                            conda activate tools_env && \
                            samtools fasta /mnt/input/$base > /mnt/input/${filename}.fasta"
            
            input_file="$sequencing_reads_main_abs/${filename}.fasta"
            
            # Extract BAM header information if it's a BAM file
            if [ "$extension" = "bam" ]; then
                apptainer exec \
                    --bind "$sequencing_reads_main_abs":/mnt/input \
                    "$straincascade_genome_assembly" \
                    /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                                conda activate tools_env && \
                                (samtools view -H /mnt/input/$base | awk '
                                BEGIN {
                                    print \"BAM File: '$bam_file'\"
                                    print \"--------------------\"
                                }
                                {
                                    if (\$1 ~ /^@(RG|SQ|PG|HD)/) {
                                        for (i=1; i<=NF; i++) {
                                            tag = substr(\$i,1,3)
                                            value = substr(\$i, 4)
                                            if (tag ~ /^(ID|SM|LB|PL|PU|SN|LN|PN|CL|VN|SO|DS):/) {
                                                print tag \" \" value
                                            }
                                        }
                                    }
                                }' > /mnt/input/${filename}_bam_info.txt) || \
                                echo \"Error: Failed to extract BAM file information\" > /mnt/input/${filename}_bam_info.txt"
            fi
            ;;
        *)
            log "$logs_dir" "Input_file_handler.log" "Error: Unsupported file extension"
            exit 1
            ;;
        esac
    fi
else
    log "$logs_dir" "Input_file_handler.log" "Error: No input file identified"
    exit 1
fi

# Output original_input, input_file, and bam_file, separated by tabs
echo -e "$original_input\t$input_file\t$bam_file"