#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_assembly_polishing.sh
# Description: Unified assembly polishing module combining long-read (Arrow/Racon/Medaka) and short-read (Polypolish) polishing

set -euo pipefail

# Function to display usage information
show_usage() {
    cat << EOF
Usage: $0 <script_dir> <logs_dir> <log_name> <utils_file> <apptainer_images_dir> <input_file>
          <bam_file> <output_dir> <sample_name> <sequencing_type> <threads>
          <genome_assembly_dir> <reproducibility_mode> <results_integration_dir> <version>
          <short_reads_r1> <short_reads_r2>
EOF
    exit 1
}

# Validate input parameters
[[ $# -eq 17 ]] || show_usage

# Constants from command line arguments
readonly SCRIPT_DIR="$1"
readonly LOGS_DIR="$2"
readonly LOG_NAME="$3"
readonly UTILS_FILE="$4"
readonly APPTAINER_DIR="$5"
readonly INPUT_FILE="$6"
readonly BAM_FILE="$7"
readonly OUTPUT_DIR="$8"
readonly SAMPLE_NAME="$9"
readonly SEQUENCING_TYPE="${10}"
readonly INITIAL_THREADS="${11}"
readonly GENOME_ASSEMBLY_DIR="${12}"
readonly REPRODUCIBILITY_MODE="${13}"
readonly RESULTS_INTEGRATION_DIR="${14}"
readonly VERSION="${15}"
readonly SHORT_READS_R1="${16}"
readonly SHORT_READS_R2="${17}"

# Set threads based on reproducibility mode
if [[ "${REPRODUCIBILITY_MODE}" == "deterministic" ]]; then
    readonly THREADS=1
else
    readonly THREADS="${INITIAL_THREADS}"
fi

# Polishing configuration
readonly POLISHING_ROUNDS=2

# Determine long-read polishing tool
determine_longread_polisher() {
    if [[ "$SEQUENCING_TYPE" == *"pacbio"* ]]; then
        if [[ -n "$BAM_FILE" && "$BAM_FILE" != "not_available" && -f "$BAM_FILE" ]]; then
            echo "arrow"
        else
            echo "racon"
        fi
    elif [[ "$SEQUENCING_TYPE" == *"nano"* ]]; then
        echo "medaka"
    else
        echo "none"
    fi
}

readonly LONGREAD_POLISHER=$(determine_longread_polisher)

# Derived constants
readonly POLISHING_OUTPUT_DIR="$OUTPUT_DIR/assembly_polishing_results"
readonly LONGREAD_OUTPUT_DIR="$POLISHING_OUTPUT_DIR/longread_polishing"
readonly SHORTREAD_OUTPUT_DIR="$POLISHING_OUTPUT_DIR/shortread_polishing"
readonly FIXED_BAM_FILE="$LONGREAD_OUTPUT_DIR/fixed_$(basename "${BAM_FILE:-dummy.bam}")"

# Source utility functions
source "$UTILS_FILE"

log "$LOGS_DIR" "$LOG_NAME" "Starting unified assembly polishing module"
log "$LOGS_DIR" "$LOG_NAME" "Sequencing type: $SEQUENCING_TYPE"
log "$LOGS_DIR" "$LOG_NAME" "Long-read polisher: $LONGREAD_POLISHER"
log "$LOGS_DIR" "$LOG_NAME" "Threads: $THREADS (mode: $REPRODUCIBILITY_MODE)"

# Determine if short-read polishing is available
has_short_reads() {
    [[ -n "$SHORT_READS_R1" && "$SHORT_READS_R1" != "not_provided" && -f "$SHORT_READS_R1" ]] && \
    [[ -n "$SHORT_READS_R2" && "$SHORT_READS_R2" != "not_provided" && -f "$SHORT_READS_R2" ]]
}

if has_short_reads; then
    log "$LOGS_DIR" "$LOG_NAME" "Short-read polishing: Polypolish (R1: $(basename "$SHORT_READS_R1"), R2: $(basename "$SHORT_READS_R2"))"
else
    log "$LOGS_DIR" "$LOG_NAME" "Short-read polishing: Skipped (no paired-end short reads provided)"
fi

# Create required directories
create_directory "$POLISHING_OUTPUT_DIR"
create_directory "$LONGREAD_OUTPUT_DIR"
if has_short_reads; then
    create_directory "$SHORTREAD_OUTPUT_DIR"
fi

# Find required SIF file and assembly file
readonly STRAINCASCADE_QC_SIF=$(find_apptainer_sif_file "$APPTAINER_DIR" 'straincascade_assembly_qc_refinement*.sif')
analysis_assembly_file=$(find_analysis_assembly_file "$GENOME_ASSEMBLY_DIR")

# Validate prerequisites
if [[ -z "$analysis_assembly_file" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "Error: No assembly files found. Skipping polishing."
    exit 0
fi

log "$LOGS_DIR" "$LOG_NAME" "Found assembly to polish: $(basename "$analysis_assembly_file")"

# Validate input files
if [[ ! -f "$INPUT_FILE" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "Error: Input reads file not found: $INPUT_FILE. Skipping polishing."
    exit 0
fi

# Entropy source binding for deterministic reproducibility only
# Using /dev/zero (infinite stream) instead of a finite file to prevent
# exhaustion-related hangs during long-running polishing stages
ENTROPY_ARGS=()
if [[ "${REPRODUCIBILITY_MODE}" == "deterministic" ]]; then
    ENTROPY_ARGS=(--bind /dev/zero:/dev/random --bind /dev/zero:/dev/urandom)
    log "$LOGS_DIR" "$LOG_NAME" "Deterministic mode: binding /dev/zero as entropy source"
fi

# Track current assembly through polishing stages
current_assembly="$analysis_assembly_file"

# =============================================================================
# LONG-READ POLISHING STAGE
# =============================================================================

if [[ "$LONGREAD_POLISHER" != "none" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "Starting long-read polishing with $LONGREAD_POLISHER ($POLISHING_ROUNDS rounds)"
    
    # Arrow-specific BAM preparation
    if [[ "$LONGREAD_POLISHER" == "arrow" ]]; then
        log "$LOGS_DIR" "$LOG_NAME" "Preparing BAM file for Arrow polishing"
        
        bam_basename=$(basename "$BAM_FILE")
        
        output=$(apptainer exec \
            --bind "$(dirname "$BAM_FILE")":/mnt/input \
            --bind "$LONGREAD_OUTPUT_DIR":/mnt/output \
            ${ENTROPY_ARGS[@]+"${ENTROPY_ARGS[@]}"} \
            "$STRAINCASCADE_QC_SIF" \
            /bin/bash -c "
            new_rg_id=\$(printf '%08x' \$RANDOM) && \
            source /opt/conda/etc/profile.d/conda.sh && \
            conda activate tools_env && \
            samtools addreplacerg \
                -r \"ID:\${new_rg_id}\" \
                -r \"PL:PACBIO\" \
                -r \"PU:$SAMPLE_NAME\" \
                -r \"SM:$SAMPLE_NAME\" \
                -o '/mnt/output/fixed_${bam_basename}' \
                '/mnt/input/${bam_basename}' && \
            samtools index '/mnt/output/fixed_${bam_basename}' && \
            samtools view -H '/mnt/output/fixed_${bam_basename}'
            " 2>&1) || true
        
        if echo "$output" | grep -q "This does not appear to be a valid PacBio BAM file"; then
            log "$LOGS_DIR" "$LOG_NAME" "BAM file not compatible with Arrow. Falling back to Racon polishing."
            LONGREAD_POLISHER_ACTUAL="racon"
        elif [[ ! -f "$FIXED_BAM_FILE" ]]; then
            log "$LOGS_DIR" "$LOG_NAME" "Warning: BAM preparation failed. Falling back to Racon polishing."
            LONGREAD_POLISHER_ACTUAL="racon"
        else
            LONGREAD_POLISHER_ACTUAL="arrow"
        fi
    else
        LONGREAD_POLISHER_ACTUAL="$LONGREAD_POLISHER"
    fi
    
    # Perform polishing rounds
    for ((round=1; round<=POLISHING_ROUNDS; round++)); do
        round_dir="$LONGREAD_OUTPUT_DIR/round_${round}"
        consensus_file="$round_dir/consensus.fasta"
        create_directory "$round_dir"
        
        log "$LOGS_DIR" "$LOG_NAME" "Long-read polishing round $round/$POLISHING_ROUNDS with $LONGREAD_POLISHER_ACTUAL"
        
        case "$LONGREAD_POLISHER_ACTUAL" in
            "arrow")
                apptainer exec \
                    --bind "$(dirname "$current_assembly")":/mnt/input_assembly \
                    --bind "$round_dir":/mnt/output \
                    --bind "$(dirname "$FIXED_BAM_FILE")":/mnt/bam_file_dir \
                    ${ENTROPY_ARGS[@]+"${ENTROPY_ARGS[@]}"} \
                    "$STRAINCASCADE_QC_SIF" \
                    bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                            conda activate genomicconsensus_env && \
                            variantCaller --algorithm=arrow \
                                -j $THREADS \
                                -r /mnt/input_assembly/$(basename "$current_assembly") \
                                -o /mnt/output/consensus.fasta \
                                /mnt/bam_file_dir/$(basename "$FIXED_BAM_FILE") && \
                            [[ -f /mnt/output/consensus.fasta ]]" || {
                    log "$LOGS_DIR" "$LOG_NAME" "Error: Arrow polishing round $round failed. Skipping long-read polishing."
                    break
                }
                ;;
            
            "racon")
                # Racon requires alignment via minimap2
                # Use appropriate minimap2 preset based on PacBio data type
                if [[ "$SEQUENCING_TYPE" == "pacbio-hifi" ]]; then
                    racon_minimap2_preset="map-hifi"
                else
                    racon_minimap2_preset="map-pb"
                fi
                log "$LOGS_DIR" "$LOG_NAME" "Aligning reads with minimap2 (preset: $racon_minimap2_preset) for Racon"
                
                apptainer exec \
                    --bind "$(dirname "$INPUT_FILE")":/mnt/input_reads \
                    --bind "$(dirname "$current_assembly")":/mnt/input_assembly \
                    --bind "$round_dir":/mnt/output \
                    ${ENTROPY_ARGS[@]+"${ENTROPY_ARGS[@]}"} \
                    "$STRAINCASCADE_QC_SIF" \
                    bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                            conda activate racon_env && \
                            minimap2 -x $racon_minimap2_preset -t $THREADS \
                                /mnt/input_assembly/$(basename "$current_assembly") \
                                /mnt/input_reads/$(basename "$INPUT_FILE") > /mnt/output/overlaps.paf && \
                            racon -t $THREADS \
                                /mnt/input_reads/$(basename "$INPUT_FILE") \
                                /mnt/output/overlaps.paf \
                                /mnt/input_assembly/$(basename "$current_assembly") > /mnt/output/consensus.fasta && \
                            rm -f /mnt/output/overlaps.paf && \
                            [[ -f /mnt/output/consensus.fasta ]]" || {
                    log "$LOGS_DIR" "$LOG_NAME" "Error: Racon polishing round $round failed. Skipping long-read polishing."
                    break
                }
                ;;
            
            "medaka")
                # Medaka 2.0+ uses auto model selection by default, no --model needed
                # The -m parameter is optional and auto-detection usually works best
                log "$LOGS_DIR" "$LOG_NAME" "Using Medaka with automatic model detection"
                
                # Try with decreasing batch sizes
                medaka_success=false
                for batch_size in "" "100" "50"; do
                    batch_param=""
                    [[ -n "$batch_size" ]] && batch_param="-b $batch_size"
                    
                    log "$LOGS_DIR" "$LOG_NAME" "Attempting Medaka polishing${batch_size:+ with batch size $batch_size}"
                    
                    if apptainer exec \
                        --bind "$(dirname "$INPUT_FILE")":/mnt/input_file \
                        --bind "$(dirname "$current_assembly")":/mnt/input_assembly \
                        --bind "$round_dir":/mnt/output \
                        ${ENTROPY_ARGS[@]+"${ENTROPY_ARGS[@]}"} \
                        "$STRAINCASCADE_QC_SIF" \
                        bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                                conda activate medaka_env && \
                                export TF_FORCE_GPU_ALLOW_GROWTH=true && \
                                medaka_consensus \
                                    -i /mnt/input_file/$(basename "$INPUT_FILE") \
                                    -d /mnt/input_assembly/$(basename "$current_assembly") \
                                    -o /mnt/output \
                                    $batch_param \
                                    -t $THREADS && \
                                [[ -f /mnt/output/consensus.fasta ]]"; then
                        medaka_success=true
                        break
                    else
                        log "$LOGS_DIR" "$LOG_NAME" "Medaka failed with batch size ${batch_size:-default}. Retrying..."
                        rm -f "$round_dir"/*
                    fi
                done
                
                if [[ "$medaka_success" != "true" ]]; then
                    log "$LOGS_DIR" "$LOG_NAME" "Error: Medaka polishing round $round failed with all batch sizes. Skipping long-read polishing."
                    break
                fi
                ;;
        esac
        
        # Verify consensus was produced
        if [[ -f "$consensus_file" && -s "$consensus_file" ]]; then
            current_assembly="$consensus_file"
            log "$LOGS_DIR" "$LOG_NAME" "Long-read polishing round $round completed successfully"
        else
            log "$LOGS_DIR" "$LOG_NAME" "Warning: No consensus produced in round $round. Stopping long-read polishing."
            break
        fi
    done
    
    # Update assembly file with long-read polished version
    if [[ "$current_assembly" != "$analysis_assembly_file" ]]; then
        log "$LOGS_DIR" "$LOG_NAME" "Updating assembly with long-read polished version"
        cp "$current_assembly" "$analysis_assembly_file"
        # Update reference to point to the actual file
        analysis_assembly_file="$analysis_assembly_file"
    fi
    
    # Cleanup Arrow temporary files
    if [[ "$LONGREAD_POLISHER_ACTUAL" == "arrow" && -f "$FIXED_BAM_FILE" ]]; then
        log "$LOGS_DIR" "$LOG_NAME" "Cleaning up Arrow temporary files"
        rm -f "$FIXED_BAM_FILE" "${FIXED_BAM_FILE}.bai" "${FIXED_BAM_FILE}.pbi"
    fi
else
    log "$LOGS_DIR" "$LOG_NAME" "Long-read polishing: Skipped (no compatible sequencing type detected)"
fi

# =============================================================================
# SHORT-READ POLISHING STAGE (Polypolish)
# =============================================================================

if has_short_reads; then
    log "$LOGS_DIR" "$LOG_NAME" "Starting short-read polishing with Polypolish"
    
    # Refresh assembly reference after long-read polishing
    current_assembly=$(find_analysis_assembly_file "$GENOME_ASSEMBLY_DIR")
    if [[ -z "$current_assembly" ]]; then
        log "$LOGS_DIR" "$LOG_NAME" "Error: Assembly file not found after long-read polishing. Skipping short-read polishing."
    else
        readonly ASSEMBLY_BASENAME=$(basename "$current_assembly" .fasta)
        readonly POLISHED_ASSEMBLY="${ASSEMBLY_BASENAME}_polished.fasta"
        readonly ALIGNMENTS_R1="alignments_R1.sam"
        readonly ALIGNMENTS_R2="alignments_R2.sam"
        readonly FILTERED_R1="filtered_R1.sam"
        readonly FILTERED_R2="filtered_R2.sam"
        
        # Step 1: Index the assembly with BWA
        log "$LOGS_DIR" "$LOG_NAME" "Polypolish Step 1/4: Indexing assembly with BWA"
        apptainer exec \
            --bind "$GENOME_ASSEMBLY_DIR":/mnt/assembly \
            --bind "$SHORTREAD_OUTPUT_DIR":/mnt/output \
            "$STRAINCASCADE_QC_SIF" \
            bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                     conda activate tools_env && \
                     cp /mnt/assembly/$(basename "$current_assembly") /mnt/output/ && \
                     bwa index /mnt/output/$(basename "$current_assembly")" || {
            log "$LOGS_DIR" "$LOG_NAME" "Error: BWA indexing failed. Skipping short-read polishing."
            exit 0
        }
        
        # Step 2: Align short reads (R1 and R2 in a single container invocation)
        log "$LOGS_DIR" "$LOG_NAME" "Polypolish Step 2/4: Aligning short reads with BWA mem"
        
        apptainer exec \
            --bind "$(dirname "$SHORT_READS_R1")":/mnt/short_r1 \
            --bind "$(dirname "$SHORT_READS_R2")":/mnt/short_r2 \
            --bind "$SHORTREAD_OUTPUT_DIR":/mnt/output \
            "$STRAINCASCADE_QC_SIF" \
            bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                     conda activate tools_env && \
                     bwa mem -t $THREADS -a /mnt/output/$(basename "$current_assembly") \
                     /mnt/short_r1/$(basename "$SHORT_READS_R1") > /mnt/output/$ALIGNMENTS_R1 && \
                     bwa mem -t $THREADS -a /mnt/output/$(basename "$current_assembly") \
                     /mnt/short_r2/$(basename "$SHORT_READS_R2") > /mnt/output/$ALIGNMENTS_R2" || {
            log "$LOGS_DIR" "$LOG_NAME" "Error: BWA alignment failed. Skipping short-read polishing."
            exit 0
        }
        
        log "$LOGS_DIR" "$LOG_NAME" "BWA alignment completed for both read pairs"
        
        # Step 3: Filter alignments with polypolish filter
        log "$LOGS_DIR" "$LOG_NAME" "Polypolish Step 3/4: Filtering alignments"
        apptainer exec \
            --bind "$SHORTREAD_OUTPUT_DIR":/mnt/output \
            "$STRAINCASCADE_QC_SIF" \
            bash -c "polypolish filter --in1 /mnt/output/$ALIGNMENTS_R1 --in2 /mnt/output/$ALIGNMENTS_R2 \
                     --out1 /mnt/output/$FILTERED_R1 --out2 /mnt/output/$FILTERED_R2" || {
            log "$LOGS_DIR" "$LOG_NAME" "Error: Polypolish filter failed. Skipping short-read polishing."
            exit 0
        }
        
        log "$LOGS_DIR" "$LOG_NAME" "Alignment filtering completed"
        
        # Step 4: Run Polypolish polishing
        log "$LOGS_DIR" "$LOG_NAME" "Polypolish Step 4/4: Running Polypolish polishing"
        apptainer exec \
            --bind "$SHORTREAD_OUTPUT_DIR":/mnt/output \
            "$STRAINCASCADE_QC_SIF" \
            bash -c "polypolish polish /mnt/output/$(basename "$current_assembly") \
                     /mnt/output/$FILTERED_R1 /mnt/output/$FILTERED_R2 > /mnt/output/$POLISHED_ASSEMBLY" || {
            log "$LOGS_DIR" "$LOG_NAME" "Error: Polypolish polishing failed. Skipping short-read polishing."
            exit 0
        }
        
        # Check if polished assembly was created and is non-empty
        if [[ ! -s "$SHORTREAD_OUTPUT_DIR/$POLISHED_ASSEMBLY" ]]; then
            log "$LOGS_DIR" "$LOG_NAME" "Error: Polypolish produced an empty assembly. Keeping current assembly."
        else
            # Normalize contig headers: Polypolish adds " polypolish" suffix to contig names
            # Remove this suffix to maintain consistent naming throughout the pipeline
            log "$LOGS_DIR" "$LOG_NAME" "Normalizing contig headers (removing polypolish suffix)"
            sed -i.bak 's/ polypolish$//' "$SHORTREAD_OUTPUT_DIR/$POLISHED_ASSEMBLY"
            rm -f "$SHORTREAD_OUTPUT_DIR/${POLISHED_ASSEMBLY}.bak"
            
            # Calculate polishing statistics
            readonly ORIGINAL_SIZE=$(grep -v ">" "$current_assembly" | tr -d '\n' | wc -c)
            readonly POLISHED_SIZE=$(grep -v ">" "$SHORTREAD_OUTPUT_DIR/$POLISHED_ASSEMBLY" | tr -d '\n' | wc -c)
            log "$LOGS_DIR" "$LOG_NAME" "Assembly size before short-read polishing: $ORIGINAL_SIZE bp"
            log "$LOGS_DIR" "$LOG_NAME" "Assembly size after short-read polishing: $POLISHED_SIZE bp"
            
            # Replace the current assembly with polished version (no backup in main dir)
            # The polishing subdirectory retains the history if needed
            readonly FINAL_POLISHED_NAME="$(basename "${current_assembly%.fasta}").fasta"
            cp "$SHORTREAD_OUTPUT_DIR/$POLISHED_ASSEMBLY" "$GENOME_ASSEMBLY_DIR/$FINAL_POLISHED_NAME"
            log "$LOGS_DIR" "$LOG_NAME" "Polished assembly saved as: $FINAL_POLISHED_NAME"
            
            # Also copy to results integration directory
            cp "$SHORTREAD_OUTPUT_DIR/$POLISHED_ASSEMBLY" "$RESULTS_INTEGRATION_DIR/${SAMPLE_NAME}_v${VERSION}_final_polished.fasta"
            log "$LOGS_DIR" "$LOG_NAME" "Polished assembly copied to results integration directory"
            
            log "$LOGS_DIR" "$LOG_NAME" "Polypolish short-read polishing completed successfully"
        fi
        
        # Clean up intermediate SAM files
        log "$LOGS_DIR" "$LOG_NAME" "Cleaning up intermediate alignment files"
        rm -f "$SHORTREAD_OUTPUT_DIR/$ALIGNMENTS_R1" "$SHORTREAD_OUTPUT_DIR/$ALIGNMENTS_R2"
        rm -f "$SHORTREAD_OUTPUT_DIR/$FILTERED_R1" "$SHORTREAD_OUTPUT_DIR/$FILTERED_R2"
        rm -f "$SHORTREAD_OUTPUT_DIR/$(basename "$current_assembly")".* # Remove BWA index files
        rm -f "$SHORTREAD_OUTPUT_DIR/$(basename "$current_assembly")"   # Remove copied assembly
    fi
fi

# =============================================================================
# POST-POLISHING EVALUATION (QUAST statistics update)
# =============================================================================

# Check if any polishing was performed by looking for a polished assembly
polished_assembly_file=$(find_analysis_assembly_file "$GENOME_ASSEMBLY_DIR")

# Track if polishing occurred (assembly was modified)
polishing_performed=false
if [[ "$LONGREAD_POLISHER" != "none" ]] || has_short_reads; then
    polishing_performed=true
fi

if [[ "$polishing_performed" == "true" && -n "$polished_assembly_file" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "Starting post-polishing evaluation with QUAST"
    
    readonly POST_POLISH_EVAL_DIR="$POLISHING_OUTPUT_DIR/post_polish_evaluation"
    readonly QS_FILES_DIR="$RESULTS_INTEGRATION_DIR/qs_files"
    readonly R_SCRIPT_DIR="$SCRIPT_DIR/R_scripts"
    create_directory "$POST_POLISH_EVAL_DIR"
    create_directory "$QS_FILES_DIR"
    
    # Find R SIF file for processing
    r_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'r_4.4.1*.sif')
    
    # Backup existing chosen_assembly.tsv if it exists
    readonly CHOSEN_ASSEMBLY_TSV="$GENOME_ASSEMBLY_DIR/chosen_assembly.tsv"
    if [[ -f "$CHOSEN_ASSEMBLY_TSV" ]]; then
        readonly PRE_POLISH_BACKUP="$GENOME_ASSEMBLY_DIR/pre_polish_chosen_assembly.tsv"
        cp "$CHOSEN_ASSEMBLY_TSV" "$PRE_POLISH_BACKUP"
        log "$LOGS_DIR" "$LOG_NAME" "Backed up chosen_assembly.tsv to pre_polish_chosen_assembly.tsv"
    fi
    
    # Run QUAST on the polished assembly
    log "$LOGS_DIR" "$LOG_NAME" "Running QUAST on polished assembly: $(basename "$polished_assembly_file")"
    
    readonly QUAST_OUTPUT_DIR="$POST_POLISH_EVAL_DIR/quast_output"
    create_directory "$QUAST_OUTPUT_DIR"
    
    apptainer exec \
        --bind "$GENOME_ASSEMBLY_DIR":/data \
        --bind "$QUAST_OUTPUT_DIR":/output \
        "$STRAINCASCADE_QC_SIF" \
        bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                 conda activate quast_env && \
                 quast.py '/data/$(basename "$polished_assembly_file")' \
                     -o /output \
                     -t $THREADS \
                     --min-contig 0" || {
        log "$LOGS_DIR" "$LOG_NAME" "Warning: QUAST post-polish evaluation failed. Skipping statistics update."
    }
    
    # Update chosen_assembly.tsv with polished assembly stats
    # The format must match exactly what the Python selection scripts produce:
    # - Line 1: Header row from QUAST transposed_report.tsv
    # - Line 2: Data row with the polished assembly stats
    # - Optionally appended: selection note (after newline)
    if [[ -f "$QUAST_OUTPUT_DIR/transposed_report.tsv" ]]; then
        # Copy QUAST transposed_report.tsv as chosen_assembly.tsv (maintains exact format)
        cp "$QUAST_OUTPUT_DIR/transposed_report.tsv" "$CHOSEN_ASSEMBLY_TSV"
        
        # Append the selection note (same format as eval3)
        echo -e "\n$(basename "$polished_assembly_file") is the final polished assembly for all downstream analysis" >> "$CHOSEN_ASSEMBLY_TSV"
        
        log "$LOGS_DIR" "$LOG_NAME" "Updated chosen_assembly.tsv with post-polishing statistics"
        
        # Also save polishing metadata to a separate file for reference
        {
            echo "Polishing Summary"
            echo "================"
            echo "Long-read polisher: ${LONGREAD_POLISHER}"
            echo "Short-read polisher: $(has_short_reads && echo "Polypolish" || echo "none")"
            echo "Polished assembly: $(basename "$polished_assembly_file")"
            echo "Pre-polish stats: pre_polish_chosen_assembly.tsv"
        } > "$GENOME_ASSEMBLY_DIR/polishing_summary.txt"
        
        # Replace pre-polish QUAST output with post-polish stats (final state is what matters)
        rm -rf "$GENOME_ASSEMBLY_DIR/detailed_quast_output_of_chosen_assembly"
        cp -r "$QUAST_OUTPUT_DIR" "$GENOME_ASSEMBLY_DIR/detailed_quast_output_of_chosen_assembly"
        
        # Process with R for qs_files
        if [[ -f "$r_sif" && -d "$R_SCRIPT_DIR" ]]; then
            log "$LOGS_DIR" "$LOG_NAME" "Processing polished assembly statistics with R"
            
            # Process QUAST results
            apptainer exec \
                --bind "$R_SCRIPT_DIR":/mnt/r_script_dir \
                --bind "$GENOME_ASSEMBLY_DIR":/mnt/input \
                --bind "$QS_FILES_DIR":/mnt/output \
                "$r_sif" \
                Rscript "/mnt/r_script_dir/R_process_quast.R" \
                --output_dir "/mnt/output" \
                --tsv "/mnt/input/chosen_assembly.tsv" \
                --version "$VERSION" || {
                log "$LOGS_DIR" "$LOG_NAME" "Warning: R processing of QUAST results failed"
            }
            
            # Process selected assembly
            apptainer exec \
                --bind "$R_SCRIPT_DIR":/mnt/r_script_dir \
                --bind "$GENOME_ASSEMBLY_DIR":/mnt/input \
                --bind "$QS_FILES_DIR":/mnt/output \
                "$r_sif" \
                Rscript "/mnt/r_script_dir/R_process_selected_assembly.R" \
                --output_dir "/mnt/output" \
                --fasta "/mnt/input/$(basename "$polished_assembly_file")" \
                --version "$VERSION" || {
                log "$LOGS_DIR" "$LOG_NAME" "Warning: R processing of selected assembly failed"
            }
            
            log "$LOGS_DIR" "$LOG_NAME" "R processing of polished assembly completed"
        fi
    else
        log "$LOGS_DIR" "$LOG_NAME" "Warning: QUAST report not found. chosen_assembly.tsv not updated."
    fi
    
    log "$LOGS_DIR" "$LOG_NAME" "Post-polishing evaluation completed"
else
    log "$LOGS_DIR" "$LOG_NAME" "No polishing was performed. Skipping post-polishing evaluation."
fi

# =============================================================================
# STEP: NORMALIZE FINAL ASSEMBLY CONTIG NAMES
# =============================================================================
# This is the critical step that ensures all downstream analysis tools use
# consistent, standardized contig names (contig_1, contig_2, etc. sorted by size).
# This MUST happen AFTER polishing but BEFORE any analysis tools run.

log "$LOGS_DIR" "$LOG_NAME" "Starting final assembly contig name normalization"

# Find the final assembly (may or may not have been polished)
final_assembly_for_normalization=$(find_analysis_assembly_file "$GENOME_ASSEMBLY_DIR")

if [[ -n "$final_assembly_for_normalization" && -f "$final_assembly_for_normalization" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "Normalizing assembly: $(basename "$final_assembly_for_normalization")"
    
    # Ensure QS_FILES_DIR is set (may not be if polishing was skipped)
    if [[ -z "${QS_FILES_DIR:-}" ]]; then
        QS_FILES_DIR="$RESULTS_INTEGRATION_DIR/qs_files"
        create_directory "$QS_FILES_DIR"
    fi
    
    # Ensure R_SCRIPT_DIR is set
    if [[ -z "${R_SCRIPT_DIR:-}" ]]; then
        R_SCRIPT_DIR="$SCRIPT_DIR/R_scripts"
    fi
    
    # Ensure r_sif is set
    if [[ -z "${r_sif:-}" ]]; then
        r_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'r_4.4.1*.sif')
    fi
    
    # Check if Circlator results exist and need updating
    CIRCLATOR_QS=$(find "$QS_FILES_DIR" -name 'circlator_results.qs' -print -quit 2>/dev/null || true)
    circlator_arg=""
    if [[ -n "$CIRCLATOR_QS" && -f "$CIRCLATOR_QS" ]]; then
        circlator_arg="--circlator_qs /mnt/qs_files/$(basename "$CIRCLATOR_QS")"
        log "$LOGS_DIR" "$LOG_NAME" "Will update Circlator results with normalized contig names"
    fi
    
    # Run normalization
    # This will:
    # 1. Sort contigs by length (descending)
    # 2. Rename to contig_1, contig_2, etc.
    # 3. Overwrite the FASTA file in place
    # 4. Create contig_name_mapping.tsv for traceability
    # 5. Update Circlator results if they exist
    apptainer exec \
        --bind "${R_SCRIPT_DIR}:/mnt/r_script_dir" \
        --bind "$(dirname "$final_assembly_for_normalization")":/mnt/assembly \
        --bind "${QS_FILES_DIR}:/mnt/qs_files" \
        "$r_sif" \
        Rscript "/mnt/r_script_dir/R_normalize_final_assembly.R" \
        --fasta "/mnt/assembly/$(basename "$final_assembly_for_normalization")" \
        --output_dir "/mnt/qs_files" \
        $circlator_arg || {
        log "$LOGS_DIR" "$LOG_NAME" "Error: Assembly normalization failed. This may cause issues with downstream tools."
        # Don't exit - try to continue with original names
    }
    
    # Verify normalization succeeded by checking for mapping file
    if [[ -f "${QS_FILES_DIR}/contig_name_mapping.tsv" ]]; then
        contig_count=$(tail -n +2 "${QS_FILES_DIR}/contig_name_mapping.tsv" | wc -l)
        log "$LOGS_DIR" "$LOG_NAME" "Assembly normalized successfully: $contig_count contigs renamed (contig_1 to contig_$contig_count)"
        log "$LOGS_DIR" "$LOG_NAME" "Contig name mapping saved to: ${QS_FILES_DIR}/contig_name_mapping.tsv"
        
        # Also copy mapping file to genome assembly directory for easy access
        cp "${QS_FILES_DIR}/contig_name_mapping.tsv" "${GENOME_ASSEMBLY_DIR}/contig_name_mapping.tsv"
    else
        log "$LOGS_DIR" "$LOG_NAME" "Warning: Contig name mapping file not created. Normalization may have failed."
    fi
    
    # Re-process selected assembly with normalized names
    log "$LOGS_DIR" "$LOG_NAME" "Re-processing selected assembly with normalized contig names"
    apptainer exec \
        --bind "${R_SCRIPT_DIR}:/mnt/r_script_dir" \
        --bind "$(dirname "$final_assembly_for_normalization")":/mnt/assembly \
        --bind "${QS_FILES_DIR}:/mnt/output" \
        "$r_sif" \
        Rscript "/mnt/r_script_dir/R_process_selected_assembly.R" \
        --fasta "/mnt/assembly/$(basename "$final_assembly_for_normalization")" \
        --output_dir "/mnt/output" \
        --version "$VERSION" || {
        log "$LOGS_DIR" "$LOG_NAME" "Warning: R processing of normalized assembly failed"
    }
    
    log "$LOGS_DIR" "$LOG_NAME" "Assembly normalization completed"
else
    log "$LOGS_DIR" "$LOG_NAME" "Warning: No final assembly found for normalization"
fi

# =============================================================================
# FINAL CLEANUP: Remove intermediate files from main assembly directory
# =============================================================================

log "$LOGS_DIR" "$LOG_NAME" "Cleaning up intermediate assembly files"

# Move intermediate files to a dedicated archive subdirectory
readonly INTERMEDIATE_ARCHIVE_DIR="$GENOME_ASSEMBLY_DIR/assembly_intermediates"
create_directory "$INTERMEDIATE_ARCHIVE_DIR"

# Identify and move non-final assembly files
# Keep only: *_final.fasta and essential metadata files
while IFS= read -r -d '' file; do
    filename=$(basename "$file")
    
    # Skip if it's the final assembly (contains _final but not _pre_)
    if [[ "$filename" == *"_final.fasta" && "$filename" != *"_pre_"* ]]; then
        continue
    fi
    
    # Skip non-assembly files (coverage, reports, etc.)
    if [[ "$filename" != *.fasta ]]; then
        continue
    fi
    
    # Move intermediate assemblies to archive
    log "$LOGS_DIR" "$LOG_NAME" "Archiving intermediate: $filename"
    mv "$file" "$INTERMEDIATE_ARCHIVE_DIR/"
done < <(find "$GENOME_ASSEMBLY_DIR" -maxdepth 1 -type f -name "*.fasta" -print0)

# Log final state
final_assembly_count=$(find "$GENOME_ASSEMBLY_DIR" -maxdepth 1 -type f -name "*_final.fasta" | wc -l)
archived_count=$(find "$INTERMEDIATE_ARCHIVE_DIR" -maxdepth 1 -type f -name "*.fasta" 2>/dev/null | wc -l)
log "$LOGS_DIR" "$LOG_NAME" "Final assemblies in main directory: $final_assembly_count"
log "$LOGS_DIR" "$LOG_NAME" "Intermediate assemblies archived: $archived_count"

log "$LOGS_DIR" "$LOG_NAME" "Unified assembly polishing module completed"
