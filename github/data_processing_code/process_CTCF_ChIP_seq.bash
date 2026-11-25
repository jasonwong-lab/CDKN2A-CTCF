#!/usr/bin/env bash
set -euo pipefail

root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
source "$root/data_processing_code/metadata.bash"
manifest=${1:-"$root/data_processing_code/samples.tsv"}
data_dir="$root/data/CTCF_ChIP"
reference="$root/reference/hg38_bwa/hg38.fa"
alignment="$data_dir/alignment"
markdup="$data_dir/markdup"
bigwigs="$data_dir/bigwigs"

load_samples CTCF_ChIP "$manifest" "$data_dir"
[[ -f $reference ]] || { echo "Missing reference FASTA: $reference" >&2; exit 1; }
input_found=0
for id in "${SAMPLE_IDS[@]}"; do
    [[ $id != CTCF_Input ]] || input_found=1
done
(( input_found )) || { echo "Manifest must contain a CTCF_Input sample" >&2; exit 1; }
mkdir -p "$alignment" "$markdup" "$bigwigs/chr9"

# Concatenating gzip members preserves lane order without decompressing them.
tmp1=""
tmp2=""
cleanup() {
    [[ -z $tmp1 ]] || rm -f -- "$tmp1"
    [[ -z $tmp2 ]] || rm -f -- "$tmp2"
}
trap cleanup EXIT

for i in "${!SAMPLE_IDS[@]}"; do
    id=${SAMPLE_IDS[$i]}
    IFS=',' read -r -a r1 <<< "${SAMPLE_READ1[$i]}"
    IFS=',' read -r -a r2 <<< "${SAMPLE_READ2[$i]}"
    read1=${r1[0]}
    read2=${r2[0]}
    if (( ${#r1[@]} > 1 )); then
        tmp1=$(mktemp "$alignment/.read1.XXXXXX.gz")
        cat -- "${r1[@]}" > "$tmp1"
        read1=$tmp1
        if [[ $read2 != - ]]; then
            tmp2=$(mktemp "$alignment/.read2.XXXXXX.gz")
            cat -- "${r2[@]}" > "$tmp2"
            read2=$tmp2
        fi
    fi
    reads=("$read1")
    [[ $read2 == - ]] || reads+=("$read2")
    echo "Aligning $id"
    conda run --no-capture-output -n bwa_env bwa-mem2 mem -t 12 "$reference" "${reads[@]}" |
        conda run --no-capture-output -n samtools_env samtools sort -@ 10 -o "$alignment/$id.bam" -
    conda run --no-capture-output -n samtools_env samtools index -@ 10 "$alignment/$id.bam"
    cleanup
    tmp1=""
    tmp2=""

    markdup_input="$alignment/$id.bam"
    if [[ $id == CTCF_Input ]]; then
        markdup_input="$alignment/${id}_withRG.bam"
        conda run --no-capture-output -n samtools_env samtools addreplacerg \
            -r '@RG\tID:CTCF_Input\tSM:CTCF_Input\tPL:BGI\tLB:Input_Lib1' \
            -o "$markdup_input" "$alignment/$id.bam" -@ 12
        conda run --no-capture-output -n samtools_env samtools index -@ 12 "$markdup_input"
    fi
    conda run --no-capture-output -n picard_env picard MarkDuplicates \
        I="$markdup_input" O="$markdup/$id.markdup.bam" M="$markdup/$id.markdup.metrics.txt" \
        ASSUME_SORTED=true MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=50000 \
        MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 SORTING_COLLECTION_SIZE_RATIO=0.25 \
        TAG_DUPLICATE_SET_MEMBERS=false REMOVE_SEQUENCING_DUPLICATES=false \
        TAGGING_POLICY=DontTag REMOVE_DUPLICATES=false \
        DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES \
        PROGRAM_RECORD_ID=MarkDuplicates PROGRAM_GROUP_NAME=MarkDuplicates \
        VALIDATION_STRINGENCY=STRICT
    if [[ $id == CTCF_Input ]]; then
        conda run --no-capture-output -n samtools_env samtools view -b -F 1024 -@ 12 "$markdup/$id.markdup.bam" |
            conda run --no-capture-output -n samtools_env samtools sort -@ 12 -O BAM -o "$markdup/$id.remove_dup.bam" -
        conda run --no-capture-output -n samtools_env samtools index -@ 12 "$markdup/$id.remove_dup.bam"
    fi
done

for region in genome chr9; do
    out_dir=$bigwigs
    region_args=()
    if [[ $region == chr9 ]]; then
        out_dir="$bigwigs/chr9"
        region_args=(--region chr9)
    fi
    for id in "${SAMPLE_IDS[@]}"; do
        conda run --no-capture-output -n deeptools_env bamCoverage \
            -b "$alignment/$id.bam" -o "$out_dir/$id.rpgc.bw" \
            --binSize 2 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 \
            --centerReads --ignoreDuplicates --minMappingQuality 30 \
            "${region_args[@]}" --numberOfProcessors 12
    done
    for id in "${SAMPLE_IDS[@]}"; do
        [[ $id == CTCF_Input ]] && continue
        conda run --no-capture-output -n deeptools_env bigwigCompare \
            -b1 "$out_dir/$id.rpgc.bw" -b2 "$out_dir/CTCF_Input.rpgc.bw" \
            -o "$out_dir/$id.ChIP_vs_Input.log2ratio.bw" --operation log2 \
            --pseudocount 1 --binSize 2 --numberOfProcessors 12
        if [[ $region == chr9 ]]; then
            conda run --no-capture-output -n deeptools_env bigwigCompare \
                -b1 "$out_dir/$id.rpgc.bw" -b2 "$out_dir/CTCF_Input.rpgc.bw" \
                -o "$out_dir/$id.ChIP_vs_Input.ratio.bw" --operation ratio \
                --pseudocount 1 --binSize 2 --numberOfProcessors 12
        fi
    done
done
