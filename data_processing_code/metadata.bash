#!/usr/bin/env bash

# Paths in the manifest are relative to the
# assay's download directory, not the current working directory.
load_samples() {
    local assay_filter=$1 manifest=$2 data_dir=$3 require_fastq=${4:-true}
    local line assay sample read1 read2 extra path resolved header
    local -a paths
    local -A seen=()
    SAMPLE_IDS=()
    SAMPLE_READ1=()
    SAMPLE_READ2=()

    [[ -f $manifest ]] || { echo "Missing sample manifest: $manifest" >&2; return 1; }
    IFS= read -r header < "$manifest" || { echo "Empty sample manifest: $manifest" >&2; return 1; }
    header=${header%$'\r'}
    [[ $header == $'assay\tsample\tread1\tread2' ]] || {
        echo "Expected manifest header: assay<TAB>sample<TAB>read1<TAB>read2" >&2
        return 1
    }

    while IFS= read -r line || [[ -n $line ]]; do
        line=${line%$'\r'}
        [[ -z $line || $line == \#* || $line == "$header" ]] && continue
        [[ $line != $'\t'* && $line != *$'\t' && $line != *$'\t\t'* ]] || {
            echo "Manifest fields must not be empty: $line" >&2; return 1;
        }
        IFS=$'\t' read -r assay sample read1 read2 extra <<< "$line"
        [[ $assay == CTCF_ChIP || $assay == RNA_seq ]] || {
            echo "Invalid assay in manifest: $line" >&2; return 1;
        }
        [[ $assay == "$assay_filter" ]] || continue
        [[ $sample =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ && $read1 && $read2 && ! $extra ]] || {
            echo "Invalid manifest row: $line" >&2; return 1;
        }
        [[ ! ${seen[$sample]+present} ]] || { echo "Duplicate sample: $sample" >&2; return 1; }
        seen[$sample]=1
        [[ $assay != RNA_seq || $read2 != - ]] || {
            echo "RNA-seq sample $sample requires paired reads" >&2; return 1;
        }
        [[ $read1 != - && ( $sample != CTCF_Input || $read2 != - ) ]] || {
            echo "Missing required reads for $sample" >&2; return 1;
        }
        local -a resolved_reads=()
        for resolved in "$read1" "$read2"; do
            if [[ $resolved == - ]]; then
                resolved_reads+=("-")
                continue
            fi
            [[ $resolved != ,* && $resolved != *, && $resolved != *,,* ]] || {
                echo "Invalid FASTQ list for $sample: $resolved" >&2; return 1;
            }
            IFS=',' read -r -a paths <<< "$resolved"
            local joined=""
            for path in "${paths[@]}"; do
                [[ $path =~ ^[A-Za-z0-9_./-]+\.gz$ && $path != /* && /$path/ != */../* ]] || {
                    echo "Invalid gzipped FASTQ path for $sample: $path" >&2; return 1;
                }
                [[ $require_fastq == false || -f $data_dir/$path ]] || {
                    echo "Missing gzipped FASTQ for $sample: $data_dir/$path" >&2; return 1;
                }
                joined+="${joined:+,}$data_dir/$path"
            done
            resolved_reads+=("$joined")
        done
        if [[ $read2 != - ]]; then
            local -a mates1 mates2
            IFS=',' read -r -a mates1 <<< "$read1"
            IFS=',' read -r -a mates2 <<< "$read2"
            [[ ${#mates1[@]} -eq ${#mates2[@]} ]] || {
                echo "Unequal lane counts for $sample" >&2; return 1;
            }
        fi
        SAMPLE_IDS+=("$sample")
        SAMPLE_READ1+=("${resolved_reads[0]}")
        SAMPLE_READ2+=("${resolved_reads[1]}")
    done < "$manifest"
    [[ ${#SAMPLE_IDS[@]} -gt 0 ]] || {
        echo "No $assay samples in $manifest" >&2; return 1;
    }
}
