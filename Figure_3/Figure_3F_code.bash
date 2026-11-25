#!/usr/bin/env bash
set -euo pipefail

root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
source "$root/data_processing_code/metadata.bash"
manifest=${1:-"$root/data_processing_code/samples.tsv"}
data_dir="$root/data/CTCF_ChIP"
markdup="$data_dir/markdup"
reference="$root/reference/hg38_bwa/hg38.fa"
parser="$root/Figure_3/CDKN2B_parse_bam_readcount_ChIP.r"
output="$root/Figure_3/results"

load_samples CTCF_ChIP "$manifest" "$data_dir" false
[[ -f $reference && -f $reference.fai ]] || {
    echo "Missing reference FASTA or FASTA index: $reference" >&2; exit 1;
}
[[ -f $parser ]] || { echo "Missing R parser: $parser" >&2; exit 1; }

for sample in "${SAMPLE_IDS[@]}"; do
    if [[ $sample == CTCF_Input ]]; then
        bam="$markdup/$sample.remove_dup.bam"
        [[ -f $bam && ( -f $bam.bai || -f ${bam%.bam}.bai ) ]] || {
            echo "Missing duplicate-removed input BAM or index: $bam" >&2; exit 1;
        }
    else
        bam="$markdup/$sample.markdup.bam"
        [[ -f $bam ]] || { echo "Missing marked-duplicate BAM: $bam" >&2; exit 1; }
    fi
done

mkdir -p "$output/deduplicated"
for sample in "${SAMPLE_IDS[@]}"; do
    [[ $sample == CTCF_Input ]] && continue
    bam="$output/deduplicated/$sample.remove_dup.bam"
    conda run --no-capture-output -n samtools_env samtools view -b -F 1024 -@ 12 "$markdup/$sample.markdup.bam" |
        conda run --no-capture-output -n samtools_env samtools sort -@ 12 -o "$bam" -
    conda run --no-capture-output -n samtools_env samtools index -@ 12 "$bam"
done

for analysis in region position; do
    if [[ $analysis == region ]]; then
        region=chr9:22008758-22008777
    else
        region=chr9:22008768-22008768
    fi
    out_dir="$output/$analysis"
    mkdir -p "$out_dir"
    summary="$out_dir/bam_readcount_summary.tsv"
    : > "$summary"

    for sample in "${SAMPLE_IDS[@]}"; do
        if [[ $sample == CTCF_Input ]]; then
            bam="$markdup/$sample.remove_dup.bam"
        else
            bam="$output/deduplicated/$sample.remove_dup.bam"
        fi
        readcount="$out_dir/$sample.readcount"
        parsed="$out_dir/$sample.txt"
        echo "Counting $sample at $region"
        conda run --no-capture-output -n bam-readcount-env bam-readcount \
            -f "$reference" -q 20 -b 20 "$bam" "$region" > "$readcount"
        conda run --no-capture-output -n r_env Rscript "$parser" "$readcount" "$parsed"
        awk -F '\t' -v sample="$sample" '
            $1 == "chr9" && $2 == 22008768 { print $0 "\t" sample; found++ }
            END {
                if (found != 1) {
                    print "Expected one chr9:22008768 row for " sample ", found " (found + 0) > "/dev/stderr"
                    exit 1
                }
            }
        ' "$parsed" >> "$summary"
    done
    {
        printf 'chr\tpos\tref_base\tTotal_count\tA\tC\tG\tT\tSample\n'
        cat "$summary"
    } > "$out_dir/bam_readcount_summary+header.tsv"
done
