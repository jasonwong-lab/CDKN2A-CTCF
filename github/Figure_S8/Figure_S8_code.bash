#!/usr/bin/env bash
set -euo pipefail

root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
source "$root/data_processing_code/metadata.bash"
manifest=${1:-"$root/data_processing_code/samples.tsv"}
data_dir="$root/data/CTCF_ChIP"
bigwigs="$data_dir/bigwigs"
motifs="$root/annotations/Katainen_CBS_motifs_hg38.bed"
output="$root/Figure_S8/results"

load_samples CTCF_ChIP "$manifest" "$data_dir" false
[[ -f $motifs ]] || { echo "Missing motif annotation: $motifs" >&2; exit 1; }
samples=()
for sample in "${SAMPLE_IDS[@]}"; do
    [[ $sample == CTCF_Input ]] && continue
    bigwig="$bigwigs/$sample.ChIP_vs_Input.log2ratio.bw"
    [[ -f $bigwig ]] || { echo "Missing genome-wide log2-ratio bigWig: $bigwig" >&2; exit 1; }
    samples+=("$sample")
done
(( ${#samples[@]} > 0 )) || { echo "No non-input CTCF samples in $manifest" >&2; exit 1; }
mkdir -p "$output"

for sample in "${samples[@]}"; do
    id="$sample.ChIP_vs_Input"
    conda run --no-capture-output -n deeptools_env computeMatrix reference-point \
        --referencePoint center \
        -S "$bigwigs/$id.log2ratio.bw" \
        -R "$motifs" \
        -b 2000 \
        -a 2000 \
        --binSize 10 \
        --skipZeros \
        -o "$output/$id.Katainen_motif_center.mat.gz" \
        --numberOfProcessors 10
done

for sample in "${samples[@]}"; do
    id="$sample.ChIP_vs_Input"
    conda run --no-capture-output -n deeptools_env plotHeatmap \
        -m "$output/$id.Katainen_motif_center.mat.gz" \
        -out "$output/$id.CTCF_motif_center_heatmap.pdf" \
        --refPointLabel "CTCF center" \
        --regionsLabel "CTCF motif" \
        --heatmapHeight 8 --heatmapWidth 6
done
