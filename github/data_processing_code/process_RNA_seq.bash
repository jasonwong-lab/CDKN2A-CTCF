#!/usr/bin/env bash
set -euo pipefail

root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
source "$root/data_processing_code/metadata.bash"
manifest=${1:-"$root/data_processing_code/samples.tsv"}
data_dir="$root/data/RNA_seq"
genome="$root/reference/hg38_STAR"
gtf="$root/reference/hg38.ensGene.gtf"
alignment="$data_dir/alignment"

load_samples RNA_seq "$manifest" "$data_dir"
[[ -f $genome/Genome && -f $genome/SA && -f $genome/SAindex ]] || {
    echo "Missing STAR genome index in $genome" >&2; exit 1;
}
[[ -f $gtf ]] || { echo "Missing original RNA-seq annotation: $gtf" >&2; exit 1; }
mkdir -p "$alignment"

for i in "${!SAMPLE_IDS[@]}"; do
    id=${SAMPLE_IDS[$i]}
    echo "Aligning $id"
    conda run --no-capture-output -n STAR_env STAR --runThreadN 12 \
        --genomeDir "$genome" \
        --readFilesIn "${SAMPLE_READ1[$i]}" "${SAMPLE_READ2[$i]}" \
        --readFilesCommand zcat --sjdbGTFfile "$gtf" --outFilterType BySJout \
        --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "$alignment/$id"
done
