#!/usr/bin/env bash
set -euo pipefail

root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
source "$root/data_processing_code/metadata.bash"
manifest=${1:-"$root/data_processing_code/samples.tsv"}
metadata=${2:-"$root/Figure_5/sample_metadata.tsv"}
reference="$root/reference/hg38_bwa/hg38.fa"
rna_dir="$root/data/RNA_seq"
chip_dir="$root/data/CTCF_ChIP"
input_bam="$chip_dir/alignment/CTCF_Input.bam"
parser="$root/Figure_3/CDKN2B_parse_bam_readcount_ChIP.r"
output="$root/Figure_5/results"

load_samples RNA_seq "$manifest" "$rna_dir" false
rna_samples=("${SAMPLE_IDS[@]}")
load_samples CTCF_ChIP "$manifest" "$chip_dir" false
input_found=0
for sample in "${SAMPLE_IDS[@]}"; do
    [[ $sample != CTCF_Input ]] || input_found=1
done
(( input_found )) || { echo "Manifest must contain CTCF_Input" >&2; exit 1; }
[[ -f $reference && -f $reference.fai ]] || {
    echo "Missing reference FASTA or FASTA index: $reference" >&2; exit 1;
}
awk -F '\t' '$1 == "chr9" && $2 >= 22003224 { found = 1 } END { exit !found }' "$reference.fai" || {
    echo "Reference index does not contain chr9 through position 22003224: $reference.fai" >&2; exit 1;
}
[[ -f $parser ]] || { echo "Missing R parser: $parser" >&2; exit 1; }
[[ -f $metadata ]] || { echo "Missing RNA sample metadata: $metadata" >&2; exit 1; }

declare -A known_samples=() experiments=() treatments=()
for sample in "${rna_samples[@]}"; do
    known_samples[$sample]=1
done
IFS= read -r header < "$metadata" || { echo "Empty metadata: $metadata" >&2; exit 1; }
header=${header%$'\r'}
[[ $header == $'sample\texperiment\ttreatment' ]] || {
    echo "Expected metadata header: sample<TAB>experiment<TAB>treatment" >&2; exit 1;
}
while IFS= read -r line || [[ -n $line ]]; do
    line=${line%$'\r'}
    [[ -z $line || $line == \#* || $line == "$header" ]] && continue
    [[ $line != $'\t'* && $line != *$'\t' && $line != *$'\t\t'* ]] || {
        echo "Metadata fields must not be empty: $line" >&2; exit 1;
    }
    IFS=$'\t' read -r sample experiment treatment extra <<< "$line"
    [[ $sample =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ &&
       $experiment =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ &&
       $treatment =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ && ! $extra ]] || {
        echo "Invalid metadata row: $line" >&2; exit 1;
    }
    [[ ${known_samples[$sample]+present} ]] || {
        echo "Metadata sample is not in the RNA manifest: $sample" >&2; exit 1;
    }
    [[ ! ${experiments[$sample]+present} ]] || {
        echo "Duplicate RNA metadata sample: $sample" >&2; exit 1;
    }
    experiments[$sample]=$experiment
    treatments[$sample]=$treatment
done < "$metadata"

tgf_count=0
for sample in "${rna_samples[@]}"; do
    [[ ${experiments[$sample]+present} ]] || {
        echo "Missing RNA metadata for sample: $sample" >&2; exit 1;
    }
    [[ ${treatments[$sample]} != TGF ]] || tgf_count=$((tgf_count + 1))
    bam="$rna_dir/alignment/${sample}Aligned.sortedByCoord.out.bam"
    [[ -f $bam ]] || { echo "Missing RNA BAM: $bam" >&2; exit 1; }
done
(( tgf_count > 0 )) || { echo "No RNA samples have treatment TGF" >&2; exit 1; }
[[ -f $input_bam ]] || { echo "Missing CTCF input BAM: $input_bam" >&2; exit 1; }

ensure_index() {
    local bam=$1 index=""
    for candidate in "$bam.bai" "${bam%.bam}.bai" "$bam.csi" "${bam%.bam}.csi"; do
        if [[ -f $candidate ]]; then
            index=$candidate
            break
        fi
    done
    if [[ -z $index ]]; then
        echo "Indexing $bam"
        conda run --no-capture-output -n samtools_env samtools index -@ 12 "$bam"
    elif [[ $bam -nt $index ]]; then
        echo "BAM index is older than its BAM; rebuild it: $index" >&2
        return 1
    fi
}

for sample in "${rna_samples[@]}"; do
    ensure_index "$rna_dir/alignment/${sample}Aligned.sortedByCoord.out.bam"
done
ensure_index "$input_bam"

count_sample() {
    local locus=$1 position=$2 allele_column=$3 sample_type=$4 sample=$5
    local experiment=$6 treatment=$7 bam=$8 summary=$9
    local out_dir="$output/$locus/$sample_type"
    local readcount="$out_dir/$sample.readcount"
    local parsed="$out_dir/$sample.txt"
    mkdir -p "$out_dir"
    echo "Counting $sample at chr9:$position"
    conda run --no-capture-output -n bam-readcount-env bam-readcount \
        -f "$reference" -q 20 -b 20 "$bam" "chr9:$position-$position" > "$readcount"
    if [[ ! -s $readcount ]]; then
        echo "No readcount reported for $sample at chr9:$position" >&2
        printf 'chr\tpos\tref_base\tTotal_count\tA\tC\tG\tT\n' > "$parsed"
        printf 'chr9\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t%s\t%s\t%s\t%s\tno_readcount\n' \
            "$position" "$sample_type" "$sample" "$experiment" "$treatment" >> "$summary"
        return
    fi
    conda run --no-capture-output -n r_env Rscript "$parser" "$readcount" "$parsed"
    awk -F '\t' -v pos="$position" -v allele="$allele_column" \
        -v sample_type="$sample_type" -v sample="$sample" \
        -v experiment="$experiment" -v treatment="$treatment" '
        BEGIN { OFS = "\t" }
        NR == 2 && $1 == "chr9" && $2 == pos {
            vaf = $4 > 0 ? $allele / $4 : "NA"
            status = $4 > 0 ? "ok" : "zero_depth"
            row = $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS vaf OFS sample_type OFS sample OFS experiment OFS treatment OFS status
            found++
        }
        END {
            if (NR != 2 || found != 1) {
                print "Unexpected parsed readcount for " sample " at chr9:" pos > "/dev/stderr"
                exit 1
            }
            print row
        }
    ' "$parsed" >> "$summary"
}

mkdir -p "$output"
for locus in rs11515 rs3217992; do
    mkdir -p "$output/$locus"
    summary="$output/$locus/${locus}_results+header.tsv"
    printf 'chr\thg38_pos\tref_base\tread_count\tA_count\tC_count\tG_count\tT_count\tVAF\tSample_type\tSample\tExperiment\tTreatment\tStatus\n' > "$summary"
    if [[ $locus == rs11515 ]]; then
        position=21968200
        allele_column=7
    else
        position=22003224
        allele_column=8
    fi
    for sample in "${rna_samples[@]}"; do
        [[ $locus != rs3217992 || ${treatments[$sample]} == TGF ]] || continue
        count_sample "$locus" "$position" "$allele_column" RNA "$sample" \
            "${experiments[$sample]}" "${treatments[$sample]}" \
            "$rna_dir/alignment/${sample}Aligned.sortedByCoord.out.bam" "$summary"
    done
    if [[ $locus == rs11515 ]]; then
        count_sample "$locus" "$position" "$allele_column" DNA CTCF_Input \
            ChIP_Input NA "$input_bam" "$summary"
    fi
done
