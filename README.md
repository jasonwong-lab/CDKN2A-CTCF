# CDKN2B sequencing analysis

## Data and reference files

FASTQ files are available through dbGaP under **phs001177.v3**. Obtain the
necessary access and download:

- CTCF ChIP-seq, IgG control (if available), and input FASTQs to `./data/CTCF_ChIP/`.
- RNA-seq FASTQs to `./data/RNA_seq/`.

The processing scripts expect gzip-compressed FASTQs. 

The hg38 reference files and indexes are not included. Before processing,
download an hg38 genome FASTA and build the indexes below (or obtain compatible
prebuilt indexes). The scripts expect:

- `./reference/hg38_bwa/hg38.fa` and its matching BWA-MEM2 index files.
- A complete STAR genome index in `./reference/hg38_STAR/`.
- `./reference/hg38.ensGene.gtf`, the annotation used to generate that STAR
  index and originally used for RNA alignment.

For example, the following commands download UCSC's hg38 FASTA and ensGene GTF
and build the indexes. Run them from the repository root after setting up the
conda environments listed below:

```bash
mkdir -p reference/hg38_bwa reference/hg38_STAR
curl -fL -o reference/hg38_bwa/hg38.fa.gz \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gzip -dc reference/hg38_bwa/hg38.fa.gz > reference/hg38_bwa/hg38.fa
curl -fL -o reference/hg38.ensGene.gtf.gz \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ensGene.gtf.gz
gzip -dc reference/hg38.ensGene.gtf.gz > reference/hg38.ensGene.gtf

conda run --no-capture-output -n bwa_env bwa-mem2 index reference/hg38_bwa/hg38.fa
conda run --no-capture-output -n samtools_env samtools faidx reference/hg38_bwa/hg38.fa
conda run --no-capture-output -n STAR_env STAR --runMode genomeGenerate \
  --runThreadN 12 --genomeDir reference/hg38_STAR \
  --genomeFastaFiles reference/hg38_bwa/hg38.fa \
  --sjdbGTFfile reference/hg38.ensGene.gtf --sjdbOverhang 100
```

Ensure FASTA, GTF, and BED chromosome names match (including `chr9`). The supplied
`annotations/MANE_hg38.gtf` is not a substitute for ensGene without evaluating
the effect on the analysis. 

## Sample manifest

Edit `./data_processing_code/samples.tsv`: replace the commented examples with
tab-separated rows, keeping the header. Each row has four fields:

| Field | Meaning |
| --- | --- |
| `assay` | `CTCF_ChIP` or `RNA_seq` |
| `sample` | Output label, unique within the assay; use letters, digits, `_`, `-`, or `.` and start with a letter or digit. |
| `read1` | Gzipped FASTQ filename relative to the assay's download directory. |
| `read2` | Matching mate filename, or `-` for single-end CTCF ChIP-seq. RNA-seq and `CTCF_Input` must be paired-end. |

Use the labels expected by downstream analyses, including `CTCF_Input`,
`CTCF_ChIP_2690`, `CTCF_ChIP_2691`, `CTCF_ChIP_2696`, and `IgG` where
applicable; the downloaded filenames do not need to match those labels.
Check the dbGaP sample metadata before assigning labels. The examples are
placeholders, not a mapping of actual dbGaP filenames.

For multiple lanes of one sample, list filenames in each read field separated
by commas, without spaces. Paired-end lane lists must have equal lengths and
matching order. Paths may contain subdirectories, letters, digits, `_`, `-`,
`.` and `/`, but not spaces or commas in individual filenames. Use a checkout
path without spaces for STAR. The scripts validate the manifest and input
files before executing.

## Conda software

The processing scripts use `conda run` and these environment names. The
versions below are installed in the original local environments; recreate
them with compatible packages if needed. No environment lockfiles are included.

| Environment | Relevant software |
| --- | --- |
| `bwa_env` | BWA-MEM2 2.2.1 |
| `samtools_env` | samtools 1.13 |
| `picard_env` | Picard 3.0.0, OpenJDK 11.0.13 |
| `deeptools_env` | deepTools 3.5.6 |
| `STAR_env` | STAR 2.7.10b |
| `bam-readcount-env` | bam-readcount 0.8 |
| `r_env` | R 4.3.1 |
| `JRE_env` | OpenJDK 22.0.1 (used by the original Picard invocation; the updated processing script uses `picard_env`) |

The historical scripts refer to `star-env`, `bam-readcount_env`, and `R_env`;
the corresponding local environment names above differ in spelling or case.
The original Picard command referenced version 3.3.0, while the currently
installed `picard_env` contains 3.0.0.

## Processing

From the repository root, after filling in the manifest and installing the
reference files and environments:

```bash
bash data_processing_code/process_CTCF_ChIP_seq.bash
bash data_processing_code/process_RNA_seq.bash
```

Either script accepts an alternative manifest path as its first argument.
The CTCF script writes sorted BAMs to `data/CTCF_ChIP/alignment/`, duplicate
marking results to `data/CTCF_ChIP/markdup/`, and genome-wide and chr9
bigWigs to `data/CTCF_ChIP/bigwigs/` and its `chr9/` subdirectory. As in the
original workflow, duplicates are removed from the input BAM, while coverage
uses the original aligned BAMs with deepTools' `--ignoreDuplicates` option.
The RNA script writes coordinate-sorted STAR BAMs and logs to
`data/RNA_seq/alignment/`. Re-running may overwrite outputs; use a clean
output directory or back up results first. Processing requires substantial
disk space, memory, and up to 12 threads per command.

## Figures

Run processing first, then the relevant figure scripts, executing the steps
within each script in the order written. The figure directories are separate
analyses:

- `Figure_3/Figure_3F_code.bash`: after CTCF processing, run
  `bash Figure_3/Figure_3F_code.bash` from the repository root. The script
  uses the manifest's CTCF sample labels and the hg38 FASTA and `.fai` index.
  It uses the duplicate-removed input BAM and derives duplicate-removed BAMs
  for other samples from the marked-duplicate BAMs in `data/CTCF_ChIP/markdup/`.
  It computes and parses read counts for chr9:22008758-22008777 and
  chr9:22008768, then summarizes chr9:22008768. Results are in
  `Figure_3/results/region/` and `Figure_3/results/position/`; derived BAMs
  are in `Figure_3/results/deduplicated/`. An alternative manifest can be
  supplied as the first argument. FASTQs need not still be present. All
  listed samples need the relevant BAMs, and a sample with no reported count
  at chr9:22008768 causes an error. Re-running may overwrite results.
- `Figure_5/Figure_5A_5B_code.bash`: fill in `Figure_5/sample_metadata.tsv`
  with one tab-separated row per RNA sample from the main manifest: `sample`,
  `experiment`, and `treatment`. Labels must start with a letter or digit and
  contain only letters, digits, `_`, `-`, or `.`. Use exactly `TGF` for samples
  to include in the rs3217992 analysis; other treatment labels, such as
  `vehicle`, are included only in rs11515. At least one TGF sample is required.
  Then run:

  ```bash
  bash Figure_5/Figure_5A_5B_code.bash
  ```

  The first and second optional arguments specify alternative main manifest and
  RNA metadata files, respectively. After CTCF and RNA processing, the script
  counts rs11515 (chr9:21968200) in all RNA samples and the aligned,
  **not duplicate-removed** `CTCF_Input` BAM, and rs3217992 (chr9:22003224)
  in TGF RNA samples. It uses mapping and base quality thresholds of 20,
  reports G/depth for rs11515 and T/depth for rs3217992, and reuses the Figure
  3 read-count parser. It creates missing BAM indexes; stale indexes must be
  rebuilt before running. FASTQs need not still be present.

  Summaries, raw counts, and parsed counts are written under
  `Figure_5/results/rs11515/` and `Figure_5/results/rs3217992/`. A status
  column identifies missing readcounts or zero depth; the corresponding VAF is
  `NA`. Re-running may overwrite outputs. The historical Figure 5
  `*_readcount.bash` scripts and `CDKN2B_parse_bam_readcount_v2.r` are not used
  by this workflow and still contain nonportable paths or parsing issues.
- `Figure_S8/Figure_S8_code.bash`: after CTCF processing, run
  `bash Figure_S8/Figure_S8_code.bash` from the repository root. It uses the
  manifest's non-input CTCF sample labels to select genome-wide log2-ratio
  bigWigs, computes motif-centered matrices, then makes heatmaps using
  `annotations/Katainen_CBS_motifs_hg38.bed`. Matrices and PDFs are written to
  `Figure_S8/results/`. An alternative manifest can be supplied as the first
  argument. FASTQs need not still be present, but the listed bigWigs must exist.
  Re-running may overwrite these results.

The historical Figure 5 scripts contain commands that remove intermediate
files or append to existing summaries; inspect them before running them.
