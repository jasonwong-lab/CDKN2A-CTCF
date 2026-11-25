
## read count hg38 RNA_rs3217992 chr9:22,003,224-22,003,224

 conda activate bam-readcount_env

cd /home/jbarbour/CDKN2B/RNA_seq/NIH_Fibroblast_RNA_seq_August_2020/bam_readcount


for f in /home/jbarbour/CDKN2B/RNA_seq/NIH_Fibroblast_RNA_seq_August_2020/alignment/*bam ;
do echo $f ;

    id=$(echo ${f%Aligned.sortedByCoord.out.bam})
    id2=$(echo ${id#*/alignment/})

	bam-readcount -f /home/jbarbour/reference/hg38.fa $f chr9:22003224-22003224 -q 20 -b 20 > $id2.readcount &
	
done
conda activate R_env

Rscript /home/jbarbour/CDKN2B/code/CDKN2B_parse_bam_readcount.r


##

 conda activate bam-readcount_env

cd /home/jbarbour/CDKN2B/RNA_seq/NIH_Fibroblast_RNA_seq_March_2021/bam_readcount


for f in /home/jbarbour/CDKN2B/RNA_seq/NIH_Fibroblast_RNA_seq_March_2021/alignment/*bam ;
do echo $f ;

    id=$(echo ${f%Aligned.sortedByCoord.out.bam})
    id2=$(echo ${id#*/alignment/})
    echo $id2
	bam-readcount -f /home/jbarbour/reference/hg38.fa $f chr9:22003224-22003224 -q 20 -b 20 > $id2.readcount &
	
done

conda activate R_env

Rscript /home/jbarbour/CDKN2B/code/CDKN2B_parse_bam_readcount.r

##

 conda activate bam-readcount_env

cd /home/jbarbour/CDKN2B/RNA_seq/NIH_Fibroblast_RNA_seq_May_2021/bam_readcount

for f in /home/jbarbour/CDKN2B/RNA_seq/NIH_Fibroblast_RNA_seq_May_2021/alignment/*bam ;
do echo $f ;

    id=$(echo ${f%Aligned.sortedByCoord.out.bam})
    id2=$(echo ${id#*/alignment/})
    echo $id2
	bam-readcount -f /home/jbarbour/reference/hg38.fa $f chr9:22003224-22003224 -q 20 -b 20 > $id2.readcount &
	
done
conda activate R_env

Rscript /home/jbarbour/CDKN2B/code/CDKN2B_parse_bam_readcount.r



############ Compile results

cd /home/jbarbour/CDKN2B/RNA_seq/RNA_rs3217992_readcount
rm *results.txt
for f in /home/jbarbour/CDKN2B/RNA_seq/NIH_Fibroblast_RNA_seq_August_2020/bam_readcount/*txt ;
do 

    id=$(echo ${f%.txt})
    id2=$(echo ${id#*/bam_readcount/})
	exp=$(echo ${f%/bam_readcount*})
	exp2=$(echo ${exp#*/RNA_seq/})
	echo $id
	echo $exp
	echo $exp2
	
	awk -v ID=$id2 -v EXP=$exp2 'NR == 2 {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"($8/$4)"\t""RNA""\t"ID"\t"EXP}' $f >> chr9_22008768_results.txt
	
done

for f in /home/jbarbour/CDKN2B/RNA_seq/NIH_Fibroblast_RNA_seq_March_2021/bam_readcount/*txt ;
do 

    id=$(echo ${f%.txt})
    id2=$(echo ${id#*/bam_readcount/})
	exp=$(echo ${f%/bam_readcount*})
	exp2=$(echo ${exp#*/RNA_seq/})
	echo $id
	echo $exp
	echo $exp2
	
	awk -v ID=$id2 -v EXP=$exp2 'NR == 2 {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"($8/$4)"\t""RNA""\t"ID"\t"EXP}' $f >> chr9_22008768_results.txt
	
done

for f in /home/jbarbour/CDKN2B/RNA_seq/NIH_Fibroblast_RNA_seq_May_2021/bam_readcount/*txt ;
do 

    id=$(echo ${f%.txt})
    id2=$(echo ${id#*/bam_readcount/})
	exp=$(echo ${f%/bam_readcount*})
	exp2=$(echo ${exp#*/RNA_seq/})
	echo $id
	echo $exp
	echo $exp2
	
	awk -v ID=$id2 -v EXP=$exp2 'NR == 2 {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"($8/$4)"\t""RNA""\t"ID"\t"EXP}' $f >> chr9_22008768_results.txt
	
done

 grep TGF chr9_22008768_results.txt| tail -n+2 |\
 { printf 'chr\thg38_pos\tref_base\tread_count\tA_count\tC_count\tG_count\tT_count\tVAF\tSample_type\tSample\tExperiment\n'; cat - ; } > rs3217992_chr9_22003224_results+header.txt




















