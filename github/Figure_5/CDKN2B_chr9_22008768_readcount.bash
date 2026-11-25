
## read count hg38 chr9:22008768

 conda activate bam-readcount_env

##

cd /home/jbarbour/CDKN2B/RNA_seq/NIH_Fibroblast_RNA_seq_August_2020/bam_readcount


for f in /home/jbarbour/CDKN2B/RNA_seq/NIH_Fibroblast_RNA_seq_August_2020/alignment/*bam ;
do echo $f ;

    id=$(echo ${f%Aligned.sortedByCoord.out.bam})
    id2=$(echo ${id#*/alignment/})

	bam-readcount -f /home/jbarbour/reference/hg38.fa $f chr9:22008768-22008768 -q 20 -b 20 > $id2.readcount &
	
done
conda activate R_env

Rscript /home/jbarbour/CDKN2B/code/CDKN2B_parse_bam_readcount_v2.r

##

 conda activate bam-readcount_env

cd /home/jbarbour/CDKN2B/RNA_seq/NIH_Fibroblast_RNA_seq_March_2021/bam_readcount


for f in /home/jbarbour/CDKN2B/RNA_seq/NIH_Fibroblast_RNA_seq_March_2021/alignment/*bam ;
do echo $f ;

    id=$(echo ${f%Aligned.sortedByCoord.out.bam})
    id2=$(echo ${id#*/alignment/})
    echo $id2
	bam-readcount -f /home/jbarbour/reference/hg38.fa $f chr9:22008768-22008768 -q 20 -b 20 > $id2.readcount &
	
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
	bam-readcount -f /home/jbarbour/reference/hg38.fa $f chr9:22008768-22008768 -q 20 -b 20 > $id2.readcount &
	
done
conda activate R_env

Rscript /home/jbarbour/CDKN2B/code/CDKN2B_parse_bam_readcount.r


cd /home/jbarbour/CDKN2B/RNA_seq/DNA_chr9_22008768_readcount

 conda activate bam-readcount_env

for f in /home/jbarbour/CDKN2B/CHIP_seq/alignment/*Input*bam ;

do echo $f ;

    id=$(echo ${f%.bam})
    id2=$(echo ${id#*/alignment/})
    echo $id2
	bam-readcount -f /home/jbarbour/reference/hg38.fa $f chr9:22008768-22008768 -q 20 -b 20 > $id2.readcount &
	
done


conda activate R_env

Rscript /home/jbarbour/CDKN2B/code/CDKN2B_parse_bam_readcount_v2.r


############ Compile results

cd /home/jbarbour/CDKN2B/RNA_seq/RNA_chr9_22008768_readcount_results
rm chr9_22008768_results.txt
for f in /home/jbarbour/CDKN2B/RNA_seq/NIH_Fibroblast_RNA_seq_August_2020/bam_readcount/*txt ;
do 

    id=$(echo ${f%.txt})
    id2=$(echo ${id#*/bam_readcount/})
	exp=$(echo ${f%/bam_readcount*})
	exp2=$(echo ${exp#*/RNA_seq/})
	echo $id
	echo $exp
	echo $exp2
	
	awk -v ID=$id2 -v EXP=$exp2 'NR == 2 {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"($7/$4)"\t""RNA""\t"ID"\t"EXP}' $f >> chr9_22008768_results.txt
	
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
	
	awk -v ID=$id2 -v EXP=$exp2 'NR == 2 {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"($7/$4)"\t""RNA""\t"ID"\t"EXP}' $f >> chr9_22008768_results.txt
	
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
	
	awk -v ID=$id2 -v EXP=$exp2 'NR == 2 {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"($7/$4)"\t""RNA""\t"ID"\t"EXP}' $f >> chr9_22008768_results.txt
	
done



for f in /home/jbarbour/CDKN2B/RNA_seq/DNA_chr9_22008768_readcount/*txt ;
do 

    id=$(echo ${f%.txt})
    id2=$(echo ${id#*/DNA_rs11515_readcount/})
	exp=$(echo ${f%/bam_readcount*})
	exp2=$(echo ${exp#*/RNA_seq/})
	echo $id
	echo $exp
	echo $exp2
	
	awk -v ID=$id2 -v EXP=$exp2 'NR == 2 {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"($7/$4)"\t""DNA""\t"ID"\t""ChIP_Input"}' $f >> chr9_22008768_results.txt
	
done


{ printf 'chr\thg38_pos\tref_base\tread_count\tA_count\tC_count\tG_count\tT_count\tVAF\tSample_type\tSample\tExperiment\n'; cat chr9_22008768_results.txt ; } > chr9_22008768_results+header.txt




















