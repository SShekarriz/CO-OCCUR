#!/bin/bash

####################################################################################
DATE=`date '+%Y-%m-%d-%H:%M'`
exec > >(tee -i CO-OCCUR_${DATE}.log)
exec 2>&1
####################################################################################

ref_lib=$PWD/genomes
genes="folP_glmM.fa"
folP="folP.fa"
bin=$PWD/bin
res=$PWD/results
output=$PWD/outputs
mkdir -p $res $output

ref_db=$res/genomes.fa
chrom=$res/chrom.sizes
Cnames=$res/genome_contig.names

out=$chrom
if [ -f $out ]; then
	echo "The file '$out' exists."
else
	echo "The file '$out' is not found."
	cat $ref_lib/*fasta > $ref_db
	makeblastdb -in $ref_db -dbtype nucl
	samtools faidx $ref_db
	grep "^>" *.fasta | sed 's/ .*//g' > $Cnames
	cut -f 1,2 ${ref_db}.fai > $chrom
fi


blastout=$res/${folP%.fa}-against-ref.blast
out=${blastout}.5K.fasta
if [ -f $out ]; then
        echo "The file '$out' exists."
else
        echo "The file '$out' is not found."

	blastn -db $ref_db -query $folP -outfmt 6 -out $blastout
	cut -f 2,9,10 $blastout | awk '{if($2 > $3) {t = $2; $2 = $3; $3 = t; print;} else if($2 < $3) print; }' OFS='\t' |awk '{a=$2-1;print $1,a,$3;}' OFS='\t'|bedtools sort > ${blastout}.bed
	bedtools slop -i ${blastout}.bed -g $chrom -b 5000 > ${blastout}.5K.bed
	bedtools getfasta -fi $ref_db -bed ${blastout}.5K.bed -fo ${blastout}.5K.fasta
fi


cluster=${blastout}.5K.fasta
blastout=$res/${genes%.fa}-against-cluster.blast
distan=$output/${genes%.fa}-distance.txt
out=$distan
if [ -f $out ]; then
        echo "The file '$out' exists."
else
        echo "The file '$out' is not found."
	makeblastdb -in $cluster -dbtype nucl
	blastn -db $cluster -query $genes -outfmt 6 -out $blastout
	Rscript $bin/distancer.R $blastout $distan
fi


echo "################################"
echo "## Running Prokka on contigs ###"
echo "################################"
#activate conda env for prokka, i have it in gtdbk
source /home/shekas3/anaconda3/bin/activate py3

Prok_out=$res/folP-against-ref.blast.5K.prokka
Prok=$Prok_out/out.gff
if [ -f $Prok ]; then
	echo "The file '$Prok' exists."
else
        echo "The file '$Prok' is not found."
        echo "Running Prokka on the contig file now"
	seqkit replace -p '.+' -r 'cluster_{nr}' $cluster > ${cluster%.fasta}_edited.fasta
	grep -e "^>" $cluster > head1.txt
	grep -e "^>" ${cluster%.fasta}_edited.fasta > head2.txt
	paste -d"," head1.txt head2.txt > ${cluster%.fasta}.cluster-names
	rm head*.txt
	prokka --cpus 15 --outdir $Prok_out --prefix out ${cluster%.fasta}_edited.fasta
fi


source /home/shekas3/anaconda3/bin/activate base
annotate=$output/cluster_annotations.txt
out=$annotate
if [ -f $out ]; then
        echo "The file '$out' exists."
else
        echo "The file '$out' is not found."
	Rscript $bin/cluster.annotator.R $Prok ${cluster%.fasta}.cluster-names $Cnames $annotate  
fi



echo "################### DONE ####################################"
date -u
echo "#############################################################"


