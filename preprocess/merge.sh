# PATHS

echo "Copying selected samples to the current directory"
while read samples
do
find /space/sequences/ -name $samples -exec cp -urp -t . {} \;
done < samples

echo "gunzip"
gunzip *.gz

echo "sampling libs"
maxseqs=120000
nlines=$((maxseqs*4))
head -q -n $nlines *R1* > r1.fastq
head -q -n $nlines *R2* > r2.fastq


echo "trimming seqs to 240 bp"
cat *R1* | cut -c 1-240 > r1.fastq
cat *R2* | cut -c 1-240 > r2.fastq


echo "pandaseq"
psoverlap=15
psminlen=245
psmaxlen=270

pandaseq -f r1.fastq -r r2.fastq -o $psoverlap -l $psminlen -L $psmaxlen > merged.fasta 2> temp.panda.txt
tail temp.panda.txt > pandaseq.stats.txt

echo "Generating a sample id file"
head -n 1 *R1* | sed -n '1~3p' | cut -f2 -d " " | cut -f1 -d "_" | cat -n > sample.name.txt
head -n 1 *R1* | sed -n '2~3p' | sed 's/\@/>/'  | cat -n > sample.header.txt
join sample.header.txt sample.name.txt | sed 's/ /,/' | cut -f2-3 -d , | rev | sed 's/ /,/' | rev | sed 's/,/\t/' > sampleid.txt

#data unique=3

unique=3
errorprob=0
subsample=100000
scripts="/space/users/malber06/miscperlscripts/"

perl $scripts/pandaseq.to.qiime.pl -i merged.fasta -s sampleid.txt -u $unique -p $errorprob -m $subsample > filter-qiime.log 

