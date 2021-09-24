
source ../../setGlobalVars.sh
echo $pipeline_root

########download refFlat files (Refseq transcripts genome coordinates information, download these files regularly)
for geno in canFam3; do # hg19 mm9 rn6 canFam3
 echo __________ $geno __________
 dir=$pipeline_root/ReferenceDB/ucsc/refFlat/$geno
 mkdir -p $dir
 cd $dir
 wget http://hgdownload.cse.ucsc.edu/goldenPath/$geno/database/refFlat.txt.gz
 gunzip refFlat.txt.gz
done


######download Refseq annotations file from UCSC (if the refFlat file is not complete)
for geno in rn6 canFam3; do # rheMac10 rn6 canFam3
 echo __________ $geno __________
 dir=$pipeline_root/ReferenceDB/ucsc/ncbiRefSeq/$geno
 mkdir -p $dir
 cd $dir
 wget http://hgdownload.soe.ucsc.edu/goldenPath/$geno/database/ncbiRefSeq.txt.gz
 gunzip ncbiRefSeq.txt.gz
done


######### download the genome sequence files (fasta format) (these files only need to download once)
for geno in hg19 mm9; do   # hg19 mm9 rheMac10
 echo __________ download the genome sequence of $geno __________
 dir=$pipeline_root/ReferenceDB/ucsc/genomes/$geno
 mkdir -p $dir
 cd $dir
 rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/$geno/chromosomes/ ./
done

for geno in rheMac10; do   # rheMac10
 echo __________ download the genome sequence of $geno __________
 dir=$pipeline_root/ReferenceDB/ucsc/genomes/$geno
 mkdir -p $dir
 cd $dir
 wget http://hgdownload.soe.ucsc.edu/goldenPath/rheMac10/bigZips/$geno.fa.gz
done



#unzip and combine fasta files:
for geno in hg19 mm9; do   # hg19 mm9
 echo __________ unzip and combine genome sequence for $geno __________
 dir=$pipeline_root/ReferenceDB/ucsc/genomes/$geno
 cd $dir
 gunzip *gz
 cat chr*.fa >$geno.fa
done



########## download ensGene (Ensembl transcripts genome coordinates information, download these files regularly)
for geno in rn6 canFam3; do #rheMac10 rn6 canFam3
 echo __________ $geno __________
 dir=$pipeline_root/ReferenceDB/ucsc/ensGene/$geno
 cd $root_dir/data/ucsc/ensGene/
 mkdir -p $dir
 cd $dir
 wget http://hgdownload.cse.ucsc.edu/goldenPath/$geno/database/ensGene.txt.gz
 gunzip *gz
done

