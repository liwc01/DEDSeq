
source ../../setGlobalVars.sh
cd $pipeline_root/Software/seq_ana/maxent
gene_model="refflat"; geno=rn6;  #hg19 mm9 rn6
gene_model="ensGene"; geno=canFam3;  #rheMac2 canFam3

io_dir=$pipeline_root/ReferenceDB/Splicing
mkdir -p $io_dir/02.ss_score/
perl score5.pl $io_dir/01.splicingSite/$geno.$gene_model.ss5.fa >$io_dir/02.ss_score/$geno.$gene_model.ss5.fa.out &
perl score3.pl $io_dir/01.splicingSite/$geno.$gene_model.ss3.fa >$io_dir/02.ss_score/$geno.$gene_model.ss3.fa.out &


##extract ss flanking sequence:
cd $pipeline_root/SequenceAnalysis 
gene_model="refflat"; geno=hg19;  #hg19 mm9 rn6
gene_model="ensGene"; geno=rheMac2;  #rheMac2
study_name=$geno.$gene_model

#extract 100 bp upstream and downstream of splice site
perl 01.extract_gene_subRegion_seq.pl -d $study_name -g $geno -b 1 -s $io_dir/01.splicingSite/$geno.$gene_model.intron.ano.txt \
  -i "contig strand ss5" -e ss5 -n ss5_PM10 -c m9p10 \
  -o $io_dir/SS_flank_seq/$study_name/ss5_PM10.fa 
