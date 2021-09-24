
##this pipeline start with fastq files, do mapping (using STAR), gene expression quantification, differential gene expression analysis and differential splicing analysis
#programs needed to be installed in linux:
# perl, R
# STAR
# samtools
# bedtools
# R library: DESeq2, ggplot2, edgeR, BiocParallel, fmsb, gplots, openxlsx, DEXSeq, knitr, pander, rmarkdown.


study_name=project_name #project name
pipeline_root=/directory/to/Pipeline/ #this is root directory of the pipeline/scripts (where the setGlobalVars.sh is located)
project_root=/directory/to/$study_name/ #project folder (all the input and output files related to this project)
mapping_data_root=/directory/to/mapping/$study_name/ #use this directory to save large files (eg. fastq files and mapped bam files). expect a fastq folder with all fastq files under this directory

source $pipeline_root/setGlobalVars.sh
if_PE=1; #whether this is paired-end sequencing ()
if_stranded=1; #whether the RNA library is stranded (contain gene strand information)
if_read1_antiSense=1; #whether the first read is antisense to gene. This is usually true for TruSeq Stranded mRNA Library (second read is sense strand)
if_do_splicing_analysis=1; #whether to do splicing analysis

nCPUs=8
nNodes=1
DEXSeq_nCPUs=3
DEG_method=DESeq2 #DESeq2 fisher
Gblock_anaMethod=DEXSeq #DEXSeq FET
rnum_idname=refseqid
GexFoldchange_cut=2; GexP_cut=0.05; GexHigher_avg_RPKM_cutoff=1

geno=hg19
samples=(siA1 siA2 siB1 siB2 siCtl1 siCtl2) #all sample names
sample_repl_l="siA:siA1 siA2;siB:siB1 siB2;siCtl:siCtl1 siCtl2;" #sample replicates information, follow the same format
test_samples=(siA siB) #test sample groups
ref_samples=(siCtl siCtl) #control sample groups

##PSI analysis setting:
PSIuse_p_type="Pfisher"; PSI_P_cut="0.001"; PSI_cal_rnumMin=20; delta_PSI_cut=10; PSIana_setting=delta_PSI${delta_PSI_cut}_${PSIuse_p_type}$PSI_P_cut


####R.1, mapping
##using STAR to map

cd $mapping_data_root
mkdir STAR_map
fq_dir=$mapping_data_root/fastq/ #the directory of the fastq files 
genomeDir=$pipeline_root/ReferenceDB/map_index/STAR/$geno #the directory of mapping index file
chr_sizes_f=$genomeDir/chrNameLength.txt


if [ $if_PE -eq 1 ]; then
   echo "#treat pair-end read..."
   ## unzip (in some server, STAR will report error if not using unzipped files)
   p=0
   for sample1 in ${samples[*]}; do
     fq1=$fq_dir/${sample1}.1.fq.gz
     fq2=$fq_dir/${sample1}.2.fq.gz
     gunzip $fq1 &
     gunzip $fq2 &
     p=$(($p+2)); echo $p; if [ $p -ge $nCPUs ]; then  p=0 ;  wait; fi
   done
   wait
   
   for sample1 in ${samples[*]}; do
     fq1=$fq_dir/${sample1}.1.fq
     fq2=$fq_dir/${sample1}.2.fq
     out_dir=STAR_map/$sample1/
     mkdir $out_dir
     STAR --runMode alignReads -f --genomeDir $genomeDir --readFilesIn $fq1 $fq2 --runThreadN 5 --outSAMtype BAM SortedByCoordinate \
       --genomeLoad LoadAndKeep --limitBAMsortRAM  30000000000 \
       --outSAMstrandField intronMotif   --outSAMattributes NH HI NM MD jM jI XS --outFileNamePrefix $out_dir
   done
   #    --readFilesCommand "gunzip -c" \

else #single end reads
   echo "#treat single-end read..."
   p=0
   for sample1 in ${samples[*]}; do
     fq1=$fq_dir/${sample1}.fq.gz
     gunzip $fq1 &
     p=$(($p+1)); echo $p; if [ $p -ge $nCPUs ]; then  p=0 ;  wait; fi
   done
   wait
   
   for sample1 in ${samples[*]}; do
     fq1=$fq_dir/${sample1}.fq
     out_dir=STAR_map/$sample1/
     mkdir $out_dir
     STAR --runMode alignReads -f --genomeDir $genomeDir --readFilesIn $fq1 --runThreadN 5 --outSAMtype BAM SortedByCoordinate \
       --genomeLoad LoadAndKeep --limitBAMsortRAM  30000000000 \
       --outSAMstrandField intronMotif   --outSAMattributes NH HI NM MD jM jI XS --outFileNamePrefix $out_dir
   done
   #    --readFilesCommand "gunzip -c" \
fi


STAR  --genomeDir $genomeDir --genomeLoad Remove  #release memory
rm -fr _STARtmp
rm Log.out
rm Log.progress.out
rm Aligned.out.sam


cd STAR_map
out_readnum_f=Input.read.num.log.out
echo -e "Sample\tRaw read" >$out_readnum_f
grep "Number of input reads" */Log.final.out  >>$out_readnum_f
sed -i -e 's/\/Log.final.out:\s*Number of input reads |//g' $out_readnum_f


##create index for IGV visualization
cd $mapping_data_root/STAR_map
p=0
for sample1 in ${samples[*]}; do
  echo $sample1
  ln $sample1/Aligned.sortedByCoord.out.bam $sample1.bam
  samtools index $sample1.bam &
  p=$(($p+1)); echo $p; if [ $p -ge $nCPUs ]; then  p=0 ;  wait; fi
done
wait





####R.2, analyze RNA-seq data
####R2.1  find uniquely mapped reads with good quality from sam or bam file and output table format files for downstream analysis
cd $pipeline_root/RNAseq

p=0
for sample in ${samples[*]}; do
  perl 01.Cal_geno_junc_readsnum_frSAM.pl -n $study_name -s "$sample"  -d "$mapping_data_root/STAR_map/" -o "$project_root/01.ReadTable/" \
    -e ".bam" -t STAR -q 10 -m 5 -u 0  &
  p=$(($p+1)); echo $p; if [ $p -ge $nCPUs ]; then  p=0 ;  wait; fi
done
wait


cd $pipeline_root/RNAseq/
perl 01a.cal_readnum.pl -o $project_root/01.ReadTable



if [ $if_PE -eq 1 ]; then
  ####R2.2  extract pair-end reads mapping info
  cd $pipeline_root/RNAseq
  
  p=0
  for sample in ${samples[*]}; do
    perl 02.Cal_PEreadsCov_frSam.pl -s $study_name -S "$sample"  -i "$mapping_data_root/STAR_map/" -e ".bam" \
      -q 10  -d 0 -t STAR -o "$project_root/02.PE.ReadTable/" &
    p=$(($p+1)); echo $p; if [ $p -ge $nCPUs ]; then  p=0 ;  wait; fi
  done
  wait
  
  #combine PE read stats
  cd $pipeline_root/RNAseq
  Rscript 02a.combinePE.readstats.R -in_dir "$project_root/02.PE.ReadTable/"
  
  
  ####R2.3a calculate gene expression using 14cal_gene_cds_readsnum.pl
  cd $pipeline_root/RNAseq
  p=0
  for i in `seq 1 ${#samples[*]}`; do
    indx=$(($i-1))
    perl 04.cal_gene_cds_readnum.pl -s $study_name -g $geno -i $rnum_idname -p refseqcds_PE. \
      -e CDS -d "$if_stranded" -r "$if_read1_antiSense" -a "${samples[indx]}" \
      -D $project_root/02.PE.ReadTable/ -E reads.cov.tbl -o $project_root/04.GeneReadNum/ &
    p=$(($p+1));  if [ "$p" -ge $nCPUs ]; then p=0; wait;  fi
  done
  wait
  
  cd $pipeline_root/RNAseq/
  perl 04a.cal_readMapStats.pl -o $project_root/04.GeneReadNum/ -p refseqcds_PE
fi



####R2.3b treat all reads as single end (SE) read, calculate junction read map to gene info (as well as gene expression based on CDS region of refseq transcripts)
##the output of this step is not used for PE reads gene expression analysis. bjut it will be used for splcing analysis for both PE and SE reads
cd $pipeline_root/RNAseq
p=0
for sample in ${samples[*]}; do
  perl 04.cal_gene_cds_readnum.pl -s $study_name -g $geno -i $rnum_idname -p refseqcds_SE. -e CDS -d "$if_stranded" -r "$if_read1_antiSense" -a "$sample" \
    -D $project_root/01.ReadTable/ -j 1 -E "JunReadNum GenoMapNum" -o $project_root/04.GeneReadNum/ &
    p=$(($p+1)); echo $p; if [ $p -ge $nCPUs ]; then  p=0 ;  wait; fi
done
wait

cd $pipeline_root/RNAseq/
perl 04a.cal_readMapStats.pl -o $project_root/04.GeneReadNum/ -p refseqcds_SE



####R2.4 calculate gene expression regulation 
cd $pipeline_root/RNAseq
if [ $if_PE -eq 1 ]; then
  gene_rnum_f=$project_root/04.GeneReadNum/refseqcds_PE.ReadNum.tbl
else
  gene_rnum_f=$project_root/04.GeneReadNum/refseqcds_SE.ReadNum.tbl
fi

out_prefix=""

/bin/Rscript 05.cal_gexch.template.R -p_cut $GexP_cut -pval_fun "$DEG_method" -foldchange_cut $GexFoldchange_cut -IfCal_RPM 0 -higher_avg_RPKM_cutoff $GexHigher_avg_RPKM_cutoff \
  -gene_model refseqcds  -geno $geno -study_name $study_name -gene_rnum_f $gene_rnum_f \
  -output_root  $project_root/05.Gene_DE/ \
  -nCores $nCPUs \
  -all_sample_name "${samples[*]}" \
  -gene_len_var cds_len -test_samples "${test_samples[*]}" -ref_samples "${ref_samples[*]}" \
  -g_ano_f "../ReferenceDB/gene/02transcript_gene_ano/$geno.refflat.desc.txt" -rnum_idname "$rnum_idname" \
  -comb_sample_l "$sample_repl_l" \
  -out_prefix "$out_prefix" \
  -what2do  all

  # -what2do  all
  # -what2do  formatGeneTb
  # -what2do  clust_rawSample_exp




####R2.5, create bigwig files for visualization:
PROG_DIR=$pipeline_root/Software/UCSC
genofasta=$pipeline_root/ReferenceDB/ucsc/genomes/$geno/$geno.fa
fext='.bam' #string after sample name in input sam file name

function sam2bw {
 totalReadNum=`samtools view -bc -q 10 $sample$fext`
 #totalReadNum=`wc -l $sample$fext | sed s/[[:blank:]].*//`
 echo "step2, ____ $sample for file $sample$fext, TotalReadNum=$totalReadNum ____"
 #echo step3, ____ $sample, convert to bam ____
 #samtools view -S -b -T $genofasta -o $sample$fext.bam  $sample$fext 
 echo "step3, ____ $sample filter original bam file (keep uniquely mapped reads only) ____"
 samtools view -b -q 10 -o $sample.filtered.bam $sample$fext

 echo step4,  ____ $sample, sort bam file ____
 #samtools sort  -o $sample$fext.bam.sorted.bam $sample$fext.bam  ##!!!no .bam in output name
 echo step5, ____ $sample, convert to bedgraph ____
 genomeCoverageBed -bg -split -ibam $sample.filtered.bam -g $chr_sizes_f > $sample.bedgraph
 #genomeCoverageBed -bg -split -ibam $sample$fext.bam.2.bam -strand + -g $chr_sizes_f > $sample.p.bedgraph #-split 
 #genomeCoverageBed -bg -split -ibam $sample$fext.bam.2.bam -strand '-' -g $chr_sizes_f > $sample.m.bedgraph #-split 
 echo step6, ____ $sample, normolize bedgraph counts ____
 norm_bedgraph.pl -t $totalReadNum -i "$sample.bedgraph" 
 echo step8, ____ $sample, convert to bw ____
 $PROG_DIR/bedGraphToBigWig $sample.bedgraph.normolized $chr_sizes_f $sample.bw 
 #$PROG_DIR/bedGraphToBigWig $sample.p.bedgraph.normolized $chr_sizes_f $sample.p.bw 
 #$PROG_DIR/bedGraphToBigWig $sample.m.bedgraph.normolized $chr_sizes_f $sample.m.bw 
}


p=0
cd $mapping_data_root/STAR_map/
for sample in  ${samples[*]}; do 
  sam2bw  &
  p=$(($p+1)); echo $p; if [ $p -ge $nCPUs ]; then  p=0 ;  wait; fi
done
wait

cd $mapping_data_root/STAR_map/
rm *.bedgraph
rm *.bedgraph.normolized
rm *.bam.unique
rm *.filtered.bam



##R2.6, generate RNA-seq report (read number and gene expression analysis)
cd $pipeline_root/RNAseq
out_prefix=""

/bin/Rscript 06.combine_read_stats.R -study_name "$study_name" -geno "$geno" -project_root "$project_root" -mapping_data_root "$mapping_data_root" \
  -if_PE "$if_PE" -if_do_splicing_analysis "$if_do_splicing_analysis" -if_stranded "$if_stranded" \
  -all_sample_name "${samples[*]}" -p_cut 0.05 -pval_fun "$DEG_method" -foldchange_cut 1.5 -gene_model refseqcds \
  -test_samples "${test_samples[*]}" -ref_samples "${ref_samples[*]}" -out_prefix "$out_prefix" 

cd $project_root/Report
/bin/Rscript -e " rmarkdown::render('RNAseq_report.Rmd', output_format='html_document') " 




if [ $if_do_splicing_analysis != 1 ]; then 
  echo "No splicing analysis was requested. Exit!"
  exit 0; 
fi
echo "Start doing splicing analysis..."




#####Splicing analysis pipeline
####S.1 combine junction read mapping info 
cd $pipeline_root/RNAseq_Splicing
Rscript 01.comb_junc_map_info.R -study_name $study_name -geno $geno -samples "${samples[*]}" -indir $project_root/04.GeneReadNum/refseqcds_SE. \
  -out_cbf $project_root/RNAseq_Splicing/01.comb_junc_map_info/combine.junc2gene.tbl \
  -id_header $rnum_idname



####S.2, cut genes to blocks (exon and intron part) and do expression analysis
cd $pipeline_root/RNAseq_Splicing
perl 02.cut_gene3blocks_basedOnSplicing.pl  -g $geno -R 4 -S 2 -F 0.05 -i $rnum_idname \
  -j $project_root/RNAseq_Splicing/01.comb_junc_map_info/combine.junc2gene.tbl \
  -o $project_root/RNAseq_Splicing/02.geneBlocks

##create bed and bigbed file
perl 02b.Gblock2bigBed.pl -g $geno -a $project_root/RNAseq_Splicing/02.geneBlocks/geneBlock.flat.tbl.ano



####S.3, extract splice site flanking region sequence
cd $pipeline_root/SequenceAnalysis

#extract 100 bp upstream and downstream of splice site (this input file has more junctions, some of them are from gene model only, not from reads)
perl 01.extract_gene_subRegion_seq.pl -d $study_name -g $geno -b 1 -s $project_root/RNAseq_Splicing/02.geneBlocks/AllJunc.tbl \
  -i "contig strand juncpos5" -e juncpos5 -n ss5_PM100 -c m99p100 \
  -o $project_root/RNAseq_Splicing/SS_flank_seq/ss5_PM100.fa &
 
perl 01.extract_gene_subRegion_seq.pl -d $study_name -g $geno -b 1 -s $project_root/RNAseq_Splicing/02.geneBlocks/AllJunc.tbl \
  -i "contig strand juncpos3" -e juncpos3 -n ss3_PM100 -c m100p99 \
  -o $project_root/RNAseq_Splicing/SS_flank_seq/ss3_PM100.fa &
wait


####S.4, calculate gene block expression (exonic and intronic parts) 
cd $pipeline_root/RNAseq
p=0
for sample in ${samples[*]}; do
  perl 04.cal_gene_cds_readnum.pl -s $study_name -g $geno -i Gblock_id -p Gblock. -e transcript -d "$if_stranded" -r "$if_read1_antiSense" -a "$sample" \
    -D $project_root/01.ReadTable/ -j 0 -E "JunReadNum GenoMapNum" -o $project_root/RNAseq_Splicing/geneBlockReadNum/ \
    -A $project_root/RNAseq_Splicing/02.geneBlocks/geneBlock.flat.tbl  &
    p=$(($p+1)); echo $p; if [ $p -ge $nCPUs ]; then  p=0 ;  wait; fi
done
wait
cd $pipeline_root/RNAseq/
perl 04a.cal_readMapStats.pl -o $project_root/RNAseq_Splicing/geneBlockReadNum/ -p Gblock



####S.5 do splicing analysis using DEXSeq 
cd $pipeline_root/RNAseq_Splicing


ana_unit=geneBlock 
#note: need relatively large memory, do not use too many workers if memory is small
Rscript 05.ana_splicing_DEXseq.R -study_name $study_name -sample_repl_l "$sample_repl_l" \
   -test_samples "${test_samples[*]}" -ref_samples "${ref_samples[*]}" -all_samples "${samples[*]}" \
   -workers $DEXSeq_nCPUs \
   -geno $geno \
   -gene_id_header $rnum_idname \
   -in_cbf $project_root/RNAseq_Splicing/geneBlockReadNum/Gblock.ReadNum.tbl \
   -juncAno_f  $project_root/RNAseq_Splicing/02.geneBlocks/geneBlock.flat.tbl.ano \
   -as_ano_root $project_root/RNAseq_Splicing/02.geneBlocks/ASano. \
   -geneAno_f "../ReferenceDB/gene/02transcript_gene_ano/$geno.refflat.desc.txt" \
   -ana_unit $ana_unit \
   -workers $DEXSeq_nCPUs -calP_method $Gblock_anaMethod -out_root $project_root/RNAseq_Splicing/05.$Gblock_anaMethod -what2do all -use_p_type Padj -P_cut 0.001 -Log2ratio_cut 1
 


####S.6, study all exon splicing using PSI
cd $pipeline_root/RNAseq_Splicing

Rscript 06.ana_exon_PSI.R -study_name $study_name -sample_repl_l "$sample_repl_l" \
   -test_samples "${test_samples[*]}" -ref_samples "${ref_samples[*]}" -all_samples "${samples[*]}" \
   -gblock_f "$project_root/RNAseq_Splicing/05.$Gblock_anaMethod/geneBlock/combine_d.tbl" \
   -read_num_f $project_root/RNAseq_Splicing/01.comb_junc_map_info/combine.junc2gene.tbl \
   -out_root $project_root/RNAseq_Splicing/06.exon_PSI/ \
   -workers $nCPUs \
   -use_p_type "Pfisher" -P_cut "0.001" -PSI_cal_rnumMin 20 -delta_PSI_cut 10 \
   -what2do all


##study all annotated exons
Rscript 06.ana_exon_PSI.R -study_name $study_name -sample_repl_l "$sample_repl_l" \
   -test_samples "${test_samples[*]}" -ref_samples "${ref_samples[*]}" -all_samples "${samples[*]}" \
   -read_num_f $project_root/RNAseq_Splicing/01.comb_junc_map_info/combine.junc2gene.tbl \
   -out_root $project_root/RNAseq_Splicing/06.exon_PSI/ \
   -workers $nCPUs \
   -use_p_type "Pfisher" -P_cut "0.001" -PSI_cal_rnumMin 20 -delta_PSI_cut 10 \
   -study_exon_type  annotated_exon -what2do all


