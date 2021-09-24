##find uniquely mapped reads with good quality from sam or bam file and output table format files for downstream analysis:

#calculate reads mapped number in genome position or junction from "tophat" output SAM format file:
#junction read strand infomation come from XS tag in tophat output file
#non-junction reads strand infomation come from flag string (0 is + strand, 16 is - strand)
#07/13/2011 add function to calculate total junction and non-junction read number for every reference sequence (could be useful for mapping results not using genome sequence as index, such as mRNA)
#07/29/2011 add function to write a sam file with all uniqely mapped reads with mismatch less than a cutoff (set $if_out_unique_reads_sam)
#03/29/12 change exon read number output: previously only including mapping start position, now including both mapping start and end positions
#7/2/2015 can open bam files using samtools if the file extension is .bam
#7/14/2015 modified junction read identification part. now one junction-read may support multi-junctions (see function from_cigar2juncinfo)
#7/14/2015 modified code to deal with cigar string with deletion or insertion
#7/23/2015 reverse mate2 strand for pair end sequencing ($if_reverse_read_strand_mate2 default is 1)
#4/13/2018 add ifOutJuncReas2Exon, #set 1 to allow output junction read start and end position to exon output file (no 5'ss and 3'ss information will be put in exon output file)
#10/24/2019 add function to output unmapped reads in fastq format (good for reads mappers that were not able to output unmapped reads)

use Getopt::Std; 
use Cwd;

require "../SharedCodes/perl.fun.inc.pl";
require "01.sam.parse.inc.pl";

#common
$debug=0;
$OutSamFrom="STAR"; ####program generating sam file: bwasw, tophat, bowtie, STAR
#$mismatch_cut =4; #maximum allowed mismatch number
$min_mismatch_in100bp=0; #5.5/100 = about 2 mismatch/36bp
$clip_cut="none"; #maximun allowed total clip region at either end of the reads (eg. cigar= 3S..M, 12S..M)
$uniqueness_mapq_minSco=10; #for judge uniqueness of read mapping 
	#in bowtie (not bowtie2) or tophat1, mapg=255 indicating unique mapping; 
	#in bowtie2, or a tophat version using bowtie2, mapg=  -10 log10 p (p is an estimate of the probability that the alignment does not correspond to the read's true point of origin), set 10 (P<0.1)
	#set "" to disable using this to filter, instead, will use occurance of read ID in sam file 
$if_out_unique_reads_sam=1; #if output a sam file with all uniqely mapped reads with mismatch less than a cutoff
$out_unmapped_reads_fq_f="";

@sam_header_arr=split(",","qname,flag,rname,posi,mapq,cigar,mrnm,mpos,isize,seq,qual,tags");

##parameter from command line
#eg: perl 6cal_geno_junc_readsnum_frSAM.pl -n dsRBD_RIP -s "GFP p10A p20 PACT"  -d "../othProject/dsRBD_RIP/3tophat_Map/" -e "/accepted_hits.sam" -t tophat -q 255 -m 5.5 -u 1 -p SE

getopt("sndetqmcupoRIbU",\%args);
if ($sample_name_str=$args{s}){ #sample name string
	$sample_name_str=~s/^\s+|\s+$//g;
	@sample_names_arr=split(/;|\s+/,$sample_name_str);
	$study_name=$args{n} if $args{n};
	$indir=$args{d} if $args{d};
	$samfilename=$args{e} if $args{e};
	$OutSamFrom=$args{t} if $args{t};
	$uniqueness_mapq_minSco=$args{q} if $args{q} ne "";
	$min_mismatch_in100bp=$args{m} if $args{m}  ne "";
	$clip_cut=$args{c} if $args{c}  ne "";
	$if_out_unique_reads_sam=$args{u} if $args{u} ne "";
	$out_unmapped_reads_fq_f=$args{U} if $args{U} ne "";
	$outdir=$args{o}?$args{o}:"$indir/../01.ReadTable";
	$if_reverse_read_strand_mate2=$args{R} ne "" ? $args{R} : 1;
	$ifOutJuncReas2Exon=$args{I} ne "" ? $args{I} : 0; #set 1 to allow output junction read start and end position to exon output file (no 5'ss and 3'ss information will be put in exon output file)
	$debug=$args{b} ne "" ? $args{b} : 0;
}


##########################RUN

print "\n".join("",map("$_=$$_\n", ("study_name","sample_name_str","indir","samfilename","OutSamFrom","uniqueness_mapq_minSco","min_mismatch_in100bp","if_out_unique_reads_sam","outdir","debug")));

foreach $samplename(@sample_names_arr){
	$sam_file="$indir$samplename$samfilename";
	my $openSamCmd=$sam_file;
	if($samfilename=~/.bam$/){
		$openSamCmd="samtools view $sam_file | ";
	}
	print "\n\$sam_file=$openSamCmd\n";
	create_dir_ifNotExist("$outdir/$samplename/");
	if($if_out_unique_reads_sam){
		open(OUTSAM_F,">$sam_file.unique") || die "error write $sam_file.unique\n";
	}
	if($out_unmapped_reads_fq_f){
		create_dir_ifNotExist("$out_unmapped_reads_fq_f");
		open(OUT_UNMAPPED_FQ_F,">$out_unmapped_reads_fq_f") || die "error write $out_unmapped_reads_fq_f\n";
	}
	$readsid_tmp_file="$outdir/$samplename/readsIDs.all.tmp";
	
	$log_file="$outdir/$samplename.run.log";
	open (LOGF,">$log_file") || die "error write $log_file\n";
	%all_non_unique_qname=();
	%counts=();
	

	
	#step2: for unique mapped reads, count geno and junc mapped numbers
	#for non-unique mapped reads, output the map positions
	print LOGF return_time()."|read SAM file again, count reads map numbers\n";
	%exon_readsnum_hash=();
	%junc_readsnum_hash=();
	%ununique_readsinfo_hash=();
	
	open (GENO_UNUNI_MAP,">$outdir/$samplename.GenoMapUNunique") || die "error write $outdir/$samplename.GenoMapUNunique\n";
	print GENO_UNUNI_MAP "rname	reads_strand	posi	posto	readsnum\n";
	open (JUNC_UNUNI_MAP,">$outdir/$samplename.JuncMapUNunique") || die "error write $outdir/$samplename.JuncMapUNunique\n";
	print JUNC_UNUNI_MAP "rname	junc_strand	juncpos5	juncpos3	exonreg5len	readsnum\n";
	open(SAM_F,$openSamCmd) || die "error open $openSamCmd\n";
	$rowi=0;
	while(<SAM_F>){
		next if /^@/; #comment line: @HD     VN:1.0  SO:sorted
		s/\s+$//;
		($qname,$flag,$rname,$posi,$mapq,$cigar,$mrnm,$mpos,$isize,$seq,$qual,@tags_arr)=split(/\t/);
		$rowi++;
		last if($rowi>100000 && $debug); ##debug
		if($rowi % 1000000==0 || eof){
			print LOGF return_time()."|read SAM file $samplename line $rowi\n";
		}
		$counts{"mapping.rows.$samplename"}++;
		if ($rname=~/\*/){#unmapped
			$counts{"reads.$samplename.1.unmapped"}++;
			if($out_unmapped_reads_fq_f){
				print OUT_UNMAPPED_FQ_F '@'."$qname\n$seq\n+\n$qual\n";
			}
			next ; 
		}
		#uniqueness
		if($uniqueness_mapq_minSco && $mapq < $uniqueness_mapq_minSco){ #judge uniqueness by $mapq
			$all_non_unique_qname{$qname}=1;
		}
		if($flag & 0x100){ # secondary alignment
			next;
		}
		#mis-match or clip
		$mismatch_cut=$min_mismatch_in100bp*(&cigar2cover_len($cigar)/100);
		if(! judge_reads_quality($OutSamFrom,$mismatch_cut,$clip_cut,$cigar,@tags_arr)){
			if(!$all_non_unique_qname{$qname}){ #unique but bad quality
				$counts{"reads.$samplename.3.unique.clip>$clip_cut.mismatch>$min_mismatch_in100bp/100bp"}++;
			}
			next;
		}
		
		#read mapping type: exonic or junction
		($reads_cigarcode,$readmap_type)=judge_readstype_fr_ciger($cigar);
		$reads_strand=from_flag2strand($flag);
		if($if_reverse_read_strand_mate2){
			if($flag & 0x0080){
				$reads_strand=$reads_strand=~/-/?"+":"-";
			}
		}
		

		if(!$all_non_unique_qname{$qname}){ #unique
			$counts{"reads.$samplename.4.unique.readstype=$readmap_type"}++;
			print OUTSAM_F "$_\n" if $if_out_unique_reads_sam;
		}
		my $posto=$posi+cigar2cover_len($cigar)-1;
		if($readmap_type eq "exon"){
			if($all_non_unique_qname{$qname}){
				$ununique_readsinfo_hash{$readmap_type}{"$rname	$reads_strand	$posi	$posto"}++;
			}else{ #unique
				$exon_readsnum_hash{"$rname	$reads_strand	$posi	$posto	e"}++;
				$ref2rnum_h{$rname}{$samplename}{"$readmap_type"}++;
			}
		}elsif($readmap_type eq "junc"){
			$tags=join(" ",@tags_arr);
			my ($juncpos5s,$juncpos3s,$junc_strand,$exonreg5lens,$exonreg3lens)=from_cigar2juncinfo($posi,$cigar,$tags);
			my @juncpos5s_arr=split(/;/,$juncpos5s);
			my @juncpos3s_arr=split(/;/,$juncpos3s);
			my @exonreg5lens_arr=split(/;/,$exonreg5lens);
			my @exonreg3lens_arr=split(/;/,$exonreg3lens);
			foreach my $junc_j(1..(scalar @juncpos5s_arr )){
				my $juncpos5=$juncpos5s_arr[$junc_j-1];
				my $juncpos3=$juncpos3s_arr[$junc_j-1];
				my $exonreg5len=$exonreg5lens_arr[$junc_j-1];
				if($all_non_unique_qname{$qname}){
					$ununique_readsinfo_hash{$readmap_type}{"$rname	$junc_strand	$juncpos5	$juncpos3	$exonreg5len"}++;
				}else{ #unique
					$junc_readsnum_hash{"$rname	$junc_strand	$juncpos5	$juncpos3"}{"num"}++; #each junction read may count >1 times here
					$junc_readsnum_hash{"$rname	$junc_strand	$juncpos5	$juncpos3"}{"poss"}.="$exonreg5len,"; 
				}
			}
			if(!$all_non_unique_qname{$qname}){ #each junction read will count once here
				$ref2rnum_h{$rname}{$samplename}{"$readmap_type"}++;
				if($ifOutJuncReas2Exon){
					$exon_readsnum_hash{"$rname	$reads_strand	$posi	$posto	j"}++;
				}
			}
		}#end elsif($readmap_type eq "junc")
	}
	close SAM_F;
	close OUTSAM_F if $if_out_unique_reads_sam;
	close OUT_UNMAPPED_FQ_F if $out_unmapped_reads_fq_f;
	
	##2.2 output non-uniquely mapped reads info
	print LOGF return_time()."|for non-unique mapped reads, output the map positions\n";
	print GENO_UNUNI_MAP join("\n", map(
		"$_	".$ununique_readsinfo_hash{"exon"}{$_},  sort keys %{$ununique_readsinfo_hash{"exon"}}
	) );
	close GENO_UNUNI_MAP;
	
	print JUNC_UNUNI_MAP join("\n", map(
		"$_	".$ununique_readsinfo_hash{"junc"}{$_},  sort keys %{$ununique_readsinfo_hash{"junc"}}
	) );
	close JUNC_UNUNI_MAP;
	undef %ununique_readsinfo_hash;
	
	##2.3 write counts
	$counts{"reads.$samplename.2.non-unique"}=scalar keys %all_non_unique_qname;
	undef %all_non_unique_qname;
	print LOGF return_time()."|write counts information\n";
	open (READSCOUNTS_OUT, ">$outdir/$samplename/readsIDs.counts.txt") || die "error write readsIDs.counts.txt\n";
	print "counts:\n". join("", map("$_	$counts{$_}\n", sort keys %counts)). "\n";
	print LOGF "counts:\n". join("", map("$_	$counts{$_}\n", sort keys %counts)). "\n";
	print READSCOUNTS_OUT join("", map("$_	$counts{$_}\n", sort keys %counts)). "\n";
	close READSCOUNTS_OUT;

	##2.4 output uniquely mapped reads numbers
	print LOGF return_time()."|for unique mapped reads, count genome and junction mapped numbers\n";
	#write geno map readsnum	
	open (GENOMAP_NUM,">$outdir/$samplename.GenoMapNum") || die "error write $outdir/$samplename.GenoMapNum\n";
	print GENOMAP_NUM "rname	reads_strand	posfr	posto	read_type	readsnum\n";
	print LOGF return_time()."|write genome reads number to  $outdir/$samplename.GenoMapNum\n";
	foreach $geno_posinfo(sort keys %exon_readsnum_hash){
		print GENOMAP_NUM "$geno_posinfo	".$exon_readsnum_hash{$geno_posinfo}."\n";
	}
	close GENOMAP_NUM;
	undef %exon_readsnum_hash;
	
	#write junc map readsnum	
	open (JUNCMAP_NUM,">$outdir/$samplename.JunReadNum") || die "error write $outdir/$samplename.JunReadNum\n";
	print JUNCMAP_NUM "rname	junc_strand	juncpos5	juncpos3	readsnum	exonreg5lens\n";
	print LOGF return_time()."|write junc reads number to $outdir/$samplename.GenoMapNum\n";
	foreach $geno_posinfo(sort keys %junc_readsnum_hash){
		print JUNCMAP_NUM	"$geno_posinfo	".$junc_readsnum_hash{$geno_posinfo}{"num"}."	".$junc_readsnum_hash{$geno_posinfo}{"poss"}."\n";
	}
	close JUNCMAP_NUM;
	undef %junc_readsnum_hash;
	
	#2.5 write rname readnum for 1 sample
	open (RNAME_RNUM,">$outdir/$samplename.refReadNum") || die "error write $outdir/$samplename.refReadNum\n";
	print RNAME_RNUM "rname	ERnum_$samplename	JRnum_$samplename\n";
	foreach my $rname(sort keys %ref2rnum_h){ #$ref2rnum_h{$rname}{$samplename}{"$readmap_type"}
		print RNAME_RNUM "$rname	". join("	",map($ref2rnum_h{$rname}{$samplename}{$_},("exon","junc")) ). "\n";
	}
	close RNAME_RNUM;
	close LOGF;

}#end foreach $samplename

#write rname readnum for all samples
open (RNAME_RNUM,">$outdir/allsp.refReadNum") || die "error write $outdir/allsp.refReadNum\n";
print RNAME_RNUM "rname	". join("	",map("ERnum_$_	JRnum_$_", @sample_names_arr) )."\n";
foreach my $rname(sort keys %ref2rnum_h){ #$ref2rnum_h{$rname}{$samplename}{"$readmap_type"}
	print RNAME_RNUM "$rname	". join("	",map($ref2rnum_h{$rname}{$_}{"exon"}."	".$ref2rnum_h{$rname}{$_}{"junc"}, @sample_names_arr) ). "\n";
}
close RNAME_RNUM;


