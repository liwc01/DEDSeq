#1, from .sam file, calculate reads coverage and then 
##10/29/2014 only output read mapped in proper pair (# in output is fragement)
#7/2/2015 can open bam files using samtools if the file extension is .bam
#4/5/2016 can output a table of read numbers for individual contigs (good for map to transcript analysis)
#2/1/2021 can output fastq files for unmapped reads (paired-end)

require "../SharedCodes/perl.fun.inc.pl";
use Getopt::Std; 
getopt("dEMC2tBsioeSpqU",\%args);


#######common
$debug=$args{d} ne ""? $args{d}:0 ; 

@sam_header_arr=split(",","qname,flag,rname,posi,mapq,cigar,mrnm,mpos,isize,seq,qual,tags");
$pairend_gapmax="50000"; ### this length are usually largest intron length in the species


$use_exist_reads_cov_file=0; #set 1 to allow using old reads.cov.tbl file (if exist) to same time.
$use_exist_reads_cov_file=$args{E} if $args{E};

$mismatch_cut=$args{M} ne ""? $args{M}: "none"; #minimal number of mismatch for both mates
$clip_cut=$args{C} ne ""? $args{C}: "none"; #minimal number of soft clip
$if_count_sencondary_map=$args{2} ne ""? $args{2}: 0; #whether or not count sencondary mapping (set to 1 if allow non-uniquely mapped reads to be counted multi-times)
$OutSamFrom=$args{t}? $args{t} : "STAR";
$if_best_score=$args{B}? $args{B} : 1;

$study_name=$args{s} if $args{s};
$sam_dir=$args{i} if $args{i};
$readsCov_outdir=$args{o}?$args{o}: "$sam_dir/../02.PE.ReadTable";
$samfilename=$args{e} if $args{e}; # .bam
if($args{S}){@sample_names_arr=split(" |,",$args{S}); }
$outPrefix=$args{p}?$args{p}:"";
$uniqueness_mapq_minSco=$args{q} ne "" ? $args{q} : 10;
$out_unmapped_reads_fq_f=$args{U} if $args{U} ne ""; #will add .1.fq and .2.fq to the file name for read 1 and 2 respectively


##################################
create_dir_ifNotExist("$readsCov_outdir/$outPrefix");

my $openSamCmd;
foreach $samplename(@sample_names_arr){
	
	$sam_file="$sam_dir$samplename$samfilename";
	if($samfilename=~/.bam$/){
		$openSamCmd="samtools view $sam_file | ";
	}else{
		$openSamCmd=" $sam_file "
	}

	$reads_cov_f="$readsCov_outdir/$outPrefix$samplename.reads.cov.tbl";
	$refname_rnum_f="$readsCov_outdir/$outPrefix$samplename.rname.readnum.tbl";
	
	if ($use_exist_reads_cov_file && -e $reads_cov_f){
		print " file $reads_cov_f already exist; \nSo skip read SAM file $sam_file\n";
	}else{
		
		$log_file="$readsCov_outdir/$outPrefix$samplename.run.log";
		open (LOGF,">$log_file") || die "error write $log_file\n";
		
		if($out_unmapped_reads_fq_f){
			create_dir_ifNotExist("$out_unmapped_reads_fq_f");
			open(OUT_UNMAPPED_FQ_F1,">$out_unmapped_reads_fq_f.1.fq") || die "error write $out_unmapped_reads_fq_f.1.fq\n";
			open(OUT_UNMAPPED_FQ_F2,">$out_unmapped_reads_fq_f.2.fq") || die "error write $out_unmapped_reads_fq_f.2.fq\n";
			my %unmapped_reads_h=();
		}
		
		###step2: for unique mapped reads, output pair-end combined and simplified reads mapping data
		my %readstype2num_hash=();
		print " open \$sam_file=$openSamCmd\n write $reads_cov_f\n";
		open(SAM_F,$openSamCmd) || die "error open $openSamCmd\n";
		open (MAPOUT,">$reads_cov_f") || die "error write $reads_cov_f\n";
		open (RNAME_RNUM_OUT,">$refname_rnum_f") || die "error write $refname_rnum_f\n";
		print MAPOUT "chromosome	mate1strand	posfr	posto	readsnum\n";
		
		my %rname2rnum_h=();
		$rowi=0;
		while(<SAM_F>){
			next if /^@/; #comment line: @HD     VN:1.0  SO:sorted
			s/\s+$//;
			($qname_raw,$flag,$rname,$posi,$mapq,$cigar,$mrnm,$mpos,$isize,$seq,$qual,@tags_arr)=split(/\t/);
			$qname=$qname_raw;
			$qname=~s/\.[12]$//; 
			$rowi++;
			if($rowi % 1000000==0){
				print  return_time()."|read SAM file $samplename line $rowi\n";
			}
			if(eof){
				print  return_time()."|read SAM file $samplename line $rowi\n";
				print LOGF "#".return_time()."|read SAM file $samplename line $rowi\n";
			}
			
			last if($rowi>100000 && $debug); 
			$counts{$samplename}{"0.reads.row_inSam"}++;
			if($flag & 0x0004){ #segment not mapped
				if ($flag & 0x8){ #next segment in the template unmapped
					my $read1or2=($flag & 0x0080) ? 2: 1;
					$unmapped_reads_h{$qname}{$read1or2}='@'."$qname $read1or2\n$seq\n+\n$qual\n";
					$counts{$samplename}{"1.reads.bothMateUnMapped.mate$read1or2"}++;
				}
				next;
			} 
			
			my $read_map_primary_orNot;
			if($flag & 0x0100){
				$read_map_primary_orNot="secondaryAli";
			}else{
				$read_map_primary_orNot="primaryAli";
				
			}
			#print "flag=$flag $read_map_primary_orNot; mapq=$mapq \n";
			if( !$if_count_sencondary_map && $read_map_primary_orNot ne "primaryAli" ){ 
				next; #not primary alignment
			}
			$counts{$samplename}{"1.reads.$read_map_primary_orNot"}++;

			#uniqueness
			if($uniqueness_mapq_minSco && $mapq < $uniqueness_mapq_minSco){ #judge uniqueness by $mapq
				#print "\$uniqueness_mapq_minSco=$uniqueness_mapq_minSco; mapq=$mapq skipped!\n";
				next;
			}

			$counts{$samplename}{"2.reads.$read_map_primary_orNot.unique"}++;
			if(!($flag & 0x0002) || ($flag & 0x0008)){next;} #if not mapped in a proper pair
			$counts{$samplename}{"3.reads.$read_map_primary_orNot.properMapped"}++;
			if($flag & 0x0080){next;} #skip second mate in a pair
			$counts{$samplename}{"4.reads.$read_map_primary_orNot.mate1"}++;
			next if (!$isize || abs($isize)>=$pairend_gapmax);
			$counts{$samplename}{"5.reads.$read_map_primary_orNot.isize<$pairend_gapmax"}++;
			
			my $tag_str=join("	",@tags_arr);
			if($if_best_score){
				if($tag_str=~/XS:i:(\d+)/i){
					$xs=$1;
					if($tag_str=~/AS:i:(\d+)/i){
						$as=$1;
						next if $as<$xs;
					}
				}
				$counts{$samplename}{"6.reads.$read_map_primary_orNot.BestAliSco"}++;
			}

			if(! judge_reads_quality($OutSamFrom,$mismatch_cut,$clip_cut,$cigar,@tags_arr)){
				next;
			}
			$counts{$samplename}{"7.reads.$read_map_primary_orNot.mismatch<=$mismatch_cut.clip<=$clip_cut"}++;

			$mate1Strand=from_flag2strand($flag);
			if($mate1Strand eq "+"){
				$frag1_str="$rname	$mate1Strand	$posi	".($posi+$isize-1)."	1\n";
			}else{
				$frag1_str="$rname	$mate1Strand	$mpos	".($mpos-$isize-1)."	1\n";
			}
			print MAPOUT $frag1_str;
			$readstype2num_hash{"PE	".abs($isize) }+=1;
			$rname2rnum_h{$rname}++;

		}
		close SAM_F;
		close MAPOUT;		
		
		##write RNAME_RNUM_OUT
		print RNAME_RNUM_OUT "rname	num_$samplename\n";
		foreach my $rname (sort keys %rname2rnum_h){
			print RNAME_RNUM_OUT "$rname	".$rname2rnum_h{$rname}."\n";
		}
		close RNAME_RNUM_OUT;

		##output PE_readsType_outf  
		open (READTYPE_NUM_OUT,">$readsCov_outdir/$outPrefix$samplename.reads_type.num.txt") || die "error write $readsCov_outdir/$outPrefix$samplename.reads_type.num.txt\n";
		print READTYPE_NUM_OUT "reads_map_type	reads_cov_len	number\n"; 
		print READTYPE_NUM_OUT join("", map("$_	$readstype2num_hash{$_}\n", sort keys %readstype2num_hash)). "\n";
		close READTYPE_NUM_OUT;
		
		if($out_unmapped_reads_fq_f){
			for my $qname (sort keys %unmapped_reads_h ){
				if($unmapped_reads_h{$qname}{1} && $unmapped_reads_h{$qname}{2}){
					print OUT_UNMAPPED_FQ_F1 $unmapped_reads_h{$qname}{1};
					print OUT_UNMAPPED_FQ_F2 $unmapped_reads_h{$qname}{2};
					$counts{$samplename}{"1.reads.bothMateUnMapped.output"}++;
				}
				# $unmapped_reads_h{$qname}{$read1or2}='@'."$qname\n$seq\n+\n$qual\n";
			}
			close OUT_UNMAPPED_FQ_F1;
			close OUT_UNMAPPED_FQ_F2;
		}

		####output numbers:
		my %all_readGrps_h=();
		foreach my $sample1( keys %counts){
			map($all_readGrps_h{$_}=1, keys %{$counts{$sample1}});
		}
		my @all_readGrps=sort keys %all_readGrps_h;
		$print_txt="Sample	".join("	",@all_readGrps)."\n";
		foreach my $sample1( keys %counts){
			$print_txt.="$sample1	". join("	",map($counts{$sample1}{$_}, @all_readGrps))."\n";
		}

		print $print_txt;
		print LOGF $print_txt;
		close LOGF;

	}#end else

}

