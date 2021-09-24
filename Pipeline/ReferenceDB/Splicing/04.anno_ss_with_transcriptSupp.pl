##annotate splice site with transcript support number (eg. Refseq, UCSC known genes, Ensembl, mRNA, EST)
#need to download files from ucsc genome browser
#need relative large memory to run (~17G for human)

require "../../SharedCodes/perl.fun.inc.pl";


@use_model_arr=("refFlat","ensGene","mrna","est"); $geno="rheMac2";
@use_model_arr=("refFlat","knownGene","ensGene","mrna","est"); $geno="mm9"; #hg19 mm9
@use_model_arr=("refFlat","ensGene"); $geno="canFam3"; #rheMac10 rn6 canFam3
$date_db="2102";
$contig_file_dir="../../../data/ucsc/genomes/$geno";
$contig_allseq_f="";

$gene_str_root="../ucsc";
$min_intron_len=20;
$debug=0;
my %gstr_f_h=(
	"refFlat"=>"$gene_str_root/refFlat/$geno/refFlat.txt",
	"knownGene"=>"$gene_str_root/knownGene/${geno}/knownGene.txt",
	"ensGene"=>"$gene_str_root/ensGene/${geno}/ensGene.txt",
	"mrna"=>"$gene_str_root/mRNA/${geno}_$date_db/all_mrna.txt",
	"est"=>"$gene_str_root/est/${geno}_$date_db/all_est.txt"
);
my %Gmodl2header_h=(
	"refFlat"=>"gene_symbol|name|contig|strand|transc_start|transc_end|cds_start|cds_end|exon_num|exon_starts|exon_ends",
	"knownGene"=>"name|contig|strand|transc_start|transc_end|cds_start|cds_end|exon_num|exon_starts|exon_ends|proteinID|alignID",
	"ensGene"=>"bin|name|contig|strand|transc_start|transc_end|cds_start|cds_end|exon_num|exon_starts|exon_ends|score|name2|cdsStartStat|cdsEndStat|exonFrames",
	"mrna"=>"bin|matches|misMatches|repMatches|nCount|qNumInsert|qBaseInsert|tNumInsert|tBaseInsert|strand|name|qSize|qStart|qEnd|contig|tSize|tStart|tEnd|exon_num|blockSizes|qStarts|exon_starts",
	"est"=>"bin|matches|misMatches|repMatches|nCount|qNumInsert|qBaseInsert|tNumInsert|tBaseInsert|strand|name|qSize|qStart|qEnd|contig|tSize|tStart|tEnd|exon_num|blockSizes|qStarts|exon_starts"
);
#need contig, strand, exon_starts, exon_ends(or blockSizes) for all gene models

$out_root="04.ss_supp/${geno}_$date_db/supp";
$log_f="$out_root.log";

##RUN
print "#####debug=$debug#########\n";
create_dir_ifNotExist($out_root);
open (LOGF, ">$log_f") || die "error $log_f\n";

##1, read transcript structure and load hash %ss2supp_h
my %ss2supp_h;
foreach my $geneModel(@use_model_arr){
	my $geneModel_f=$gstr_f_h{$geneModel};
	@headername_arr=split(/\|/, $Gmodl2header_h{$geneModel});
	open (GENE_ANO, $geneModel_f) || die "error open $geneModel_f\n";
	print "#open $geneModel_f\n";
	print "#headers=".join(" ",@headername_arr).";\n";
	print LOGF "#open $geneModel_f\n";
	print LOGF "#headers=".join(" ",@headername_arr).";\n";
	while(<GENE_ANO>){
		chomp;
		last if ($debug && $counts{"input.gene	$geneModel	1.row"}>1000);
		$counts{"input.gene	$geneModel	1.row"}++;
		$exon_ends="";
		@temp_arr = split(/\t/);
		for($i=0;$i<@headername_arr;$i++){
			${$headername_arr[$i]}=$temp_arr[$i];
			#print "$headername_arr[$i]=$temp_arr[$i];\n"
		}
		next if $exon_num<2;
		next if !($strand eq "+" || $strand eq "-");
		$counts{"input.gene	$geneModel	2.spliced"}++;
		my @exon_starts_arr=split(/,/,$exon_starts);
		
		my @exon_ends_arr=();
		if(!$exon_ends && $exon_starts && $blockSizes){ #for mrna and est file format, 
			my @blockSizes_arr=split(/,/,$blockSizes);
			foreach my $ei(0..(scalar @exon_starts_arr-1)){
				push (@exon_ends_arr, $exon_starts_arr[$ei]+$blockSizes_arr[$ei]);
			}
		}else{
			@exon_ends_arr=split(/,/,$exon_ends);
		}
		if ( scalar @exon_ends_arr != scalar @exon_starts_arr ){die "error exon_ends not defined or error for $geneModel $name;";}

		map($exon_starts_arr[$_]++,  0..(scalar @exon_starts_arr-1) ); #convert to 1-based coordinates
		if($strand eq "-"){
			my @tmparr=@exon_starts_arr;
			@exon_starts_arr=reverse(@exon_ends_arr);
			@exon_ends_arr=reverse(@tmparr);
		}
		#print "scalar \@exon_starts_arr=".(scalar @exon_starts_arr).";\n";
		my $gene_strand="";
		my $if_reversedStrand="";
		foreach my $exon_i( 1.. (scalar @exon_starts_arr-1)){
			#print "exon_starts_arr=".join(";",@exon_starts_arr)."\nexon_ends_arr=".join(";",@exon_ends_arr)."\n";
			my $ss3=$exon_starts_arr[$exon_i]; #upstream intron ss
			my $ss5=$exon_ends_arr[$exon_i-1];
			die "error $geneModel $name ss5=$ss5\n" if !$ss5;
			die "error $geneModel $name ss3=$ss3\n" if !$ss3;
			next if(abs($ss3-$ss5)<$min_intron_len && $geneModel=~/est|mrna/);
			if($geneModel=~/est|mrna/ ){
				if($if_reversedStrand){
					($ss5,$ss3)=($ss3,$ss5);
				}
				if( !$strand || sum(map($ss2supp_h{"ss5"}{"$contig	$strand	$ss5"}{$_}, (@use_model_arr)))<1 && 
					sum(map($ss2supp_h{"ss3"}{"$contig	$strand	$ss3"}{$_}, (@use_model_arr)))<1 ){
					#judge strandness of mRNA or est based on splice site sequences
					my ($exonPos1,$exonPos2)=sort {$a<=>$b} ($ss5,$ss3);
					my $ss2nt1=($exonPos1+1).":".($exonPos1+2);
					my $ss2nt2=($exonPos2-2).":".($exonPos2-1);
					my $spliceSite_nts=extract_seq_fr_contig($contig_file_dir, $contig, '+', "$ss2nt1|$ss2nt2", $contig_allseq_f);
					if( $spliceSite_nts=~/GTAG|GCAG|ATAC/ ){
						$if_reversedStrand=1 if ($strand eq "-");
						next if ($gene_strand && $gene_strand eq "-");
						$strand="+";$gene_strand="+";
						($ss5,$ss3)=($exonPos1,$exonPos2);
					}elsif(reverse_complement($spliceSite_nts)=~/GTAG|GCAG|ATAC/ ){
						$if_reversedStrand=1 if ($strand eq "+");
						next if ($gene_strand && $gene_strand eq "+");
						$strand="-";$gene_strand="-";
						($ss5,$ss3)=($exonPos2,$exonPos1);
					}else{
						$strand="";
					}
				}
			}
			next if ($strand eq "");
			$ss2supp_h{"junc"}{"$contig	$strand	$ss5	$ss3"}{$geneModel}++;
			$counts{"input.intron	$geneModel"}++;
			$ss2supp_h{"ss5"}{"$contig	$strand	$ss5"}{$geneModel}++;
			$ss2supp_h{"ss3"}{"$contig	$strand	$ss3"}{$geneModel}++;
		}
	}

}

##2, write output 

foreach my $output_item("junc","ss5","ss3"){
	my $out_supp_file="$out_root.$output_item.tbl";
	print "#write $out_supp_file\n";
	print LOGF "#write $out_supp_file\n";
	open (SS_SUPP_OUT, ">$out_supp_file") || die "error write $out_supp_file\n";
	print SS_SUPP_OUT "contig	strand	". ( $output_item eq "junc" ? "ss5	ss3":$output_item )."	".
		join("	", map( "supp_$_",  @use_model_arr ) )."\n";
	foreach my $ss_key (sort keys %{$ss2supp_h{$output_item}}){
		next if sum(map($ss2supp_h{$output_item}{$ss_key}{$_}, (@use_model_arr)) )<1 ;
		$counts{"output	$output_item"}++;
		print SS_SUPP_OUT "$ss_key	".
			join("	",map( $ss2supp_h{$output_item}{$ss_key}{$_},  @use_model_arr )). 
			"\n";
	}
	close SS_SUPP_OUT;
}


print "counts:\n". join("", map("$_	$counts{$_}\n", sort keys %counts)). "\n";
print LOGF "counts:\n". join("", map("$_	$counts{$_}\n", sort keys %counts)). "\n";
close LOGF;
