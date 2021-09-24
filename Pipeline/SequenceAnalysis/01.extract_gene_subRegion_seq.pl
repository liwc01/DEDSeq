###extract gene sub-region sequences from genome
#all sequence will be output, including redundent ones with different transcript ID
#7/8/2015 avoild output duplicate sequences
#7/8/2015 can also output a table with id and sequence
#1/29/2020 can set both $refPos_name, $refPos_name2 and $extractSeq_Fr,$extractSeq_To (eg. to extract snRNA sequences with flanking +/- 100 bp )

use Getopt::Std; 
getopt("rdgsiencftoCRbETS",\%args);
require "../SharedCodes/perl.fun.inc.pl";


### below are several examples:
my $geno="hg19";
@fa_id_headers=("name"); #using transcript id as id in fasta file
@fa_id_headers=(); #using default (">$gene_symbol	$name $$refPos_name	$contig|$strand|$regions_str")


##extract trancript start upstream region: TSS.m2000m1 TrEnd.p101p2100
###utr3; utr5; cds; exons; introns;  exon1cds(CDS in first exon); exon2utr5 (utr5 in non-first exon)
$geneStru_f="22geneStru/$geno.refFlat.stru.tbl";
$out_region_name="cds"; #header name in $geneStru_f
$refPos_name="";
$out_seq_f="23gene_subRegSeq/$geno/refFlat.$out_region_name.fa";
if($out_region_name=~/utr3|exon/){ 
	$filter_header='if_last_utr3'; $filterSel_val='1';
}else{
	$filter_header='';
}

##common
$geneStru_f="22geneStru/$geno.refFlat.stru.tbl";
$out_seq_f="23gene_subRegSeq/$geno/refFlat.$out_region_name.fa";
$if_out_seqtb=0;

# example of command line:
#  01.extract_gene_subRegion_seq.pl -d mm_stemCell -g mm9 -s ../../analyze/polyA/G_pA_fr_reads/20pAdistr/mm_stemCell/gene_tb.tbl \
#   -i club_id -e ini_transcript_start -n TSS.m4kp4k -c m4000p4000 

$db_name=$args{d} ? $args{d}: "";
$geno=$args{g} if ($args{g}); 
$contig_file_dir="../ReferenceDB/ucsc/genomes/$geno";
$geneStru_f=$args{s}?$args{s}:""; 
@fa_id_headers=$args{i}?split(/\s+/, $args{i}): (); 
$refPos_name=$args{e}?$args{e}:""; #first reference position head name
$refPos_name2=$args{E}?$args{E}:""; #second reference position head name
$out_region_name=$args{n}?$args{n}:"";
$region_header=$args{R}?$args{R}:$out_region_name;
($extractSeq_Fr,$extractSeq_To)=("","");
if($args{c}){ ($extractSeq_Fr,$extractSeq_To)=convert_regName2_coordinates($args{c}); }
$filter_header=$args{f}?$args{f}:"";
$filterSel_val=$args{t}?$args{t}:"";
$contig_header=$args{T}?$args{T}:"contig";
$strand_header=$args{S}?$args{S}:"strand";
$out_seq_f=$args{o}?$args{o}:"01.gene_subRegSeq/$db_name/$out_region_name.fa";
$contig_f=$args{C} ? $args{C} : "";
$if_out_seqtb=$args{b} if $args{b} ne "";



#################### RUN #################
create_dir_ifNotExist($out_seq_f);

#1, open gene structure file, get sequence
open (GENE_STR, $geneStru_f) || die "error open $geneStru_f";
open (OUTSEQ, ">$out_seq_f") || die " error write $out_seq_f\n";
open (OUTSEQ_TBL, ">$out_seq_f.tbl") if $if_out_seqtb;
print " open $geneStru_f\n write $out_seq_f\n";
my %treated_ids=();
while(<GENE_STR>){
	chomp;
	if(/gene_symbol/ || ($refPos_name && /\b$refPos_name\b/ )  ){ #gene_symbol     contig  strand name  transc_start  transc_end  cds_start cds_end exon_num  utr5 utr3 cds exons  introns gene_span transcript_len  utr5_len        cds_len utr3_len        intron_len      If_discard_smaller_utr3 if_last_utr3
		s/\bchromosome\b/contig/i;
		@headername_arr=split(/\t/);
		print "header line: $_;\n";
	}else{ #GUCY1A3 chr4    +       NM_001130687    156807312       156862835       156837470       156862798       9       156807312-156807567|156807948-156808021|156837358-156837469  156862799-156862835     156837470-156837724|156844534-156844595
		$counts{"num.transcript.input"}++;
		if($counts{"num.transcript.input"} % 5000 ==0){
			print " treat transcript ".$counts{"num.transcript.input"}."\n";
		}
		@temp_arr = split(/\t/);
		for($i=0;$i<@headername_arr;$i++){
			${$headername_arr[$i]}=$temp_arr[$i];
		}
		if($filter_header){
			next if $$filter_header!~/$filterSel_val/;
			$counts{"num.$filter_header=~/$filterSel_val/"}++;
		}
		my $strand_s;
		if($$strand_header =~ /\-/){
			$$strand_header="-"; $strand_s=-1;
		}else{
			$$strand_header="+"; $strand_s=1;
		}

		if($refPos_name && $$refPos_name ne ""){
			if($refPos_name2 && $$refPos_name2 ne ""){
				if($extractSeq_Fr ne ""){
					$regions_str=($$refPos_name + $extractSeq_Fr*$strand_s).":".($$refPos_name2 + $extractSeq_To*$strand_s);
				}else{
					$regions_str="$$refPos_name:$$refPos_name2";
				}
				
			}elsif($extractSeq_Fr ne ""){
				$regions_str=get_genomic_region($$refPos_name, $strand_s, $extractSeq_Fr,$extractSeq_To);
			}
		}elsif($out_region_name eq 'exon1cds'){
			if($utr5=~/\|/ || $exons=~/$utr5/ || !$cds){ 
				$regions_str="";
			}else{
				$regions_str=$cds;
				if($$strand_header eq "+"){ $regions_str=~s/\|.*//; }else{ $regions_str=~s/.*\|//; }
			}
		}elsif($out_region_name eq 'exon2utr5'){
			if($utr5=~/\|/){ 
				$regions_str=$utr5;
				if($$strand_header eq "+"){ $regions_str=~s/[^\\]*\|//; }else{ $regions_str=~s/\|[^\\]*//; }
			}else{
				$regions_str="";
			}
			
		}else{
			$regions_str=$$region_header;
		}
		my $seq_id="";
		if(scalar @fa_id_headers==0){
			$seq_id="$gene_symbol	$name $$refPos_name	$$contig_header|$$strand_header|$regions_str";
		}else{
			$seq_id=join(":",map($$_,@fa_id_headers));
		}
		next if ( $treated_ids{$seq_id} );
		$treated_ids{$seq_id}=1;

		#print "$$strand_header utr5=$utr5 regions_str=$regions_str cds=$cds\n";
		if(!$regions_str){
			$counts{"num.regions_str_notDefined"}++;
		}elsif($regions_str){ # && 
			$counts{"num.seq.output"}++;
			if(!$all_out_regions{"$$contig_header|$$strand_header|$regions_str"}) {$counts{"num.seq.output.unique"}++;}
			my $contig_f_used="$contig_file_dir/$geno.fa";
			if(! -e $contig_f_used){$contig_f_used=$contig_file_dir/$$contig_header.fa;}
			my $subseq=extract_seq_fr_contig($contig_file_dir,$$contig_header,$$strand_header,$regions_str, $contig_f?$contig_f:$contig_f_used);
			print OUTSEQ_TBL "$seq_id\t$subseq\n" if $if_out_seqtb;
			$subseq=~s/(\S{50})/$1\n/gi;
			print OUTSEQ ">$seq_id\n$subseq\n"; 
			$all_out_regions{"$$contig_header|$$strand_header|$regions_str"}=1;
		}
	}
}
close OUTSEQ;
close GENE_STR;
close OUTSEQ_TBL if $if_out_seqtb;

print "counts:\n". join("", map("$_	$counts{$_}\n", sort keys %counts)). "\n";
