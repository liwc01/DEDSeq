###get the flanking sequence for splicing site (5' and 3') and output bed format and fasta sequences:
#in UCSC table browser: for 5'ss, add 3nt (exon, total 3) at 5'end and 6nt in the 3' end
#for 3'ss, add 20nt (intron, total 20) at 5'end and 3nt(exon, total 3) in the 3' end

##9/11/2015 change coordinates to 1 based exonic coordinates; Also can use perl to extract splice site flanking sequences
##9/11/2015 output intron annotation file instead of bed file; output unique splice site fasta file
##11/29/2016 output exon annotation file

$ss5_exonLen=3; $ss5_intronLen=6; 
$ss3_exonLen=3; $ss3_intronLen=20; 

#species
$geno="hg19";  #hg19, mm9, rn6, rheMac2 canFam3
$gene_model="ensGene";$id_header="ensid"; # ensGene
$gene_model="refflat";$id_header="refseqid"; #refflat 
$transcript_desc_f="../gene/02transcript_gene_ano/$geno.$gene_model.desc.txt";

$out_intron_ano_file="01.splicingSite/$geno.$gene_model.intron.ano.txt";
$out_exon_ano_file="01.splicingSite/$geno.$gene_model.exon.ano.txt";
$out_ss5_fafile="01.splicingSite/$geno.$gene_model.ss5.fa";
$out_ss3_fafile="01.splicingSite/$geno.$gene_model.ss3.fa";
$log_f="$out_intron_ano_file.log";

$contig_file_dir="../ucsc/genomes/$geno";
$contig_file="";
$contig_file="../ucsc/genomes/$geno/$geno.fa" if ($geno=~/rn6|canFam3/); #for rn6, canFam3 ...
require "../SharedCodes/perl.fun.inc.pl";



###RUN
create_dir_ifNotExist($log_f);
open (LOGF,">$log_f")  || die "error write $log_f\n";
print  "#open $transcript_desc_f\n#write $out_intron_ano_file, $out_exon_ano_file, $out_ss5_fafile, $out_ss3_fafile\n";
print  LOGF "#open $transcript_desc_f\n#write $out_intron_ano_file, $out_exon_ano_file, $out_ss5_fafile, $out_ss3_fafile\n";

open (OUT_INTRON,">$out_intron_ano_file")  || die "error write $out_intron_ano_file\n";
open (OUT_EXON,">$out_exon_ano_file")  || die "error write $out_exon_ano_file\n";
open (OUTSS5FA,">$out_ss5_fafile")  || die "error write $out_ss5_fafile\n";
open (OUTSS3FA,">$out_ss3_fafile")  || die "error write $out_ss3_fafile\n";
open (GENE_DEFI,$transcript_desc_f) || die "error $transcript_desc_f\n";
print OUT_INTRON "gene_symbol	genef_row	$id_header	contig	strand	introni	intron_revnum	ss5	ss3	ss5_seq	ss3_seq\n";
print OUT_EXON "gene_symbol	genef_row	$id_header	contig	strand	exonID	exon_revID	exon_num	ss5	ss3	ss5_seq	ss3_seq\n";
$genef_row=0;
my %ss5Seq_hash=();
my %ss3Seq_hash=();

print "headers= ... \n";
while(<GENE_DEFI>){
	chomp;

	if(/gene_symbol/i){ #gene_symbol refseqid contig strand transc_start transc_end cds_start cds_end exon_num exon_starts exon_ends gene_id alias gene_desc
		#bin     ensid   contig  strand  transc_start    transc_end      cds_start       cds_end exon_num        exon_starts     exon_ends       score   ensGeneId       cdsStartStat    cdsEndStat      exonFrames      gene_symbol     gene_Biotype    gene_desc
		@headername_arr=split(/\t/);
		$genef_row=0;
		print join(" ",@headername_arr)."\n";
	}else{ #DLX1 NM_178120 chr2 + 172658453 172662647 172658651 172661231 3 172658453,172659627,172660976, 172658964,172659827,172662647, 1745 - distal-less homeobox 1
		$genef_row++;
		@temp_arr = split(/\t/);
		for($i=0;$i<@headername_arr;$i++){
			${$headername_arr[$i]}=$temp_arr[$i];
		}
		$counts{"input.transcript"}++;
		next if ($exon_num==1);
		@exon_starts_arr=split(/,/,$exon_starts);
		@exon_ends_arr=split(/,/,$exon_ends);
		my @intron_arr=();
		my @exon_arr=();
		foreach $i(0..(@exon_starts_arr-1)){
			push (@exon_arr,($exon_starts_arr[$i]+1,$exon_ends_arr[$i])); #convert to 1-based exon coordinates 
		}
		@intron_arr=@exon_arr[1..(@exon_arr-2)]; #remove TSS and TES
		if($strand eq "-"){
			@intron_arr=reverse(@intron_arr);
			@exon_arr=reverse(@exon_arr);
		}
		$strand_s=$strand eq "-"?-1:1;
		foreach my $introni(1..($exon_num-1)){
			$counts{"input.intron"}++;
			$intron_revnum=$introni-$exon_num;
			my $ss5=shift(@intron_arr);
			my $ss3=shift(@intron_arr);
			my $ss5_region_str=get_genomic_region($ss5, $strand_s, -$ss5_exonLen+1,$ss5_intronLen);
			my $ss3_region_str=get_genomic_region($ss3, $strand_s,-$ss3_intronLen, $ss3_exonLen-1);
			my $ss5_region_seq=$ss5Seq_hash{"$contig:$strand:$ss5"} ? $ss5Seq_hash{"$contig:$strand:$ss5"} : extract_seq_fr_contig($contig_file_dir,$contig,$strand,$ss5_region_str,$contig_file) ;
			my $ss3_region_seq=$ss3Seq_hash{"$contig:$strand:$ss3"} ? $ss3Seq_hash{"$contig:$strand:$ss3"} : extract_seq_fr_contig($contig_file_dir,$contig,$strand,$ss3_region_str,$contig_file) ;
			print OUT_INTRON "$gene_symbol	$genef_row	$$id_header	$contig	$strand	$introni	$intron_revnum	$ss5	$ss3	$ss5_region_seq	$ss3_region_seq\n";
			$ss5Seq_hash{"$contig:$strand:$ss5"}=$ss5_region_seq;
			$ss3Seq_hash{"$contig:$strand:$ss3"}=$ss3_region_seq;
		}
		foreach my $exonID (1..$exon_num){
			$exon_revID=$exonID-$exon_num -1;
			my $ss3=shift(@exon_arr);
			my $ss5=shift(@exon_arr);
			my $ss5_region_seq=$ss5Seq_hash{"$contig:$strand:$ss5"};
			my $ss3_region_seq=$ss3Seq_hash{"$contig:$strand:$ss3"};
			print OUT_EXON "$gene_symbol	$genef_row	$$id_header	$contig	$strand	$exonID	$exon_revID	$exon_num	$ss5	$ss3	$ss5_region_seq	$ss3_region_seq\n";
		}
	}
}
close OUT_INTRON;
close OUT_EXON;
foreach my $ss5id(sort keys %ss5Seq_hash){
	if(! ($ss5Seq_hash{$ss5id} =~ /[^ACGTacgt]/ ) ){ #discard ss with N
		print OUTSS5FA ">$ss5id\n".$ss5Seq_hash{$ss5id}."\n";
		$counts{"output.ss5"}++;
	}
}
foreach my $ss3id(sort keys %ss3Seq_hash){
	if(! ($ss3Seq_hash{$ss3id} =~ /[^ACGTacgt]/ ) ){ #discard ss with N
		print OUTSS3FA ">$ss3id\n".$ss3Seq_hash{$ss3id}."\n";
		$counts{"output.ss3"}++;
	}
}
close OUTSS5FA;
close OUTSS3FA;

print  "#counts:\n". join("", map("$_	$counts{$_}\n", sort keys %counts)). "\n";
print LOGF "#counts:\n". join("", map("$_	$counts{$_}\n", sort keys %counts)). "\n";

close LOGF;
