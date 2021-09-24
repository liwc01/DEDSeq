#convert gene block to bed format then to bigbed format used for visualization

require "../SharedCodes/perl.fun.inc.pl";
use Getopt::Std; 



getopt("gabB",\%args);
$geno=$args{g} ? $args{g}:"hg19";
$GeneBlock_anof=$args{a} ? $args{a} : "02.geneBlocks/geneBlock.flat.tbl.ano";
$out_bed_f=$args{b} ? $args{b} : "$GeneBlock_anof.bed";
$out_bigbed_f=$args{B} ? $args{B} : "$GeneBlock_anof.bb";


##RUN
open (ANO,$GeneBlock_anof) || die "error open $GeneBlock_anof\n";
open (BED,">$out_bed_f") || die "error write $out_bed_f\n";
while(<ANO>){
	chomp;
	if(/region_ano/i){ #eg: contig  strand  refseqid   row  Gblock_id   start_pos   end_pos startPosType  startPosAno  endPosType  endPosAno reg_len region_ano
		@headername_arr=split(/\t/);
		print $_."\n";
	}else{
		@temp_arr = split(/\t/);
		for($i=0;$i<@headername_arr;$i++){
			${$headername_arr[$i]}=$temp_arr[$i];
		}
		my($leftPos,$rightPos)=sort {$a<=>$b} ($start_pos, $end_pos);
		$leftPos--;
		my $name="$region_ano:$startPosAno->$endPosAno($reg_len)";
		my $thickStart=$leftPos; #exonic region will show thicker
		my $thickEnd=$region_ano=~/exon/i ? $rightPos: $leftPos; #exonic region will show thicker
		print BED "$contig	$leftPos	$rightPos	$name	1	$strand	$thickStart	$thickEnd\n";
	}

}

close ANO;
close BED;

#convert to bigbed
my $cmd1="sort -k 1,1d -k 2,2n $out_bed_f > $out_bed_f.sorted"; print "run: $cmd1;\n";
system($cmd1);
my $chrSize_file="../Software/UCSC/$geno.chrom.sizes";
my $cmd2="bedToBigBed $out_bed_f.sorted $chrSize_file $out_bed_f.sorted.bb"; print "run: $cmd2;\n";
system($cmd2);


