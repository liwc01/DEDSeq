#calculate read number for different types:

use Getopt::Std; 
getopt("o",\%args);
$readTableOut_dir=$args{o} if $args{o};


$inputfs="$readTableOut_dir/*/readsIDs.counts.txt";
$out_f="$readTableOut_dir/readnum.tbl";
die "\$readTableOut_dir (-o XXXX) not defined\n" if !$readTableOut_dir;

open (IN,"cat $inputfs |");

my %all_readTypes=();
while(<IN>){
	chomp;
	next if !$_;
	next if /event\./;
	next if /mapping\.rows/;
	
	my ($item,$num)=split(/\t/);
	$item=~s/num\.//;
	$item=~s/reads\.//;
	#$item=~s/([^\.]+)\.//;
	#$sample=$1; $map_type=$item;
	$item=~s/(\d+\.unmapped|\d+\.non-unique|\d+\.unique[^\s]+)//;
	$map_type=$1; $sample=$item;


	if($map_type=~s/readstype=(.+)//){
		my $cigar=$1;
		my $type=$cigar=~/N|junc/ ? "Junction": "ExonBody";
		$map_type.=$type;
	}	
	$counts{"$sample"}{$map_type}+=$num;
	$all_readTypes{$map_type}++;
	#if($map_type=~/non_unique/ || $map_type=~/non\-unique/){
	#	$counts{"$sample"}{"non-unique"}=$num;
	#}elsif($map_type=~s/\.uni\.readstype=//){
	#	my($sample,$cigar)=split(/\./,$sample);
	#	my $type=$cigar=~/N/ ? "Junction": "ExonBody";
	#	$counts{"$sample"}{"unique.$type"}+=$num;
	#	$counts{"$sample"}{"unique"}+=$num;
	#}
}
close IN;

@readTypes=("non-unique","unique","unique.ExonBody","unique.Junction");
@readTypes=sort keys %all_readTypes;
open (OUT_TB,">$out_f") || die "error write $out_f\n";
my $out_str="sample	".join("	",@readTypes)."\n";
foreach my $sample(sort keys %counts){
	my $sample_name=$sample;
	$sample_name=~s/\.$//;
	$out_str.= "$sample_name	".join("\t", map($counts{$sample}{$_}, @readTypes))."\n";
}
print $out_str;
print OUT_TB $out_str;
close OUT_TB;


#eg: format in input file
# MsLVMR.num.alievent.readstype=M25N      1
# MsLVMR.num.event.readmap        861094
# MsLVMR.num.reads.non_unique.mapped      122434
# MsLVMR.num.reads.uni.readstype=16M      17
# MsLVMR.num.reads.uni.readstype=M25N     1
# MsLVMR.num.reads.unique.mapped  94435
