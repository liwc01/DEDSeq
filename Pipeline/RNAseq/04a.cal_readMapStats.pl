#calculate read number for different mapping regions:
%readsTypes2sco_h=( #higher score represent higher priority (applicable to RNA-seq); + - represent strand relative to reference gene only useful for directional reads 
	"CDS_S"=>18,
	"CDS"=>17,
	"3UTR_S"=>16,
	"3UTR"=>15,
	"5UTR_S"=>14,
	"5UTR"=>13,
	"exon_S"=>12.9,
	"exon"=>12.8,
	"intron_S"=>12,
	"intron"=>11,	
	"UTR3e_S"=>10,
	"UTR3e"=>9,
	"intron_A"=>8,
	"UTR5e_A"=>7,
	"UTR5e"=>6,
	"UTR5e_S"=>5,
	"UTR3e_A"=>4,
	"3UTR_A"=>3,
	"CDS_A"=>2,
	"exon_A"=>1.9,
	"5UTR_A"=>1,
	"intergenic"=>0.1,
);

#perl 1cal_readnum.pl -o HuR_CLIP

use Getopt::Std; 
getopt("op",\%args);
$outdir=$args{o} if $args{o};
$prefix=$args{p} if $args{p};
$out_f="$outdir/${prefix}ReadTypeCnt.tbl";
die "\$outdir (-o XXXX) not defined\n" if !$outdir;

opendir (INDIR,"$outdir/") || die "error open dir $outdir\n";
my %all_readTypes=();
while( ($filename = readdir(INDIR))){
	next if ($filename!~/$prefix.*\.log/);
	$sample=$filename;
	$sample=~s/$prefix\.*//;
	$sample=~s/\.ReadNum\.tbl\.temp\.log//;
	print("open $filename sample=$sample\n");
	open(IN, "$outdir/$filename");
	my $if_read_counts=0;
	while(<IN>){
		s/\s+$//;
		next if !$_;
		next if /^ /;
		if(/^counts/){
			$if_read_counts=1;
			next;
		}
		next if !$if_read_counts;
		my ($map_type,$num)=split(/\t/);
		$map_type=~s/reads\.*//;
		next if !$readsTypes2sco_h{$map_type} && $prefix=~/refseqcds_PE/ ;
		$counts{"$sample"}{$map_type}=$num;
		$all_readTypes{$map_type}++;
	}
	close IN;
}
closedir(INDIR);



@readTypes=sort { $readsTypes2sco_h{$b} <=> $readsTypes2sco_h{$a} } keys %all_readTypes;
print join("	",@readTypes)."\n";
open (OUT_TB,'>', "$out_f") || die "error write $out_f\n";
my $out_str="sample	".join("	",@readTypes)."\n";
foreach my $sample(sort keys %counts){
	$out_str.= "$sample	".join("\t", map($counts{$sample}{$_}, @readTypes))."\n";
}
print $out_str;
print OUT_TB $out_str;
close OUT_TB;

