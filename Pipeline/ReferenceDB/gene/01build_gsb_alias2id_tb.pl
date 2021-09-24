#use gene information downloaded from NCBI
#build a table with all gene alias to gene id, gene official symbol table

use Getopt::Std; 
$geneinfo_file="../ncbi/gene/geneinfo/All_Mammalia.gene_info";

$taxid=9606; $organism="human"; $geneinfo_file="../ncbi/gene/geneinfo/Homo_sapiens.gene_info";
$taxid=10090; $organism="mouse"; $geneinfo_file="../ncbi/gene/geneinfo/Mus_musculus.gene_info";
$taxid=9544; $organism="monkey"; $geneinfo_file="../ncbi/gene/geneinfo/All_Mammalia.gene_info";
$taxid=10116; $organism="rat";  $geneinfo_file="../ncbi/gene/geneinfo/All_Mammalia.gene_info";
$taxid=9615; $organism="dog";  $geneinfo_file="../ncbi/gene/geneinfo/All_Mammalia.gene_info";

#example:
#perl 01build_gsb_alias2id_tb.pl -t 9606 -o human -i ../ncbi/gene/geneinfo/Homo_sapiens.gene_info
#perl 01build_gsb_alias2id_tb.pl -t 10090 -o mouse -i ../ncbi/gene/geneinfo/Mus_musculus.gene_info
#perl 01build_gsb_alias2id_tb.pl -t 10116 -o rat -i ../ncbi/gene/geneinfo/Rattus_norvegicus.gene_info
#perl 01build_gsb_alias2id_tb.pl -t 9544 -o monkey -i ../ncbi/gene/geneinfo/All_Mammalia.gene_info
#perl 01build_gsb_alias2id_tb.pl -t 9615 -o dog -i ../ncbi/gene/geneinfo/Canis_familiaris.gene_info

getopt("toi",\%args);
$taxid=$args{t} if ($args{t}); 
$organism=$args{o} if ($args{o}); 
$geneinfo_file=$args{i} if ($args{i}); 

$out_file="01gene_alias2id/$organism.$taxid.tbl";


##1, load gene info, output
system("mkdir 01gene_alias2id");
open (OUT, ">$out_file") || die "error write $out_file\n";
print OUT "alias	gene_id	gene_symbol	gene_desc\n";
open (GENEINFO,$geneinfo_file) || die "error open $geneinfo_file\n";
print " open $geneinfo_file\n output $out_file\n";
@gene_alias_arr=();
while(<GENEINFO>){
	next if ($_!~/^$taxid/);
	chomp;
	($taxid,$gene_id,$gene_symbol,$temp4,$alias,$temp6,$temp7,$temp8,$gene_desc)=split(/\t/); #9606 1 A1BG - A1B|ABG|DKFZp686F0970|GAB|HYST2477
	next if ($gene_symbol eq "-");
	$counts{"num.gene_symbol"}++;
	print OUT "$gene_symbol	$gene_id	$gene_symbol	$gene_desc\n";
	foreach $alias1(split(/\|/,$alias)){
		next if ($alias1 eq "-");
		push (@gene_alias_arr, "$alias1	$gene_id	$gene_symbol	$gene_desc\n"); 
		$counts{"num.gene.alias1"}++;
	}
}
close GENEINFO;

foreach $rowi(@gene_alias_arr){
	print OUT $rowi;
}
close OUT;

#print numbers
foreach $count_key(sort keys %counts){
	print "$count_key	$counts{$count_key}\n";
}


