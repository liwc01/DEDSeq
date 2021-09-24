#extract info for one species


$gene2accession_f="gene2accession";

foreach $tax_id(9606, 10090, 10116, 9031, 9544){ #hs 9606, mm 10090, rn 10116, galGal 9031, rheMac 9544  #
	$outf="$gene2accession_f.$tax_id";
	open (GENE2ACC,$gene2accession_f) || die "error open gene2accession file: $gene2accession_f\n";
	open (OUT,">$outf") || die "error open gene2accession file: $gene2accession_f\n";
	while(<GENE2ACC>){
		next if /^#/;
		my ($tax_id2,$gid,$status,$RNA_nucleotide_accession,@oths)=split(/\t/); #3702 814629  REVIEWED NM_126167.1 18379136 ...
		next if ($tax_id2 ne $tax_id);
		print OUT $_;
		next if ($RNA_nucleotide_accession eq "-");
		$RNA_nucleotide_accession=~s/\..*$//; #remove version, ## AAX38516.1
		my $gsb=$gid2gsb{$gid};
		if($gsb){
			$acc2gid_gsb{$RNA_nucleotide_accession}="$gid:$gsb";
			$cnt{'# of acc2gid_gsb taken'}++;
		}
		#if($RNA_nucleotide_accession =~ /NM_033324|NM_011239/){ #debug
		#	print "$tax_id2,$gid,$status,$RNA_nucleotide_accession; $gsb; $acc2gid_gsb{$RNA_nucleotide_accession}\n";
		#}
	}
	close GENE2ACC;
	close OUT;  
}
