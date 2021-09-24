##############################functions#################
sub judge_readstype_fr_ciger{
	my $ciger_str=shift; #like 75M , 34M3634N41M , 6H44M, 25M25H
	if($ciger_str eq "*"){
		return(($ciger_str,"unmap"));
	}elsif($ciger_str=~/\d+M/ && $ciger_str!~/N/){ #should be match with no intron, allow S or H(clipping at the ends), IDP
		return(($ciger_str,"exon"));
	}else{ #like 34M3634N41M
		$match_len=0;
		while($ciger_str=~/(\d+)M/g){
			$match_len+=$1;
		}
		$ciger_str=~s/[\dM]//g; #remove all numbers and letter M
		if($ciger_str=~/N/){
			return (("M$match_len$ciger_str","junc"));
		}else{
			return (("M$match_len$ciger_str","oth"));
		}
	}
}

sub from_flag2strand { ###
	my $flag_str=shift;
	return ($flag_str & 16)?"-":"+";
}


sub from_cigar2juncinfo{ #88759659,	39M227N36M,	'NM:i:0 XS:A:+ NS:i:0' 
	my($reads5pos,$ciger_str,$tags_str)=@_; #$reads5pos is 1-based coordinate
	$ciger_str=~s/\d+[SH]//g; #remove clipping like 3H25M1603N22M
	if($ciger_str=~/I|D/){ #insertion or deletion
		$ciger_str=modify_cigar_rm_I_D($ciger_str);
	}
	my $juncpos1s="";
	my $juncpos2s="";
	my $exonreglen1s="";
	my $exonreglen2s="";
	while($ciger_str=~/(\d+)M(\d+)N(\d+)M/){
		my $ori_Str=$&;
		my $convertTo_str=($3)."M";
		my ($exonreglen1,$intronlen,$exonreglen2)=($1,$2,$3);
		my $juncpos1=$reads5pos+$exonreglen1-1; ###position on exon
		my $juncpos2=$juncpos1+$intronlen+1; ###position on exon
		$juncpos1s=$juncpos1s.($juncpos1s?";":"").$juncpos1;
		$juncpos2s=$juncpos2s.($juncpos2s?";":"").$juncpos2;
		$exonreglen1s=$exonreglen1s.($exonreglen1s?";":"").$exonreglen1;
		$exonreglen2s=$exonreglen2s.($exonreglen2s?";":"").$exonreglen2;
		
		$ciger_str=~ s/$ori_Str/$convertTo_str/;
		$reads5pos=$juncpos2;
	}
	
	$tags_str=~/XS:A:(\S)/; ###different for different SAM output!!! (tophat will have this tag to indicate junction strand)
	my $junc_strand=$1;
	if($junc_strand eq "-"){
		($juncpos1s,$juncpos2s)=($juncpos2s,$juncpos1s);
		($exonreglen1s,$exonreglen2s)=($exonreglen2s,$exonreglen1s);
	}
	return ($juncpos1s,$juncpos2s,$junc_strand,$exonreglen1s,$exonreglen2s);
}

sub modify_cigar_rm_I_D{ #modify cigar string: remove insertion and deletion in cigar string
	my $in_cigar=shift;
	while($in_cigar=~/(\d+)M(\d+)I(\d+)M/){
		my $ori_Str=$&;
		my $convertTo_str=($1+$3)."M";
		$in_cigar=~ s/$ori_Str/$convertTo_str/;
	}
	while($in_cigar=~/(\d+)M(\d+)D(\d+)M/){
		my $ori_Str=$&;
		my $convertTo_str=($1+$3+$2)."M";
		$in_cigar=~ s/$ori_Str/$convertTo_str/;
	}
	$in_cigar;
}

sub from_tag2mismatch{
	my $tags_str=shift;
	my $mismatch_num="";
	if($tags_str=~/NM:i:(\S+)/i){ # NM:i:0
		$mismatch_num=$1;
	}
	return $mismatch_num;
}

1;
