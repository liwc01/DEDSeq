##function to add splice site sequence in the table
add_ss_site2data=function(input_d=in_cbd2, addSeqNames=c("ss5_seq","ss3_seq"), ss5_header="juncpos5", ss3_header="juncpos3", 
	ss5_pos_shift=rep(0,nrow(input_d)), ss3_pos_shift=rep(0,nrow(input_d)), ss5range=c(-9,10), ss3range=c(-49,10),
	ss5seq_f=paste(out_root,"/../SS_flank_seq/ss5_PM100.fa.tbl",sep=""),
	ss3seq_f=paste(out_root,"/../SS_flank_seq/ss3_PM100.fa.tbl",sep="")
	){ 
	
	ss5_seq_l=list(seqf=ss5seq_f, extractRegion=list(toupper=c(100+ss5range[1],100), tolower=c(101,100+ss5range[2]))  )
	ss3_seq_l=list(seqf=ss3seq_f, extractRegion=list(tolower=c(100+ss3range[1],100), toupper=c(101,100+ss3range[2]))  )

	addSSinfo_l=list(
		ss5_seq=list(ssSeq_l=ss5_seq_l, ssHeader=ss5_header, ssPos_shift=ss5_pos_shift), 
		ss3_seq=list(ssSeq_l=ss3_seq_l, ssHeader=ss3_header, ssPos_shift=ss3_pos_shift)
	)

	for(seq_name in addSeqNames){
		ssSeq_l=addSSinfo_l[[seq_name]]$ssSeq_l
		ssHeader=addSSinfo_l[[seq_name]]$ssHeader
		ssPos_shift=addSSinfo_l[[seq_name]]$ssPos_shift
		seqf=ssSeq_l$seqf; 
		print(paste("open",seqf))
		extractRegion=ssSeq_l$extractRegion
		seq_d=read.table(seqf, header=F, sep="\t", stringsAsFactors=F, quote="")
		seq_m=NULL
		inputSeqIds=paste(input_d$contig, input_d$strand, format(input_d[,ssHeader]+ssPos_shift,scientific=F,trim=T), sep=":")
		for(caseFunc in names(extractRegion)){
			fromPos=extractRegion[[caseFunc]][1]
			toPos=extractRegion[[caseFunc]][2]
			seq_v=do.call(caseFunc, list(x=substring(seq_d[,2], fromPos, toPos)) )
			seq_m=cbind(seq_m, seq_v)
		}
		seq_d$extractSeq=apply(seq_m,1,paste,collapse="")
		input_d[,seq_name]=seq_d$extractSeq[match(inputSeqIds, seq_d[,1])]
	}
	input_d
}

add_ss_supp2data=function(input_d=in_cbd2,  ss5_header="juncpos5", ss3_header="juncpos3", addSuppNames=c("ss5_supp","ss3_supp"),
	ss5_pos_shift=rep(0,nrow(input_d)), ss3_pos_shift=rep(0,nrow(input_d)), 
	anno_version="1511", geno="hg19", contig_header="contig", strand_header="strand",
	gene_models_l=list(mRNA_EST=c("mrna","est"), KnG_Ens=c("knownGene","ensGene"), Refseq="refFlat" ),
	ano_folder="../ReferenceDB/Splicing/04.ss_supp/"
	){ 
	
	ss5_supp_f=paste(ano_folder,geno,"_",anno_version,"/supp.ss5.tbl",sep="")
	ss3_supp_f=paste(ano_folder,geno,"_",anno_version,"/supp.ss3.tbl",sep="")

	addSSinfo_l=list(
		ss5_supp=list(ssSuppf=ss5_supp_f, ssHeader=ss5_header, ssPos_shift=ss5_pos_shift, ssSuppf_ssHeader="ss5"), 
		ss3_supp=list(ssSuppf=ss3_supp_f, ssHeader=ss3_header, ssPos_shift=ss3_pos_shift, ssSuppf_ssHeader="ss3")
	)

	for(supp_name in addSuppNames){
		print(supp_name)
		ssSuppf=addSSinfo_l[[supp_name]]$ssSuppf
		ssHeader=addSSinfo_l[[supp_name]]$ssHeader
		ssPos_shift=addSSinfo_l[[supp_name]]$ssPos_shift
		ssSuppf_ssHeader=addSSinfo_l[[supp_name]]$ssSuppf_ssHeader
		print(paste("open",ssSuppf))
		ssSupp_d=read.table(ssSuppf, header=T, sep="\t", stringsAsFactors=F, quote="")
		ssSupp_d$id=paste(ssSupp_d$contig, ssSupp_d$strand, format(ssSupp_d[,ssSuppf_ssHeader],scientific=F,trim=T), sep=":" )
		ssSupp_d[,supp_name]="None"
		for(supp_cat in names(gene_models_l) ){
			supp_names=paste("supp_",gene_models_l[[supp_cat]], sep="")
			supp_names=intersect(supp_names, names(ssSupp_d))
			if(length(supp_names)>0){
				ssSupp_d[rowSums(ssSupp_d[supp_names],na.rm=T)>0,supp_name] = supp_cat
			}
		}
		table(ssSupp_d[,supp_name])
		inputSeqIds=paste(input_d[,contig_header], input_d[,strand_header], format(input_d[,ssHeader]+ssPos_shift,scientific=F,trim=T), sep=":")
		input_d[,supp_name]=ssSupp_d[match(inputSeqIds, ssSupp_d$id), supp_name]
		input_d[is.na(input_d[,supp_name]),supp_name]="None"
		table(input_d[,supp_name])
	}
	input_d
}



add_ss_ConsSco=function(input_d=in_cbd2,  ss5_header="juncpos5", ss3_header="juncpos3", addConsNames=c("ss5_ConsSco","ss3_ConsSco"),
	ss5_pos_shift=rep(0,nrow(input_d)), ss3_pos_shift=rep(0,nrow(input_d))
	){
	
	ss5_consSco_f=paste("../../club/29SubRegion_consSco/",study_name,"/ss5.m50p100.mammal.sco.tbl",sep="")
	ss3_consSco_f=paste("../../club/29SubRegion_consSco/",study_name,"/ss3.m100p50.mammal.sco.tbl",sep="")

	addSSinfo_l=list(
		ss5_ConsSco=list(ss_consSco_f=ss5_consSco_f, ssHeader=ss5_header, ssPos_shift=ss5_pos_shift, ssConsf_ScoName="Sco_ss5.m50p100.mammal.mp10"), 
		ss3_ConsSco=list(ss_consSco_f=ss3_consSco_f, ssHeader=ss3_header, ssPos_shift=ss3_pos_shift, ssConsf_ScoName="Sco_ss3.m100p50.mammal.mp10")
	)

	for(addConsName in addConsNames){
		print(addConsName)
		ss_consSco_f=addSSinfo_l[[addConsName]]$ss_consSco_f
		ssHeader=addSSinfo_l[[addConsName]]$ssHeader
		ssPos_shift=addSSinfo_l[[addConsName]]$ssPos_shift
		ssConsf_ScoName=addSSinfo_l[[addConsName]]$ssConsf_ScoName
		print(paste("open",ss_consSco_f))
		ss_consSco_d=read.table(ss_consSco_f, header=T, sep="\t", stringsAsFactors=F, quote="")
		inputSeqIds=paste(input_d$contig, input_d$strand, format(input_d[,ssHeader]+ssPos_shift,scientific=F,trim=T), sep=":")
		input_d[,addConsName]=ss_consSco_d[match(inputSeqIds, ss_consSco_d$ref_id), ssConsf_ScoName]
	}
	input_d
}

add_ss_liftover=function(input_d=in_cbd2){
	print(paste("read",second_spl_f) )
	second_spl_d=read.table(second_spl_f, header=T, sep="\t", quote="", comment.char="", stringsAsFactors=F)
	ss5liftover_d=read.table(ss5liftover_f, header=F, sep="\t", quote="", comment.char="", stringsAsFactors=F)
	ss3liftover_d=read.table(ss3liftover_f, header=F, sep="\t", quote="", comment.char="", stringsAsFactors=F)
	names(ss5liftover_d)[4]="from_id"
	names(ss3liftover_d)[4]="from_id"
	ss5liftover_d$to_id=paste(ss5liftover_d[,1],ss5liftover_d[,6], ss5liftover_d[,3], sep=":")
	ss3liftover_d$to_id=paste(ss3liftover_d[,1],ss3liftover_d[,6], ss3liftover_d[,3], sep=":")
	strand_signs=ifelse(input_d$strand=="-",-1,1)
	if_start_ss5=!is.na(input_d$startPosType) & input_d$startPosType=="ss5"
	if_start_ss3=!is.na(input_d$startPosType) & input_d$startPosType=="ss3"
	if_end_ss5=!is.na(input_d$endPosType) & input_d$endPosType=="ss5"
	if_end_ss3=!is.na(input_d$endPosType) & input_d$endPosType=="ss3"
	strand_signs_2=ifelse(second_spl_d$strand=="-",-1,1)
	if_start_ss5_2=!is.na(second_spl_d$startPosType) & second_spl_d$startPosType=="ss5"
	if_start_ss3_2=!is.na(second_spl_d$startPosType) & second_spl_d$startPosType=="ss3"
	if_end_ss5_2=!is.na(second_spl_d$endPosType) & second_spl_d$endPosType=="ss5"
	if_end_ss3_2=!is.na(second_spl_d$endPosType) & second_spl_d$endPosType=="ss3"
	update_names=paste(second_spl_name,second_spl_update_names,sep=".")
	input_d=input_d[,!grepl(second_spl_name,names(input_d))]
	input_d[,c(paste(second_spl_name, c(".start_pos",".end_pos"),sep=""),update_names)]=NA

	end_ss5_ids_from=paste(input_d$contig, input_d$strand, input_d$end_pos, sep=":")[if_end_ss5]
	end_ss5_ids_to=ss5liftover_d$to_id[match(end_ss5_ids_from, ss5liftover_d$from_id)]
	tmp_d2=second_spl_d[if_end_ss5_2,]
	end_ss5_ids2=paste(tmp_d2$contig, tmp_d2$strand, tmp_d2$end_pos, sep=":")
	if_need_update=rowSums(!is.na(input_d[,update_names]) )==0; sum(if_need_update)
	input_d[if_end_ss5,paste(second_spl_name, ".end_pos",sep="")]=end_ss5_ids_to
	new_tb=tmp_d2[match(end_ss5_ids_to, end_ss5_ids2), second_spl_update_names]
	if(any(if_need_update[if_end_ss5])) {
		input_d[if_end_ss5 & if_need_update, update_names]= new_tb[  if_need_update[if_end_ss5], ]
	}
	
	start_ss3_ids_from=paste(input_d$contig, input_d$strand, input_d$start_pos, sep=":")[if_start_ss3]
	start_ss3_ids_to=ss3liftover_d$to_id[match(start_ss3_ids_from, ss3liftover_d$from_id)]
	tmp_d2=second_spl_d[if_start_ss3_2,]
	start_ss3_ids2=paste(tmp_d2$contig, tmp_d2$strand, tmp_d2$start_pos, sep=":")
	if_need_update=rowSums(!is.na(input_d[,update_names]) )==0; sum(if_need_update)
	input_d[if_start_ss3,paste(second_spl_name, ".start_pos",sep="")]=start_ss3_ids_to
	new_tb=tmp_d2[match(start_ss3_ids_to, start_ss3_ids2), second_spl_update_names]
	if(any(if_need_update[if_start_ss3])) {
		input_d[if_start_ss3 & if_need_update, update_names]= new_tb[  if_need_update[if_start_ss3], ]
	}

	start_ss5_ids_from=paste(input_d$contig, input_d$strand, input_d$start_pos-strand_signs, sep=":")[if_start_ss5]
	start_ss5_ids_to=ss5liftover_d$to_id[match(start_ss5_ids_from, ss5liftover_d$from_id)]
	tmp_d2=second_spl_d[if_start_ss5_2,]
	start_ss5_ids2=paste(tmp_d2$contig, tmp_d2$strand, tmp_d2$start_pos-strand_signs_2[if_start_ss5_2], sep=":")
	if_need_update=rowSums(!is.na(input_d[,update_names]) )==0; sum(if_need_update)
	input_d[if_start_ss5,paste(second_spl_name, ".start_pos",sep="")]=start_ss5_ids_to
	new_tb=tmp_d2[match(start_ss5_ids_to, start_ss5_ids2), second_spl_update_names]
	if(any(if_need_update[if_start_ss5])) {
		input_d[if_start_ss5 & if_need_update, update_names]= new_tb[  if_need_update[if_start_ss5], ]
	}

	end_ss3_ids_from=paste(input_d$contig, input_d$strand, input_d$end_pos+strand_signs, sep=":")[if_end_ss3]
	end_ss3_ids_to=ss3liftover_d$to_id[match(end_ss3_ids_from, ss3liftover_d$from_id)]
	tmp_d2=second_spl_d[if_end_ss3_2,]
	end_ss3_ids2=paste(tmp_d2$contig, tmp_d2$strand, tmp_d2$end_pos+strand_signs_2[if_end_ss3_2], sep=":")
	if_need_update=rowSums(!is.na(input_d[,update_names]) )==0; sum(if_need_update)
	input_d[if_end_ss3,paste(second_spl_name, ".end_pos",sep="")]=end_ss3_ids_to
	new_tb=tmp_d2[match(end_ss3_ids_to, end_ss3_ids2), second_spl_update_names]
	if(any(if_need_update[if_end_ss3])) {
		input_d[if_end_ss3 & if_need_update, update_names]= new_tb[  if_need_update[if_end_ss3], ]
	}
	
	rm("second_spl_d","tmp_d2","ss5liftover_d","ss3liftover_d","new_tb")
	input_d
}


cal_regu_type<-function(data, p_name, change_name, p_cut, change_cut, ReguTypes=c("UP","DN","NC","na") ){
	reguType_v=rep(ReguTypes[4],nrow(data))
	if_valued=!is.na(data[,p_name]) & !is.na(data[,change_name])
	if_sig= data[,p_name]<p_cut
	if_up=if_valued & if_sig  & data[,change_name]>change_cut
	if_dn=if_valued & if_sig  & data[,change_name]< (-change_cut)
	reguType_v[if_valued]=ReguTypes[3]
	reguType_v[if_up]=ReguTypes[1]
	reguType_v[if_dn]=ReguTypes[2]
	reguType_v
}

