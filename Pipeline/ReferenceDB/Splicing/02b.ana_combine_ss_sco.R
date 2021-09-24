
##combine exon splice site, neighbor ss, sequences and scores, etc

source("../../SharedCodes/Rfunc.inc.R")
gene_model="refflat";id_header="refseqid"; geno="hg19";  #hg19 mm9 rn6
gene_model="ensGene";id_header="ensid"; geno="canFam3";  #rheMac2 canFam3

ss5_f=paste("02.ss_score/",geno,".",gene_model,".ss5.fa.out",sep="")
ss3_f=paste("02.ss_score/",geno,".",gene_model,".ss3.fa.out",sep="")
ss5_seq_f2=paste("SS_flank_seq/",geno,".",gene_model,"/ss5_PM10.fa.tbl",sep="")
trans_ano_f=paste("../gene/02transcript_gene_ano/",geno,".",gene_model,".desc.nonRedundant.txt",sep="") 
geneBlock="exon" #exon, intron
IntronExon_ano_f=paste("10splicingSite/",geno,".",gene_model,".",geneBlock,".ano.txt",sep="")


#####09/11/2015 PTC
ss5_d=read.table(ss5_f, header=F,sep="\t")
ss3_d=read.table(ss3_f, header=F,sep="\t")
ss5_seq_d2=read.table(ss5_seq_f2, header=F,sep="\t")
names(ss5_d)=c("seq","sco")
names(ss3_d)=c("seq","sco")
names(ss5_seq_d2)=c("id","seq")
trans_ano_d=read.table(trans_ano_f, header=T, sep="\t", quote="", stringsAsFactors=F, comment.char="")

ss_out_redun_f=paste("02.ss_score/",geno,".",gene_model,".redundant.",geneBlock,".ss.sco.txt",sep="")
ss_out_f=paste("02.ss_score/",geno,".",gene_model,".",geneBlock,".ss.sco.txt",sep="")

IntronExon_ano_d=read.table(IntronExon_ano_f, header=T, sep="\t", quote="", stringsAsFactors=F, comment.char="")
IntronExon_ano_d$ss5_sco=ss5_d$sco[match(IntronExon_ano_d$ss5_seq, ss5_d$seq)]; table(is.na(IntronExon_ano_d$ss5_sco))
IntronExon_ano_d$ss3_sco=ss3_d$sco[match(IntronExon_ano_d$ss3_seq, ss3_d$seq)]; table(is.na(IntronExon_ano_d$ss3_sco))
write.table(IntronExon_ano_d, file=ss_out_redun_f, sep="\t", col.names=T, row.names=F, quote=F, na="")

##remove redundancy
IntronExon_ano_d$IE_ID=paste(IntronExon_ano_d$contig,IntronExon_ano_d$strand,IntronExon_ano_d$ss5,IntronExon_ano_d$ss3,sep=":")

if_representative_transcripts=IntronExon_ano_d[,id_header] %in% trans_ano_d[,id_header]
IntronExon_ano_d1=IntronExon_ano_d[if_representative_transcripts, ]
IntronExon_ano_d2=IntronExon_ano_d[!if_representative_transcripts, ]
IntronExon_ano_d1=IntronExon_ano_d1[!duplicated(IntronExon_ano_d1$IE_ID),]
IntronExon_ano_d2=IntronExon_ano_d2[!(IntronExon_ano_d2$IE_ID %in% IntronExon_ano_d1$IE_ID), ]
IntronExon_ano_d2=IntronExon_ano_d2[!duplicated(IntronExon_ano_d2$IE_ID),]
IntronExon_ano_d3=rbind(IntronExon_ano_d1,IntronExon_ano_d2)
setdiff(IntronExon_ano_d$IE_ID, IntronExon_ano_d3$IE_ID) #should have no extra
table(duplicated(IntronExon_ano_d3$IE_ID)) #should have no duplicated
IntronExon_ano_d3=IntronExon_ano_d3[order(IntronExon_ano_d3$gene_symbol),]
nrow(IntronExon_ano_d);
nrow(IntronExon_ano_d3);
IntronExon_ano_d3$ss5_id=paste(IntronExon_ano_d3$contig,IntronExon_ano_d3$strand,IntronExon_ano_d3$ss5,sep=":"); table(duplicated(IntronExon_ano_d3$ss5_id))
IntronExon_ano_d3$ss3_id=paste(IntronExon_ano_d3$contig,IntronExon_ano_d3$strand,IntronExon_ano_d3$ss3,sep=":"); table(duplicated(IntronExon_ano_d3$ss3_id))

#change 5'ss sequence to -4 ~ +6
IntronExon_ano_d3$ss5_seq=ss5_seq_d2$seq[match(IntronExon_ano_d3$ss5_id, ss5_seq_d2$id)]
IntronExon_ano_d3$ss5_seq=substr(IntronExon_ano_d3$ss5_seq, 7,16)
#change intron sequence to lowercase
substr(IntronExon_ano_d3$ss5_seq, 5,10) = tolower(substr(IntronExon_ano_d3$ss5_seq, 5,10))
substr(IntronExon_ano_d3$ss3_seq, 1,20) = tolower(substr(IntronExon_ano_d3$ss3_seq, 1,20))
#add intron or exon size
IntronExon_ano_d3$size=abs(IntronExon_ano_d3$ss5-IntronExon_ano_d3$ss3)+1

#change value to percentile (optional)
change2percentile_l=list(ss5_sco='ss5_id', ss3_sco='ss3_id')
for(var in names(change2percentile_l)){
	id_name=change2percentile_l[[var]]
	all_vals=IntronExon_ano_d3[!duplicated(IntronExon_ano_d3[,id_name]) & !is.na(IntronExon_ano_d3[,var]), var]
	fun=ecdf(all_vals)
	IntronExon_ano_d3[paste(var,"percentile",sep=".")]=round(fun(IntronExon_ano_d3[,var])*100,2)
}
##output
IntronExon_ano_d3=re_format_tb(IntronExon_ano_d3, delete_headers=c("genef_row","IE_ID","ss5_id","ss3_id") )
write.table(IntronExon_ano_d3, file=ss_out_f, sep="\t", col.names=T, row.names=F, quote=F,na="")




##add score percentile info (optional)
ss_f=paste("02.ss_score/",geno,".",gene_model,".ss.sco.tbl",sep="")
ss_d=read.table(ss_f, header=T, sep="\t", stringsAsFactors=F, comment.char="", quote="")

#get the ecdf function for variables:
ecdf_fs=list()

vars=c("ss5sco","ss3sco")
for(var in vars){
	ifsel=!duplicated(paste(ss_d$contig, ss_d$strand, ss_d[,var])) & !is.na(ss_d[,var])
	ss_d2=ss_d[ifsel,]
	ecdf_fs[[var]]=ecdf(ss_d2[,var])
}
#intron and exon size:(separate to first middle and last)
intron_types=c("f","m","l","s")
ss_d$intron_type="m"
ss_d$intron_type[ss_d$introni==1]="f"
ss_d$intron_type[ss_d$intron_revnum==-1]="l"
ss_d$intron_type[ss_d$introni==1 & ss_d$intron_revnum==-1]="s"
table(ss_d$intron_type)[intron_types]
ifsel=!duplicated(paste(ss_d$contig, ss_d$strand, ss_d$ss5, ss_d$ss3))
ss_d2=ss_d[ifsel & ss_d$intron_type %in% c("m","l"), ]
ecdf_fs[["exon_len.m"]]=ecdf(ss_d2$ups_exon_len)
for(intron_type in intron_types){ #intron length
	ss_d2=ss_d[ifsel & ss_d$intron_type==intron_type, ]
	for(var in c("intron_len") ){
		sco_name=paste(var,intron_type,sep=".")
		ecdf_fs[[sco_name]]=ecdf(ss_d2[,var])
	}
}

#change value to quantile:
names(ecdf_fs)
val_names=c("ss5sco","ss3sco","ups_ss3sco","dns_ss5sco","intron_len","ups_exon_len","dns_exon_len")
val_ecdf.fun_names=c("ss5sco","ss3sco","ss3sco","ss5sco","intron_len.f","exon_len.m","exon_len.m")
val_qnames=paste(val_names,  ifelse(val_names==val_ecdf.fun_names, "",val_ecdf.fun_names), "q",sep="."); val_qnames
for(i in 1:length(val_names) ){
	ss_d[,val_qnames[i]]=ecdf_fs[[val_ecdf.fun_names[i]]] ( ss_d[,val_names[i]] )
}
##output file:
out_ss_f=paste(ss_f,".percentile.tbl",sep="")
write.table(ss_d, file=out_ss_f, sep="\t", col.names=T, row.names=F, quote=F)
