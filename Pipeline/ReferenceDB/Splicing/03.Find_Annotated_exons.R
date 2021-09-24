#extract GAgt type 5'ss exons. build a database of junctions supporting its inclusion and exclusions.
#coordinates of 5'ss or 3'ss are based on 1-based exonic coordinates
#2018/3/19: can also output a file of all annotated exons based on Refseq

source("../../SharedCodes/Rfunc.inc.R") 
geno="mm9"; #mm9 hg19
gene_model="refflat";
transcript_id_name="refseqid"
transcript_ss_sco_f=paste("02.ss_score/",geno,".",gene_model,".redundant.ss.sco.txt",sep="") 


out_subset_f=paste("03.Annotated_exons/",geno,".",gene_model,".GAgt.exons.txt",sep="")
out_AllAnno_f=paste("03.Annotated_exons/",geno,".",gene_model,".Annotated.exons.txt",sep="")
mkdir_if_not_exist(out_subset_f)
ss5_seq_f=paste("SS_flank_seq/",geno,".",gene_model,"/ss5_PM10.fa.tbl",sep="")

transcript_ss_sco_d=read.table(transcript_ss_sco_f, header=T, sep="\t", stringsAsFactors=F)
#add upstream intron position and score
if_non_1st_intron=transcript_ss_sco_d$introni>1; table(if_non_1st_intron)
ss_info_names=c("ss5","ss3","ss5_sco","ss3_sco")
add_ss_info_names=paste("ups",ss_info_names,sep="_")
transcript_ss_sco_d[add_ss_info_names]=NA
target_rows=(1:nrow(transcript_ss_sco_d))[if_non_1st_intron]
src_rows=target_rows-1
transcript_ss_sco_d[target_rows,add_ss_info_names]=transcript_ss_sco_d[src_rows,ss_info_names]

if_GAgt=substr(transcript_ss_sco_d$ss5_seq, 2,5)=="GAGT"; table(if_GAgt)

#transcript_ss_sco_d=transcript_ss_sco_d[if_GAgt,]
transcript_ss_sco_d$junc_id=paste(transcript_ss_sco_d$contig,transcript_ss_sco_d$strand,transcript_ss_sco_d$ups_ss5, transcript_ss_sco_d$ups_ss3,transcript_ss_sco_d$ss5,transcript_ss_sco_d$ss3, sep=":")
transcript_ss_sco_d=transcript_ss_sco_d[!duplicated(transcript_ss_sco_d$junc_id),]; nrow(transcript_ss_sco_d)
transcript_ss_sco_d=transcript_ss_sco_d[transcript_ss_sco_d$introni>1,]; nrow(transcript_ss_sco_d) #exclude first intron

transcript_ss_sco_d$ups_intron=paste(transcript_ss_sco_d$contig, transcript_ss_sco_d$strand, transcript_ss_sco_d$ups_ss5, transcript_ss_sco_d$ups_ss3,sep=":")
transcript_ss_sco_d$dns_intron=paste(transcript_ss_sco_d$contig, transcript_ss_sco_d$strand, transcript_ss_sco_d$ss5, transcript_ss_sco_d$ss3,sep=":")
transcript_ss_sco_d$skip_intron=paste(transcript_ss_sco_d$contig, transcript_ss_sco_d$strand, transcript_ss_sco_d$ups_ss5, transcript_ss_sco_d$ss3,sep=":")
transcript_ss_sco_d$reg_len=abs(transcript_ss_sco_d$ups_ss3- transcript_ss_sco_d$ss5)+1 #exon length

##update 5'ss sequence
ss5_seq_d=read.table(ss5_seq_f, header=F, sep="\t", stringsAsFactors=F)
transcript_ss_sco_d$ss5_seq =ss5_seq_d[match(paste(transcript_ss_sco_d$contig, transcript_ss_sco_d$strand, transcript_ss_sco_d$ss5, sep=":"), ss5_seq_d[,1]), 2]; table(is.na(transcript_ss_sco_d$ss5_seq))

##output
write.table(transcript_ss_sco_d, file=out_AllAnno_f, col.names=T, row.names=F, sep="\t", quote=F)
#write.table(transcript_ss_sco_d[if_GAgt,], file=out_subset_f, col.names=T, row.names=F, sep="\t", quote=F)
