
#combine PE read stats
args<- commandArgs(TRUE)
print(args)

args_v<-NA
if(length(args)>1){
	args_v <- args[c(FALSE, TRUE)]
	names(args_v) <- sub("\\-+","", args[c(TRUE, FALSE)])
	print(args_v)
}
in_dir="/HPCTMP_NOBKUP/wl314/analyze/Projects/TJ/Feng_p14p16_RD_RNAseq/02.PE.ReadTable"
in_dir<- ifelse(is.na(args_v["in_dir"]), "", args_v["in_dir"])
out_f=ifelse(is.na(args_v["out_f"]), paste(in_dir,"/PEread.stats.txt",sep=""), args_v["out_f"])
inf_pattern=ifelse(is.na(args_v["inf_pattern"]), ".run.log", args_v["inf_pattern"])
infs=list.files(path=in_dir, pattern=inf_pattern); print(infs)

save.image(paste0(out_f,".R.image"))


dl=list()
all_headers=NULL
samples=NULL
for(i in 1:length(infs)){
	inf=paste(in_dir,infs[i],sep="/")
	print(inf)
	sample=sub(inf_pattern,"", infs[i]); sample
	d=read.table(inf, comment.char="#", stringsAsFactors=F, quote="", sep="\t", header=T, check.names=F)
	all_headers=unique(c(all_headers,names(d)))
	dl[[i]]=d
	samples=c(samples,sample)
}
samples

comb_d=data.frame(Sample=samples)
add_headers=setdiff(all_headers,c("Sample",""))
comb_d[add_headers]=0
for(i in 1:length(samples)){
	d=dl[[i]]
	update_headers=intersect(add_headers,names(d))
	comb_d[i,update_headers]=d[update_headers]
}
comb_d

save.image(paste0(out_f,".R.image"))
write.table(comb_d, file=out_f, col.names=T, row.names=F, sep="\t",quote=F)
