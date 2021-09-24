#source(paste(root_dir,"analyze/Rfunc.inc.R",sep=""))

geno_len <- list(hg18=3107677273, hg19=3137161264, mm8=2664455088, mm9=2725765481, rn4=2834127293)


tax_ids        <- c(9606,9606,9606,9606,9606,                  10090,10090,10090,10090,10090
	,10116,10116,10116,10116,10116,               9031,9031,9031,9031
	,9544,9544,9544,9544
	,9615,9615,9615)
names(tax_ids) <- c("hs","human","Homo_sapiens","hg18","hg19", "mm","mouse","Mus_musculus","mm8","mm9"
	,"rn","rat","Rattus_norvegicus","rn4","rn6", "galGal","galGal3","chicken","Gallus_gallus"
	,"monkey","rheMac2","rheMac10","Macaca_mulatta"
	,"dog","canFam3","Canis_lupus_familiaris")
spe_common_names <- c("human","mouse","rat","chicken","monkey","dog")
names(spe_common_names) <- tax_ids[spe_common_names]

spe_latin_names <- c("Homo_sapiens","Mus_musculus","Rattus_norvegicus","Macaca_mulatta","Canis_lupus_familiaris")
names(spe_latin_names) <- tax_ids[spe_latin_names]

loadSpe2TaxIDs<-function(inf="/Home/cocochen/data/ncbi/taxonomy/ucsc_genos.taxid.txt"){
	tax_d=read.table(inf,header=T,sep='\t')
	tax_ids=tax_d$taxid
	names(tax_ids)=as.character(tax_d$name)
	tax_ids
}

mkdir_if_not_exist<-function(out_root){
	out_root_dir=sub("[^/]*$","",out_root)
	system(paste("mkdir -p \'",out_root_dir,"\'", sep=""))
}

string2list<-function(s){
	l<-list()
	for(v1 in unlist(strsplit(s,";"))){ #eg: "REGU:UP DN;XXXX:XX XXX
		v2<-unlist(strsplit(v1, ":")) #PA:P0 P1 P2 P3
		l[[v2[1]]] <- unlist(strsplit(v2[2], " "))
	}
	return(l)
}


load_pAsco_data <- function(species="human"){
	if(species=="human"){
		pA_score_file <- paste(root_dir,"analyze/polyA/0_club_like/7pA_svm_res/Hs.site5.pASVM.tbl",sep="")
	}else if(species=="mouse"){
		pA_score_file <- paste(root_dir,"analyze/polyA/0_club_like/7pA_svm_res/Mm.site5.pASVM.tbl",sep="")
	}
	pA_score_data <- read.table(pA_score_file,header=T,sep="\t")
	cis_names <<- unlist(strsplit("AUE1,AUE2,AUE3,AUE4,CUE1,CUE2,CDE1,CDE2,CDE3,CDE4,ADE1,ADE2,ADE3,ADE4,ADE5", ","))
	counts["num.row.pAscore.all"] <<- nrow(pA_score_data)
	pA_score_data
}

change_values<-function(inv, from_patt, to_txt){ #change values in a vector
	out<-as.character(unlist(inv))
	for(i in 1:length(from_patt)){
		out<-gsub(from_patt[i], to_txt[i], out)
	}
	return(out)
}


load_gene_NTcontent <- function(spe="human", mode="refflat"){
	#spe= human mouse
	#mode=refflat club
	geno <- c("hg18","mm8"); names(geno) <- c("human","mouse")
	gene_ntnum_file <- paste(root_dir,"analyze/club/5gene_NTcontent/",geno[spe],".",mode,".gene.NTnum.ano",sep="")
	if(mode=="refflat"){
		gene_ano_f <- paste(root_dir,"analyze/club/3add_gene_ano/",geno[spe],".",mode,".desc",sep="")
	}else if(mode=="club"){
		gene_ano_f <- paste(root_dir,"analyze/club/0remove_club_redundance/4gene_exon_ano/",geno[spe],"_",mode,"_DOEuni.splice.gene",sep="")
	}
	
	gene_ntnum_data <- read.table(gene_ntnum_file,header=T,sep="\t")
	all_regnames <<- sub("A.","", grep('A.',names(gene_ntnum_data),value=T))
	for(regname1 in all_regnames){
		gene_ntnum_data[paste("pGC",regname1,sep=".")] <- rowSums(gene_ntnum_data[,paste(c("C","G"),regname1,sep=".")])/rowSums(gene_ntnum_data[paste(c("A","T","C","G"),regname1,sep=".")])
		gene_ntnum_data[paste("pG",regname1,sep=".")] <- (gene_ntnum_data[,paste(c("G"),regname1,sep=".")])/rowSums(gene_ntnum_data[paste(c("A","T","C","G"),regname1,sep=".")])
		gene_ntnum_data[paste("pC",regname1,sep=".")] <- (gene_ntnum_data[,paste(c("C"),regname1,sep=".")])/rowSums(gene_ntnum_data[paste(c("A","T","C","G"),regname1,sep=".")])
	}
	gene_ano_D <- read.table(gene_ano_f, header=T, sep="\t", quote="", comment.char="")
	cbind(gene_ano_D[setdiff(names(gene_ano_D),c("exon_starts","exon_ends"))], gene_ntnum_data)
}


sep2_SAPA <- function(inids,idtype="gene_symbol",species="human",pA_ano_file=""){
	#idtype = gene_symbol, gene_id, club_id, initiator
	#species= "mouse", human
	spe_latin <- spe_latin_names[as.character(tax_ids[species])]
	if (pA_ano_file==""){
		if(species=="human"){
			pA_ano_file <- paste(root_dir,"analyze/polyA/0_club_like/5b_pA_unique_ano/hg18.all.pA.unique.ano",sep="")
		}else if(species=="mouse"){
			pA_ano_file <- paste(root_dir,"analyze/polyA/0_club_like/5b_pA_unique_ano/mm8.all.pA.unique.ano",sep="")
		}		
	}

	pA_ano_data <- read.table(pA_ano_file,header=T,sep="\t")
	pA_ano_data["club_id"] <- sub("\\.[0-9]*$","", pA_ano_data$pAid)
	if(idtype=="gene_id"){
		pA_ano_data["gene_id"] <- geneidmapping(pA_ano_data$gene_symbol, inspe=spe_latin, outspe=spe_latin)
	}
	outtype <- rep("unkn", length(inids))
	outtype[inids %in% pA_ano_data[pA_ano_data$pA_type %in% "3S",idtype]] <- "SPA"
	outtype[inids %in% pA_ano_data[!(pA_ano_data$pA_type %in% "3S"),idtype]] <- "APA"
	outtype[is.na(inids)] <- NA
	print(table(outtype))
	outtype
}



gid2_gene_desc <- function(gids, spe="human", gene_info_file=NA){
	geneinfo_data=open_geneInfo_data(spe, gene_info_file=gene_info_file)
	outtb <- geneinfo_data[match(gids, geneinfo_data$gene_id), c("gene_symbol","gene_desc")]
}
gsb2_gene_desc <- function(gsbs, spe="human", gene_info_file=NA){
	geneinfo_data=open_geneInfo_data(spe, gene_info_file=gene_info_file)
	outtb <- geneinfo_data[match(gsbs, geneinfo_data$alias), c("gene_id","gene_symbol","gene_desc")]
}
open_geneInfo_data<-function(spe="human", gene_info_file=NA){
	tax_id <- tax_ids[spe]
	spe_latin_name=spe_latin_names[as.character(tax_id)]
	if(is.na(gene_info_file)){ gene_info_file <- paste("01gene_alias2id/",spe_common_names[as.character(tax_id)],".",tax_id,".tbl",sep="") }
	if(! file.exists(gene_info_file)){
		gene_info_file <- paste("../ncbi/gene/geneinfo/All_Mammalia.gene_info",sep="")
		print(gene_info_file)
		geneinfo_data <- read.table(file=gene_info_file,header=F,quote="",comment.char="", sep="\t", skip=1,stringsAsFactors=F)
		names(geneinfo_data)[c(1,2,3,5,9)] <- c("tax_id","gene_id","gene_symbol","alias","gene_desc")
		geneinfo_data <- geneinfo_data[,c(1,2,3,5,9)]
		if(!is.na(tax_id)){
			geneinfo_data <- geneinfo_data[geneinfo_data$tax_id == tax_id,]
		}
		geneinfo_data[nrow(geneinfo_data)+1,] <- rep("-",ncol(geneinfo_data)) 
	}else{
		print(gene_info_file)
		geneinfo_data <- read.table(file=gene_info_file,header=T,quote="",comment.char="", sep="\t", stringsAsFactors=F)
	}
}

#id mapping
geneidmapping <- function(invec, intype="alias", inspe="Mus_musculus", outtype="gene_id", outspe="Homo_sapiens", outspetype="gid", root_dir="/drive2/wli/",
	gsb2gid_file = paste(root_dir,"analyze/Pipeline/ReferenceDB/gene/01gene_alias2id/",spe_common_names[as.character(tax_ids[inspe])],".",tax_ids[inspe],".tbl",sep=""),
	homologene_file =paste(root_dir,"analyze/club/8HomoloGene_tb/HomoloGene.build68.tbl",sep="")
	){
	#type should be alias, gene_symbol or gene_id (header in the gsb2gid_file)
	#spe should be Homo_sapiens Mus_musculus ... (header in the homgene.64.tbl hid     gid_Mus_musculus        gsb_Mus_musculus        gid_Rattus_norvegicus   gsb_Rattus_norvegicus   gid_Magnaporthe_grisea...)
	
	if(intype==outtype){
		outvec1 <- invec
	}else{
		gsb2gid_data <- read.table(gsb2gid_file, header=T,sep="\t",quote="",comment.char="",stringsAsFactors=F)
		outvec1 <- gsb2gid_data[[outtype]][match(toupper(invec), toupper(gsb2gid_data[[intype]]) )]
		outvec1[is.na(invec)] <- NA
		print(table(is.na(invec)))
		print(table(is.na(outvec1)))
	}
	if(inspe==outspe){
		return(outvec1)
	}else{ #tranfer species
		print(homologene_file)
		homologene_data <- read.table(homologene_file, header=T,sep="\t",quote="",comment.char="")
		s_inout_type <- sub("geneid|gene_id","gid",c(outtype,outspetype))
		s_inout_type <- sub("alias|genesymbol|gene_symbol","gsb",s_inout_type)
		print(s_inout_type)
		inheader <- paste(s_inout_type[1],inspe,sep="_")
		outheader <- paste(s_inout_type[2],outspe,sep="_")
		outvec2 <- homologene_data[[outheader]][match( toupper(outvec1), toupper(homologene_data[[inheader]]) )]
		print("table(is.na(outvec2))")
		print(table(is.na(outvec2)))
		return(outvec2)
	}
}

##refseq_id to gene id
refseqid2gid <- function(refseq_id, species="human"){
	db_f <- paste("/Home/cocochen/analyze/club/15gene2refseq_tb/taxid.",tax_ids[species],".tbl",sep="")
	db_d <- read.table(db_f, header=T, sep="\t")
	db_d$refseq_id <- sub("\\..*$","",db_d$refseq)
	db_d$gid[match(refseq_id, db_d$refseq_id)]
}
out_cor_param_str<-function(x,y,cor_method="pea"){
	if(cor_method=="both"){
		paste(out_cor_param_str(x,y,"s"), ";",out_cor_param_str(x,y,"p"),sep="")
	}else{
		corres0 <- cor.test(x, y, method=cor_method)
		corres0_str=paste(  cor_method, ":r",round(corres0$estimate,3)," P", format(corres0$p.value, scientific=T, digits=1), sep="")
	}
}
out_cor_param<-function(x,y,cor_method="pea"){
	if(cor_method=="both"){
		c(out_cor_param(x,y,"s"), out_cor_param(x,y,"p") )
	}else{
		corres0 <- cor.test(x, y, method=cor_method)
		c(  round(corres0$estimate,3), corres0$p.value ) # output r and p in a vector
	}
}

draw_scatter_plot<-function(x,y, if_highlight_v=rep(T,length(x)), xlab="", ylab="", highlight_col="black", 
	if_highlight_v2=rep(T,length(x)), highlight_col2="red", cor_method="pea", title="", 
	verline=0, hline=0, lm_line=2, pch=".",pch2=19, cex0=1,cex1=1,cex2=1, bg_col="gray", vh_line_col="gray", bg_line_col=bg_col,  ... ){
	
	if_plotable<- abs(x)!=Inf & abs(y)!=Inf & !is.na(x) & !is.na(y); print(sum(if_plotable))
	print(table(!is.na(if_plotable)))
	print(sum(if_highlight_v))
	print(sum(if_highlight_v2))
	if(any(!if_plotable)){
		x=x[if_plotable]
		y=y[if_plotable]
		if_highlight_v=if_highlight_v[if_plotable]
		if_highlight_v2=if_highlight_v2[if_plotable]
	}
	if(!is.na(cor_method)){
		corres0_str=out_cor_param_str(x,y,cor_method=cor_method)
		if(sum(if_highlight_v)>4 & sum(if_highlight_v)<length(if_highlight_v) ){
			corres_str=out_cor_param_str(x[if_highlight_v], y[if_highlight_v],cor_method=cor_method)
		}else{
			corres_str=NULL
		}
		if(sum(if_highlight_v2)>4 & sum(if_highlight_v2)<length(if_highlight_v2) ){
			corres2_str=out_cor_param_str(x[if_highlight_v2], y[if_highlight_v2],cor_method=cor_method)
		}else{
			corres2_str=NULL
		}
	}else{
		corres0_str=NULL; corres_str=NULL; corres2_str=NULL; 
	}
	title <- paste(title,  paste("#all",bg_col, length(x)), corres0_str, 
		ifelse(any(if_highlight_v) & !all(if_highlight_v), paste(highlight_col, sum(if_highlight_v)), ''),   corres_str, "\n",
		ifelse(any(if_highlight_v2) & !all(if_highlight_v2), paste(highlight_col2, sum(if_highlight_v2)), ''),  corres2_str,  sep=" ")
	print(title)
	plot(x,y, pch=pch, col=bg_col, main= title, xlab=xlab, ylab=ylab,cex=cex0, ...)
	if( sum(if_highlight_v)>0 & sum(if_highlight_v)<length(if_highlight_v) ){
		points(x[if_highlight_v], y[if_highlight_v], pch=pch2, col=highlight_col,cex=cex1)
	}
	if( sum(if_highlight_v2)>0 & sum(if_highlight_v2)<length(if_highlight_v2) ){
		points(x[ if_highlight_v2] , y[ if_highlight_v2], pch=pch2, col=highlight_col2,cex=cex2)
	}
	if(0 %in% lm_line & length(x)>4){
		abline( lm(y~x), lwd=2, col=bg_line_col)
	}
	if(1 %in% lm_line & sum(if_highlight_v, na.rm=T)>4){
		abline( lm(y[if_highlight_v]~x[if_highlight_v]), lwd=2, col=highlight_col)
	}
	if(2 %in% lm_line & sum(if_highlight_v2, na.rm=T)>4){
		abline( lm(y[ if_highlight_v2]~x[ if_highlight_v2]), lwd=2, col=highlight_col2)
	}
	if(any(!is.na(verline)) ){ abline(v=verline, col=vh_line_col, lty=2, lwd=2)	}
	if(any(!is.na(hline)) ){ abline(h=hline, col=vh_line_col, lty=2, lwd=2)	}
}




draw_grped_boxplot <- function(datasheet,grp_varname,grp_values,plot_var_names,plot_var_snames=plot_var_names, title="", xlabel="", ylabel="", 
	hline_at=NULL, ifoutline=T, mycolor=rainbow(length(grp_values)),  if_legend=F,
	ymax=quantile(datasheet[plot_var_names][abs(datasheet[plot_var_names])!=Inf],na.rm=T, probs=ifelse(ifoutline,1, 0.99)),
	ymin=quantile(datasheet[plot_var_names][abs(datasheet[plot_var_names])!=Inf],na.rm=T, probs=ifelse(ifoutline,0, 0.01)), 
	vlines=F, colorBy_vars=F, join_by_grps=F,
	test_sample_pairs="first_last", test_method="ks.test", spaceBwGrp=0.5, las=3, ...	 ){
		
	if(colorBy_vars){
		names(mycolor) <- as.character(plot_var_names)
	}else{
		names(mycolor) <- as.character(grp_values)
	}
	
	test1reg<-function(regname1){
		comp_grps_pval(data=datasheet, x=datasheet[,regname1], grpname=grp_varname,grpvals=grp_values,
			test_sample_pairs=test_sample_pairs, test_method=test_method)
	}
	pvals <- sapply(plot_var_names,test1reg)
	title <- paste(c(title, paste(grp_values,"=", table(datasheet[[grp_varname]])[grp_values], sep=""), 
		"\n",test_method,".P=", unlist(pvals)), collapse=" " )
	
	if(join_by_grps){
		xmax=(spaceBwGrp+length(grp_values))*length(plot_var_names)+1-spaceBwGrp*2
		group1s_x=1 +(1:length(plot_var_names)-1)*(spaceBwGrp+length(grp_values))
	}else{ #defaulf color by group
		xmax=(spaceBwGrp+length(plot_var_names))*length(grp_values)+1-spaceBwGrp*2
		group1s_x=1 + (1:length(grp_values)-1)*(spaceBwGrp+length(plot_var_names))
	}
	print(xmax)
	print(group1s_x)
	plot(0, type='n', xlim=c(spaceBwGrp,xmax), ylim=c(ymin,ymax),xlab=xlabel,ylab=ylabel, xaxt='n', main=title, ...)
	for (i in 1:length(grp_values)){
		grp_value=grp_values[i]
		if(join_by_grps){
			at_coors=group1s_x+i-1
		}else{
			at_coors=group1s_x[i]+1:length(plot_var_names) -1
		}
		
		if(colorBy_vars){
			col=mycolor
		}else{
			col=mycolor[as.character(grp_value)]
		}
		boxplot(datasheet[datasheet[[grp_varname]]==grp_value,plot_var_names], add=T, at=at_coors, 
			col=col, names=plot_var_snames, ylab="", xlab="",xaxt='n',  outline=ifoutline, las=las, lwd=1, ...)
	}

	if(join_by_grps){
		axis(1, at=group1s_x, labels=plot_var_snames, tick=F, las=las)
	}else{
		axis(1, at=group1s_x, labels=grp_values, tick=F, las=las)
	}
	
	if(!is.null(hline_at)){
		abline(h=hline_at,col="gray",lty=1,lwd=1)
	}
	if(vlines){
		abline(v=group1s_x-1 ,col="gray",lty=1, lwd=1)
	}
	if(if_legend){
		plot(0, xlim=0:1,ylim=0:1, type='n')
		if(colorBy_vars){
			legend(x=0, y=1, legend=plot_var_snames,fil=mycolor)
		}else{
			legend(x=0, y=1, legend=grp_values,fil=mycolor)
		}
	}
}





draw_2grped_boxplot <- function(datasheet,grp_varname,grp_values,grp2name,grp2vals,plot_var_name,plot_var_snames, 
	title="", xlabel=grp2name, ylabel=plot_var_name, hline_at=NULL, ifoutline=T, mycolor=rainbow(length(grp_values)), 
	iflegend=F, boxwex=0.5, 
	ymax=quantile(datasheet[,plot_var_name][abs(datasheet[,plot_var_name])!=Inf],probs=0.98,na.rm=T),
	ymin=quantile(datasheet[,plot_var_name][abs(datasheet[,plot_var_name])!=Inf],probs=0.02,na.rm=T), ylogbase=NULL,
	method="boxplot" ){
	#in boxplot: grp2vals is the bigger group in x axis, grp_varname is small group used in legend
	if(!is.null(ylogbase)){
		datasheet[[plot_var_name]] <- log(datasheet[[plot_var_name]], base=ylogbase)
		ymin <- log(ymin, base=ylogbase)
		ymax <- log(ymax, base=ylogbase)
		title <- paste("log",ylogbase," ",title,sep="" )
	}
	test1reg<-function(grp2val){
		return_ks_pval( datasheet[datasheet[[grp_varname]]==grp_values[1] & datasheet[[grp2name]]==grp2val,plot_var_name], 
			datasheet[datasheet[[grp_varname]]==grp_values[length(grp_values)] & datasheet[[grp2name]]==grp2val,plot_var_name] )
	}
	return_ks_pval<-function(vec1,vec2){
		if(length(vec1[!is.na(vec1)])<3 |  length(vec2[!is.na(vec2)])<3){return(NA)}
		format(ks.test(vec1,vec2)$p.value, scientific=T, digits=1)
	}
	kspvalues <- sapply(grp2vals,test1reg)
	datasheet <- datasheet[!is.na(datasheet[,plot_var_name]),]
	numtb <- table(datasheet[[grp_varname]],datasheet[[grp2name]])[grp_values,grp2vals]
	#print(numtb)
	title=paste(c(title, plot_var_name,"~", grp2name, "~", grp_varname, "\n#=", as.numeric(numtb), "\nKS.P=", kspvalues), collapse=" ")
	print (title)
	x_shift <- ( (1:length(grp_values))-mean(1:length(grp_values)) ) / (length(grp_values)-(1+length(grp_values))/2) / 4
	names(x_shift)=grp_values
	names(mycolor) <- grp_values
	
	tmptb <- datasheet[datasheet[[grp_varname]]==grp_values[2] & datasheet[[grp2name]]==grp2vals[1],]
	boxplot(tmptb[[plot_var_name]], xlim=c(0.5,length(grp2vals)+0.5), ylim=c(ymin,ymax),xlab=xlabel,ylab=ylabel,at=1:length(plot_var_name)+x_shift[2], main=title, boxwex =boxwex,col=mycolor[2], names=paste("    ",plot_var_snames), ylab="", las=3, outline=ifoutline)
	for (grp_value in grp_values){
		for(i in 1:length(grp2vals)){
			tmptb <- datasheet[datasheet[[grp_varname]]==grp_value & datasheet[[grp2name]]==grp2vals[i],]
			boxplot(tmptb[[plot_var_name]], add=T, at=i+x_shift[grp_value], boxwex = boxwex, col=mycolor[grp_value], names=rep("",length(plot_var_name)), ylab="", xlab="", outline=ifoutline)
		}
	}
	axis(side=1,label=grp2vals,at=1:length(grp2vals))
	#legend(x=0.5, y=ymax,legend=grp_values,fil=mycolor)		
	text(kspvalues,x=1:length(grp2vals),y=ymax+(ymax-ymin)*0.1,col="red")
	if(!is.null(hline_at)){
		abline(h=hline_at,col="grey",lty=4,lwd=3)
	}
	if(iflegend){
		plot(0,xlim=0:1,ylim=0:1,main=grp_varname)
		legend(0,1, grp_values,fill=mycolor)
	}
}

draw_ecdf_density <- function(data, vars, grpname, grpvals, mycols=rainbow(length(grpvals)), desc="", 
	method="ecdf", xlog=F,logbase=10, iflegend=F, xmin=NULL, xmax=NULL, combine_vars=F, combine_var_name="all_vars",
	test_method="ks.test", test_sample_pairs="first_last", mu=0,  importedPvals=NULL, xlab=NA,  ltys=rep(1,length(grpvals)), ...){
	#1/23/2013 if set grpname=NA, the function will use each var as a group, grpname='var'
	names(mycols) <- grpvals
	names(ltys) <- grpvals
	print(nrow(data))
	if(!is.na(grpname)){
		data<- data[!is.na(data[,grpname]) & data[,grpname] %in% grpvals, ]; print(nrow(data))
	}else{ #each var is a group
		grpname='var'
		grpvals=vars
		combine_vars=T
	}
	if(xlog){
		data[,vars] <- log(data[,vars], base=logbase)
		logstr <- paste("log", logbase,sep="")
	}else{logstr=NULL}
	
	if(combine_vars){
		if(grpname=='var'){
			data <- data.frame( unlist(data[,vars]), rep(vars,each=nrow(data)) )
		}else{
			data <- data.frame( unlist(data[,vars]), rep(data[,grpname],length(vars)) )
		}
		names(data) <- c(combine_var_name,grpname)
		vars <- combine_var_name
	}
	for(var in vars){
		x <- data[[var]]
		x[abs(x)==Inf] <- NA
		xlims= c(ifelse(is.null(xmin),quantile(x,na.rm=T,probs=0.01),xmin), ifelse(is.null(xmax), quantile(x,na.rm=T,probs=0.99), xmax))
		if(method=="ecdf"){
			ylims <- c(0,1)
		}else{
			ylims <- c(0, max (sapply(grpvals, function(v){max( density(x[data[[grpname]]==v],na.rm=T)$y )} )))
		}
		
		if(test_sample_pairs=="first_last"){
			test_grp_matrix=rbind(grpvals[1], grpvals[length(grpvals)])
		}else if(test_sample_pairs=="first_oths"){ #first one vs others
			test_grp_matrix=rbind( rep(grpvals[1], length(grpvals)-1), grpvals[-1] )
		}else if(test_sample_pairs=="adjacent"){
			test_grp_matrix=rbind( grpvals[-length(grpvals)], grpvals[-1] )
		}else if(test_sample_pairs=="All2constant"){ #all value to a constant value (eg. 0)
			test_grp_matrix=rbind( grpvals, "mu" )
		}
		if(!is.null(importedPvals)){
			pval=importedPvals
		}else{
			pval=apply(test_grp_matrix, 2, function(grp_2vals){
				grp1_xs<- x[!is.na(x) & data[[grpname]]==grp_2vals[1]]
				if(! (grp_2vals[2]=="mu")){
					grpn_xs<- x[!is.na(x) & data[[grpname]]==grp_2vals[2]]
					if(length(grp1_xs)>=3 & length(grpn_xs)>=3){
						pval1 <- format(do.call(test_method, list(x=grp1_xs, y=grpn_xs))$p.value, scientific=T, digits=2)
					}else{
						pval1 <- NA
					}
				}else{ ##mu value defined, for All2constant comparison
					if(length(grp1_xs)>=3 ){
						pval1 <- format(do.call(test_method, list(x=grp1_xs, mu=mu))$p.value, scientific=T, digits=2)
					}else{
						pval1 <- NA
					}
				}
				
				pval1
			} )			
		}
		
		title <- paste( c(desc,   grpname,
			"\n",  paste(grpvals,"=", table(data[[grpname]])[grpvals], sep=""), 
			"\n",test_method,"p=",pval 
			), collapse=" ")
		
		print(title)
		plot( 1, xlim=xlims, ylim=ylims, type="n", main=title, xlab=ifelse(is.na(xlab),paste(logstr,var),xlab), ylab=method, ...)
		for(grpval in grpvals){
			if(sum(data[[grpname]]==grpval & !is.na(data[,var]))<5){next}
			if(method=="ecdf"){
				lines(ecdf(data[data[[grpname]]==grpval,var]),  col=mycols[grpval], lty=ltys[grpval], ...)
			}else if(method=="density"){
				lines(density(data[data[[grpname]]==grpval,var], na.rm=T), col=mycols[grpval], lty=ltys[grpval],...)
			}
		}
	}
	
	if(iflegend){
		plot(c(0,1),c(0,1),type="n")
		legend(0,1, legend=grpvals, col=mycols, lwd=2)
	}
}

draw_grped_barplots <- function(data, c_grpname, r_grpname, c_grpvals=NULL, r_grpvals=NULL, perc_byrow=T, row_cols=rainbow(length(r_grpvals)), desc="", iflegend=F, xlim=NULL, ylim=NULL){
	names(row_cols) <- r_grpvals
	data<- data[!is.na(data[,r_grpname]) & !is.na(data[,c_grpname]), ]
	if(is.null(r_grpvals)){
		r_grpvals <- as.character(unique(data[,r_grpname]))
		row_cols=rainbow(length(r_grpvals))
	}else{
		data<- data[data[,r_grpname] %in% r_grpvals, ]
	}
	if(is.null(c_grpvals)){
		c_grpvals <- as.character(sort(unique(data[,c_grpname])))
	}
	
	
	
	numtb <- table(data[,r_grpname], data[,c_grpname])[r_grpvals, c_grpvals]
	numtb [is.na(numtb)] <- 0
	print(numtb)
	if(perc_byrow){
		perc_tb <- numtb/rowSums(numtb)
		tlnum <- table(data[,r_grpname])[r_grpvals]
		desc <- paste(desc, c("# of", r_grpname, paste(r_grpvals,"=", tlnum, ",", sep="")), collapse=" ")
		beside<- T
	}else{
		perc_tb <- t(t(numtb)/colSums(numtb))
		tlnum <- table(data[,c_grpname])[c_grpvals]
		desc <- paste(desc, c("#|", c_grpname,r_grpname, paste(c_grpvals,"=", tlnum, ",", sep="")), collapse=" ")
		beside<- F
	}
	barplot(perc_tb,col=row_cols,xlab=c_grpname,ylab="Percent in group",beside=beside,border=T,main=desc,xlim=xlim, ylim=ylim)
	
	if(iflegend){
		plot(c(0,1),c(0,1),type="n",main=r_grpname)
		legend(0,1, legend=r_grpvals, fill =row_cols)
	}	
}


##find cross of two density plot
find_density_cross<-function(a,b){
	a<- a[!is.na(a) & abs(a)!=Inf]
	b<- b[!is.na(b) & abs(b)!=Inf]
	Fa<-ecdf(a)
	Fb<-ecdf(b)
	all_val <- c(a,b)
	#all_x <- seq(min(all_val),max(all_val),length.out=1000)
	all_val[ which.max(abs(Fa(all_val)-Fb(all_val))) ]
}
##eg
#a<- rnorm(1000)
#b<- rnorm(1000,mean=3)
#plot(density(a))
#lines(density(b))
#abline(v=find_density_cross(a,b))



draw_line_profile<- function(datasheet,grp_varname,grp_values,grp_colors,plotcolus,xVals,plottitle="",xlabel="Position",ylabel="",halfwin=0,
	ymax=max(datasheet[,plotcolus],na.rm=T),
	ymin=min(datasheet[,plotcolus],na.rm=T), blackl10p=10, iflegend=T){
	if(halfwin>0){ #window size
		D1 <- datasheet[,plotcolus]
		D2 <- D1; 
		for(shifti in setdiff(-halfwin:halfwin,0) ){
			colus <- (halfwin+1):(length(plotcolus)-halfwin)
			D2[,colus] <- D2[,colus] + D1[,colus+shifti]
		}
		D2 <- D2/(2*halfwin+1)
		datasheet[,plotcolus] <- D2
		#print(datasheet[1,])
		plotcolus <- plotcolus[(halfwin+1):(length(plotcolus)-halfwin)]
		xVals <- xVals[(halfwin+1):(length(xVals)-halfwin)]
	}
	tmptb <- table(datasheet[[grp_varname]])
	plottitle <- paste(c(plottitle, "Winsize=", (2*halfwin+1), paste(names(tmptb), "=", tmptb, sep="") ), collapse=" ")
	print(plottitle)
	plot(xVals,colMeans(datasheet[datasheet[[grp_varname]]==grp_values[1],plotcolus]), ylim=c(ymin,ymax), col=grp_colors[1],ylab=ylabel,xlab=xlabel,main=plottitle)
	for(i in 1:length(grp_values)){
		datasheet1 <- datasheet[datasheet[[grp_varname]]==grp_values[i],]
		ymeans <- colMeans(datasheet1[,plotcolus],na.rm=T)
		ysds <- apply(datasheet1[,plotcolus],2,sd,na.rm=T)
		ysems <- ysds/(nrow(datasheet1)^0.5)
		points(xVals,ymeans, col=grp_colors[i])
		#s.e.m
		for(j in 1:length(xVals)){
			lines(x=rep(xVals[j],2), y=c(ymeans[j]+ysems[j]*c(1,-1)), col=grp_colors[i])
		}
	}
	if(iflegend){
		legend(x=0, y=ymax-(ymax-ymin)*0.06,legend=grp_values,fil=grp_colors) ### remove comment #
	}
	###pvalue color
	do_col_ttest <- function(column_i){
		pval <- t.test(datasheet[datasheet[[grp_varname]]==grp_values[1],column_i], datasheet[datasheet[[grp_varname]]==grp_values[length(grp_values)],column_i],na.rm=T)$p.value
		-log10(pval)
	}
	ttestpvals <- sapply(plotcolus,do_col_ttest)
	# color
	#blackl10p <- max(ttestpvals,na.rm=T)*0.8
	white_log10p <- -log10(0.05)
	#ttestpvalcols <- sapply(ttestpvals,froml10p2Col,blackl10p=blackl10p)
	#ttestpvalcols <- sapply(ttestpvals,froml10p2Col,blackl10p=10)
	ttestpvalcols <- froml10p2Col(ttestpvals, blackl10p)
	
	for(j in 1:length(xVals)){
		lines(x=rep(xVals[j],2), y=c(ymax,ymax-(ymax-ymin)*0.05), col=ttestpvalcols[j],lwd=2)
	}	
	
	##draw color bar:
	draw_l10col_bar(white_log10p, blackl10p, plottitle)
	#barcols <- froml10p2Col(c( rep(NA,10), seq(from=white_log10p,to=blackl10p, length.out=50), rep(blackl10p,10) ), blackl10p)
	#barplot(rep(0.1,71),ylim=c(0,1),xlim=c(1,90),main=plottitle, col=barcols, border = NA, space=0,axes=F, xlab="-log10(P value)",mar=c(4,4,10,1))
	#pvec <- c(white_log10p,2:blackl10p)
	#pquantile <- ecdf(seq(from=white_log10p,to=blackl10p,length.out=100))(pvec)
	#axis(side=1, at=c(1,quantile(10:60,probs=pquantile)),las=1,labels=c(0, pvec), tick=F, tck=0.2,lwd=1,mgp=c(0,0,0))
}

draw_l10col_bar <- function(white_log10p=-log10(0.05), blackl10p=3, plottitle="",col_model="rainbow", col_scale=2/3){
	barcols <- froml10p2Col(c( rep(NA,10), seq(from=white_log10p,to=blackl10p, length.out=50), rep(blackl10p,10) ), blackl10p, col_model=col_model, col_scale=col_scale)
	barplot(rep(0.1,71),ylim=c(0,1),xlim=c(1,90),main=plottitle, col=barcols, border = NA, space=0,axes=F, xlab="-log10(P value)",mar=c(4,4,10,1))
	pvec <- c(white_log10p,2:blackl10p)
	pquantile <- ecdf(seq(from=white_log10p,to=blackl10p,length.out=100))(pvec)
	axis(side=1, at=c(1,quantile(10:60,probs=pquantile)),las=1,labels=c(0, pvec), tick=F, tck=0.2,lwd=1,mgp=c(0,0,0))
}

froml10p2Col <- function(l10p,blackl10p=3,whitel10p=-log10(0.05),col_model="rainbow", col_scale=2/3){
	#if(is.na(l10p)){
	#	gray(whitel10p)
	#}else{
		l10p[is.na(l10p)]<-whitel10p
		degreeOfblack <- (blackl10p-l10p)/(blackl10p-whitel10p)
		degreeOfblack[degreeOfblack<0]<-0
		degreeOfblack[degreeOfblack>1]<-1
		if(col_model=="rainbow"){
			retcol <- unlist(sapply(degreeOfblack, function(d){rainbow(1,start=d*col_scale)}))
		}else{
			retcol <- unlist(sapply(degreeOfblack,gray))
		}
	#}
	retcol[is.na(l10p) | l10p<=whitel10p] <- gray(1)
	retcol
}


sel_samp_with_similar_background <- function(datasheet, bg_var, grp_var, grp1, grp2){ 
	rownames(datasheet)<- 1:nrow(datasheet)
	datasheet <- datasheet[datasheet[[grp_var]] %in% c(grp1, grp2),]
	datasheet <- datasheet[order(datasheet[[bg_var]]),]
	datasheet["bi_grp"]<- rep("",nrow(datasheet))
	datasheet$bi_grp[2:nrow(datasheet)] <- paste(datasheet[1:(nrow(datasheet)-1),grp_var], datasheet[2:nrow(datasheet),grp_var])
	datasheet <- datasheet[datasheet$bi_grp %in% c(paste(grp1,grp2),paste(grp2,grp1)),]
	as.numeric(rownames(datasheet)) #return rows in the original datasheet
}

##old one, delete later!
fromValue_quantile2Class <- function(value_vec,quantile_seps,class_names){ #length(quantile_seps)=length(class_names)+1
	quantile_vec <- quantile(value_vec,probs=(0:100)/100,na.rm=T)
	return_classes <- rep(NA,length(value_vec))
	for(i in 1:length(class_names)){
		return_classes[!is.na(value_vec) & value_vec>=quantile_vec[quantile_seps[i]] & value_vec<quantile_vec[quantile_seps[i+1]]]<-class_names[i]
		if(i==length(class_names)){
			return_classes[!is.na(value_vec) & value_vec==quantile_vec[quantile_seps[i+1]]] <- class_names[i]
		}
	}
	return_classes
}
##new one
fromValue_quantile2Class <- function(value_vec,quantile_seps,class_names){ #length(quantile_seps)=length(class_names)+1
	quantile_vec <- quantile(value_vec,probs=quantile_seps,na.rm=T)
	return_classes <- rep(NA,length(value_vec))
	for(i in 1:length(class_names)){
		return_classes[ !is.na(value_vec) & value_vec>=quantile_vec[i] & value_vec<quantile_vec[i+1] ]<-class_names[i]
		if(i==length(class_names)){
			return_classes[ !is.na(value_vec) & value_vec==quantile_vec[i+1] ] <- class_names[i]
		}
	}
	return_classes
}

sep_val2block <- function(invalues,blocknum){
	block_ids <- rep(NA,length(invalues))
	rank_input <- rank(invalues, ties.method="first")
	quantile_val <- quantile(rank_input,probs=(0:blocknum)/blocknum) #quantile_val[1]=min value
	for (quid in 1:blocknum){
		block_ids[!is.na(invalues) & rank_input>=quantile_val[quid] & rank_input<=quantile_val[quid+1]] <- quid
	}
	return(block_ids)
}


fromValue_break2Class <- function(value_vec,break_seps,class_names){ #length(break_seps)=length(class_names)-1
	return_classes <- rep(NA,length(value_vec))
	return_classes[!is.na(value_vec) & value_vec<break_seps[1]] <- class_names[1]
	for(i in 1:length(break_seps)){
		return_classes[!is.na(value_vec) & value_vec>=break_seps[i]] <- class_names[i+1]
	}
	return_classes
}


cut_val_basedOn_val_and_quantile<-function(vals, inival_cuts, block_num, returnMinBinVal=T){ #useful to split discreate data eg. binding site number in 3'UTR
	val_cuts=quantile( vals[vals>inival_cuts[length(inival_cuts)]], probs=0:(block_num-1)/block_num, na.rm=T )
	val_cuts=c(inival_cuts, val_cuts, Inf)
	val_cuts=val_cuts[!duplicated(val_cuts)]
	if(returnMinBinVal){
		valbins = cut(vals, breaks=val_cuts, right=F, labels=val_cuts[-length(val_cuts)] )
		valbins=as.numeric(as.character(valbins)) #return a numeric value reflect the min value in each bin
	}else{
		valbins = cut(vals, breaks=val_cuts, right=F, labels=paste("bin",1:(length(val_cuts)-1),sep="_") )
	}
	valbins
	
}


draw_density_map <- function(rankx,ranky,xname,yname,grid_num=100, plot_title=NULL, 
	low_col_perc=0.05,high_col_perc=0.95, labx=xname,laby=yname,colbar_x_lab="", 
	mycolor=rev(rainbow(100,start=0,end=4/6)), byratio=T, mid_val=NULL, 
	bandwidths="default", out.breaks=NULL,ifaxes=T, draw_color_bar=T, onlyCalculateDensity=F){
	#isbothvalued <- apply(!is.na(cbind(rankx,ranky)),1,all)
	isbothvalued <- rowSums(is.na(cbind(rankx,ranky)))==0
	rankx <- rankx[isbothvalued]
	ranky <- ranky[isbothvalued]
	corres <- cor.test(rankx,ranky)
	tnum <- length(rankx)
	
	plot_title <- paste(plot_title, "#",tnum, "p=",format(corres$p.value,scientific=T,digits=2), "rho=", round(corres$estimate,3))
	print (plot_title)
	library("MASS")
	if(bandwidths=="default"){
		f1 <- kde2d(rankx, ranky,               n=grid_num, lims=c(min(rankx[abs(rankx)!=Inf]),max(rankx[abs(rankx)!=Inf]), min(ranky[abs(ranky)!=Inf]),max(ranky[abs(ranky)!=Inf])))
	}else if( is.null(bandwidths) ){
		rankx_bins=sep_val2block(rankx, grid_num)
		ranky_bins=sep_val2block(ranky, grid_num)
		z=table(rankx_bins,ranky_bins) #a matrix
		f1= list(x=1:nrow(z), y=1:ncol(z), z=z)
	}else{
		f1 <- kde2d(rankx, ranky, h=bandwidths, n=grid_num, lims=c(min(rankx[abs(rankx)!=Inf]),max(rankx[abs(rankx)!=Inf]), min(ranky[abs(ranky)!=Inf]),max(ranky[abs(ranky)!=Inf])))
	}
	if(!onlyCalculateDensity){
		draw_density_map_fr_kde2dfun(f1, low_col_perc=low_col_perc,high_col_perc=high_col_perc,main=plot_title,labx=labx,laby=laby,colbar_x_lab=colbar_x_lab, mycolor=mycolor, byratio=byratio, mid_val=mid_val, out.breaks=out.breaks, ifaxes=ifaxes, draw_color_bar=draw_color_bar)
	}
	f1$z
}


draw_density_map_fr_kde2dfun <- function(kde2dfun,low_col_perc=0.05,high_col_perc=0.95,main="",labx="",laby="",colbar_x_lab="", mycolor=rev(rainbow(100,start=0,end=4/6)), 
	byratio=T, mid_val=NULL,out.breaks=NULL, ifaxes=T, draw_color_bar=T, ifgrid=F ){
	if(is.matrix(kde2dfun)){
		kde2dfun <- list(x=1:nrow(kde2dfun), y=1:ncol(kde2dfun), z=kde2dfun)
	}
	if(byratio){
		kde2dfun$z <- log2( kde2dfun$z/median(kde2dfun$z,na.rm=T) )
	}
	if(is.null(mid_val)){
		q_probs <- c(0,seq(from=low_col_perc,to=high_col_perc,length.out=99),1)
	}else{
		mid_q <- ecdf(kde2dfun$z)(mid_val)
		q_probs <- c(0,seq(from=low_col_perc,to=mid_q,length.out=50), seq(from=mid_q,to=high_col_perc,length.out=50)[-1], 1)
	}
	if(is.null(out.breaks)){
		out.breaks <- quantile(kde2dfun$z, prob=q_probs, na.rm=T)
	}
	
	#density plot
	image(kde2dfun, col=mycolor, breaks=out.breaks, frame.plot=F, xlab=labx,ylab=laby,main=main,axes=ifaxes)
	if(ifgrid){
		box()
		grid(nx=nrow(kde2dfun$z), ny=ncol(kde2dfun$z), col="black", lty=1)
	}
	#color scale 
	if(draw_color_bar){
		y_val<- rep(1,length(mycolor))
		colbar_x_val <- quantile(out.breaks[-c(1,length(out.breaks))], prob=seq(0,1,length.out=5), na.rm=T)
		#colbar_x_val[c(1,5)]<-out.breaks[c(2,length(out.breaks)-1)]
		colbar_x_val <- round(colbar_x_val,3)
		if(!is.null(mid_val)){
			colbar_x_val[3]<-mid_val
		}
		colbar_x_lab <- paste(colbar_x_lab, ifelse(byratio, "Log2 ratio value/median", "value"))
		barplot(y_val,col=mycolor,border = NA, space=0,axes=F,xlab=colbar_x_lab)
		if(ifaxes){axis(side=1, at=seq(0,1,0.25)*length(mycolor),las=1,labels=colbar_x_val, tck=0.2,lwd=2,mgp=c(0,0,0))}
	}
}

return_ks_spval <-function(v1,v2,method="ks.test"){
	median1 <- quantile(v1, probs=0.5, na.rm=T)
	median2 <- quantile(v2, probs=0.5, na.rm=T)
	c(sum(!is.na(v1)), median1, sign(median1-median2)*(-log10(do.call(method, list(x=v1,y=v2))$p.value)) )
}
return_2grp_test_SS <-function(v1,v2, method="ks.test"){
	median1 <- quantile(v1, probs=0.5, na.rm=T)
	median2 <- quantile(v2, probs=0.5, na.rm=T)
	sign(median1-median2)*(-log10(do.call(method, list(x=v1,y=v2))$p.value)) 
}
cal_multi_grp_test_SS<-function(datalist, vars=names(datalist), method="ks.test"){
	var_num=length(vars)
	SS_tb=matrix(rep(NA,var_num^2),  nrow=var_num )
	for(i in 1:(var_num-1)){
		for( j in (i+1):var_num ){
			SS_tb[i,j]=return_2grp_test_SS(datalist[[vars[i]]], datalist[[vars[j]]], method=method )
	}
	}
	rownames(SS_tb)=vars
	colnames(SS_tb)=vars
	SS_tb
}
	
dim2_val_tb <- function(grpvecx,grpnumx,grpvecy,grpnumy, z){
	grpx <- fromValue_quantile2Class(grpvecx,seq(0,1,length.out=grpnumx+1),1:grpnumx)
	grpy <- fromValue_quantile2Class(grpvecy,seq(0,1,length.out=grpnumy+1),1:grpnumy)
	
	res_arr <- apply(tapply(z, list(grpx, grpy), return_ks_spval, z), c(1,2), unlist)
	list(nums=res_arr[1,,], medians=res_arr[2,,], sl10pvals=res_arr[3,,])
}


#############polyA related functions

load_pAsco_data <- function(species="human"){
	if(species=="human"){
		pA_score_file <- paste(root_dir,"analyze/polyA/0_club_like/7pA_svm_res/Hs.site5.pASVM.tbl",sep="")
	}else if(species=="mouse"){
		pA_score_file <- paste(root_dir,"analyze/polyA/0_club_like/7pA_svm_res/Mm.site5.pASVM.tbl",sep="")
	}
	pA_score_data <- read.table(pA_score_file,header=T,sep="\t")
	cis_names <<- unlist(strsplit("AUE1,AUE2,AUE3,AUE4,CUE1,CUE2,CDE1,CDE2,CDE3,CDE4,ADE1,ADE2,ADE3,ADE4,ADE5", ","))
	#counts["num.row.pAscore.all"] <<- nrow(pA_score_data)
	pA_score_data
}

load_pAano_sco_data <- function(species="human"){
	pA_score_data <- load_pAsco_data(species)
	# pA anotation
	if(species=="human"){
		pA_ano_file <- paste(root_dir,"analyze/polyA/0_club_like/5b_pA_unique_ano/hg18.all.pA.unique.ano",sep="")
	}else if(species=="mouse"){
		pA_ano_file <- paste(root_dir,"analyze/polyA/0_club_like/5b_pA_unique_ano/mm8.all.pA.unique.ano",sep="")
	}
	pA_ano_data <- read.table(pA_ano_file,header=T,sep="\t")
	#counts["num.row.1pAano.all"] <- nrow(pA_ano_data)
	
	#combine pA score and pA ano
	combined_pAscore_data <- cbind(pA_ano_data,pA_score_data[match(pA_ano_data$pAid,pA_score_data$pAid),])
	combined_pAscore_data <- combined_pAscore_data[!is.na(combined_pAscore_data$SVM_sco),]
	combined_pAscore_data["club_id"]<- sub("\\.[0-9]*$","",combined_pAscore_data$pAid)
	combined_pAscore_data
}

load_lastpAsco_forgene <- function(inids, idtype="initiator", species="human"){
	pAsco_D <- load_pAano_sco_data(species)
	pAsco_D <- pAsco_D[pAsco_D$pA_type %in% c("3S","3L"),]
	pAsco_D$SVM_sco[ match(inids, pAsco_D[[idtype]]) ]
}


###
do_fisher_test<-function(count_vec){
	fisher_matrix <- matrix(count_vec,nrow=2)
	ft_result <- fisher.test(fisher_matrix)
	return( round(-sign(ft_result$estimate-1)*log10(ft_result$p.value),2) )
}



do_chisq_test <- function(count_vec){
	count_vec <- as.numeric(count_vec)
	if(all(count_vec[1:2]==0) | any(is.na(count_vec))){
		return(NA)
	}else{
		observed_nums <- count_vec[1:2]
		expected_fre <- count_vec[1:2]+count_vec[3:4]
		chi_res <- chisq.test(observed_nums, p=expected_fre, rescale.p=T, correct = T)
		round (-sign(count_vec[1]-chi_res$expected[1]) * log10(chi_res$p.value), 2)
	}
}


#############for SAM analysis (microarray)
load_sam <- function(){
	##load sam functions
	samlib_dir <- "/Home/cocochen/soft/biosoft/array/SAM/samr/R"
	sapply(list.files(samlib_dir, full.names=T), source)
}

run_sam_ana <- function(x, y, geneid=as.character(1:nrow(x)), genenames=geneid, logged2=T, delta=NULL, fdr=0.05, fc=0, if_QQplot=F, gtypes=c("up","dn","nc")){
	#x is a matrix with raw data (x must not contain NA values), y contains group info such as c(1,1,1,2,2,2), 
	
	set.seed(100)
	gene_types <- rep(NA, length(genenames)) 
	pval_out <-  rep(NA, length(genenames)) 
	ifvalued <- rowSums(is.na(x))==0
	gene_types[ifvalued] <- gtypes[3]
	data <- list(x=x[ifvalued,], y=y, geneid=geneid, genenames=genenames, logged2=logged2)
	samr.obj <- samr(data,  resp.type="Two class unpaired", nperms=100)
	pval=samr.pvalues.from.perms(samr.obj$tt, samr.obj$ttstar)
	delta.table <- samr.compute.delta.table(samr.obj, min.foldchange=fc, nvals=200)
	delta.table2 <- data.frame(delta.table)
	if(is.null(delta)){ #fdr -> delta
		idx <- which.min(abs(delta.table2$median.FDR - fdr)) #index
		delta <- delta.table2$delta[idx]
	}else{ #delta -> fdr_actural
		idx <- which.min(abs(delta.table2$delta - delta)) #indev
	}
	fdr_actural <- delta.table2$median.FDR[idx] 
	
	if(if_QQplot){
		samr.plot(samr.obj, delta, fc)
	}
	siggenes.table<-samr.compute.siggenes.table(samr.obj, delta, data, delta.table, fc) 
	st <- samr.tail.strength(samr.obj)
	ts <- format(st$ts*100, digits=3, nsmall=1)
	se <- format(st$se.ts*100, digits=3, nsmall=1)
	gene_types[ifvalued][ as.numeric(siggenes.table$genes.up[,"Row"])-1 ] <- gtypes[1]
	gene_types[ifvalued][ as.numeric(siggenes.table$genes.lo[,"Row"])-1 ] <- gtypes[2]
	pval_out[ifvalued] <- pval
	print (c("fdr=",fdr, ", fdr_actural=", fdr_actural, "#=", table(gene_types)[gtypes]), collapse=" ")
	list(out_types=gene_types, pval=pval_out, fdr=fdr_actural, numbers=table(gene_types)[gtypes])
}

####try example:
#rownum <- 2000
#x<-matrix(rnorm(rownum*6),ncol=6)
#dd<-sample(1:rownum,size=100)
#
#u<-matrix(2*rnorm(100),ncol=3,nrow=100)
#x[dd,4:6]<-x[dd,4:6]+u
#y<-c(rep(1,3),rep(2,3))
#sam_res <- run_sam_ana(x,y)
#


do_line_plot<-function(xposs, valtb, error_tb, cols=rainbow(ncol(valtb)), error_cols=cols, smooth_f=0, title='', xlab='', ylab='', 
	vline=NULL,lwd=5,xlim=NULL,ylim=NULL,if_minor_shift_x=0,add=F, ltys=rep(1,ncol(valtb)), ... ){
		#given a set of value and error bar, color, do smoothed line graph.
		sed_div_by=ifelse(smooth_f>0, (length(xposs)*smooth_f)^0.5, 1)
		sem_max=max(error_tb, na.rm=T)/sed_div_by
		bin_size=abs(xposs[2]-xposs[1])
		colums_mid=(ncol(valtb)+1)/2
		y_vals=NULL;
		for(i in 1:ncol(valtb)){
			if(smooth_f>0){
				smooth_l=lowess(xposs, valtb[,i], f=smooth_f)
				y_vals=c(y_vals,smooth_l$y)
			}else{
				y_vals=c(y_vals,valtb[,i])
			}
		}
		if(is.null(xlim)){ xlim=quantile(xposs,probs=0:1)+c(-2,2) }
		if(is.null(ylim)){ ylim=quantile(y_vals,probs=0:1,na.rm=T)+sem_max*c(-1,1) }
		if(!add){
			plot(0, xlim=xlim, ylim=ylim, type='n', ylab=ylab, xlab=xlab, main=title,... )
		}
		for(i in 1:ncol(valtb)){
			if(smooth_f>0){
				smooth_l=lowess(xposs, valtb[,i], f=smooth_f)
			}else{
				smooth_l=list(x=xposs, y=valtb[,i])
			}
			arrows(smooth_l$x+if_minor_shift_x*bin_size/5*(i-colums_mid), smooth_l$y+error_tb[,i]/sed_div_by, 
				     smooth_l$x+if_minor_shift_x*bin_size/5*(i-colums_mid), smooth_l$y-error_tb[,i]/sed_div_by, 
				     angle=90, code=3, length=0.01, col=error_cols[i], lwd=1)
			lines(smooth_l$x+if_minor_shift_x*bin_size/5*(i-colums_mid), smooth_l$y, col=cols[i], lwd=lwd,lty=ltys[i], ...)
		}
		if(!is.null(vline)){ abline(v=vline, lty=2,col="black", lwd=2) }
}


cal_running_average<-function(v, n){
	#for vector v, for each elelemnt i, calculate average for neighbor 2n+1 elements (centered by i)
	tb=NULL
	ori_index=1:length(v)
	names(v)=as.character(ori_index)
	for(i in (-n):n){
		v2=v[as.character(ori_index+i)]
		tb<-rbind(tb,v2)
	}
	as.numeric(colMeans(tb,na.rm=T))
}


sel_gene_with_similar_background<-function(data, bgvars, if_ref, binnum, selgeneNum=sum(if_ref), if_random_sel=!if_ref , use_probablity=T, return_row_ids=F){
	#given a distribution of a variable (data[,bgvar]), select genes with similar distribution from a list data[if_ref,bgvar]
	ifsel=if_ref | if_random_sel
	data1=data[ifsel,]
	if_ref=if_ref[ifsel]
	if_random_sel=if_random_sel[ifsel]
	bin_tb<-sapply(bgvars,function(bgvar){
		sep_val2block(data1[,bgvar], binnum)
	})
	data1$bin<-apply(bin_tb,1,paste,collapse=' ');
	unique_bins=unique(data1$bin)
	bin_nums <- table(factor(data1$bin[if_ref],levels=unique_bins))[unique_bins]
	randbin_nums <- table(factor(data1$bin[if_random_sel],levels=unique_bins))[unique_bins]
	bin_nums[is.na(bin_nums)]<-0
	#names(bin_nums) <- as.character(1:binnum)
	bin_nums2=bin_nums[bin_nums>0]
	if(use_probablity){
		rand_probs= bin_nums[as.character(data1$bin)[if_random_sel]] / randbin_nums[as.character(data1$bin)[if_random_sel]]
		sample_rowids <- sample( (1:nrow(data1))[if_random_sel], selgeneNum, prob=rand_probs )
	}else{
		sample_rowids = unlist( sapply(names(bin_nums2), function(binName1){
			ifrand_bin1=if_random_sel & data1$bin==as.character(binName1)
			if(any(ifrand_bin1)){
				return( sample((1:nrow(data1))[ifrand_bin1], bin_nums2[as.character(binName1)], replace=T ) )
			}else{
				print (paste("warning: when control background,", binName1, "not found in background", bgvars, bin_nums2[binName1], "found in reference"))
				return (NULL)
			}
		}))
	}
	if(return_row_ids){ #some may be duplicate
		return( (1:nrow(data))[ifsel][sample_rowids] )
	}else{
		out_if_sel=rep(F, nrow(data))
		out_if_sel[ifsel][sample_rowids]=T
		out_if_sel
	}
}


cal_tb_fisher_ps<-function(tb, row2s, col2s=NULL, comp_num=NULL, ...){ #comparerow2s and column 1 vs2, 3 vs 4 etc in tb
	if(!is.null(col2s)){
		fisher_res=fisher.test( tb[row2s,col2s],...)
		format(fisher_res$p.value * sign(ifelse(fisher_res$estimate==1,2,fisher_res$estimate)-1), scientific=T, digits=2)
	}else{
		unlist(sapply(1:comp_num, function(i){
			col2s=(i-1)*2+1:2
			fisher_res=fisher.test( tb[row2s,col2s],...)
			format(fisher_res$p.value * sign(ifelse(fisher_res$estimate==1,2,fisher_res$estimate)-1), scientific=T, digits=2)
		}))
	}
}

re_format_tb<-function(data, delete_headers=NULL, del_pattern=NULL, Inf_headers=NULL, front_headers=NULL, back_headers=NULL, ch_header_from=NULL,ch_header_to=NULL ){ #reformat table, delete columns, change Inf to 999, change headings, etc
	if(!is.null(delete_headers)){
		print (paste(c("delete:",delete_headers)) )
		data=data[,setdiff(names(data), delete_headers)]
	}
	if(!is.null(del_pattern)){
		print (paste(c("delete column pattern:",del_pattern)) )
		del_cols=grep(del_pattern,names(data))
		if(length(del_cols)>0){
			data=data[,-del_cols]
		}
	}
	if(!is.null(Inf_headers)){
		print (paste(c("Inf headers=:",Inf_headers)) )
		data[Inf_headers][data[Inf_headers]==Inf]= 999
		data[Inf_headers][data[Inf_headers]==-Inf]= -999
	}
	if(!is.null(front_headers)){
		front_headers=intersect(front_headers,names(data))
		data=data[,c(front_headers,setdiff(names(data), front_headers))]
	}
	if(!is.null(back_headers)){
		back_headers=intersect(back_headers,names(data))
		data=data[,c(setdiff(names(data), back_headers), back_headers)]
	}
	if(!is.null(ch_header_from)){
		names(data)= change_values(names(data), ch_header_from, ch_header_to)
	}
	data
}

split_read_a_huge_table<-function(inf, read_rows_eachtime=100000, apply_fun=NA){
	
	totalrow=system(paste("wc -l",inf),intern=T); totalrow=as.numeric(sub(" .*","",totalrow)); totalrow
	read_rows=0
	while(read_rows<totalrow){
		d_i=read.table(inf, header=ifelse(read_rows==0,T,F), 
			nrows=ifelse(read_rows_eachtime+read_rows>totalrow, totalrow-read_rows, read_rows_eachtime), 
			skip=read_rows, 
			sep='\t', stringsAsFactors=F, comment.char="", quote="", fill=T)
		if(read_rows>0){
			names(d_i)=names(d)
		}
		#do something for d_i
		if(!is.na(apply_fun)){
			d_i2=lapply(list(d_i), apply_fun);
			d_i=d_i2[[1]]
		}
		
		if(read_rows==0){
			d=d_i
			read_rows=1
		}else{
			d=rbind(d,d_i)
		}
		read_rows=read_rows+read_rows_eachtime
		print (paste(read_rows, nrow(d)))
	}
	nrow(d)
	d
}

from_RegStr2length<-function(s){ #eg: s="100-200|300-400"; output length 202
	v=as.numeric(unlist(strsplit(s,"\\-|\\|")))
	sum(abs(v[c(F,T)]-v[c(T,F)])+1)
}

pvalue_cutoff2type<-function(deltaVal, pvals, ss_cut, delta_valcut, outtypes=c("proximal","distal","NC")){
	if_higher_x<- deltaVal>= delta_valcut & pvals>= ss_cut
	if_higher_y<- deltaVal<= (-delta_valcut) & pvals<= (-ss_cut)
	out_type<-rep(outtypes[3], length(deltaVal)); out_type[if_higher_x]<-outtypes[1]; out_type[if_higher_y]<-outtypes[2]; 
	out_type
}
bootstrap_counts=function(v){
	pool=factor(rep(1:length(v), v), levels=1:length(v)) 
	table(sample(pool, length(pool), replace=T))[ as.character(1:length(v)) ]
}

filter_tb<-function(genetb_data, filter_names=NULL, filter_unSelVal_list=NULL, filter_val_list=NULL){
	if(!is.null(filter_names) ){
		if( length(filter_names)>0){
			print(paste("use filter: original #row=",nrow(genetb_data) ))
			if(!is.null(filter_unSelVal_list) ){ #not selected values
				names(filter_unSelVal_list) <- filter_names
				for(filter_name in filter_names){
					genetb_data <- genetb_data[!(genetb_data[[filter_name]] %in% filter_unSelVal_list[[filter_name]]), ]
					print (paste(c(filter_name,"!=", filter_unSelVal_list[[filter_name]])))
					print(paste("row# after filtring: ", nrow(genetb_data) ))
				}
			}else{ #use filter_val_list
				names(filter_val_list) <- filter_names
				for(filter_name in filter_names){
					genetb_data <- genetb_data[genetb_data[[filter_name]] %in% filter_val_list[[filter_name]],]
					print (paste(c(filter_name,"=", filter_val_list[[filter_name]])))
					print(paste("row# after filtring: ", nrow(genetb_data) ))
				}
			}
		}
	}
	genetb_data
}
combine_table_columns<-function(genetb_data, gene_grp_names=NULL, gene_grp_name){
	if("gene_grp_names" %in% ls() ){ #combine several column as gene group 
		if( length(gene_grp_names)>1){
			print(paste(c('combine several columns (', gene_grp_names, ')TO',gene_grp_name), collapse=' '))
			print(genetb_data[1,])
			genetb_data[gene_grp_name] <- apply(genetb_data[gene_grp_names],1,paste,collapse="_")
			print(table(genetb_data[gene_grp_name]))
	}}
	genetb_data
}

#temp:
#library("RColorBrewer")
# mypalette<- c(rev(brewer.pal(9,"Blues")), "white", brewer.pal(9,"Reds"))
# image(1:9,1,as.matrix(1:9),col=mypalette,xlab="Greens (sequential)",ylab="",xaxt="n",yaxt="n",bty="n")

strReverse <- function(x){
	sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
}

rev_compl_seq=function(x){
	y=strReverse(toupper(x))
	y=sapply(y, chartr, old="ACGTU",new="TGCAA")
	names(y)=x
	y
}

define_parallel_fun<-function(nCores=1){
	if (nCores > 1) {
			library("parallel")
		if (!is.loaded("mc_fork", PACKAGE = "parallel")) {
			stop("Please load first parallel package or set parameter nCores to 1...")
		}else {
			myApply <<- function(X, FUN) {
				parallel::mclapply(X, FUN, mc.cores = nCores)
			}
		}
	}else {
		myApply <<- lapply
	}
}



compare_types_among_vars<-function(out_statsf, data, type_headers1, type_simp_headers1=type_headers1, Types, 
	type_headers2=type_headers1, type_simp_headers2=type_headers2, Types2=Types, prefix="", ifout_f=!is.na(out_statsf)){
	
	if(ifout_f){write(paste("\n#",prefix,"\n#cross table: (row ~ column); "), file=out_statsf, append=T) }
	if_Symmetrical=T
	if(length(type_headers1) != length(type_headers2)){if_Symmetrical=F}
	if(!all(type_headers1 %in% type_headers2) | !all(type_headers2 %in% type_headers1)){if_Symmetrical=F}
	SS_tb=matrix(rep(NA,length(type_headers1)*length(type_headers2)), length(type_headers1)); SS_tb
	kappa_tb=matrix(rep(NA,length(type_headers1)*length(type_headers2)), length(type_headers1)); kappa_tb
	raw_num_l=list()
	raw_enrich_tb=list()
	for(i in 1:length(type_headers1)){ #rows
		for(j in 1:length(type_headers2)){ #columns
			if(j<i & if_Symmetrical){next()}
			tb <- table( factor(data[, type_headers1[i]], levels=Types), factor(data[, type_headers2[j]],levels=Types2) )[Types,Types2]
			if(all(tb[1:2,1:2]==0)){
				tb2=matrix(rep(NA,4),2)
			}else{
				tb2=chisq.test(tb[1:2,1:2])$expected
			}
			#enrich_tb=log2( tb[1:2,1:2] / tb2) #log2 ratio
			enrich_tb= (tb[1:2,1:2]/sum(tb[1:2,1:2]) - tb2/sum(tb2)) #delta percent
			fisherPval=round(do_fisher_test(tb[1:2,1:2]), 2)
			if(fisherPval==Inf){fisherPval=999}
			if(fisherPval==-Inf){fisherPval=-999}
			kappa=Kappa.test(tb[1:2,1:2])$Result$estimate
			desc=paste('desc=',prefix," ", type_headers1[i], "(r)~(c)",type_headers2[j], " SS=",fisherPval," kappa=",kappa,sep="")
			#print(desc); print( tb )
			#if(type_headers1[i] != type_headers2[j]){
				SS_tb[i,j]=fisherPval; if(if_Symmetrical){ SS_tb[j,i]=fisherPval;  }
				kappa_tb[i,j]=kappa; if(if_Symmetrical){ kappa_tb[j,i]=kappa; }
			#}
			raw_num_l[[paste(i, j)]]=tb; if(if_Symmetrical){ raw_num_l[[paste(j, i)]]=t(tb);  }
			raw_enrich_tb[[paste(i, j)]]=enrich_tb; if(if_Symmetrical){ raw_enrich_tb[[paste(j, i)]]=t(enrich_tb);  }
		}	
	}
	kappa_tb[abs(kappa_tb)==Inf]=NA
	colnames(SS_tb)=type_simp_headers2; rownames(SS_tb)=type_simp_headers1
	colnames(kappa_tb)=type_simp_headers2; rownames(kappa_tb)=type_simp_headers1
	
	all_num_tb=NULL
	all_enrich_tb=NULL
	for(i in 1:length(type_headers1)){ #rows
		all_num_tb_j=NULL
		all_enrich_tb_j=NULL
		for(j in 1:length(type_headers2)){ #cols
			all_num_tb_j=cbind(all_num_tb_j, raw_num_l[[paste(i, j)]])
			all_enrich_tb_j=cbind(all_enrich_tb_j, raw_enrich_tb[[paste(i, j)]])
		}
		all_num_tb=rbind(all_num_tb, all_num_tb_j)
		all_enrich_tb=rbind(all_enrich_tb, all_enrich_tb_j)
	}
	tb_names1=paste( rep(type_simp_headers1,each=length(Types)), Types ); print(tb_names1)
	tb_names2=paste( rep(type_simp_headers2,each=length(Types2)), Types2 ); print(tb_names2)
	print(all_num_tb)
	colnames(all_num_tb)=tb_names2
	rownames(all_num_tb)=tb_names1

	tb_names1=paste( rep(type_simp_headers1,each=2), Types[1:2] ); print(tb_names1)
	tb_names2=paste( rep(type_simp_headers2,each=2), Types2[1:2] ); print(tb_names2)
	colnames(all_enrich_tb)=tb_names2
	rownames(all_enrich_tb)=tb_names1

	#output SS_tb and kappa_tb
	
	if(ifout_f){
		result_list= list(SS_table=SS_tb, kappa_tb=kappa_tb)
		table_desc=list(SS_table="Table of -Log10 P-value of Fisher\'s Exact Test:", kappa_tb="Kappa Statistics Table:")
		write.table(cbind(rownames(all_num_tb),all_num_tb), file=out_statsf, col.names=T, row.names=F, sep="\t", quote=F, append=T)
		#write(paste("\nall_enrich_tb "), file=out_statsf, append=T) 
		#write.table(cbind(rownames(all_enrich_tb),all_enrich_tb), file=out_statsf, col.names=T, row.names=F, sep="\t", quote=F, append=T)
		for (tb_name in names(result_list) ){
			tb=result_list[[tb_name]]
			write(paste("\n#nonClustered results ",prefix, table_desc[[tb_name]]), file=out_statsf, append=T) 
			write.table(cbind(rownames(tb),tb), file=out_statsf, col.names=T, row.names=F, sep="\t", quote=F, append=T)
		}
	}else{
		result_list= list(SS_table=SS_tb, kappa_tb=kappa_tb, all_num_tb=all_num_tb)
		return(result_list)
	}
}


calcu_cor_among_vars<-function(data, vars, vars2=vars, corMeth="spe", P_max=0.05){
	cor_m=matrix(rep(NA, length(vars)*length(vars2)),length(vars))
	if_Symmetrical=T
	if(length(vars) != length(vars2)){if_Symmetrical=F}
	if(!all(vars %in% vars2) | !all(vars2 %in% vars)){if_Symmetrical=F}
	data[vars][abs(data[vars])==Inf]=NA
	data[vars2][abs(data[vars2])==Inf]=NA
	for(i in 1:length(vars)){
		for(j in 1:length(vars2)){
			if(j<i & if_Symmetrical){next()}
			if( sum( !is.na(data[,vars[i]]) & !is.na(data[,vars2[j]]) )<=3 ){
				cor=NA
			}else{
				cor_res=cor.test(data[,vars[i]], data[,vars2[j]], method=corMeth, na.rm=T )
				cor=cor_res$estimate
				if(cor_res$p.value>P_max){cor=NA}
			}
			cor_m[i,j]=cor
			if(if_Symmetrical){ cor_m[j,i]=cor;  }
		}
	}
	rownames(cor_m)=vars
	colnames(cor_m)=vars2
	cor_m
}


combine_regu_type=function(data,type_names, cb_type_name, comb_method="union", types=c("UP","DN","NC") ){#comb_method=union, subtract or consistent
	cb_types=rep(types[3], nrow(data))
	raw_type_d=data[type_names]
	raw_type_d[is.na(raw_type_d)]="NA"
	raw_type_d[!is.na(raw_type_d) & raw_type_d=="na"]="NA"
	if(comb_method=="union"){
		cb_types[ rowSums(raw_type_d==types[1])>0 ]=types[1]
		cb_types[ rowSums(raw_type_d==types[2])>0 ]=types[2]
		cb_types[ rowSums(raw_type_d==types[2])>0 & rowSums(raw_type_d==types[1])>0 ]="NA"
		cb_types[ rowSums(raw_type_d=="NA")==length(type_names)  ]="NA"
	}else if(comb_method=="subtract"){ #substract inconsistent ones
		cb_types=raw_type_d[,1]
		cb_types[ rowSums(raw_type_d==types[2])>0 & rowSums(raw_type_d==types[1])>0 ]="NA"
	}else{ #consistent
		cb_types[ rowSums(raw_type_d==types[1])==length(type_names) ]=types[1]
		cb_types[ rowSums(raw_type_d==types[2])==length(type_names) ]=types[2]
	}
	data[cb_type_name]=cb_types
	data
}


combine_regu_type_for_data<-function(data, replicates_combine_list1=replicates_combine_list, comb_method="consistent", 
	type_prefix="pAutype_", cbtype_prefix=type_prefix ){
	for(cb_to_sample in names(replicates_combine_list1)){
		cb_fr_samples=gsub("/","_",replicates_combine_list1[[cb_to_sample]])
		if(!all(paste(type_prefix,cb_fr_samples,sep="") %in% names(data))){
			print (paste(c("not found ", paste(type_prefix,cb_fr_samples,sep="")), collapse=" "))
			next
		}
		print(paste(c("combine ", type_prefix, comb_method, cb_fr_samples, "=>", cb_to_sample),collapse=" "))
		data=combine_regu_type(data, paste(type_prefix, cb_fr_samples, sep=''), 
			paste(cbtype_prefix, cb_to_sample, sep=''), comb_method=comb_method )
	}
	data
}

change_inf<-function(ind, headers){
	for(header1 in headers){
		ind[,header1][ind[,header1]==Inf]=999
		ind[,header1][ind[,header1]==-Inf]=-999
	}
	ind
}

sel_gene_by_GO=function(organism="mm",go_dir=paste(root_dir,"/analyze/go/build_go2gid/2go2gene/",sep=""), sel_GO){
	go2geneid_file=paste(go_dir,organism,".go2gene.tbl", sep="")
	go2geneid_d<-read.table(go2geneid_file, header=T, sep="\t")
	sel_gids<-go2geneid_d$Gene_ID[go2geneid_d$GO_ID %in% sel_GO]
	sel_gids
}

cal_delta_median<-function(data, val_var, type_var, types, if_bg=!(data[,type_var] %in% types), reSample_time=50 ){
	ifsel= !is.na(data[,type_var])  & !is.na(data[,val_var])
	data=data[ifsel, ]
	if_bg=if_bg[ifsel]
	
	return_tb=NULL
	for(type1 in types){
		vals=data[data[,type_var] %in% type1, val_var]
		val_median=median(vals)
		bg_vals=data[if_bg, val_var]
		bg_medians=NULL; pvals=NULL
		for(i in 1:reSample_time){
			bg_vals_resample=sample(bg_vals, length(vals), replace=T )
			bg_medians=c(bg_medians, median(bg_vals_resample))
			pvals=c(pvals, ks.test(vals, bg_vals_resample)$p.value)
		}
		return_tb=cbind(return_tb, c(val_median-median(bg_medians), median(pvals)))
	}
	rownames(return_tb)=c("delta.median","pval")
	colnames(return_tb)=types
	return_tb
}

resampling_cal_delta_q_value<-function( data, val_var, type_var, types, if_bg=!(data[,type_var] %in% types), reSample_time=5000 ){
	#an alternative of KS test, but resolve the number of cases issue
        ifsel= !is.na(data[,type_var])  & !is.na(data[,val_var])
        data=data[ifsel, ]
        if_bg=if_bg[ifsel]
	
        return_tb=NULL
        for(type1 in types){
                vals=data[data[,type_var] %in% type1, val_var]
                bg_vals=data[if_bg, val_var]
                obs_delta=median(vals)-median(bg_vals)
                pool=c(vals,bg_vals)
                bs_deltas=NULL
                bs_deltas=unlist(myApply(1:reSample_time, function(i){
                	 bg_ids=sample(1:length(pool), length(vals), replace=F )
                	 bs_delta=median(pool[-bg_ids]) - median(pool[bg_ids])
                	 bs_delta
                }))
                #for(i in 1:reSample_time){
                #	 bg_ids=sample(1:length(pool), length(vals), replace=F )
                #	 bs_delta=median(pool[-bg_ids]) - median(pool[bg_ids])
                #	 bs_deltas=c(bs_deltas,bs_delta)
                #}
                qval=sum(abs(bs_deltas)>abs(obs_delta))/reSample_time
                return_tb=cbind(return_tb, c(obs_delta, qval))
        }
        rownames(return_tb)=c("delta.median","qval")
        colnames(return_tb)=types
        return_tb
}


compare_bg_controlled_vals<-function(data, control_var, var_name, grp_var, grp_vals, ctrl_grps, bin_bg_num=10, sample_time=10){
	#for a data frame, compare var_name for each group in grp_var (eg. AUUUA # in UTR), vs. ctrl_grps (eg. no AUUUA); samping values in ctrl_grps to have similar distribution of control_var (eg. UTR length)
	val_l=list()
	pvals=NULL
	control_var_Ps=NULL
	for(cisNumBinType1 in grp_vals){
		if_cis=data[,grp_var]==cisNumBinType1
		if_bg=data[,grp_var] %in% ctrl_grps
		
		val_l[[paste(cisNumBinType1,"exp",sep=".")]]=NULL
		obs_vals=data[if_cis, var_name]
		control_var_P1s=NULL
		pvals1=NULL

		#(use the following multi-core will cause error in some image function )
		#res=myApply(1:sample_time, function(i){
		#	bg_sample_rows= sel_gene_with_similar_background(data, control_var, if_cis, bin_bg_num,  if_random_sel=if_bg, 
		#		use_probablity=F, return_row_ids=T)  
		#	control_var_P=ks.test(data[bg_sample_rows,control_var], data[if_cis,control_var])$p.value #  this should  be less significant
		#	#control_var_P1s<-c(control_var_P1s, control_var_P)
		#	expVals=data[bg_sample_rows, var_name]
		#	#val_l[[paste(cisNumBinType1,"exp",sep=".")]] <- c(val_l[[paste(cisNumBinType1,"exp",sep=".")]], expVals )
		#	#pvals1=c( pvals1, ks.test(obs_vals, expVals)$p.value )
		#	list(expVals, ks.test(obs_vals, expVals)$p.value, control_var_P)
		#})
		#val_l[[paste(cisNumBinType1,"exp",sep=".")]]=unlist(sapply(res, function(l){l[[1]]}))
		#pvals=c(pvals,   median(unlist(sapply(res,function(l){l[[2]]})), na.rm=T) )
		#control_var_Ps=c(control_var_Ps,   median(unlist(sapply(res,function(l){l[[3]]})), na.rm=T) )

		for(i in 1:sample_time){
			bg_sample_rows= sel_gene_with_similar_background(data, control_var, if_cis, bin_bg_num,  if_random_sel=if_bg, 
				use_probablity=F, return_row_ids=T)  
			control_var_P=ks.test(data[bg_sample_rows,control_var], data[if_cis,control_var])$p.value #  this should  be less significant
			control_var_P1s<-c(control_var_P1s, control_var_P)
			expVals=data[bg_sample_rows, var_name]
			val_l[[paste(cisNumBinType1,"exp",sep=".")]] <- c(val_l[[paste(cisNumBinType1,"exp",sep=".")]], expVals )
			pvals1=c( pvals1, ks.test(obs_vals, expVals)$p.value )
		}
		pvals=c(pvals,   median(pvals1, na.rm=T) )
		control_var_Ps=c(control_var_Ps,   median(control_var_P1s, na.rm=T) )
		
		val_l[[paste(cisNumBinType1,"obs",sep=".")]] = obs_vals
	}
	list(val_l, pvals, control_var_Ps)
}


do_pairs_scatter_plot<-function(data, val_headers, corMeth="spe", cor.cex=2, text_cex=1.5, pch=".", col="black"){
	panel.cor <- function(x, y, digits=2, prefix="",  ...){
	    usr <- par("usr"); on.exit(par(usr))
	    par(usr = c(0, 1, 0, 1))
	    num=sum( !is.na(x) & !is.na(y) )
	    r <- round(cor.test(x, y, method=corMeth, na.rm=T)$estimate, 3)
	    txt <- paste(prefix, r,"\n#",num, sep="")
	    text(0.5, 0.5, txt, cex = cor.cex)
	}
	text.panel<-function(x, y, labels,  ...){
	 	text(0.5, 0.5, labels, cex=text_cex, col="red")
	}
	pairs(data[,val_headers], lower.panel=panel.cor, corMeth=corMeth, 
		text.panel=text.panel, 
		pch=pch, col=col)
	
}



from_combn2str<-function(totalnum, selnum, all_char="N", sel_char="Y"){
	rmID_matrix=combn(1:totalnum, selnum)
	timeSet_str_m=matrix(rep(all_char, totalnum*ncol(rmID_matrix)), ncol=totalnum)
	for(i in 1:ncol(rmID_matrix)){
		timeSet_str_m[i, rmID_matrix[,i]]=sel_char
	}
	apply(timeSet_str_m,1,paste, collapse="")
}

comp_grps_pval<-function(data, x, grpname,grpvals,test_sample_pairs="first_last",test_method="ks.test",test_grp_matrix=NULL){
	x[abs(x)==Inf]=NA
	if(is.null(test_grp_matrix)){
		if(test_sample_pairs=="first_last"){
			test_grp_matrix=rbind(grpvals[1], grpvals[length(grpvals)])
		}else if(test_sample_pairs=="first_oths"){ #first one vs others
			test_grp_matrix=rbind( rep(grpvals[1], length(grpvals)-1), grpvals[-1] )
		}else if(test_sample_pairs=="adjacent"){
			test_grp_matrix=rbind( grpvals[-length(grpvals)], grpvals[-1] )
		}
	}
	pval=apply(test_grp_matrix, 2, function(grp_2vals){
		grp1_xs<- x[!is.na(x) & data[[grpname]]==grp_2vals[1]]
		grpn_xs<- x[!is.na(x) & data[[grpname]]==grp_2vals[2]]
		if(length(grp1_xs)>=3 & length(grpn_xs)>=3){
			pval1 <- format(do.call(test_method, list(x=grp1_xs, y=grpn_xs))$p.value, scientific=T, digits=2)
		}else{
			pval1 <- NA
		}
		pval1
	})			
	pval
}

filter_FC_using_readNumber<-function(data, readNum_names, FC_name, readNumSum_cut=20, readNumMin_cut=1){
	total_rnums=rowSums(data[,readNum_names],na.rm=T)
	min_readnums=apply(data[,readNum_names],1,min, na.rm=T)
	ifsel=!is.na(total_rnums) & total_rnums>=readNumSum_cut & !is.na(min_readnums) & min_readnums>readNumMin_cut
	print(paste("readNumSum_cut=",readNumSum_cut, "; readNumMin_cut=",readNumMin_cut, "valid",FC_name,"#=",sum(ifsel) ))
	data[!ifsel, FC_name]=NA
	data
}


myDist_fun=function(x, dist_fun="selfDefined", method="pearson"){
	if(dist_fun=="selfDefined"){ 
		cor.res=cor(as.matrix( t(x)), use="pairwise.complete.obs", method=method)
		dissimilarity =1-cor.res
		dist.res=as.dist(dissimilarity)
		dist.res
	}else{
		do.call(dist, list(x=x, method=method) ) 
	}
}

cal_ss3_score<-function(seqs, prog_dir=paste(root_dir,"analyze/AS/ss_score/",sep="") ){
	previous_dir=getwd()
	setwd(prog_dir)
	scores=NULL
	while(length(seqs)>0){
		if(length(seqs)>2000){
			cmd=paste(c("perl 2cal_23nt_ss3_score.pl 0 ",seqs[1:2000]), collapse=" ")
			seqs=seqs[-(1:2000)]
		}else{
			cmd=paste(c("perl 2cal_23nt_ss3_score.pl 0 ",seqs), collapse=" ")
			seqs=NULL
		}
		scores1=system(cmd, intern=T)
		scores=c(scores, as.numeric(scores1))
	}
	setwd(previous_dir)
	scores
}


change_plateau_val<-function(v, change_from=c(999,Inf), change_to_max_fold=1.5, changeToMin=5 ){ #change unplotable values (Inf or 999) to a fixed plateau value
	v_sign=sign(v)
	abs_v=abs(v)
	if_change=!is.na(v) & abs_v %in% change_from
	change_to_abs=max(changeToMin, change_to_max_fold*max(abs_v[!if_change],na.rm=T)  )
	v[if_change]=change_to_abs * v_sign[if_change]
	v
}
