#for RNA-seq differential expression (DE) analysis, calculate rpkm change and p-value (fisher's exact test, chi square test, DEseq2)
#3/3/2014, can allow parallal computing when nCores was set >1
#3/31/2015, add function to use gene expression normalization method in edgeR (GexNorm_method)
#9/1/2015, allow analysis using DESeq2
#3/16/2016, draw clustering of raw samples based on RPKM or RPM, identify sample outliers. (what2do=clust_rawSample_exp)
#3/17/2016, can calculate gene expression with a phenotype (eg., tumor size), what2do=calCorWithPhenotype
#1/29/2019, can set fold change to 1, if P value>setFC1_p_cut
#7/15/2019, modified to allow using raw P-value (instead of default adjusted P-value) in DESeq2, set useP_type="pvalue"
#6/7/2021, add function to allowing selection of significantly regulated genes using an higher_avg_RPKM_cutoff (between ) cutoff (eg. for the comparion two groups of samples, the higher one of the average rpkm within a group should be higher than a cutoff)
#7/21/2021, added function to output a table (topGeneList_tb) of top DEGs for each treatment

#available what2do: all (default), formatGeneTb, clust_rawSample_exp, calCorWithPhenotype, make_plots, rawSample_ScatterPlot

args<- commandArgs(TRUE)
print(args)

args_v<-NA
if(length(args)>1){
	args_v <- args[c(FALSE, TRUE)]
	names(args_v) <- sub("\\-+","", args[c(TRUE, FALSE)])
	print(args_v)
}
source("../SharedCodes/Rfunc.inc.R")


#parameters from command line 
what2do<- "all"
if(!is.na(args_v["what2do"])){what2do=unlist(strsplit(as.character(args_v["what2do"])," "))}
out_prefix<- ifelse(is.na(args_v["out_prefix"]), "", args_v["out_prefix"])
p_cut <- ifelse(is.na(args_v["p_cut"]), 0.01, as.numeric(args_v["p_cut"])); 
setFC1_p_cut <- ifelse(is.na(args_v["setFC1_p_cut"]), 0.5, as.numeric(args_v["setFC1_p_cut"]));  #set fold change to 1, if P value>setFC1_p_cut
if_correct_pval<- ifelse(is.na(args_v["if_correct_pval"]), FALSE, args_v["if_correct_pval"]=="1" ); if_correct_pval
pval_fun<- ifelse(is.na(args_v["pval_fun"]), "DESeq2", args_v["pval_fun"]); #fisher chisq none DESeq2
foldchange_cut <- ifelse(is.na(args_v["foldchange_cut"]), 1.5, as.numeric(args_v["foldchange_cut"])); foldchange_cut
IfCal_RPKM <- ifelse(is.na(args_v["IfCal_RPKM"]), TRUE, args_v["IfCal_RPKM"]=="1" ); IfCal_RPKM
IfCal_RPM <- ifelse(is.na(args_v["IfCal_RPM"]), TRUE, args_v["IfCal_RPM"]=="1" ); IfCal_RPM  ##calculate RPM, eg. for small RNA-seq; 3'READS set T, otherwise set F
gene_model<- ifelse(is.na(args_v["gene_model"]), "refseqcds", args_v["gene_model"]); gene_model ##club refseq enscds repeat LTR Rcluster pAR2club
tlRnum_using<-ifelse(is.na(args_v["tlRnum_using"]), "colSum", args_v["tlRnum_using"]); tlRnum_using #"colSum", using read number column sum (Read mapped to all gene models) as total read number; "uniqueMapped", using uniquely mapped reads number
comb_sample_l=list()
if(!is.na(args_v["comb_sample_l"])){ comb_sample_l<- string2list(args_v["comb_sample_l"]) }; comb_sample_l
nCores <- ifelse(is.na(args_v["nCores"]), 1, as.numeric(args_v["nCores"])); 
higher_avg_RPKM_cutoff <- ifelse(is.na(args_v["higher_avg_RPKM_cutoff"]), 1, as.numeric(args_v["higher_avg_RPKM_cutoff"])); higher_avg_RPKM_cutoff

geno<- ifelse(is.na(args_v["geno"]), "hg19", args_v["geno"])
study_name<- ifelse(is.na(args_v["study_name"]), NA, args_v["study_name"]) ##project name
output_root<- ifelse(is.na(args_v["output_root"]), paste("05.Gene_DE/",study_name,"/",sep=""), args_v["output_root"]) 
out_plot_root=paste(output_root, out_prefix,gene_model,".",pval_fun,".format_tb/GexPlot/",sep="");

if( !is.na(args_v["gene_rnum_f"]) ){ 	gene_rnum_f=args_v["gene_rnum_f"] } #gene read counts file
gene_rnum_f;
gene_len_var <- ifelse(is.na(args_v["gene_len_var"]), NA, args_v["gene_len_var"]) #gene length variable (column name)
test_samples<- unlist(strsplit(as.character(args_v["test_samples"])," ")) #must set a value 
ref_samples<- unlist(strsplit(as.character(args_v["ref_samples"])," ")) #must set a value 
if(!is.na(args_v["other_samples"])){ 
	other_samples<- unlist(strsplit(as.character(args_v["other_samples"])," ")) 
}
if( any( c( is.na(test_samples) , is.na(ref_samples) ) ) ){		print ("test_samples or ref_samples is NA!"); quit(); 	}
sel_gex_ids=1:length(test_samples)
if(!is.na(args_v["sel_gex_ids"])){ 
	sel_gex_ids<- as.numeric( unlist(strsplit(args_v["sel_gex_ids"]," ") ) ) 
}

GexNorm_method<-ifelse(is.na(args_v["GexNorm_method"]), "totalMillionRead", args_v["GexNorm_method"]);GexNorm_method; #totalMillionRead(RPM or RPKM), TMM (weighted trimmed mean of M-values, used on edgeR)

out_file <- paste(output_root, out_prefix, gene_model,
	ifelse( GexNorm_method=="totalMillionRead","",paste(".",GexNorm_method,sep="") ),
	".",pval_fun,".p",sep="");
if( !is.na(args_v["out_file"]) ){ out_file=args_v["out_file"] }
g_ano_f<-ifelse(is.na(args_v["g_ano_f"]), NA, args_v["g_ano_f"])
anof_idname<-ifelse(is.na(args_v["anof_idname"]), NA, args_v["anof_idname"])
rnum_idname<-ifelse(is.na(args_v["rnum_idname"]), "refseqid", args_v["rnum_idname"])
geneSym_name<-ifelse(is.na(args_v["geneSym_name"]), "gene_symbol", args_v["geneSym_name"])
geneType_anof=paste("../../club/21transcriptAno/",geno,".transcript.ano", sep="") #have gene type annotation
if(!is.na(args_v["geneType_anof"])){
	geneType_anof=ifelse(args_v["geneType_anof"] %in% c("NA","na"), NA, args_v["geneType_anof"])
}
geneType_ano_idName=ifelse(is.na(args_v["geneType_ano_idName"]), "name", args_v["geneType_ano_idName"])

#for what2do=formatGeneTb
formatGeneTbName<-ifelse(is.na(args_v["formatGeneTbName"]), "", args_v["formatGeneTbName"])
L2FC_scale_Max=ifelse(is.na(args_v["L2FC_scale_Max"]), 2, as.numeric(args_v["L2FC_scale_Max"]) )
filter_header=ifelse(is.na(args_v["filter_header"]), NA, args_v["filter_header"]) #used for filter the table when plotting the clustering plot
filter_text=ifelse(is.na(args_v["filter_text"]), NA, unlist(strsplit(as.character(args_v["filter_text"])," "))) #only plot rows with these filter_text values
max_geneList_num=ifelse(is.na(args_v["max_geneList_num"]), 100, as.numeric(args_v["max_geneList_num"]) )

comment_f=paste(out_file,"comment.tbl",sep=".")

define_parallel_fun(nCores=nCores)
print(paste("nCores=",nCores))
mkdir_if_not_exist(out_plot_root)

##common setup:
gexch_types <- c("up","dn","nc","na")
if(!("other_samples" %in% ls())){	other_samples<-NULL }
if(!is.na(args_v["all_sample_name"])){ 
	all_sample_name<- unlist(strsplit(as.character(args_v["all_sample_name"])," ")) 
}else if("sample_grps" %in% ls()){
	all_sample_name=unique(unlist(sample_grps));
}else{
	all_sample_name <- unique(c(ref_samples,test_samples,other_samples))
}
print(all_sample_name)

reads_num_names <- paste("num_",all_sample_name,sep="")
rpkm_names <- paste(ifelse(IfCal_RPM,"rpm_","rpkm_"), all_sample_name,sep="")
samplePairs=paste(test_samples,ref_samples,sep="_")
gchtp_names <- paste("GexType_",test_samples,"_",ref_samples,sep=""); gchtp_names
l2gch_names <- paste("L2FC_",test_samples,"_", ref_samples,sep="")
pval_names <- paste("SLog10P_",test_samples, "_", ref_samples, sep="")
adjPval_names <- paste("adjSLog10P_",test_samples, "_", ref_samples, sep="")
corr_p_cut <- p_cut/ifelse(if_correct_pval ,nrow(cmb_d),1) #Bonferroni corrected p value cut
isexpP_names <- paste("expP_",all_sample_name,sep="") #is expressed p value (poisson distribution)
isexp_names <- paste("exp_",all_sample_name,sep="")

if(pval_fun %in% c("DESeq2") ){ ##automatic define P-value type
	useP_type="padj"
}else{
	useP_type="pvalue"
}
if(!is.na(args_v["useP_type"])){ useP_type=args_v["useP_type"] }
if(useP_type %in% c("padj") ){ ##automatic define P-value type
	use_pval_names=adjPval_names
}else{
	use_pval_names=pval_names
}

min_TotalreadNum <- ifelse(is.na(args_v["min_TotalreadNum"]), length(reads_num_names), as.numeric(args_v["min_TotalreadNum"])); 


#for what2do=calCorWithPhenotype
# phenotypeName="tumorVol"
# phenotypeSamples=unlist(strsplit("Veh_1.1 Veh_1.2 Veh_1.3 Veh_1.7 Veh_1.8 PTC596_5mg_kg_2.1 PTC596_5mg_kg_2.4 PTC596_5mg_kg_2.9 PTC596_10mg_kg_4.1 PTC596_10mg_kg_4.3 PTC596_10mg_kg_4.9 PTC596_15mg_kg_6.4 PTC596_15mg_kg_6.8 PTC596_15mg_kg_6.10"," ")) 
# phenotypeVals=as.numeric(unlist(strsplit("1536 1797 1447 1856 2016 1744 1888 1712 961 950 1090 562 377 325"," ")) )
phenotypeName<-ifelse(is.na(args_v["phenotypeName"]), NA, args_v["phenotypeName"])
if(!is.na(args_v["phenotypeSamples"])){ phenotypeSamples=unlist(strsplit(as.character(args_v["phenotypeSamples"])," ")) }
if(!is.na(args_v["phenotypeVals"])){ phenotypeVals=as.numeric(unlist(strsplit(as.character(args_v["phenotypeVals"])," "))) }

figure_out_format=ifelse(capabilities()['png'],"png","pdf")
figure_size_fac=ifelse(figure_out_format=="png",1,1/72) #convert pixel in png to inch in pdf

################### RUN
if(any(what2do %in% c("all")) ){
	print(paste("read gene_rnum_f:",gene_rnum_f))
	gene_rnum_D <- read.table(gene_rnum_f, sep="\t", header=T, comment.char="", quote='')
	gene_rnum_D <- gene_rnum_D [,apply(!is.na(gene_rnum_D),2,any)]

	###combine with gene anotation
	if( !("g_ano_f" %in% ls()) ){g_ano_f<-NA}
	if(!is.na(g_ano_f)){
		print(paste("read g_ano_f: ",g_ano_f))
		g_ano_d <- read.table(g_ano_f, header=T, sep="\t", quote="", comment.char="")
		g_ano_d<-g_ano_d[setdiff(names(g_ano_d),c("exon_starts","exon_ends","refseqid.1"))]
		if("gene_Biotype" %in% names(g_ano_d) ){
			g_ano_d=g_ano_d[order(!is.na(g_ano_d$gene_Biotype) & g_ano_d$gene_Biotype %in% c("protein_coding"), decreasing=T ), ] #put protein coding genes on top
		}
		if( !is.na(anof_idname) & !is.na(rnum_idname) ){
			#cmb_d <- cbind(g_ano_d,gene_rnum_D[match(g_ano_d[,anof_idname], gene_rnum_D[,rnum_idname]),])
			cmb_d <- gene_rnum_D
			update_names=setdiff(names(g_ano_d),names(gene_rnum_D))
			cmb_d[,update_names]=g_ano_d[match(gene_rnum_D[,rnum_idname], g_ano_d[,anof_idname]), update_names]
		}else{
			if(nrow(g_ano_d)!=nrow(gene_rnum_D)){
				print("Error: row number in gene_rnum_D do not match g_ano_d!")
				quit()
			}else{
				cmb_d <- cbind(g_ano_d, gene_rnum_D[setdiff(names(gene_rnum_D),names(g_ano_d))])
			}
		}
	}else{
		cmb_d <- gene_rnum_D
	}
	rm ("gene_rnum_D","g_ano_d")
	save.image( paste(out_file,".img",sep="") )

	#combine reads number
	raw_num_names <- grep("num_",names(cmb_d),value=T)
	cmb_d[raw_num_names][is.na(cmb_d[raw_num_names])] <- 0
	for(i in 1:length(reads_num_names)){
		reads_num_name=reads_num_names[i]
		if( ! (reads_num_name %in% names(cmb_d)) ){
			if(all_sample_name[i] %in% names(comb_sample_l)){
					combine_frs <- paste('num_',comb_sample_l[[all_sample_name[i]]],sep='')
			}else{
				search_patt=paste("num_",".*",all_sample_name[i],sep="")
				combine_frs=grep(search_patt,names(cmb_d),value=T)
			}
			print (paste(c("combine read number for sample",reads_num_name, " from", combine_frs), collapse=" " ) )
			cmb_d[reads_num_name] <- rowSums(cmb_d[combine_frs], na.rm=T)
		}
	}
	cmb_d[reads_num_names][is.na(cmb_d[reads_num_names])] <- 0

	##remove redundancy in cmb_d
	#sort with read number; remove rows with duplicated rnum_idname
	cmb_d=cmb_d[order(rowSums(cmb_d[reads_num_names]),decreasing=T),]; nrow(cmb_d)
	if(!is.na(rnum_idname)){
		cmb_d=cmb_d[!duplicated(cmb_d[,rnum_idname]),]; nrow(cmb_d)
	}
	#remove transcripts with too small number of reads
	ifExpressed=rowSums(cmb_d[reads_num_names])>=min_TotalreadNum; sum(ifExpressed)
	cmb_d=cmb_d[ifExpressed,]; nrow(cmb_d)

	##calculate scaling factor if GexNorm_method=="TMM"
	if(GexNorm_method=="TMM"){
		library('edgeR')
		scaling_facs=calcNormFactors(cmb_d[,reads_num_names])
	}else{
		scaling_facs=rep(1,length(reads_num_names))
	}
	names(scaling_facs)=reads_num_names

	#calculate t_readsnum
	if(tlRnum_using=="colSum"){
		t_readsnum <- apply(cmb_d[reads_num_names],2,sum, na.rm=T); 
	}else if(tlRnum_using=="uniqueMapped"){
		t_readsnum <- tl_mapped_readsnum[all_sample_name]
		names(t_readsnum) <- reads_num_names
	}
	#calculate rpkm or rpm
	t_readsnum=round( t_readsnum*(scaling_facs) )
	t_readsnum
	if(gene_model=="club"){ 
		gene_len_var <- paste("_len",sep="")
		print (paste("gene_len_var=",gene_len_var),quote=F)
	}else if(gene_model=="repeat"){ 
		if("Trnum_cut" %in% ls()){ #filter low read number rows
			print (paste("original rows:", nrow(cmb_d)))
			cmb_d <- cmb_d[rowSums(cmb_d[reads_num_names],na.rm=T)>=Trnum_cut,]
			print (paste("after removing rows with total read number<",Trnum_cut, nrow(cmb_d), "rows remaining."))
		}
		cmb_d[,gene_len_var] <- abs(cmb_d$repTo-cmb_d$repFr)
	}
	if( is.null(cmb_d[1,rpkm_names[1]]) | IfCal_RPKM | IfCal_RPM){
		if(IfCal_RPM){
			cmb_d[rpkm_names] <- round(t( t(cmb_d[reads_num_names]) / (t_readsnum/1000000) ) , 3)
		}else{
			if(gene_len_var=="cds_len" & "transcript_len" %in% names(cmb_d)){
				if_change=!is.na(cmb_d$cds_len) & cmb_d$cds_len==0
				cmb_d$cds_len[if_change]=cmb_d$transcript_len[if_change]
			}
			cmb_d[rpkm_names] <- round(t( t(cmb_d[reads_num_names]/(cmb_d[,gene_len_var]/1000)) / (t_readsnum/1000000) ) , 3)
		}
	}

}

if(any(what2do %in% c("all")) ){
	if(pval_fun =="DESeq2"){
		library("DESeq2")
		library("BiocParallel")
		register(MulticoreParam(nCores))
		for(i in 1:length(test_samples)){
			samp_pair <- c(test_samples[i], ref_samples[i])
			names(samp_pair)=c("test","ctrl")

			gchtp_name <- gchtp_names[i]
			l2gch_name <- l2gch_names[i]
			print (samp_pair)
			reads_num_2names <- paste("num_",samp_pair,sep="")
			#build countData and colData used for DESeqDataSetFromMatrix
			condition=NULL
			DeSeq_samples=NULL
			for(sample_type in names(samp_pair) ){
				sample1=samp_pair[sample_type]
				if(sample1 %in% names(comb_sample_l) ){
					DeSeq_samples=c(DeSeq_samples, comb_sample_l[[sample1]])
					condition=c(condition, rep(sample_type, length(comb_sample_l[[sample1]])) )
				}else{
					DeSeq_samples=c(DeSeq_samples, sample1)
					condition=c(condition, sample_type )
				}
			}
			coldata=data.frame(condition=condition); rownames(coldata)=make.unique(DeSeq_samples)
			countData=cmb_d[paste("num_",DeSeq_samples,sep="")]
			if_no_counts=rowSums(countData)<1
			colnames(countData)=make.unique(DeSeq_samples)
			rownames(countData)=cmb_d[,rnum_idname]

			dds <- DESeqDataSetFromMatrix(countData = countData[!if_no_counts,], colData = coldata, design = ~ condition )
			dds <- DESeq(dds, parallel=TRUE)
			res <- results(dds)
			print(summary(res))

			cmb_d[,l2gch_name]=NA
			cmb_d[!if_no_counts,l2gch_name] = res$log2FoldChange
			
			if(useP_type %in% c("padj") ){
				cmb_d[,adjPval_names[i]]=NA
				cmb_d[!if_no_counts,adjPval_names[i]]= -sign(res$log2FoldChange) * log10(res$padj)
			}else{
				cmb_d[,pval_names[i]]=NA
				cmb_d[!if_no_counts,pval_names[i]]= -sign(res$log2FoldChange) * log10(res$pvalue)
			}
		}
		cmb_d[l2gch_names][ is.na(cmb_d[use_pval_names]) | abs(cmb_d[use_pval_names])< abs(log10(setFC1_p_cut))  ] = 0
		rm("dds","res")

	}else{ #fisher's exact test
		if(!is.null(test_samples)){
			#fold change:
			cmb_d[l2gch_names] <- round(log2( cmb_d[paste(ifelse(IfCal_RPM,"rpm_","rpkm_"),test_samples,sep="")] / cmb_d[paste(ifelse(IfCal_RPM,"rpm_","rpkm_"),ref_samples,sep="")] ) ,3)
			
			#fisher's exact test/chi2 test,  for gene expression change of every pair
			for(i in 1:length(test_samples)){
				samp_pair <- c(test_samples[i], ref_samples[i])
				gchtp_name <- gchtp_names[i]
				l2gch_name <- l2gch_names[i]
				print (samp_pair)
				reads_num_2names <- paste("num_",samp_pair,sep="")
				numtb <- cbind(cmb_d[reads_num_2names], t(t_readsnum[reads_num_2names]-t(cmb_d[reads_num_2names])))
				pval_name <- pval_names[i]
				numtb=as.matrix(numtb)
				#cmb_d[pval_name] <- apply(numtb,1,paste("do",pval_fun,"test",sep="_") )
			 	
			 	if_do_test=rowSums(numtb[,1:2],na.rm=T)>=7; table(if_do_test)
			 	if(pval_fun !="none"){
				 	cmb_d[pval_name]=0
				 	cmb_d[if_do_test,pval_name] <- unlist(myApply((1:nrow(numtb))[if_do_test], function(row_i){ 
				 		do.call(paste("do",pval_fun,"test",sep="_"), list(count_vec=round(numtb[row_i,])) ) 
				 	 } ))
			 	}
			}
			cmb_d[l2gch_names][ is.na(cmb_d[pval_names]) | abs(cmb_d[pval_names])< abs(log10(setFC1_p_cut))  ] = 0
		}

	}
	


	save.image( paste(out_file,".img",sep="") )
}	
		


if(any(what2do %in% c("all")) ){
	####add gene annotations to the table
	if(! "gene_Biotype" %in% names(cmb_d) ){
		if( !is.na(rnum_idname) & !is.na(geneType_anof) ){
			gene_anod<- read.table(geneType_anof, header=T, sep="\t", quote="", comment.char="")
			if(any(cmb_d[,rnum_idname] %in% gene_anod[,geneType_ano_idName])){
				cmb_d$gene_Biotype<- gene_anod[match(cmb_d[,rnum_idname],gene_anod[,geneType_ano_idName] ), "gene_Biotype"]
				cmb_d$gene_Biotype<- change_values(cmb_d$gene_Biotype, c("lincRNA","protein_coding"), c("lncRNA","mRNA"))
				table(cmb_d$gene_Biotype)
				#intersect(cmb_d$gene_symbol[cmb_d$gene_Biotype=='lncRNA'], cmb_d$gene_symbol[cmb_d$gene_Biotype=='mRNA']) #check
			}
			rm("gene_anod")
		}
		
	}

	#####output
	desc_str <- paste (
		"#foldchange_cut=",foldchange_cut,
		"\n#pval_fun=",pval_fun, "\n#rows=",nrow(cmb_d),
		"\n#p_cut=",p_cut,   "\n", 
		"GexNorm_method=",GexNorm_method,"\n",
		sep=""
	)
	desc_str

	cmb_d<-cmb_d[ setdiff(names(cmb_d),c("gene_id.1")) ]
	out_file
	write(desc_str, file=comment_f)
	tb=data.frame(sample=all_sample_name, normalized_total_read=t_readsnum, scaling_fac=scaling_facs)
	write.table(tb, file=comment_f, col.names=T, row.names=F, sep="\t", quote=F, append=T)
	write.table(cmb_d, file=out_file, col.names=T, row.names=F, sep="\t", quote=F)

	save.image( paste(out_file,".img",sep="") )


}
	

#####output a formatted table of all genes or top regulated genes (for visulization, etc)
if(any(what2do %in% c("all","formatGeneTb")) ){
	cmb_d=read.table(out_file, header=T, sep="\t", quote="", comment.char="", stringsAsFactors=F)
	
	cmb_d2=cmb_d
	# if("gene_Biotype" %in% names(cmb_d2)){
	# 	cmb_d2=cmb_d2[!is.na(cmb_d2$gene_Biotype) & cmb_d2$gene_Biotype %in% c("mRNA","lncRNA","protein_coding"), ]; nrow(cmb_d2)
	# }
	#remove gene symbol with duplicated transcript IDs
	cmb_d2=cmb_d2[order(rowSums(cmb_d2[reads_num_names]),decreasing=T),]; nrow(cmb_d2)
	cmb_d2[is.na(cmb_d2$gene_symbol),geneSym_name]= cmb_d2[is.na(cmb_d2$gene_symbol), rnum_idname]
	cmb_d2=cmb_d2[ !duplicated(cmb_d2[,geneSym_name]),]; nrow(cmb_d2) #!is.na(cmb_d2[,geneSym_name]) &
	
	#calculate average expression
	for(sample1 in names(comb_sample_l)){
		avg_rpkm_name=paste(ifelse(IfCal_RPM,"rpm_","rpkm_"),sample1,sep="")
		raw_rpkm_names=paste(ifelse(IfCal_RPM,"rpm_","rpkm_"),comb_sample_l[[sample1]],sep="")
		if(!(avg_rpkm_name %in% names(cmb_d2))){
			cmb_d2[,avg_rpkm_name]=rowMeans(cmb_d2[,raw_rpkm_names,drop=F])
		}
	}

	for(i in 1:length(test_samples)){
		gchtp_name <- gchtp_names[i]
		l2gch_name <- l2gch_names[i]
		avg_rpkm_2names=paste(ifelse(IfCal_RPM,"rpm_","rpkm_"),c(test_samples[i],ref_samples[i]),sep="")
		higher_rpkm_vals=apply(cmb_d2[,avg_rpkm_2names],1,max)
		cmb_d2[,gchtp_name]="nc"
		if_expressed=higher_rpkm_vals>higher_avg_RPKM_cutoff & !is.na(higher_rpkm_vals)
		if_not_na= !is.na(cmb_d2[,l2gch_name]) & !is.na(cmb_d2[,use_pval_names[i]])
		if_up=if_not_na & if_expressed & cmb_d2[,use_pval_names[i]] > -log10(p_cut) & cmb_d2[,l2gch_name] > log2(foldchange_cut)
		if_dn=if_not_na & if_expressed & cmb_d2[,use_pval_names[i]] < log10(p_cut) & cmb_d2[,l2gch_name] < (-log2(foldchange_cut))
		cmb_d2[if_up,gchtp_name]="up"
		cmb_d2[if_dn,gchtp_name]="dn"
		cmb_d2[!if_not_na,gchtp_name]="na"
		print(table(cmb_d2[,gchtp_name])[gexch_types])
	}

	#output a table with all genes
	cmb_d2$max_SS=apply(abs(cmb_d2[use_pval_names[sel_gex_ids]]),1,max,na.rm=T)
	cmb_d2=cmb_d2[order(cmb_d2$max_SS,decreasing=T),]
	out_d=re_format_tb(cmb_d2, delete_headers=c("cds_len",isexp_names, isexpP_names, "ifexp1","ifreg1","ifGch","TRnum","coding","transc_start","transc_end","cds_start","cds_end","name2","cdsStartStat","cdsEndStat","MsP","MGFc"),
		Inf_headers=c(use_pval_names,l2gch_names,"max_SS"), front_headers=c(geneSym_name,"gene_id",'gene_Biotype','chromosome','strand'), back_headers=c(l2gch_names[sel_gex_ids],use_pval_names[sel_gex_ids],gchtp_names[sel_gex_ids], 'gene_desc')
	)
	out_format_f= paste(output_root, out_prefix,gene_model,".",pval_fun,".format_tb/",formatGeneTbName,"allGene.txt" ,sep="");
	mkdir_if_not_exist(out_format_f)
	write.table(out_d, file=out_format_f, col.names=T,row.names=F, sep="\t", quote=F)
	#write Excel file:
	if("openxlsx" %in% installed.packages()){
		library("openxlsx")
		out_xlsx_f= paste(output_root, out_prefix,gene_model,".",pval_fun,".format_tb/",formatGeneTbName,"allGene.","regulated.P",p_cut,".Ch",foldchange_cut,".Exp",higher_avg_RPKM_cutoff,".xlsx" ,sep="");
		print(paste("write",out_xlsx_f))
		wb <- createWorkbook()
		addWorksheet(wb, sheetName = "AllGeneExp")
		freezePane(wb, sheet = 1, firstRow = TRUE, firstCol = TRUE)
		writeDataTable(wb, sheet = 1, x = out_d, colNames = T, rowNames = F,  firstColumn=T)
		LeftBorderStyle <- createStyle( border = "left")
		leftBorderCols=c(min(grep("^num_",names(out_d))), match(c( rpkm_names[1], gchtp_names[1], l2gch_names[1], use_pval_names[1],"gene_desc") ,names(out_d)) )
		for(leftBorderCol in leftBorderCols){
			addStyle(wb, 1, style = LeftBorderStyle, rows = 1:(nrow(out_d)+1), cols=leftBorderCol, gridExpand = T)
		}
		headerStyle <- createStyle(textDecoration = "bold", halign = "left", valign = "bottom", border = "Bottom", textRotation=90)
		addStyle(wb, 1, style = headerStyle, rows = 1, cols=1:ncol(out_d), gridExpand = F, stack=T)
		setColWidths(wb, sheet = 1, cols =1:ncol(out_d), widths = 5)
		shortCol_headers=intersect(c("strand", "exon_num", gchtp_names), names(out_d))
		setColWidths(wb, sheet = 1, cols = match(shortCol_headers, names(out_d) ), widths = 3)
		setColWidths(wb, sheet = 1, cols = match("gene_desc", names(out_d) ), widths = 30)
		setCol_headers=intersect(c("gene_symbol", "gene_Biotype", rnum_idname), names(out_d))
		setColWidths(wb, sheet = 1, cols = match(setCol_headers, names(out_d) ), widths = 8)
		val_style=createStyle( numFmt = "0.0") #decimal point of values
		for(val_col in match(c(rpkm_names,l2gch_names,use_pval_names), names(out_d) )){
			addStyle(wb, 1, style = val_style, rows = 1+1:(nrow(out_d)), cols=val_col, gridExpand = T, stack=T)
		}
		upStyle <- createStyle( bgFill = "pink")
		dnStyle <- createStyle( bgFill = "royalblue") 
		conditionalFormatting(wb, 1, cols=match(gchtp_names,names(out_d)), rows=1+1:(nrow(out_d)), type = "contains", rule="up", style = upStyle)
		conditionalFormatting(wb, 1, cols=match(gchtp_names,names(out_d)), rows=1+1:(nrow(out_d)), type = "contains", rule="dn", style = dnStyle)
		conditionalFormatting(wb, 1, cols=match(l2gch_names,names(out_d)), rows=1+1:(nrow(out_d)), type = "colourScale", rule=c(-3,0,3), style = c("royalblue","white","red"))
		saveWorkbook(wb, out_xlsx_f, overwrite = TRUE)
		#write.xlsx(out_d, file = out_xlsx_f, asTable = TRUE)
	}

	#convert and write a long table (convert all data columns with sample pairs to rows)
	varheader_patts=c("^L2FC_","^adjSLog10P_","^GexType_") #header pattern for sample comparison variables
	names(varheader_patts)=c("GexL2FC","adjSLog10P","GexType")
	common_headers=intersect(c("gene_symbol","gene_id","gene_Biotype","gene_desc"),names(out_d))
	longData_i=out_d[rep(1:nrow(out_d),each=length(samplePairs)), common_headers, drop=F]
	longData_i$Project=study_name
	longData_i$TestSample=test_samples
	longData_i$CtlSample=ref_samples
	avg_exp_l=list(AvgTest_rpkm=test_samples, AvgCtl_rpkm=ref_samples)
	for(avg_exp_name in names(avg_exp_l)){
		avg_exps_m=sapply(avg_exp_l[[avg_exp_name]], function(sample1){
			#if(sample1 %in% names(comb_sample_l) ){
			#	raw_rpkm_names=paste(ifelse(IfCal_RPM,"rpm_","rpkm_"),comb_sample_l[[sample1]],sep="")
			#}else{
			#	raw_rpkm_names=paste(ifelse(IfCal_RPM,"rpm_","rpkm_"),sample1,sep="")
			#}
			#round(rowMeans(out_d[,raw_rpkm_names,drop=F],na.rm=T),2)
			raw_rpkm_names=paste(ifelse(IfCal_RPM,"rpm_","rpkm_"),sample1,sep="")
			round(out_d[,raw_rpkm_names],2)
		} )		
		longData_i[,avg_exp_name]=unlist(c(t(avg_exps_m)))
	}
	for(j in 1:length(varheader_patts)){
	    var_headers=grep(varheader_patts[j], names(out_d), value=T)
	    if(length(var_headers)>=length(samplePairs)){
	        var_headers2=paste0(sub("^\\^","",varheader_patts[j]),samplePairs)
	        longData_i[,names(varheader_patts)[j]]=c(unlist(t(out_d[,var_headers2])))
	    }
	}
	numeric_headers=intersect(c("GexL2FC","adjSLog10P"), names(longData_i))
	longData_i[,numeric_headers]=round(longData_i[,numeric_headers],2)
	longData_i=re_format_tb(longData_i, back_headers=c("gene_desc"))
	out_longTb_f= paste(output_root, out_prefix,gene_model,".",pval_fun,".format_tb/",study_name,".",formatGeneTbName,"allGene.longTable.txt" ,sep="");
	write.table(longData_i, file=out_longTb_f, col.names=T,row.names=F, sep="\t", quote=F)



	#output a table with regulated genes only
	ifreg <- rowSums( !is.na(cmb_d2[gchtp_names[sel_gex_ids]]) & (cmb_d2[gchtp_names[sel_gex_ids]]=="up" | cmb_d2[gchtp_names[sel_gex_ids]]=="dn") 
		& abs(cmb_d2[l2gch_names[sel_gex_ids]])> log2(foldchange_cut)  ) >= 1 #gene are regulated in >=XX of the samples
	if_visulizable <- rowSums( abs(cmb_d2[l2gch_names[sel_gex_ids]])==Inf | is.na(cmb_d2[l2gch_names[sel_gex_ids]]) ) ==0
	sum(ifreg); sum(if_visulizable & ifreg)
	out_f_base=paste(output_root, out_prefix,gene_model,".",pval_fun,".format_tb/",formatGeneTbName,"regulated.P",p_cut,".Ch",foldchange_cut,".Exp",higher_avg_RPKM_cutoff ,sep="");
	
	out_format_f= paste(out_f_base,".txt" ,sep="");
	write.table(out_d[ifreg,], file=out_format_f, col.names=T,row.names=F, sep="\t", quote=F)

	out_stats_f=paste(out_f_base,".ReguNum.txt",sep="")
	out_stats_f2=paste(out_f_base,".ReguNum.crossTable.txt",sep="")
	write(paste("#Differential expression (DE) analysis:", " Gene model=", gene_model,", P-value using ",pval_fun,", P<",p_cut,", fold change>",foldchange_cut,", higher one of the normalized expression>",higher_avg_RPKM_cutoff ,sep="") , file=out_stats_f)
	if(pval_fun !="none" ){
		tb=as.matrix(apply(out_d[gchtp_names],2,function(v){table( factor(v,levels=gexch_types) )})[gexch_types,])
		colnames(tb)=paste(test_samples,"/",ref_samples,sep="")
		write.table(cbind(rownames(tb),tb), file=out_stats_f, col.names=T, row.names=F, sep="\t", quote=F, append=T) 
		out_img_f=paste(out_plot_root, "GeneExpChange.Number.Barplot.",figure_out_format, sep="")
		do.call(figure_out_format, list(out_img_f, width=600*figure_size_fac, height=(300+ncol(tb)*60)*figure_size_fac ) )

		par(las=1,mar=c(5,10,4,2)+0.1)
		barplot(tb[1:2,], beside=T, horiz=T, col=c("red","royalblue"), main="Number of genes regulated", xlab="Number of genes",
			sub=paste0("(P<",p_cut,", Fold>",foldchange_cut,")"), names.arg=colnames(tb), legend.text=gexch_types[1:2], axisnames=T,
			space=c(0,1) )
		dev.off()
		##write regulation gene number in comment file
		if(length(gchtp_names)>1){
			library("fmsb")
			compare_types_among_vars(out_stats_f2, out_d, type_headers1=gchtp_names, Types=gexch_types)
		}
	}

	#output a table with rows as each comparison (sample pair), columns as numbers of up or down and sorted gene list DEGs of up or down
	#exclude small RNA genes
	if_exclude=!(out_d$gene_Biotype %in% c("mRNA","lncRNA","ncIsoform")); table(if_exclude)
	tb=as.matrix(apply(out_d[!if_exclude,gchtp_names],2,function(v){table( factor(v,levels=gexch_types) )})[gexch_types,])
	topGeneList_tb=data.frame(project=rep(study_name,length(samplePairs)), sample_pair=samplePairs)
	topGeneList_tb[paste0("Num_",gexch_types[1:2])]=t(tb[1:2,])
	tmp_tb=t(matrix(unlist(myApply(1:length(samplePairs), function(i){
		if_up=!if_exclude & out_d[,gchtp_names[i]] %in% gexch_types[1]
		if_dn=!if_exclude & out_d[,gchtp_names[i]] %in% gexch_types[2]
		return_v=c(0,0,"","")
		if(any(if_up | if_dn)){
			return_v[1:2]=apply(abs(out_d[if_up | if_dn, c(use_pval_names[i], l2gch_names[i] ), drop=F]), 2, max, na.rm=T)
		}
		if(any(if_up)){
			up_list=out_d[if_up, c("gene_symbol", l2gch_names[i],use_pval_names[i])]
			up_list=up_list[order(up_list[,2]*up_list[,3],decreasing=T),]
			if(nrow(up_list)>max_geneList_num){up_list=up_list[1:max_geneList_num,]}
			return_v[3]=paste(up_list[,1], collapse=" ")
		}
		if(any(if_dn)){
			dn_list=out_d[if_dn, c("gene_symbol", l2gch_names[i], use_pval_names[i])]
			dn_list=dn_list[order(dn_list[,2]*dn_list[,3],decreasing=T),]
			if(nrow(dn_list)>max_geneList_num){dn_list=dn_list[1:max_geneList_num,]}
			return_v[4]=paste(dn_list[,1], collapse=" ")
		}
		return(return_v)
	})), nrow=4))
	topGeneList_tb[,c("max_absSS","max_absL2FC",paste0("List_",gexch_types[1:2]))] = tmp_tb
	topGeneList_tb[c("max_absSS","max_absL2FC")]=as.numeric(unlist(topGeneList_tb[c("max_absSS","max_absL2FC")]))
	out_topGeneList_f=paste(out_f_base,".topGeneList.csv",sep="")
	write.csv(topGeneList_tb, file=out_topGeneList_f,  row.names=F)

	if(!is.na(filter_header)){
		if_select_plot=out_d[,filter_header] %in% filter_text & !is.na(out_d[,filter_header])
	}else{
		if_select_plot=rep(T, nrow(out_d))
	}
	#draw heatmap:
	library("gplots")
	plot_d=out_d[ifreg & if_visulizable & if_select_plot, ]
	sel_log2ratio_names=l2gch_names[sel_gex_ids]
	sel_sample_pair_names=paste(test_samples[sel_gex_ids], ref_samples[sel_gex_ids],sep="/")
	plot_d[,sel_sample_pair_names]=plot_d[,sel_log2ratio_names]
	if(length(sel_sample_pair_names)==1){sel_sample_pair_names=rep(sel_sample_pair_names,2)}
	out_img_f=paste(out_plot_root, "GeneExpChange.heatmap.",figure_out_format, sep="")
	mycolor=colorRampPalette(c("royalblue","white", "red"))(20);
	valbreaks=c(seq(-L2FC_scale_Max,L2FC_scale_Max, length.out=21))
	do.call(figure_out_format, list(out_img_f, width=800*figure_size_fac, height=800*figure_size_fac ) )
	
	labRow=plot_d[,geneSym_name]
	if(nrow(plot_d)>100){ labRow=rep("",nrow(plot_d)) }
	hv=heatmap.2(as.matrix(plot_d[,sel_sample_pair_names]), 
		ylab=paste("Reguated genes, #=",nrow(plot_d)) , margins=c(12,10), labRow=labRow, 
		col=mycolor
		, trace="none", scale="none", density.info="none"
		, breaks=valbreaks, 
		lwid=c(2,8), lhei=c(2,8),
		dendrogram="both", Colv=T,
		cexCol=1.5 ,
		key.title="Expression", key.xlab="Log2 ratio"
	)
	dev.off()
	
}
if(any(what2do %in% c("clust_rawSample_exp","rawSample_ScatterPlot","do_samplePair_plots")) ){
	out_format_f= paste(output_root, out_prefix,gene_model,".",pval_fun,".format_tb/",formatGeneTbName,"allGene.txt" ,sep="");
	cmb_d2=read.table(out_format_f, header=T, sep="\t", quote="", comment.char="", stringsAsFactors=F)

}

if(any(what2do %in% c("all","clust_rawSample_exp")) ){
	cor_method="spearman"
	nrow(cmb_d2)
	gex_d=cmb_d2
	#if("gene_Biotype" %in% names(cmb_d2)){
	#	gex_d=cmb_d2[!is.na(cmb_d2$gene_Biotype) & cmb_d2$gene_Biotype %in% c("mRNA","lncRNA","protein_coding"), ]; nrow(cmb_d2)
	#}
	sel_sample_ids=1:length(all_sample_name);
	sel_samples=all_sample_name[sel_sample_ids]
	sel_rnum_names=paste("num_",all_sample_name[sel_sample_ids],sep="")
	sel_rpkm_names=paste("rpkm_",all_sample_name[sel_sample_ids],sep="")
	if_expressed=rowMeans(gex_d[,sel_rpkm_names])>2; sum(if_expressed)
	if_regulated_matrix=!is.na(gex_d[gchtp_names]) & (gex_d[gchtp_names]=="up" | gex_d[gchtp_names]=="dn") & abs(gex_d[l2gch_names])> log2(foldchange_cut)
	if_regulated <- rowSums( if_regulated_matrix  ) >= 1 #gene are regulated in >=XX of the samples
	table(if_expressed,if_regulated)
	
	gex_d2=gex_d[if_expressed,]
	if(nrow(gex_d2)>10){
		out_image_f= paste(out_plot_root,"GeneExp.dendroplot.AllSamples.",figure_out_format ,sep="");
		mkdir_if_not_exist(out_image_f)
		hc=hclust( myDist_fun(t( gex_d2[,rpkm_names] ), method=cor_method ) , "ave" )
		hc$labels<-sub("rpkm_|rpm_","",rpkm_names); hc$labels
		do.call(figure_out_format, list(out_image_f, width=480*figure_size_fac, height=(300+length(rpkm_names) *10)*figure_size_fac ) )
		
		par(mar=c(4,3,3,10))
		plot(as.dendrogram(hc),horiz = TRUE, main=paste("Distance=1-correlation (",cor_method,"), #",nrow(gex_d2)))
		dev.off()
	}

	l2ratios=gex_d[l2gch_names]
	l2ratios[is.na(gex_d[l2gch_names]) | !if_regulated_matrix ]=0
	if(ncol(l2ratios)>1){
		l2ratios= change_plateau_val( apply(abs(l2ratios),1,max))
	}else{
		l2ratios=abs(unlist(l2ratios))
	}
	absL10Ps= change_plateau_val ( apply(abs(gex_d[use_pval_names]),1,max)  )
	if_top_fc_rows=(1:nrow(gex_d))[if_regulated] [ rank(-abs(l2ratios[if_regulated]), ties.method="first") %in% (1:50) ]
	if_top_P_rows=(1:nrow(gex_d))[if_regulated] [ rank(-absL10Ps[if_regulated], ties.method="first") %in% (1:50) ]
	if_topRegulated=(1:nrow(gex_d)) %in% if_top_fc_rows | (1:nrow(gex_d)) %in% if_top_P_rows 

	##do heat map of normalized expression for all genes or regulated genes:
	sel_list=list(regulatedGenes=if_expressed & if_regulated, allExpressedGenes=if_expressed, topRegulated=if_topRegulated)
	library("gplots")
	for (plotSet in names(sel_list)){
		out_img_f= paste(out_plot_root,"GeneExp.heatmap.",plotSet,".",figure_out_format ,sep="");
		do.call(figure_out_format, list(out_img_f, width=(900+length(sel_rpkm_names)*20)*figure_size_fac, height=1200*figure_size_fac ) )

		print(plotSet)
		plot_d=gex_d[sel_list[[plotSet]],  ]
		plot_d[,sel_samples]=log2( plot_d[,sel_rpkm_names]+0.01 )
		plot_d[,sel_samples]=plot_d[,sel_samples]- rowMeans(plot_d[,sel_samples])
		mycolor=colorRampPalette(c("royalblue","white", "red"))(50);
		valbreaks=c(seq(-2,2, length.out=51))
		labrow=plot_d[,geneSym_name]
		if(nrow(plot_d)>100){labrow=rep("",nrow(plot_d))}
		hv=heatmap.2(as.matrix(plot_d[,sel_samples]), 
		  ylab=paste(plotSet, ", #",nrow(plot_d)) , margins=c(10,10), labRow=labrow, 
		  col=mycolor
		  , trace="none", scale="none", density.info="none"
		  , breaks=valbreaks, 
		  dendrogram="both", Colv=T,
		  lwid=c(2,8), lhei=c(2,8),
		  key.title="Expression", key.xlab="mean-centered log2 RPKM"
		)
		dev.off()
	}
}

if(any(what2do %in% c("all","rawSample_ScatterPlot")) ){
	cor_method="spearman"
	nrow(cmb_d2)
	gex_d=cmb_d2
	sel_sample_ids=1:length(all_sample_name);
	sel_samples=all_sample_name[sel_sample_ids]
	sel_rpkm_names=paste("rpkm_",all_sample_name[sel_sample_ids],sep="")
	if_expressed=rowMeans(gex_d[,sel_rpkm_names])>2; sum(if_expressed)
	gex_d2=gex_d[if_expressed,]
	plot_d=gex_d2[,sel_rpkm_names]
	plot_d[,sel_samples]=log10(plot_d[,sel_rpkm_names]+0.01)

	panel.hist <- function(x, ...){
	    usr <- par("usr"); on.exit(par(usr))
	    par(usr = c(usr[1:2], 0, 1.5) )
	    h <- hist(x, plot = FALSE)
	    breaks <- h$breaks; nB <- length(breaks)
	    y <- h$counts; y <- y/max(y)
	    rect(breaks[-nB], 0, breaks[-1], y, col = "gray", ...)
	}

	panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...){
	    usr <- par("usr"); on.exit(par(usr))
	    par(usr = c(0, 1, 0, 1))
	    cor_r <- round(cor.test(x, y, method="spearman")$estimate,2)
	    txt <- paste0(prefix, cor_r)
	    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
	    mycolors=colorRampPalette(c("blue", "black", "red"))(20);
	    mycolor=as.character(cut(cor_r, breaks=seq(-1,1,length.out=21), labels=mycolors, include.lowest=T  ))
	    text(0.5, 0.5, txt, cex = cex.cor * abs(cor_r), col=mycolor )
	}

	panel.scatter <- function(x, y, ...){
	    cor_r <- round(cor.test(x, y, method="spearman")$estimate,2)
	    mycolors=colorRampPalette(c("blue", "black", "red"))(20);
	    mycolor=as.character(cut(cor_r, breaks=seq(-1,1,length.out=21), labels=mycolors, include.lowest=T  ))
	    points(x,y, col=mycolor, ...)
	    #abline(lm(y~x), lwd=2)
	}

	if(nrow(plot_d)>10){
		out_image_f= paste(out_plot_root,"GeneExp.Scatter-plot.AllSamples.",figure_out_format ,sep="");
		mkdir_if_not_exist(out_image_f)
		do.call(figure_out_format, list(out_image_f, width=(200+length(sel_rpkm_names)*60)*figure_size_fac, height=(200+length(sel_rpkm_names) *60)*figure_size_fac ) )

		pairs(plot_d[,sel_samples],main=paste0("Scatter plot and Spearman correlation of gene expression (log10(RPKM+0.01)), n=",nrow(plot_d) ), pch=".",
			diag.panel = panel.hist, upper.panel = panel.cor, lower.panel = panel.scatter )
		dev.off()
	}
}

if(any(what2do %in% c("all","do_samplePair_plots")) ){ #do volcano plot, scatter plot, MA-plot, etc
	gex_d=cmb_d2
	##do volcano plot:
	for(i in 1:length(test_samples)){
		if_plot=!is.na(gex_d[,l2gch_names[i]]) & !is.na(gex_d[,use_pval_names[i]])
		gex_d2=gex_d[if_plot, ]
		gex_d2=gex_d2[order(abs(gex_d2[,use_pval_names[i]]), abs(gex_d2[,l2gch_names[i]]), decreasing=T ),]
		out_img_f=paste(out_plot_root, "GeneExpChange.volcano.",test_samples[i],".vs.", ref_samples[i],".",figure_out_format ,sep="")
		sample_pair=paste(test_samples[i],".vs.", ref_samples[i], sep="")
		gchtp_name <- gchtp_names[i]
		l2ratios= change_plateau_val( gex_d2[,l2gch_names[i]] )
		absL10Ps= change_plateau_val ( abs(gex_d2[,use_pval_names[i]] ) )
		ymax=max(max(absL10Ps)*1.05, 5)
		xmin=min(min(l2ratios)*1.05, -3)
		xmax=max(max(l2ratios)*1.05, 3)
		if_up=!is.na(gex_d2[,gchtp_names[i]]) & gex_d2[,gchtp_names[i]]=="up"
		if_dn=!is.na(gex_d2[,gchtp_names[i]]) & gex_d2[,gchtp_names[i]]=="dn"
		if_regu=if_up | if_dn
		if_top_fc_rows=(1:nrow(gex_d2))[if_regu] [ rank(-abs(l2ratios[if_regu]), ties.method="first") %in% (1:60) ]
		if_top_P_rows=(1:nrow(gex_d2))[if_regu] [ rank(-absL10Ps[if_regu], ties.method="first") %in% (1:60) ]
		if_label=(1:nrow(gex_d2)) %in% if_top_fc_rows | (1:nrow(gex_d2)) %in% if_top_P_rows 
		do.call(figure_out_format, list(out_img_f, width=1500*figure_size_fac, height=800*figure_size_fac ) )
		layout(matrix(1:2, byrow=T, ncol=2))
		draw_scatter_plot(l2ratios,absL10Ps, xlim=c(xmin,xmax), ylim=c(0, ymax), xlab=paste("Log2 ratio",sample_pair), ylab="Absolute -Log10 (P-value)",
			title=paste("Volcano plot:",sample_pair),lm_line=NA, cor_method=NA,  
			if_highlight_v=if_up, if_highlight_v2=if_dn, highlight_col="red", highlight_col2="blue",
			verline=c(-1,1)*log2(foldchange_cut), hline=abs(-log10(p_cut)), pch2=1 )
		draw_scatter_plot(l2ratios,absL10Ps, xlim=c(xmin,xmax), ylim=c(0, ymax), xlab=paste("Log2 ratio",sample_pair), ylab="Absolute -Log10 (P-value)",
			title=paste("Volcano plot:",sample_pair),lm_line=NA, cor_method=NA,  
			if_highlight_v=if_up, if_highlight_v2=if_dn, highlight_col="pink", highlight_col2="lightblue",
			verline=c(-1,1)*log2(foldchange_cut), hline=abs(-log10(p_cut)), pch2=20 )
		if(any(if_label)){
			text(l2ratios[if_label],absL10Ps[if_label], labels=gex_d2[if_label, geneSym_name], cex=0.8)
		}
		dev.off()
	}

	##do MA-plot:
	for(i in 1:length(test_samples)){
		if_plot=!is.na(gex_d[,l2gch_names[i]]) & !is.na(gex_d[,use_pval_names[i]])
		gex_d2=gex_d[if_plot, ]
		gex_d2=gex_d2[order(abs(gex_d2[,use_pval_names[i]]), abs(gex_d2[,l2gch_names[i]]), decreasing=T ),]
		out_img_f=paste(out_plot_root, "GeneExpChange.MA-plot.",test_samples[i],".vs.", ref_samples[i],".",figure_out_format ,sep="")
		sample_pair=paste(test_samples[i],".vs.", ref_samples[i], sep="")
		gchtp_name <- gchtp_names[i]
		l2ratios= change_plateau_val( gex_d2[,l2gch_names[i]] )
		
		raw_samples=unlist(sapply(c(test_samples[i],ref_samples[i]), function(sample1){
			if(sample1 %in% names(comb_sample_l) ){
				comb_sample_l[[sample1]]
			}else{sample1}
		} ) )
		raw_rpkm_names=paste(ifelse(IfCal_RPM,"rpm_","rpkm_"),raw_samples,sep="")
		avg_exps=log10(rowMeans(gex_d2[,raw_rpkm_names],na.rm=T))

		ymin=min(min(l2ratios)*1.05, -3)
		ymax=max(max(l2ratios)*1.05, 3)
		if_up=!is.na(gex_d2[,gchtp_names[i]]) & gex_d2[,gchtp_names[i]]=="up"
		if_dn=!is.na(gex_d2[,gchtp_names[i]]) & gex_d2[,gchtp_names[i]]=="dn"
		if_regu=if_up | if_dn
		if_top_fc_rows=(1:nrow(gex_d2))[if_regu] [ rank(-abs(l2ratios[if_regu]), ties.method="first") %in% (1:60) ]
		if_top_P_rows=(1:nrow(gex_d2))[if_regu] [ rank(-absL10Ps[if_regu], ties.method="first") %in% (1:60) ]
		if_label=(1:nrow(gex_d2)) %in% if_top_fc_rows | (1:nrow(gex_d2)) %in% if_top_P_rows 
		do.call(figure_out_format, list(out_img_f, width=1200*figure_size_fac, height=600*figure_size_fac ) )
		layout(matrix(1:2, byrow=T, ncol=2))
		draw_scatter_plot(avg_exps,l2ratios,  ylim=c(ymin, ymax), xlab=paste("Average expression (Log10 ",ifelse(IfCal_RPM,"rpm","rpkm"),")",sep=""), 
			ylab=paste("Log2 ratio of expression (",sample_pair,")", sep=""),
			title=paste("MA plot:",sample_pair),lm_line=NA, cor_method=NA,  
			if_highlight_v=if_up, if_highlight_v2=if_dn, highlight_col="red", highlight_col2="blue",
			hline=0, verline=NA, pch2=1 )
		draw_scatter_plot(avg_exps,l2ratios,  ylim=c(ymin, ymax), xlab=paste("Average expression (Log10 ",ifelse(IfCal_RPM,"rpm","rpkm"),")",sep=""), 
			ylab=paste("Log2 ratio of expression (",sample_pair,")", sep=""),
			title=paste("MA plot:",sample_pair),lm_line=NA, cor_method=NA,  
			if_highlight_v=if_up, if_highlight_v2=if_dn, highlight_col="pink", highlight_col2="lightblue",
			hline=0, verline=NA, pch2=20 )
		if(any(if_label)){
			text(avg_exps[if_label],l2ratios[if_label], labels=gex_d2[if_label, geneSym_name], cex=0.8)
		}
		dev.off()
	}

	##do scatter plot:
	for(i in 1:length(test_samples)){
		if_plot=!is.na(gex_d[,l2gch_names[i]]) & !is.na(gex_d[,use_pval_names[i]])
		gex_d2=gex_d[if_plot, ]
		gex_d2=gex_d2[order(abs(gex_d2[,use_pval_names[i]]), abs(gex_d2[,l2gch_names[i]]), decreasing=T ),]
		out_img_f=paste(out_plot_root, "GeneExpChange.Scatter-plot.",test_samples[i],".vs.", ref_samples[i],".",figure_out_format ,sep="")
		sample_pair=paste(test_samples[i],".vs.", ref_samples[i], sep="")
		gchtp_name <- gchtp_names[i]
		l2ratios= change_plateau_val( gex_d2[,l2gch_names[i]] )
		
		raw_samples=sapply(c(test_samples[i],ref_samples[i]), function(sample1){
			if(sample1 %in% names(comb_sample_l) ){
				comb_sample_l[[sample1]]
			}else{sample1}
		}, simplify=F )
		print(raw_samples)
		test_rpkm_names=paste(ifelse(IfCal_RPM,"rpm_","rpkm_"),raw_samples[[1]],sep="")
		ref_rpkm_names=paste(ifelse(IfCal_RPM,"rpm_","rpkm_"),raw_samples[[2]],sep="")
		gex_d2[test_rpkm_names][gex_d2[test_rpkm_names]==0]=NA
		gex_d2[ref_rpkm_names][gex_d2[ref_rpkm_names]==0]=NA
		test_avg_exp=rowMeans(log2(gex_d2[test_rpkm_names]), na.rm=T)
		ref_avg_exp=rowMeans(log2(gex_d2[ref_rpkm_names]), na.rm=T)
		print(c(test_rpkm_names, ref_rpkm_names))
		l2ratios=test_avg_exp-ref_avg_exp

		if_up=!is.na(gex_d2[,gchtp_names[i]]) & gex_d2[,gchtp_names[i]]=="up" & l2ratios>log2(foldchange_cut)
		if_dn=!is.na(gex_d2[,gchtp_names[i]]) & gex_d2[,gchtp_names[i]]=="dn" & l2ratios< -log2(foldchange_cut)
		if_regu=if_up | if_dn
		if_top_fc_rows=(1:nrow(gex_d2))[if_regu] [ rank(-abs(l2ratios[if_regu]), ties.method="first") %in% (1:60) ]
		if_top_P_rows=(1:nrow(gex_d2))[if_regu] [ rank(-absL10Ps[if_regu], ties.method="first") %in% (1:60) ]
		if_label=(1:nrow(gex_d2)) %in% if_top_fc_rows | (1:nrow(gex_d2)) %in% if_top_P_rows 
		do.call(figure_out_format, list(out_img_f, width=1200*figure_size_fac, height=600*figure_size_fac ) )
		layout(matrix(1:2, byrow=T, ncol=2))
		draw_scatter_plot(ref_avg_exp,test_avg_exp,   ylab=paste("Average expression (Log2 ",ifelse(IfCal_RPM,"rpm","rpkm")," ",test_samples[i],")",sep=""), 
			xlab=paste("Average expression (Log2 ",ifelse(IfCal_RPM,"rpm","rpkm"),ref_samples[i],")",sep=""),
			title=paste("Scatter plot:",sample_pair),lm_line=NA, cor_method=NA,  
			if_highlight_v=if_up, if_highlight_v2=if_dn, highlight_col="red", highlight_col2="blue",
			hline=NA, verline=NA, pch2=1 )
		draw_scatter_plot(ref_avg_exp,test_avg_exp,   ylab=paste("Average expression (Log2 ",ifelse(IfCal_RPM,"rpm","rpkm")," ",test_samples[i],")",sep=""), 
			xlab=paste("Average expression (Log2 ",ifelse(IfCal_RPM,"rpm","rpkm"),ref_samples[i],")",sep=""),
			title=paste("Scatter plot:",sample_pair),lm_line=NA, cor_method=NA,  
			if_highlight_v=if_up, if_highlight_v2=if_dn, highlight_col="pink", highlight_col2="lightblue",
			hline=NA, verline=NA, pch2=20 )
		if(any(if_label)){
			text(ref_avg_exp[if_label],test_avg_exp[if_label], labels=gex_d2[if_label, geneSym_name], cex=0.8)
		}
		dev.off()
	}

}


if(any(what2do %in% c("calCorWithPhenotype")) ){
	print(paste("read",out_file))
	cmb_d2=read.table(out_file, header=T, sep="\t", quote="", comment.char="", stringsAsFactors=F)
	names(phenotypeVals)=phenotypeSamples
	use_rpkm_names=rpkm_names[match(phenotypeSamples,all_sample_name)]
	if_cal_cor=rowSums(!is.na(cmb_d2[use_rpkm_names]) & cmb_d2[use_rpkm_names]>0)>2; table(if_cal_cor)
	tmp_d=cmb_d2[if_cal_cor,use_rpkm_names]
	cor_names=paste("cor.",c("pea","spe"),".", phenotypeName,sep="")
	cmb_d2[,cor_names]=NA

	
	tmp_tb=unlist(myApply(1:nrow(tmp_d), function(row_i){
		gex_v=unlist(tmp_d[row_i,])
		cor_pea=cor.test(gex_v, phenotypeVals,method="pearson")
		cor_spe=cor.test(gex_v,phenotypeVals,method="spearman")
		c(cor_pea$estimate,cor_spe$estimate)
	}))
	cmb_d2[if_cal_cor,cor_names]= t(matrix(tmp_tb,nrow=2))
	rm(tmp_tb, tmp_d)
	write.table(cmb_d2, file=out_file, col.names=T, row.names=F, sep="\t", quote=F)
}
