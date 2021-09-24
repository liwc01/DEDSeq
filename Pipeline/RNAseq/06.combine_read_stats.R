
args<- commandArgs(TRUE)
print(args)

args_v<-NA
if(length(args)>1){
	args_v <- args[c(FALSE, TRUE)]
	names(args_v) <- sub("\\-+","", args[c(TRUE, FALSE)])
	print(args_v)
}
source("../SharedCodes/Rfunc.inc.R")
study_name<- ifelse(is.na(args_v["study_name"]), "Project name", args_v["study_name"]) ##project name
geno<- ifelse(is.na(args_v["geno"]), "genome", args_v["geno"]) ##genome version
project_root<- ifelse(is.na(args_v["project_root"]), NA, args_v["project_root"]) ##project name
mapping_data_root<- ifelse(is.na(args_v["mapping_data_root"]), project_root, args_v["mapping_data_root"]) 
if_PE<-ifelse(is.na(args_v["if_PE"]), T, args_v["if_PE"] %in% c("1","T")) 
if_do_splicing_analysis<-ifelse(is.na(args_v["if_do_splicing_analysis"]), F, args_v["if_do_splicing_analysis"] %in% c("1","T")) 
if_stranded<-ifelse(is.na(args_v["if_stranded"]), F, args_v["if_stranded"] %in% c("1","T")) 
report_out_dir=ifelse(is.na(args_v["report_out_dir"]), paste(project_root,"/Report/",sep=""), args_v["report_out_dir"]) 
p_cut <- ifelse(is.na(args_v["p_cut"]), 0.01, as.numeric(args_v["p_cut"])); 
foldchange_cut <- ifelse(is.na(args_v["foldchange_cut"]), 1.5, as.numeric(args_v["foldchange_cut"])); foldchange_cut
higher_avg_RPKM_cutoff <- ifelse(is.na(args_v["higher_avg_RPKM_cutoff"]), 1, as.numeric(args_v["higher_avg_RPKM_cutoff"])); higher_avg_RPKM_cutoff

gene_model<- ifelse(is.na(args_v["gene_model"]), "refseqcds", args_v["gene_model"]); gene_model ##club refseq enscds 
pval_fun<- ifelse(is.na(args_v["pval_fun"]), "DESeq2", args_v["pval_fun"]); #fisher chisq none DESeq2
out_prefix<- ifelse(is.na(args_v["out_prefix"]), "", args_v["out_prefix"])
test_samples<- unlist(strsplit(as.character(args_v["test_samples"])," ")) #must set a value 
ref_samples<- unlist(strsplit(as.character(args_v["ref_samples"])," ")) #must set a value 


out_log_f=paste(report_out_dir,"report.log",sep="")
out_R_inc_f=paste(report_out_dir,"report.inc.R",sep="")
out_data_img_f=paste(report_out_dir,"report.R.img",sep="")
if(!is.na(args_v["all_sample_name"])){ 
	all_sample_name<- unlist(strsplit(as.character(args_v["all_sample_name"])," ")) 
}

input_readnum_f<- ifelse(is.na(args_v["input_readnum_f"]), paste(mapping_data_root,"/STAR_map/Input.read.num.log.out",sep=""), args_v["input_readnum_f"]) 

SEreadmap_stats_f=paste(project_root,"/01.ReadTable/readnum.tbl",sep="")
PEreadmap_stats_f=paste(project_root,"/02.PE.ReadTable/PEread.stats.txt",sep="")
SE_2gene_stats_f=paste(project_root,"/04.GeneReadNum/",gene_model,"_SEReadTypeCnt.tbl",sep="")
PE_2gene_stats_f=paste(project_root,"/04.GeneReadNum/",gene_model,"_PEReadTypeCnt.tbl",sep="")
out_Gex_dir=paste(project_root, "/05.Gene_DE/", out_prefix,gene_model,".",pval_fun,".format_tb/" ,sep="");
out_f_base=paste(out_Gex_dir,"regulated.P",p_cut,".Ch",foldchange_cut,".Exp",higher_avg_RPKM_cutoff, sep="");
gene_DE_stats_f=paste(out_f_base,".ReguNum.txt",sep="")
gene_DE_stats_f2=paste(out_f_base,".ReguNum.crossTable.txt",sep="")
out_xlsx_f=paste(out_Gex_dir,"allGene.","regulated.P",p_cut,".Ch",foldchange_cut,".Exp",higher_avg_RPKM_cutoff,".xlsx",sep="") #
Gex_plot_folder=paste(out_Gex_dir,"GexPlot/",sep="");

# Spl_DEXSeq_stats_f=paste(project_root,"/RNAseq_Splicing/05.DEXSeq/geneBlock/formatted_tb/reguType.Padj0.001.Ch1.tbl.stats.tbl",sep="")
# Spl_PSI_allExon_stats_f=paste(project_root,"/RNAseq_Splicing/06.exon_PSI/all_exon/formatted_tb/exons.delta_PSI10_Pfisher0.001.stats.tbl",sep="")
# Spl_PSI_AnoExon_stats_f=paste(project_root,"/RNAseq_Splicing/06.exon_PSI/annotated_exon/formatted_tb/exons.delta_PSI10_Pfisher0.001.stats.tbl",sep="")


readStats_f_list=list(
	input_readnum_f=input_readnum_f,
	SEreadmap_stats_f=SEreadmap_stats_f,
	PEreadmap_stats_f=PEreadmap_stats_f,
	SE_2gene_stats_f=SE_2gene_stats_f,
	PE_2gene_stats_f=PE_2gene_stats_f
)
if(!if_PE){
	readStats_f_list$PEreadmap_stats_f=""
	readStats_f_list$PE_2gene_stats_f=""
}
if(!if_do_splicing_analysis){
	readStats_f_list$Spl_DEXSeq_stats_f=""
	readStats_f_list$Spl_PSI_allExon_stats_f=""
	readStats_f_list$Spl_PSI_AnoExon_stats_f=""
}

GexStats_f_list=list(
	gene_DE_stats_f=gene_DE_stats_f,
	gene_DE_stats_f2=gene_DE_stats_f2,
	out_xlsx_f=out_xlsx_f,
	Gex_plot_folder=Gex_plot_folder
)
if(length(test_samples)==1){
	GexStats_f_list$gene_DE_stats_f2=""
}

mkdir_if_not_exist(out_data_img_f)
save.image(out_data_img_f)


#copy files to report root
copy_f_list=list(ReadStats=readStats_f_list, GexStats=GexStats_f_list)

write("", file=out_R_inc_f)
write("", file=out_log_f)

write(paste("study_name=\'",study_name,"\';",sep=""),file=out_R_inc_f,append=T)
write(paste("if_PE=",if_PE,";",sep=""),file=out_R_inc_f,append=T)
write(paste("if_do_splicing_analysis=",if_do_splicing_analysis,";",sep=""),file=out_R_inc_f,append=T)


for(copy_out_folder in names(copy_f_list)){
	if_file_exist_v=NULL
	CopyOut_dir=paste(report_out_dir, "/", copy_out_folder, "/", sep="")
	write(paste("#Read stats report output directory:",CopyOut_dir), file=out_log_f, append=T)
	file_list1=copy_f_list[[copy_out_folder]]
	mkdir_if_not_exist(CopyOut_dir)
	for(file_var_name in names(file_list1)){
		file_path=file_list1[[file_var_name]]
		if(grepl("/$",file_path)){ #copy a folder
			cmd=paste("cp ", file_path,"*  \'", CopyOut_dir, "\'", sep=""); print(cmd)
			system(cmd)
			if_file_exist_v=c(if_file_exist_v,T)
		}else{ #copy a file
			file_name=sub(".*/","",file_path)
			write(paste(file_var_name,"=\'",copy_out_folder,"/",file_name,"\';",sep=""), file=out_R_inc_f,append=T)
			if_file_exist=F
			if(file_path != ""){
				print(file_path)
				if(file.exists(file_path)){
					if_file_exist=T
					system(paste("cp \'", file_path,"\'  \'", CopyOut_dir, "\'", sep=""))
				}
			}
			if_file_exist_v=c(if_file_exist_v,if_file_exist)
		}
	}
	report_files_data=data.frame(file_var=names(file_list1), if_file_exist=if_file_exist_v, file_path=unlist(file_list1))
	write.table(report_files_data, file=out_log_f, append=T, col.names=T, row.names=F, sep="\t", quote=F)
}

system( paste("cp RNAseq_report.Rmd ",report_out_dir,  sep="") )
save.image(out_data_img_f)
