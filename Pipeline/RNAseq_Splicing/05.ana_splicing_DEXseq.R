##analyze RNA-seq junction reads by DEXSeq

#available what2do: all, save_parameters; add_ss_seq; cal_reguType_output; add_as_ano(add alternative splicing annotation); add_ss_supp;
# add_ss_ConsSco add_ss_liftover 
#2016/12/15 can deal with samples without replicates by not using DEXseq2
#2016/12/21 can deal with samples without replicates using Fisher's exact test and also calculate fold change of gene block 
#2018/05/20 add function to deal with ana_unit=junc5 or junc3 (compare one junction with all other junctions use the same 5'(junc5) or 3'(junc3) splice sites)

args<- commandArgs(TRUE)
print(args)

args_v<-NA
if(length(args)>1){
	args_v <- args[c(FALSE, TRUE)]
	names(args_v) <- sub("\\-+","", args[c(TRUE, FALSE)])
	print(args_v)
}

source("../SharedCodes/Rfunc.inc.R") 
source("05.inc.func.R")

gene_id_header="gene_symbol"
workers=4
ana_unit="ss" #junc or ss #unit of analysis in DEXSeq: can be junction (junc) ot splice site (ss)

study_name="huntington_monkey_mouse"
#samples=unlist(strsplit("N_Control_Rep1 N_Control_Rep2 N_Control_Rep3 N_Cpd_Rep1 N_Cpd_Rep2 N_Cpd_Rep3"," ")); samples
in_cbf=paste("14gene_Rnum/",study_name,"/refseqcds.combine.junc2gene.tbl",sep="")
sample_repl_l=list(N_Control=c("N_Control_Rep1","N_Control_Rep2","N_Control_Rep3"), N_Cpd=c("N_Cpd_Rep1","N_Cpd_Rep2","N_Cpd_Rep3"))
sample_compare_matrix=rbind(c("N_Cpd"),c("N_Control"))

#parameters from command line 
study_name<- ifelse(is.na(args_v["study_name"]), "", args_v["study_name"])
gene_id_header<- ifelse(is.na(args_v["gene_id_header"]), "gene_symbol", args_v["gene_id_header"])
transc_id_header<- ifelse(is.na(args_v["transc_id_header"]), "refseqid", args_v["transc_id_header"])
in_cbf<- ifelse(is.na(args_v["in_cbf"]), paste("01.comb_junc_map_info","/combine.junc2gene.tbl",sep=""), args_v["in_cbf"])

juncAno_f=ifelse(is.na(args_v["juncAno_f"]), NA, args_v["juncAno_f"])
juncAno_match_headers=c("Gblock_id","Gblock_id")
juncAno_update_headers=c("contig","strand",transc_id_header,"start_pos","end_pos","startPosType","startPosAno","endPosType","endPosAno","reg_len","region_ano")
if(!is.na(args_v["juncAno_match_headers"])){ juncAno_match_headers= unlist(strsplit(as.character(args_v["juncAno_match_headers"])," "))}
if(!is.na(args_v["juncAno_update_headers"])){ juncAno_update_headers= unlist(strsplit(as.character(args_v["juncAno_update_headers"])," "))}

as_ano_root=ifelse(is.na(args_v["as_ano_root"]), NA, args_v["as_ano_root"])
as_ano_type_l=list(se=c("exon","exon_i"), a5ss=c("exon_a5ss","exon_i"),  a3ss=c("exon_a3ss","exon_i"))
as_ano_update_headers=c("exon_i_inFrame","exon_i_stopCodon")
Gblock_anos=c("exon","exon_a5ss","exon_a3ss","exon_part","intron","intron_a5ss","intron_a3ss","intron_part")

geneAno_f=ifelse(is.na(args_v["geneAno_f"]), NA, args_v["geneAno_f"]) 
geneAno_match_headers=c(transc_id_header,transc_id_header)

if(!is.na(args_v["geneAno_match_headers"])){ geneAno_match_headers= unlist(strsplit(as.character(args_v["geneAno_match_headers"])," "))}


if(!is.na(args_v["sample_repl_l"])){ sample_repl_l<- string2list(args_v["sample_repl_l"]) }
sample_compare_matrix=rbind(unlist(strsplit(as.character(args_v["test_samples"])," ")), unlist(strsplit(as.character(args_v["ref_samples"])," ")))
workers <- ifelse(is.na(args_v["workers"]), 16, as.numeric(args_v["workers"])); 
ana_unit<- ifelse(is.na(args_v["ana_unit"]), "junc", args_v["ana_unit"])
what2do<- ifelse(is.na(args_v["what2do"]), "all", args_v["what2do"])
all_samples=unlist(strsplit(as.character(args_v["all_samples"])," "))
use_p_type=ifelse(is.na(args_v["use_p_type"]), "P", args_v["use_p_type"])  #P or Padj
P_cut=ifelse(is.na(args_v["P_cut"]), 0.05, as.numeric(args_v["P_cut"])); 
delta_PSI_cut=ifelse(is.na(args_v["delta_PSI_cut"]), 10, as.numeric(args_v["delta_PSI_cut"])); 
if_use_PSI_change=F
if(grepl("junc",ana_unit)){
	if_use_PSI_change=ifelse(is.na(args_v["if_use_PSI_change"]), T, (args_v["if_use_PSI_change"] %in% c("T","Yes","True","1","Y","TURE")) ); 
}
PSI_cal_rnumMin=5

Log2ratio_cut=ifelse(is.na(args_v["Log2ratio_cut"]), log2(1.2), as.numeric(args_v["Log2ratio_cut"])); 

geno<- ifelse(is.na(args_v["geno"]), "hg19", args_v["geno"])
anno_version<- ifelse(is.na(args_v["anno_version"]), "1511", args_v["anno_version"]) #used for add_ss_supp

#add second species data -what2do add_ss_liftover
second_spl_name="PTC906_mm"
second_spl_f="27ana_juncreads/huntington_monkey_mouse/geneBlock/combine_d.tbl" #used in add_ss_liftover
second_spl_update_names=c("reg_len","region_ano" ,"Padj_N_Cpd_N_Control","Log2R_N_Cpd_N_Control")

second_spl_name<- ifelse(is.na(args_v["second_spl_name"]), "", args_v["second_spl_name"])
second_spl_f<- ifelse(is.na(args_v["second_spl_f"]), "", args_v["second_spl_f"])
if(!is.na(args_v["second_spl_update_names"])){ second_spl_update_names= unlist(strsplit(as.character(args_v["second_spl_update_names"])," "))}
ss5liftover_f<- ifelse(is.na(args_v["ss5liftover_f"]), paste("29geneBlocks/",study_name,"/liftover.hg19ToMm9/AllJunc.tbl.juncpos5.bed.out",sep=""), args_v["ss5liftover_f"])
ss3liftover_f<- ifelse(is.na(args_v["ss3liftover_f"]), paste("29geneBlocks/",study_name,"/liftover.hg19ToMm9/AllJunc.tbl.juncpos3.bed.out",sep=""), args_v["ss3liftover_f"])


suppNum_names=paste("num_",all_samples,sep="")
PSI_names=paste("PSI_",all_samples,sep="")
sample_pairs=paste(sample_compare_matrix[1,],sample_compare_matrix[2,],sep="_")
delta_PSI_names=paste("deltaPSI_", sample_pairs,sep="")
pval_names=paste("P_", sample_pairs,sep="")
adjpval_names=paste("Padj_", sample_pairs,sep="")
ReguType_names=paste("ReguType_", sample_pairs,sep="")
ReguTypes=c("UP","DN","NC","na")
sel_comp_ids=1:ncol(sample_compare_matrix)
Log2Ratio_names=paste("Log2R_", sample_pairs,sep="")

out_root<- ifelse(is.na(args_v["out_root"]), paste("05.DEXSeq/",study_name,"/",sep=""), args_v["out_root"])

out_alltb_f=paste(out_root,"/",ana_unit,"/combine_d.tbl",sep="")
out_allimg_f=paste(out_root,"/",ana_unit,"/combine_d.img",sep="")
out_ParameterImg_f=paste(out_root,"/",ana_unit,"/parameters.img",sep="")


if(use_p_type=="P"){
	use_P_names=pval_names
}else{
	use_P_names=adjpval_names
}

ssSupp_headers=c("startSS_supp","endSS_supp")
ssSupp_types=c("Refseq","KnG_Ens","mRNA_EST","None")
if(grepl("junc",ana_unit)){
	ss_seq_names=c("ss5_seq","ss3_seq")
	ss_pos_names=c("juncpos5","juncpos3")
}else if(ana_unit=="ss"){
	ss_seq_names=c("ss_seq")
	ss_pos_names=c("ss_pos")
}else if(ana_unit=="geneBlock"){
	ss_seq_names=c("ss5_seq","ss3_seq")
	ss_pos_names=c("juncpos1","juncpos2")
	ssSupp_headers=c("startSS_supp","endSS_supp")
}

if(transc_id_header=="refseqid"){
	gene_info_headers=c("gene_symbol","gene_id","contig","strand","gene_desc")
	geneAno_update_headers=c("gene_symbol","gene_id","gene_desc")
}else if(transc_id_header=="ensid"){
	gene_info_headers=c("gene_symbol","gene_Biotype","contig","strand","gene_desc")
	geneAno_update_headers=c("gene_symbol","gene_desc")
}
if(!is.na(args_v["geneAno_update_headers"])){ geneAno_update_headers= unlist(strsplit(as.character(args_v["geneAno_update_headers"])," "))}

calP_method<- ifelse(is.na(args_v["calP_method"]), "DEXSeq", args_v["calP_method"]) #DEXSeq, Ttest (for PSI) or FET, None (do not do any test)
if(max(sapply(sample_repl_l,length))==1){
	print("Warning: DEXseq will not be used because no replicates are defined!")
	calP_method="FET" #Fisher's exact test
}


if(what2do %in% c("all","save_parameters") ){
	mkdir_if_not_exist(out_ParameterImg_f)
	save.image(out_ParameterImg_f)
}

define_parallel_fun(nCores=workers)

###RUN

if(what2do %in% c("all") ){
	if(calP_method %in% c("DEXSeq","DEXseq")){library("DEXSeq")}
	##load input read counts data
	in_cbd=read.table(in_cbf, header=T, sep="\t", quote="", comment.char="", stringsAsFactors=F)

	###0, add annotations:
	if(!is.na(juncAno_f)){
		print (paste("read",juncAno_f))
		juncAno_d=read.table(juncAno_f, header=T, sep="\t", quote="", comment.char="", stringsAsFactors=F)
		if(ana_unit=="geneBlock"){
			in_cbd[,juncAno_update_headers]=juncAno_d[match( in_cbd[,juncAno_match_headers[1]], juncAno_d[,juncAno_match_headers[2]] ), juncAno_update_headers]
		}else if(grepl("junc",ana_unit)){
			juncAno_d$site_id=paste(juncAno_d[,transc_id_header],juncAno_d$strand, juncAno_d$site_pos,sep=":")
			in_cbd_ss5_id=paste(in_cbd[,transc_id_header], in_cbd$strand, in_cbd$juncpos5, sep=":")
			in_cbd_ss3_id=paste(in_cbd[,transc_id_header], in_cbd$strand, in_cbd$juncpos3, sep=":")
			in_cbd[,c("ss5_Ano")]=juncAno_d[match(in_cbd_ss5_id, juncAno_d$site_id) ,"PosAno"]
			in_cbd[,c("ss3_Ano")]=juncAno_d[match(in_cbd_ss3_id, juncAno_d$site_id) ,"PosAno"]
		}
		
		rm("juncAno_d")
	}
	if(!is.na(geneAno_f)){
		print (paste("read",geneAno_f))
		geneAno_d=read.table(geneAno_f, header=T, sep="\t", quote="", comment.char="", stringsAsFactors=F)
		in_cbd[,geneAno_update_headers]=geneAno_d[match( in_cbd[,geneAno_match_headers[1]], geneAno_d[,geneAno_match_headers[2]] ), geneAno_update_headers]
		rm("geneAno_d")
	}

	###1, use DEXSeq to study splicing
	#citation("DEXSeq") #check if DEXSeq is installed in R
	out_gff_f=paste(out_root,"/",ana_unit,"/",ana_unit,".gff",sep="")
	gff_raw_f=paste(out_gff_f,".raw.tbl",sep="")
	out_counts_folder=paste(out_root,"/",ana_unit,"/counts/",sep="")
	out_splicing_folder=paste(out_root,"/",ana_unit,"/splicing/",sep="")
	mkdir_if_not_exist(out_splicing_folder)
	mkdir_if_not_exist(out_counts_folder)

	#1.1 create a GFF format file as input of DEXSeq:
	print(paste("create GFF file:",out_gff_f))
	if(grepl("junc",ana_unit)){
		in_cbd2=in_cbd; nrow(in_cbd2)
		in_cbd2=in_cbd2[order(in_cbd2[,gene_id_header], in_cbd2$juncpos5*ifelse(in_cbd2$strand=="-",-1,1), in_cbd2$juncpos3*ifelse(in_cbd2$strand=="-",-1,1)), ]
		in_cbd2$GeneUnit_id= unlist(tapply(1:nrow(in_cbd2), in_cbd2[,gene_id_header],function(v){ 	1:length(v) }))
		in_cbd2$unit_id=paste(in_cbd2[,gene_id_header],":", in_cbd2$GeneUnit_id,    sep="")
		table(duplicated(in_cbd2$unit_id))

		gff_d=data.frame(chr=in_cbd2$contig, source="Seq", feature="exonic_part", 
			start=ifelse(in_cbd2$strand=="-",in_cbd2$juncpos3, in_cbd2$juncpos5), 
			end=ifelse(in_cbd2$strand=="-",in_cbd2$juncpos5, in_cbd2$juncpos3), score =1, 
			strand=in_cbd2$strand, frame='.')
		gff_d$attribute=paste("gene_id \"",in_cbd2[,gene_id_header], "\"; transcripts \"",in_cbd2[,gene_id_header],"\"; exonic_part_number \"",in_cbd2$GeneUnit_id,  "\"", sep="") #column 9 of gff file
		gff_d$unit_id=in_cbd2$unit_id
		table(duplicated(gff_d$unit_id))
	}else if(ana_unit=="ss"){
		ss5_ids=paste(in_cbd[,gene_id_header],in_cbd$juncpos5,sep=":")
		ss3_ids=paste(in_cbd[,gene_id_header],in_cbd$juncpos3,sep=":")
		ss5_d=data.frame(ss_id=unique(ss5_ids))
		ss3_d=data.frame(ss_id=unique(ss3_ids))
		#add gene information
		update_headers=unique(c(gene_id_header,gene_info_headers))
		ss5_d[,update_headers]=in_cbd[match(ss5_d$ss_id, ss5_ids),update_headers]
		ss3_d[,update_headers]=in_cbd[match(ss3_d$ss_id, ss3_ids),update_headers]
		ss5_d$ss_pos=in_cbd[match(ss5_d$ss_id, ss5_ids),"juncpos5"]
		ss3_d$ss_pos=in_cbd[match(ss3_d$ss_id, ss3_ids),"juncpos3"]
		ss5_d$ss_type="ss5"
		ss3_d$ss_type="ss3"
		for(i in 1:length(all_samples)){
			print(paste("calculate splice site read counts for",all_samples[i]))
			ss5_total_rnums=tapply( in_cbd[,suppNum_names[i]], ss5_ids, sum, na.rm=T)
			ss3_total_rnums=tapply( in_cbd[,suppNum_names[i]], ss3_ids, sum, na.rm=T)
			ss5_d[,suppNum_names[i]]=ss5_total_rnums[ss5_d$ss_id]
			ss3_d[,suppNum_names[i]]=ss3_total_rnums[ss3_d$ss_id]
		}
		in_cbd2=rbind(ss5_d,ss3_d)
		in_cbd2=in_cbd2[ order(in_cbd2[,gene_id_header], in_cbd2$ss_pos*ifelse(in_cbd2$strand=="-",-1,1)), ]
		in_cbd2$GeneUnit_id= unlist(tapply(1:nrow(in_cbd2), in_cbd2[,gene_id_header],function(v){ 	1:length(v) }))
		in_cbd2$unit_id=paste(in_cbd2[,gene_id_header],":", in_cbd2$GeneUnit_id,    sep="")
		gff_d=data.frame(chr=in_cbd2$contig, source="Seq", feature="exonic_part", 
			start=in_cbd2$ss_pos, end=in_cbd2$ss_pos,  score =1, 
			strand=in_cbd2$strand, frame='.')
		gff_d$attribute=paste("gene_id \"",in_cbd2[,gene_id_header], "\"; transcripts \"",in_cbd2[,gene_id_header],"\"; exonic_part_number \"",in_cbd2$GeneUnit_id,  "\"", sep="") #column 9 of gff file
		gff_d$unit_id=in_cbd2$unit_id
	}else if(ana_unit=="geneBlock"){
		in_cbd2=in_cbd; nrow(in_cbd2)
		in_cbd2=in_cbd2[order(in_cbd2[,gene_id_header], in_cbd2$start_pos*ifelse(in_cbd2$strand=="-",-1,1), in_cbd2$end_pos*ifelse(in_cbd2$strand=="-",-1,1)), ]
		in_cbd2$GeneUnit_id= unlist(tapply(1:nrow(in_cbd2), in_cbd2[,gene_id_header],function(v){ 	1:length(v) }))
		in_cbd2$unit_id=paste(in_cbd2[,gene_id_header],":", in_cbd2$GeneUnit_id,    sep="")
		gff_d=data.frame(chr=in_cbd2$contig, source="Seq", feature="exonic_part", 
			start=ifelse(in_cbd2$strand=="-",in_cbd2$end_pos, in_cbd2$start_pos), 
			end=ifelse(in_cbd2$strand=="-",in_cbd2$start_pos, in_cbd2$end_pos),  
			score =1, 
			strand=in_cbd2$strand, frame='.')
		gff_d$attribute=paste("gene_id \"",in_cbd2[,gene_id_header], "\"; transcripts \"",in_cbd2[,gene_id_header],"\"; exonic_part_number \"",in_cbd2$GeneUnit_id,  "\"", sep="") #column 9 of gff file
	}
		
	write.table(gff_d[,1:9], file=out_gff_f, col.names=F, row.names=F, sep="\t", quote=F )
	write.table(gff_d, file=gff_raw_f, col.names=T, row.names=F, sep="\t", quote=F )
	
	rm("in_cbd") #to save memory
	rm("gff_d") #to save memory
	write.table(in_cbd2, file=out_alltb_f, col.names=T, row.names=F, sep="\t", quote=F)
	save.image(out_allimg_f)


	if(grepl("junc",ana_unit)){
		#1.2 calculate PSI (percent spliced in) based on number of read supporting one junction among all reads using either 5'ss or 3'ss of the same junction
		ss5_ids=paste(in_cbd2[,gene_id_header],in_cbd2$juncpos5,sep=":")
		ss3_ids=paste(in_cbd2[,gene_id_header],in_cbd2$juncpos3,sep=":")
		for(i in 1:length(all_samples)){
			print(paste("calculate",PSI_names[i]))
			ss5_total_rnums=tapply( in_cbd2[,suppNum_names[i]], ss5_ids, sum, na.rm=T)
			ss3_total_rnums=tapply( in_cbd2[,suppNum_names[i]], ss3_ids, sum, na.rm=T)
			
			if(ana_unit=="junc"){
				total_rnums=ss5_total_rnums[ss5_ids]+ss3_total_rnums[ss3_ids]-in_cbd2[,suppNum_names[i]]
			}else if(ana_unit=="junc5"){
				total_rnums=ss5_total_rnums[ss5_ids]
			}else if(ana_unit=="junc3"){
				total_rnums=ss3_total_rnums[ss3_ids]
			}
			in_cbd2[,PSI_names[i]]=round(in_cbd2[,suppNum_names[i]]/total_rnums*100, 2)
			in_cbd2[total_rnums<PSI_cal_rnumMin  ,PSI_names[i]]=NA
		}
		#calculate delta-PSI
		for(comp_id in sel_comp_ids){
			test_psi_names=paste("PSI_",sample_repl_l[[sample_compare_matrix[1,comp_id]]], sep="")
			ctrl_psi_names=paste("PSI_",sample_repl_l[[sample_compare_matrix[2,comp_id]]], sep="")
			in_cbd2[,delta_PSI_names[comp_id]]= rowMeans(in_cbd2[test_psi_names],na.rm=T)-rowMeans(in_cbd2[ctrl_psi_names],na.rm=T)
		}
	}


	
	if(calP_method %in% c("DEXSeq","DEXseq")){
		##1.3 output counts. table
		for(suppNum_name in suppNum_names){
			outtmp_d=data.frame(count=in_cbd2[,suppNum_name])
			rownames(outtmp_d)=in_cbd2$unit_id
			out_cnt_f=paste(out_counts_folder,suppNum_name,".txt",sep="")
			print(paste("create counts file:",out_cnt_f))
			write.table(outtmp_d, file=out_cnt_f, row.names=T, col.names=F, sep="\t", quote=F)
		}
		rm("outtmp_d")

		#1.4 build sampleTable (define test and control sample, replicates)
		for(comp_id in sel_comp_ids){
			samples_i=c(sample_compare_matrix[,comp_id])
			repNums=sapply(sample_repl_l[samples_i], length)
			sampleTable=data.frame(
				countFile=paste(out_counts_folder,"num_",unlist(sample_repl_l[samples_i]),".txt",sep=""), 
				condition=rep(c("test","control"),repNums) 
			)
			rownames(sampleTable)=make.unique(unlist(sample_repl_l[samples_i]))
			out_sampleTable_f=paste(out_splicing_folder,sample_pairs[comp_id],".sampleTable.txt",sep="")
			print(paste("create sampleTable file:",out_sampleTable_f))
			write.table(sampleTable, file=out_sampleTable_f, col.names=T, row.names=T, sep="\t", quote=F)
		}


		##1.5 compare samples using DEXSeq:
		for(comp_id in sel_comp_ids){ 
			out_sampleTable_f=paste(out_splicing_folder,sample_pairs[comp_id],".sampleTable.txt", sep="")
			out_splicing_res_f=paste(out_splicing_folder,sample_pairs[comp_id],".res.tbl", sep="")
			out_sizefac_res_f=paste(out_splicing_folder,sample_pairs[comp_id],".sizefac.tbl", sep="")
			out_dxd_Robject_f=paste(out_splicing_folder,sample_pairs[comp_id],".dxd.obj", sep="")
			
			sampleTable=read.table(out_sampleTable_f, header=T, sep="\t")
			dxd <- DEXSeqDataSetFromHTSeq( as.character(sampleTable$countFile), sampleTable, flattenedfile=out_gff_f ) #design= ~ sample + exon + condition:exon
			
			colData(dxd)
			head( counts(dxd), 5 )
			head( featureCounts(dxd), 5 )
			sampleAnnotation( dxd )
			dxd = estimateSizeFactors( dxd )
			
			sizefacs=sizeFactors(dxd)
			write.table(sampleAnnotation( dxd ), file=out_sizefac_res_f, col.names=T, row.names=F, sep="\t", quote=F)
			
			dxd = estimateDispersions( dxd, BPPARAM=MulticoreParam(workers=workers) ) #do not work for samples without replicates
			dxd = testForDEU( dxd, BPPARAM=MulticoreParam(workers=workers) )
			dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition", BPPARAM=MulticoreParam(workers=workers))

			dxr1 = DEXSeqResults( dxd )
			head(dxr1)
			table( dxr1$padj < 0.1 )
			table ( tapply( dxr1$padj < 0.1, dxr1$groupID, any ) ) #genes affected
			#dxr1_reg=dxr1[!is.na(dxr1$padj) & dxr1$padj<0.1, ]
			#output:
			print(paste("write",out_splicing_res_f))
			write.table(dxr1, file=out_splicing_res_f, col.names=T, row.names=F, sep="\t", quote=F)
			#write.table(dxr1_reg, file=paste(out_splicing_res_f,".reg.tbl",sep=""), col.names=T, row.names=F, sep="\t", quote=F)
			#save(dxd, file=out_dxd_Robject_f)
			#load(out_dxd_Robject_f)
			rm("dxd","dxr1")
		}


		##1.6 intergrate the results in DEXSeq to the original table (in_cbd2)
		unit_ids1=paste(in_cbd2[,gene_id_header],":E", in_cbd2$GeneUnit_id,    sep="")
		for(comp_id in sel_comp_ids){ 
			out_sampleTable_f=paste(out_splicing_folder,sample_pairs[comp_id],".sampleTable.txt", sep="")
			out_splicing_res_f=paste(out_splicing_folder,sample_pairs[comp_id],".res.tbl", sep="")
			dxr1=read.table(out_splicing_res_f, header=T, sep="\t", quote="", comment.char="", stringsAsFactors=F)
			dim(dxr1)
			unit_ids2=paste(dxr1$groupID, dxr1$featureID,sep=":")
			print(c(pval_names[comp_id], adjpval_names[comp_id]))
			in_cbd2[c(pval_names[comp_id], adjpval_names[comp_id])]=dxr1[match(unit_ids1,unit_ids2), c("pvalue","padj")]
			if("log2fold_control_test" %in% names(dxr1) ){
				dxr1$log2fold_test_control= -dxr1$log2fold_control_test
			}
			in_cbd2[,Log2Ratio_names[comp_id]]= round(dxr1[match(unit_ids1,unit_ids2), "log2fold_test_control"], 2)
		}
		rm("dxr1")
	}else if(calP_method=="FET"){ #use fisher's exact test to calculate P value, use log2 ratio of ratio of geneBlock read number vs. gene read number to calculate the fold change
		for(comp_id in sel_comp_ids){
			test_readnums=rowSums(in_cbd2[paste("num_",sample_repl_l[[sample_compare_matrix[1,comp_id]]], sep="")])
			ref_readnums=rowSums(in_cbd2[paste("num_", sample_repl_l[[sample_compare_matrix[2,comp_id]]], sep="")])
			test_gene_readnums=tapply(test_readnums, in_cbd2[,gene_id_header],sum,na.rm=T)
			ref_gene_readnums=tapply(ref_readnums, in_cbd2[,gene_id_header],sum,na.rm=T)
			readnum_tb=cbind(test_readnums,ref_readnums, test_gene_readnums[in_cbd2[,gene_id_header]]-test_readnums, ref_gene_readnums[in_cbd2[,gene_id_header]]-ref_readnums  )
			#fold change
			Log2FCs=log2(test_readnums+0.1) - log2(ref_readnums+0.1) - ( log2(test_gene_readnums[in_cbd2[,gene_id_header]]+0.1) - log2(ref_gene_readnums[in_cbd2[,gene_id_header]]+0.1) )
			in_cbd2[,Log2Ratio_names[comp_id]]=round(Log2FCs,2)
			#P-value based on Fisher's exact test
			in_cbd2[,pval_names[comp_id]] = unlist( myApply(1:nrow(readnum_tb), function(row_i){  fisher.test(matrix(readnum_tb[row_i,],2))$p.value } ) )
		}
		rm(readnum_tb,Log2FCs,test_readnums,ref_readnums,test_gene_readnums,ref_gene_readnums)
	}else if(calP_method=="Ttest" & if_use_PSI_change){
		for(comp_id in sel_comp_ids){
			test_psi_names=paste("PSI_",sample_repl_l[[sample_compare_matrix[1,comp_id]]], sep="")
			ctrl_psi_names=paste("PSI_",sample_repl_l[[sample_compare_matrix[2,comp_id]]], sep="")
			comp_nums1=apply(in_cbd2[,test_psi_names],1,function(v){ length(unique(v[!is.na(v)]) ) })
			comp_nums2=apply(in_cbd2[,ctrl_psi_names],1,function(v){ length(unique(v[!is.na(v)]) ) })
			if_do_ttest=comp_nums1>1 & comp_nums2>1; table(if_do_ttest)
			in_cbd2[,pval_names[comp_id]]=NA
			in_cbd2[if_do_ttest,pval_names[comp_id]]= unlist( myApply( (1:nrow(in_cbd2))[if_do_ttest], function(row_i){  
				t.test( in_cbd2[row_i,test_psi_names], in_cbd2[row_i,ctrl_psi_names] )$p.value
			} ) )
		}

	}
	if( !(calP_method=="None") ){
		write.table(in_cbd2, file=out_alltb_f, col.names=T, row.names=F, sep="\t", quote=F)
		save.image(out_allimg_f)
	}
}

if(what2do %in% c("add_ss_seq","cal_reguType_output","add_as_ano","add_ss_supp","add_ss_ConsSco","add_ss_liftover") ){ #these function need to reload in_cbd2
	print(paste("re-open",out_alltb_f))
	in_cbd2=read.table(out_alltb_f, header=T, sep="\t", quote="", comment.char="", stringsAsFactors=F)
}

if(what2do %in% c("all","add_ss_seq") ){
	##	add splice site sequence to the table (in_cbd2)
	if(grepl("junc",ana_unit)){
		in_cbd2=add_ss_site2data(input_d=in_cbd2)
	}else if(ana_unit=="ss"){
		ifss5=in_cbd2$ss_type=="ss5"
		ifss3=in_cbd2$ss_type=="ss3"
		in_cbd2$ss_seq=""
		in_cbd2$ss_seq[ifss5]=add_ss_site2data(input_d=in_cbd2[ifss5,], addSeqNames=c("ss5_seq"), ss5_header="ss_pos")$ss5_seq
		in_cbd2$ss_seq[ifss3]=add_ss_site2data(input_d=in_cbd2[ifss3,], addSeqNames=c("ss3_seq"), ss3_header="ss_pos")$ss3_seq
	}else if(ana_unit=="geneBlock"){
		strand_signs=ifelse(in_cbd2$strand=="-",-1,1)
		if_start_ss5=!is.na(in_cbd2$startPosType) & in_cbd2$startPosType=="ss5"
		if_start_ss3=!is.na(in_cbd2$startPosType) & in_cbd2$startPosType=="ss3"
		if_end_ss5=!is.na(in_cbd2$endPosType) & in_cbd2$endPosType=="ss5"
		if_end_ss3=!is.na(in_cbd2$endPosType) & in_cbd2$endPosType=="ss3"
		in_cbd2$startSS_seq=""; 
		in_cbd2$startSS_seq[if_start_ss5]=add_ss_site2data( input_d=in_cbd2[if_start_ss5,], addSeqNames=c("ss5_seq"), ss5_header="start_pos", ss5_pos_shift=-strand_signs[if_start_ss5] )$ss5_seq
		in_cbd2$startSS_seq[if_start_ss3]=add_ss_site2data(input_d=in_cbd2[if_start_ss3,], addSeqNames=c("ss3_seq"), ss3_header="start_pos")$ss3_seq
		in_cbd2$endSS_seq=""
		in_cbd2$endSS_seq[if_end_ss5]=add_ss_site2data(input_d=in_cbd2[if_end_ss5,], addSeqNames=c("ss5_seq"), ss5_header="end_pos")$ss5_seq
		in_cbd2$endSS_seq[if_end_ss3]=add_ss_site2data(input_d=in_cbd2[if_end_ss3,], addSeqNames=c("ss3_seq"), ss3_header="end_pos", ss3_pos_shift=strand_signs[if_end_ss3])$ss3_seq

	}
	##output in_cbd2
	write.table(in_cbd2, file=out_alltb_f, col.names=T, row.names=F, sep="\t", quote=F)
	save.image(out_allimg_f)
}

as_ano_type_l=list(se=c("exon","exon_i"), a5ss=c("exon_a5ss","exon_i"),  a3ss=c("exon_a3ss","exon_i"))
as_ano_update_headers=c("exon_i_inFrame","exon_i_stopCodon")

if(what2do %in% c("all","add_as_ano") ){
	if(!is.na(as_ano_root) & ana_unit=="geneBlock" ){
		ind_ids=paste(in_cbd2[,transc_id_header],":",format(in_cbd2$start_pos,scientific=F,trim=T), "-",format(in_cbd2$end_pos,scientific=F,trim=T), sep="")
		for(as_type in names(as_ano_type_l) ){
			as_ano_f=paste(as_ano_root,as_type,".txt", sep="")
			print(paste("open",as_ano_f))
			update_region_type=as_ano_type_l[[as_type]][1]
			update_match_header=as_ano_type_l[[as_type]][2]
			as_ano_d=read.table(as_ano_f, header=T, sep="\t", quote="", stringsAsFactors=F, comment.char="")
			if("allJunc_ReadNum" %in% names(as_ano_d)){
				as_ano_d=as_ano_d[order(as_ano_d$allJunc_ReadNum, decreasing=T),]
			}
			ano_d_ids=paste(as_ano_d[,transc_id_header],":", as_ano_d[,update_match_header],sep="")
			as_ano_update_headers2=intersect(as_ano_update_headers, names(as_ano_d))
			if(length(as_ano_update_headers2)>0){
				new_headers=setdiff(as_ano_update_headers2, names(in_cbd2))
				if(length(new_headers)>0){
					in_cbd2[new_headers]=NA
				}
				if_update=!is.na(in_cbd2$region_ano) & in_cbd2$region_ano %in% update_region_type
				in_cbd2[if_update,as_ano_update_headers2]=as_ano_d[match(ind_ids[if_update], ano_d_ids),as_ano_update_headers2]
			}
		}
		rm(as_ano_d,ano_d_ids)
		write.table(in_cbd2, file=out_alltb_f, col.names=T, row.names=F, sep="\t", quote=F)
		save.image(out_allimg_f)
	}

}


##add splice site support information (eg. Refseq, Ensembl, known gene, mRNA, est)
if(what2do %in% c("all","add_ss_supp") ){
	if(grepl("junc",ana_unit)){
		# in_cbd2=add_ss_site2data(input_d=in_cbd2)
	}else if(ana_unit=="ss"){
		# ifss5=in_cbd2$ss_type=="ss5"
		# ifss3=in_cbd2$ss_type=="ss3"
		# in_cbd2$ss_seq=""
		# in_cbd2$ss_seq[ifss5]=add_ss_site2data(input_d=in_cbd2[ifss5,], addSuppNames=c("ss5_supp"), ss5_header="ss_pos")$ss5_supp
		# in_cbd2$ss_seq[ifss3]=add_ss_site2data(input_d=in_cbd2[ifss3,], addSuppNames=c("ss3_supp"), ss3_header="ss_pos")$ss3_supp
	}else if(ana_unit=="geneBlock"){
		strand_signs=ifelse(in_cbd2$strand=="-",-1,1)
		if_start_ss5=!is.na(in_cbd2$startPosType) & in_cbd2$startPosType=="ss5"
		if_start_ss3=!is.na(in_cbd2$startPosType) & in_cbd2$startPosType=="ss3"
		if_end_ss5=!is.na(in_cbd2$endPosType) & in_cbd2$endPosType=="ss5"
		if_end_ss3=!is.na(in_cbd2$endPosType) & in_cbd2$endPosType=="ss3"
		in_cbd2$startSS_supp="na";
		in_cbd2$startSS_supp[if_start_ss5]=add_ss_supp2data( input_d=in_cbd2[if_start_ss5,], addSuppNames=c("ss5_supp"), ss5_header="start_pos", ss5_pos_shift=-strand_signs[if_start_ss5], anno_version=anno_version, geno=geno )$ss5_supp
		in_cbd2$startSS_supp[if_start_ss3]=add_ss_supp2data(input_d=in_cbd2[if_start_ss3,], addSuppNames=c("ss3_supp"), ss3_header="start_pos", anno_version=anno_version, geno=geno)$ss3_supp
		in_cbd2$endSS_supp="na"
		in_cbd2$endSS_supp[if_end_ss5]=add_ss_supp2data(input_d=in_cbd2[if_end_ss5,], addSuppNames=c("ss5_supp"), ss5_header="end_pos", anno_version=anno_version, geno=geno)$ss5_supp
		in_cbd2$endSS_supp[if_end_ss3]=add_ss_supp2data(input_d=in_cbd2[if_end_ss3,], addSuppNames=c("ss3_supp"), ss3_header="end_pos", ss3_pos_shift=strand_signs[if_end_ss3], anno_version=anno_version, geno=geno)$ss3_supp
		##output in_cbd2
		write.table(in_cbd2, file=out_alltb_f, col.names=T, row.names=F, sep="\t", quote=F)
		save.image(out_allimg_f)
	}
}

##add splice site conservation (PhastCons score) information 
if(what2do %in% c("add_ss_ConsSco") ){
	if(grepl("junc",ana_unit)){
		# in_cbd2=add_ss_site2data(input_d=in_cbd2)
	}else if(ana_unit=="ss"){
		# ifss5=in_cbd2$ss_type=="ss5"
		# ifss3=in_cbd2$ss_type=="ss3"
		# in_cbd2$ss_seq=""
		# in_cbd2$ss_seq[ifss5]=add_ss_site2data(input_d=in_cbd2[ifss5,], addConsNames=c("ss5_ConsSco"), ss5_header="ss_pos")$ss5_ConsSco
		# in_cbd2$ss_seq[ifss3]=add_ss_site2data(input_d=in_cbd2[ifss3,], addConsNames=c("ss3_ConsSco"), ss3_header="ss_pos")$ss3_ConsSco
	}else if(ana_unit=="geneBlock"){
		strand_signs=ifelse(in_cbd2$strand=="-",-1,1)
		if_start_ss5=!is.na(in_cbd2$startPosType) & in_cbd2$startPosType=="ss5"
		if_start_ss3=!is.na(in_cbd2$startPosType) & in_cbd2$startPosType=="ss3"
		if_end_ss5=!is.na(in_cbd2$endPosType) & in_cbd2$endPosType=="ss5"
		if_end_ss3=!is.na(in_cbd2$endPosType) & in_cbd2$endPosType=="ss3"
		in_cbd2$startSS_ConsSco="na";
		in_cbd2$startSS_ConsSco[if_start_ss5]=add_ss_ConsSco( input_d=in_cbd2[if_start_ss5,], addConsNames=c("ss5_ConsSco"), ss5_header="start_pos", ss5_pos_shift=-strand_signs[if_start_ss5] )$ss5_ConsSco
		in_cbd2$startSS_ConsSco[if_start_ss3]=add_ss_ConsSco(input_d=in_cbd2[if_start_ss3,], addConsNames=c("ss3_ConsSco"), ss3_header="start_pos")$ss3_ConsSco
		in_cbd2$endSS_ConsSco="na"
		in_cbd2$endSS_ConsSco[if_end_ss5]=add_ss_ConsSco(input_d=in_cbd2[if_end_ss5,], addConsNames=c("ss5_ConsSco"), ss5_header="end_pos")$ss5_ConsSco
		in_cbd2$endSS_ConsSco[if_end_ss3]=add_ss_ConsSco(input_d=in_cbd2[if_end_ss3,], addConsNames=c("ss3_ConsSco"), ss3_header="end_pos", ss3_pos_shift=strand_signs[if_end_ss3])$ss3_ConsSco
		##output in_cbd2
		write.table(in_cbd2, file=out_alltb_f, col.names=T, row.names=F, sep="\t", quote=F)
		save.image(out_allimg_f)
	}
}

##add regulation data based on splice site liftovered to another species
if(what2do %in% c("add_ss_liftover") ){
	if(ana_unit=="geneBlock"){
		in_cbd2=add_ss_liftover(in_cbd2)
		##output in_cbd2
		write.table(in_cbd2, file=out_alltb_f, col.names=T, row.names=F, sep="\t", quote=F)
		save.image(out_allimg_f)
	}
}


if(what2do %in% c("all","cal_reguType_output") & !(calP_method=="None") ){	
	#calculate minimal p value
	if(ana_unit=="ss"){
		in_cbd2$ss_seq_id=paste(in_cbd2$contig, in_cbd2$strand, in_cbd2$ss_pos, sep=":")  #for cis elements analysis
	}
	##calculate regulation type
	if(use_p_type=="P"){
		delete_P_names=adjpval_names
	}else{
		delete_P_names=pval_names
	}
	if(length(use_P_names)>1){
		in_cbd2$min_Padj=apply(in_cbd2[use_P_names],1,min,na.rm=T)
	}
	for(i in 1:ncol(sample_compare_matrix)){
		if(grepl("junc",ana_unit)){
			if(if_use_PSI_change){
				in_cbd2[,ReguType_names[i]]=cal_regu_type(data=in_cbd2, p_name=use_P_names[i], change_name=delta_PSI_names[i], p_cut=P_cut, change_cut=delta_PSI_cut )
			}else{
				print(paste(ReguType_names[i],":", use_P_names[i],"<",P_cut,"&", Log2Ratio_names[i], ">",Log2ratio_cut))
				in_cbd2[,ReguType_names[i]]=cal_regu_type(data=in_cbd2, p_name=use_P_names[i], change_name=Log2Ratio_names[i], p_cut=P_cut, change_cut=Log2ratio_cut )
			}
			
		}else if(ana_unit %in% c("ss","geneBlock") ){
			in_cbd2[,ReguType_names[i]]=cal_regu_type(data=in_cbd2, p_name=use_P_names[i], change_name=Log2Ratio_names[i], p_cut=P_cut, change_cut=Log2ratio_cut )
		}
	}
	#add ss5 and ss3 position (exonic +1 based position for all region types)
	if(ana_unit=="geneBlock"){
		if_start_ss5=!is.na(in_cbd2$startPosType) & in_cbd2$startPosType=="ss5"
		if_start_ss3=!is.na(in_cbd2$startPosType) & in_cbd2$startPosType=="ss3"
		if_end_ss5=!is.na(in_cbd2$endPosType) & in_cbd2$endPosType=="ss5"
		if_end_ss3=!is.na(in_cbd2$endPosType) & in_cbd2$endPosType=="ss3"
		strand_signs=ifelse(in_cbd2$strand=="-",-1,1)

		in_cbd2$ss5_pos=NA; in_cbd2$ss3_pos=NA; 
		in_cbd2$ss5_pos[if_end_ss5]=in_cbd2$end_pos[if_end_ss5]
		in_cbd2$ss5_pos[if_start_ss5]=in_cbd2$start_pos[if_start_ss5]-strand_signs[if_start_ss5]
		in_cbd2$ss3_pos[if_start_ss3]=in_cbd2$start_pos[if_start_ss3]
		in_cbd2$ss3_pos[if_end_ss3]=in_cbd2$end_pos[if_end_ss3]+strand_signs[if_end_ss3]
	}


	##output a table with all regulation information (include those without regulation, can be used for cis-elements analysis or GO analysis, etc)
	out_reg_f=paste(out_root,"/",ana_unit,"/formatted_tb/reguType.",use_p_type,P_cut,".Ch",
		ifelse(grepl("junc",ana_unit) & if_use_PSI_change, delta_PSI_cut,round(Log2ratio_cut,2)),
		".tbl",sep="")
	mkdir_if_not_exist(out_reg_f)
	out_format_d=re_format_tb(in_cbd2, 
		delete_headers=c("GeneUnit_id","row","gene_desc","min_P","ss_id","cds_len","transcript_len",suppNum_names,PSI_names,delta_PSI_names,Log2Ratio_names,pval_names,adjpval_names), 
		front_headers=c(setdiff(gene_info_headers,"gene_desc"), ss_pos_names,ss_seq_names), back_headers=c(ReguType_names,"gene_desc"))
	if( ana_unit=="geneBlock" ){
		out_format_d$Gblock_id2=paste(out_format_d$contig, out_format_d$strand, format(out_format_d$start_pos,scientific=F,trim=T), format(out_format_d$end_pos,scientific=F,trim=T),  sep=":")
	}
	write.table(out_format_d, file=out_reg_f, col.names=T, row.names=F, sep="\t", quote=F)

	#stats of regulation type
	tb=as.matrix( apply(in_cbd2[ReguType_names], 2, function(v){ table(factor(v,levels=ReguTypes))} ) [ReguTypes,] , nrow=length(ReguTypes) )
	colnames(tb)=ReguType_names
	print(tb)
	out_stats_f=paste(out_reg_f,".stats.tbl",sep="")
	write(paste("Regulation"),  file=out_stats_f)
	write.table(cbind(ReguTypes,tb), file=out_stats_f, col.names=T, row.names=F, sep="\t", quote=F, append=T)
	if( ana_unit=="geneBlock" ){
		#cross table between region_ano and regulation type ReguTypes
		region_anos=unique(in_cbd2$region_ano) ; region_anos=region_anos[!is.na(region_anos)]
		for(region_ano in region_anos){
			print(region_ano)
			ifsel=!is.na(in_cbd2$region_ano) & in_cbd2$region_ano==region_ano
			tb=as.matrix( apply(as.matrix(in_cbd2[ifsel, ReguType_names]), 2, function(v){ table(factor(v,levels=ReguTypes))} ) [ReguTypes,] , nrow=length(ReguTypes) )
			colnames(tb)=ReguType_names
			write(paste("\nRegion annotation=",region_ano),  file=out_stats_f, append=T)
			write.table(cbind(ReguTypes,tb), file=out_stats_f, col.names=T, row.names=F, sep="\t", quote=F, append=T)
		}
		if(any(ssSupp_headers %in% names(in_cbd2) ) ){
			ssSupp_headers1=intersect(ssSupp_headers, names(in_cbd2))
			for(region_ano in region_anos){
				print(region_ano)
				ifsel=!is.na(in_cbd2$region_ano) & in_cbd2$region_ano==region_ano
				sum_tb=NULL
				headers1=NULL; headers2=NULL; headers3=NULL
				for(ssSupp_header in ssSupp_headers1){
					for(ReguType_name in ReguType_names){
						tb=table( factor(in_cbd2[ifsel, ssSupp_header],levels=ssSupp_types),  factor(in_cbd2[ifsel, ReguType_name],levels=ReguTypes))[ssSupp_types,ReguTypes]
						headers3=c(headers3,colnames(tb))
						headers2=c(headers2,c(ReguType_name, rep("", ncol(tb)-1) ) )
						sum_tb=cbind(sum_tb,tb)
					}
					headers1=c(headers1, c(ssSupp_header, rep("", length(ReguType_names)*length(ReguTypes)-1 )  ) )
				}
				sum_tb2=rbind(headers1, headers2, headers3, sum_tb)
				sum_tb2=cbind(c("","","", ssSupp_types),sum_tb2)
				write(paste("\nRegion annotation=",region_ano),  file=out_stats_f, append=T)
				write.table(sum_tb2, file=out_stats_f, col.names=F, row.names=F, sep="\t", quote=F, append=T)
			}
		}

	}
	#cross match between multi samples:
	if(length(ReguType_names)>1){
		library("fmsb")
		compare_types_among_vars(out_stats_f, in_cbd2, type_headers1=ReguType_names, Types=ReguTypes)
	}


	##output a table with regulated events only:
	out_reg_f2=paste(out_reg_f,".regu.tbl",sep="")
	ifsel=rowSums ( !is.na(in_cbd2[ReguType_names]) & (in_cbd2[ReguType_names]=="UP" | in_cbd2[ReguType_names]=="DN") )>0 ; table(ifsel)
	out_format_d=in_cbd2[ifsel, ]
	out_format_d=re_format_tb(out_format_d, 
		delete_headers=c("GeneUnit_id","unit_id","row","ss_id","ss_seq_id","cds_len","transcript_len",delete_P_names), 
		front_headers=c( setdiff(gene_info_headers,"gene_desc"), ss_pos_names,ss_seq_names), 
		back_headers=c(Log2Ratio_names, use_P_names,ReguType_names,"gene_desc"))
	write.table(out_format_d, file=out_reg_f2, col.names=T, row.names=F, sep="\t", quote=F)
	if(length(ReguType_names)>1){
		cor_tb=calcu_cor_among_vars(in_cbd2[ifsel, ], Log2Ratio_names, corMeth="pea" )
		write(paste("\ncorrelation of regulation (Pearson correlation coefficient)"),  file=out_stats_f, append=T)
		write.table(cbind(rownames(cor_tb),cor_tb), file=out_stats_f, col.names=T, row.names=F, sep="\t", quote=F, append=T)
	}
	rm("out_format_d")
}


