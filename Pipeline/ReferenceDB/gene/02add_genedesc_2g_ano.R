#add more anotation to gene: eg. "alias" and "gene_desc" based on ncbi/gene/geneinfo/All_Mammalia.gene_info


source(paste("../../SharedCodes/Rfunc.inc.R",sep=""))


species <- "Canis_lupus_familiaris" #in HomoGene table header: Rattus_norvegicus, Homo_sapiens, Mus_musculus Macaca_mulatta Canis_lupus_familiaris
tax_id <- tax_ids[species]

#######for refflat transcripts
	geno <- "mm9" # mm9  hg19
	anot_folder <- paste("../ucsc/refFlat/",geno,"/",sep="") 
	ano_file <- "refFlat.txt"
	anno_headers=c("gene_symbol","refseqid","contig","strand","transc_start","transc_end","cds_start","cds_end","exon_num","exon_starts","exon_ends")
	
	geno <- "canFam3" # rheMac10 rn6 canFam3
	anot_folder <- paste("../ucsc/ncbiRefSeq/",geno,"/",sep="") 
	ano_file <- "ncbiRefSeq.txt"
	anno_headers=unlist(strsplit("bin refseqid contig strand transc_start transc_end cds_start cds_end exon_num exon_starts exon_ends score gene_symbol cdsStartStat cdsEndStat exonFrames"," "))

	outfile <- paste("02transcript_gene_ano/",geno,".refflat.desc.txt",sep="")
	out_nonredundant_f <- paste("02transcript_gene_ano/",geno,".refflat.desc.nonRedundant.txt",sep="")
	
	##
	old_ano_data <- read.table(file=paste(anot_folder,ano_file,sep=""), header=F, sep="\t") ; nrow(old_ano_data)
	names(old_ano_data) <- anno_headers

	#old_ano_data["gene_id"] <- geneidmapping(old_ano_data$gene_symbol,intype="alias",outtype="geneid",inspe=species, outspe=species)
	old_ano_data[c("gene_id","gene_desc")]=gsb2_gene_desc( old_ano_data$gene_symbol, spe=species) [,c("gene_id","gene_desc")]
	table(is.na(old_ano_data$gene_desc))
	table(is.na(old_ano_data$gene_id))
	##add gene_Biotype based on Refseq ID, gene symbol, transcript length etc
	old_ano_data$gene_Biotype=""
	old_ano_data$gene_Biotype[grepl("^NM_|^XM_|^YP_", old_ano_data$refseqid)]="mRNA"
	if_noncoding=grepl("^NR_|^XR_", old_ano_data$refseqid)
	old_ano_data$gene_Biotype[if_noncoding & old_ano_data$exon_num==1 & abs(old_ano_data$transc_start-old_ano_data$transc_end)<200]="small_noncoding"
	old_ano_data$gene_Biotype[if_noncoding & old_ano_data$gene_Biotype !="small_noncoding" ]="lncRNA"
	coding_gids=old_ano_data$gene_id[old_ano_data$gene_Biotype=="mRNA"]; length(unique(coding_gids))
	old_ano_data$gene_Biotype[ !is.na(old_ano_data$gene_id) & old_ano_data$gene_id %in% coding_gids & old_ano_data$gene_Biotype=="lncRNA" ]="ncIsoform"
	old_ano_data$gene_Biotype[old_ano_data$gene_Biotype=="small_noncoding" & grepl("^MIR|^Mir",old_ano_data$gene_symbol)]="miRNA"
	old_ano_data$gene_Biotype[old_ano_data$gene_Biotype=="small_noncoding" & grepl("^SNO|^Sno",old_ano_data$gene_symbol)]="snoRNA"
	table(old_ano_data$gene_Biotype)



############################RUN

####output
table(is.na(old_ano_data$gene_desc))
mkdir_if_not_exist(outfile)
write.table(old_ano_data,file=outfile,col.names=T,row.names=F,quote=F,sep="\t")

#remove redundant transcript (based on length)
gene_ano_d1=old_ano_data; nrow(gene_ano_d1)
gene_ano_d1$ori_rowID=1:nrow(gene_ano_d1) #original row ID
gene_ano_d1$gene_span=abs(gene_ano_d1$transc_start-gene_ano_d1$transc_end)
##generate gene cluster ID (genes overlap with each other), this is to keep the genes map to >1 non-contiguous regions of the same chromosome (eg. SMN1 and SMN2)
gene_ano_d1=gene_ano_d1[order(gene_ano_d1$contig, gene_ano_d1$strand,gene_ano_d1$transc_start,gene_ano_d1$transc_end), ]
total_rows=nrow(gene_ano_d1)
prev_rows=1:(total_rows-1)
next_rows=2:total_rows
pre_ends=gene_ano_d1$transc_end[prev_rows]
next_starts=gene_ano_d1$transc_start[next_rows]
if_same_cluster=gene_ano_d1$contig[prev_rows]==gene_ano_d1$contig[next_rows] & gene_ano_d1$strand[prev_rows]==gene_ano_d1$strand[next_rows] & next_starts<pre_ends
if_new_cluster=c(T,!if_same_cluster)
gene_ano_d1$clu_ID=cumsum(if_new_cluster)

gene_ano_d1=gene_ano_d1[order(gene_ano_d1$gene_symbol, as.numeric(gene_ano_d1$gene_Biotype=="mRNA"), gene_ano_d1$exon_num,gene_ano_d1$gene_span, decreasing=T), ]
gene_ano_d1=gene_ano_d1[!duplicated( paste(gene_ano_d1$gene_symbol, gene_ano_d1$clu_ID)), ]; nrow(gene_ano_d1)
table(is.na(gene_ano_d1$gene_symbol))
table(old_ano_data$contig)
table(gene_ano_d1$contig)
write.table(gene_ano_d1,file=out_nonredundant_f,col.names=T,row.names=F,quote=F,sep="\t")
