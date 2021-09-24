

##download Refflat files from UCSC and gene information files from NCBI and prepare reference gene/transcript files for RNA-seq data analysis

##download from ucsc:
./ucsc/download_ucsc.sh

##download from ncbi:
./ncbi/download_ncbi.sh

##build gene alias to gene official symbol and gene ID table
cd ./gene
perl 01build_gsb_alias2id_tb.pl -t 9606 -o human -i ../ncbi/gene/geneinfo/Homo_sapiens.gene_info
perl 01build_gsb_alias2id_tb.pl -t 10090 -o mouse -i ../ncbi/gene/geneinfo/Mus_musculus.gene_info
perl 01build_gsb_alias2id_tb.pl -t 9544 -o monkey -i ../ncbi/gene/geneinfo/All_Mammalia.gene_info

##add gene description, ID, gene biotype to gene table
02add_genedesc_2g_ano.R  # manualy run this R script


