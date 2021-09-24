##download files from ncbi ftp site

source ../../setGlobalVars.sh
echo $pipeline_root

##download gene information (NCBI gene symbol, gene id, description etc)
downloadDir=$pipeline_root/ReferenceDB/ncbi/gene/geneinfo
mkdir -p $downloadDir
cd $downloadDir
wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/All_Mammalia.gene_info.gz 
wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz
wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz

wget https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/All_Mammalia.gene_info.gz
wget https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Rattus_norvegicus.gene_info.gz
wget https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Canis_familiaris.gene_info.gz
gunzip *.gz

