#!/bin/bash
# download and process databse file 
VEP_version="105"
echo "download reference file"
cd database

# data to download: reference, 1000G, dbsnp, vep, star-fusion, funcotator, pep, cdna

# download reference gene fasta file
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip hg19.ga.gz

# 1000G, dbsnp, mills download and unzip
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg19/1000G_phase1.indels.hg19.sites.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg19/1000G_phase1.indels.hg19.sites.vcf.idx.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg19/dbsnp_138.hg19.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg19/dbsnp_138.hg19.vcf.idx.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx.gz

# download vep reference file and unzip
wget http://ftp.ensembl.org/pub/grch37/release-${VEP_version}/variation/vep/homo_sapiens_vep_${VEP_version}_GRCh37.tar.gz
tar xvzf homo_sapiens_vep_${VEP_version}_GRCh37.tar.gz

# download dataset for star-fusion and unzip
wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/__genome_libs_StarFv1.10/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play.tar.gz
tar xvzf GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play.tar.gz

# download pep&cnda and unzip
wget http://ftp.ensembl.org/pub/grch37/release-${VEP_version}/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.cdna.all.fa.gz
gunzip Homo_sapiens.GRCh37.cdna.all.fa.gz
wget http://ftp.ensembl.org/pub/grch37/release-${VEP_version}/fasta/homo_sapiens/pep/Homo_sapiens.GRCh37.pep.all.fa.gz
gunzip Homo_sapiens.GRCh37.pep.all.fa.gz

# download funcotator dataset
wget https://console.cloud.google.com/storage/browser/_details/broad-public-datasets/funcotator/funcotator_dataSources.v1.6.20190124s.tar.gz
tar xzvf funcotator_dataSources.v1.6.20190124s.tar.gz

# download refseq annotation for hg19
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.refGene.gtf.gz
