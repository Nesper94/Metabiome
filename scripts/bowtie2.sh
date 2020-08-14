#!/bin/sh

echo "Performing bowtie2 alignment. You should have a conda  environment in order to run this script:"

##Variables##:

threads=$1
host=$2
#name_env=$3
env_path="$CONDA_PREFIX/bin"
original_dir=$(pwd)

echo "Here lies the packages from your environment: $env_path"
cd $env_path

###############Installing required packages###################:

if [ -e bowtie2* ];then
	echo "bowtie2 installed"
else
	echo "Installing Bowtie2" &&  conda install bowtie2 --yes
fi
: '
if [ -e entrez-direct* ]; then
	echo "entrez-direct installed"
else
	echo "Installing entrez-direct: " && conda install entrez-direct --yes
fi
'
#################Downloading Human and PhiX reference genomes#################:
: '
for i in {1..10};do
	echo "PhiX and human genomes:" && esearch -db nucleotide -query "NC_001422.1" | efetch -format fasta > NC_001422.1.fasta
	esearch -db nucleotide -query "GCA_000001405.28"  | efetch -format fasta > GCA_000001405.28.fasta
	echo "Confirming downloaded genomes: "
	if  [ -e NC_001422.1.fasta ];then
			echo "PhiX genome was downloaded"
			if  [ -e GCA_000001405.28.fasta ];then
				echo "Human genome was downloaded"
				break
			fi
	fi;done
'
cd $original_dir

###Concatenating genomes to be aligned#####:

#cat $host NC_001422.1.fasta GCA_000001405.28.fasta > Mixed.fasta
cat $host NC_001422.fasta GCF_000001405.39_GRCh38.p13_genomic > Mixed.fasta # Se eliminó genoma humano porque no estaba todavía en el cluster

#########Building genome index and bowtie alignment###################:
echo "Building genome index:"
bowtie2-build Mixed.fasta Mix --threads $threads && echo 'Indexing genomes to filter out...'


############################PE alignment##############################:
mkdir aligned
	for i in $(ls *f-paired.fq.gz -1);do
		#cd aligned
		#touch $(echo $i | sed 's/f-paired.fq.gz/PE_aligned.sam/')
		#touch $(echo $i | sed 's/f-paired.fq.gz/SE_aligned.sam/')
		#cd ..
		echo "PE alignment: "
		bowtie2 -x Mix -1 $i -2 $(echo $i | sed 's/f-paired/r-paired/') \
		--un-conc-gz aligned/$(echo $i | sed 's/f-paired.fq.gz/PE_aligned.fq.gz/') \
		-q -p $threads --met-file $(echo $i | sed 's/f-paired.fq.gz/PE_summary.txt/');done

############################SE alignment###############################:
	for i in $(ls *unpaired.fq.gz -1);do
		#cd aligned
		#touch $(echo $i | sed 's/unpaired.fastq.gz/unpaired_aligned.sam/')
		#cd ..
		echo "SE alignemnt: "
		bowtie2 -x Mix -U $i \
		--un-gz aligned/$(echo $i | sed 's/unpaired.fastq.gz/unpaired.fq.gz/') \
		-q -p $threads --met-file $(echo $i | sed 's/unpaired.fastq.gz/unpaired_summary.txt/');done
