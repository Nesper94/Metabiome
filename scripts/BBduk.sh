

echo "Performing kmer 16S rDNA kmer picking strategy with bbDuk: "

##Variables##:

threads=$1
env_path="$CONDA_PREFIX"
dir="$(pwd)"

echo "Here lies the packages from your environment: $env_path"
cd $env_path/bin

###############Installing required packages###################:

if [ -e bbmap* ];then
	echo "bbmap is installed"
else
	echo "Installing bbmap" &&  conda install bbmap --yes
fi

cd dir
##Matching reads to the 16S RNA SSU from SILVA Database, based on a Kmer approach by BBduk:
#For PE reads:
for i in $(ls *f-paired.fq.gz -1);do
  bbduk.sh in=$i in2=$(echo $i | sed 's/f-paired/r-paired/') ref="" outm=$(echo $i | sed 's/f-paired/f-paired-16S/') \
  outm2=$(echo $i | sed 's/f-paired/r-paired-16S/') outs=$(echo $i | sed 's/f-paired/unpaired-16S/') \
  stats=$(echo $i | sed 's/f-paired/statistics_16S/') ordered=T ;done


for i in $(ls *unpaired.fq.gz -1);do
  bbduk.sh in=$i ref="" outm=$(echo $i | sed 's/f-paired/unpaired-16S/')  stats=$(echo $i | sed 's/f-paired/statistics_16S/') \
  ordered=T;done
