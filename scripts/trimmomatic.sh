#!/bin/bash
# Este script instala y ejecuta Trimmomatic usando como input una carpeta con
# archivos en formato FastQ con los reads en bruto.
# 2020-07-10

echo 'Welcome to the best pipeline ever (?)
'

script_dir=$(readlink -f $(dirname "$0")) #Lines to get path to main directory
main_dir=$(cd $script_dir ; cd .. ; pwd)
dir=$(readlink -f $1)

if [ ! -d "$dir" ]; then
   echo "Error: You need to specify the directory containing reads in FastQ"
   exit 1
fi

echo "Your read directory is $dir"

cd $main_dir

if [ ! -d "log-files" ]; then
    echo 'Creating log-files/ directory...'; mkdir log-files
fi

# Install Trimmomatic

for f in Trimmomatic*; do
    [ -d "$f" ] && echo "Trimmomatic is already installed." ||
    (echo 'Installing Trimmomatic...'
    (wget -nd -P software/ http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip || exit 1)
    cd software/
    unzip Trimmomatic-0.39.zip
    cd $main_dir)
    break
done

trim_loc=$main_dir/software/Trimmomatic-0.39/trimmomatic-0.39.jar

if [ ! -d "results/trimmed-reads" ]; then
    echo "Creating directory 'results/trimmed-reads/'"
	mkdir -p results/trimmed-reads
fi

read -p "Write arguments for Trimmomatic: " trimopt
# Las opciones utilizadas con los reads de prueba el 2020-07-13
# fueron: MINLEN:140 TRAILING:25 HEADCROP:20
# Es importante recordar que en Trimmomatic el orden de las opciones indica
# su orden de ejecución.

cd $dir

for i in $(ls *R1*.fastq -1);
do java -jar $trim_loc PE $i $(echo $i | sed 's/R1/R2/') \
$main_dir/results/trimmed-reads/$(echo $i | sed 's/R1.*/f-paired.fq.gz/') $main_dir/results/trimmed-reads/$(echo $i | sed 's/R1.*/f-unpaired.fq.gz/') \
$main_dir/results/trimmed-reads/$(echo $i | sed 's/R1.*/r-paired.fq.gz/') $main_dir/results/trimmed-reads/$(echo $i | sed 's/R1.*/r-unpaired.fq.gz/') \
$trimopt 2>&1 | tee -a ../log-files/trimmomatic.log; done # Cambio

exit
