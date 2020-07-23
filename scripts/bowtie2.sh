#!/bin/bash
# 2020-07-22

script_dir=$(readlink -f $(dirname "$0"))
echo "$script_dir"
main_dir=$(cd $script_dir ; cd .. ; pwd)
echo "$main_dir"

cd $main_dir

echo 'Installing Bowtie2...'

for f in software/bowtie2*; do
    [ -d "$f" ] && echo "Bowtie2 is already installed." ||
    wget -nd -P software/ https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.1/bowtie2-2.4.1-linux-x86_64.zip || exit 1
    cd software
    echo 'Extracting bowtie2-2.4.1-linux-x86_64.zip...'
    unzip bowtie2-2.4.1-linux-x86_64.zip || exit 1
    break
done

#Create BT2_HOME environment variable
export BT2_HOME=$(readlink -f bowtie2-2.4.1-linux-x86_64)

echo "$BT2_HOME"
