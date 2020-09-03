#!/bin/bash
# Metagenomic assembly using metaSPAdes from SPAdes-3.12.0
# Written by: Phagomica_club

#Ejecuta el comando de actualización para actualizar los repositorios de paquetes
# y obtener la información más reciente de los paqutes disponibles(para qué el -y?)

set -e # Para que el script termine si encuentra un error en un comando.

function usage() {
    echo "Usage: $0 -i <input directory> -o <output directory> [-t <threads>] [-k <kmers>]"
    echo "Output directory will be created if it doesn't exists."
}

if [[ "$#" == 0 ]]; then
    echo "No arguments given."
    usage
    exit 1
fi

while [[ -n "$1" ]]; do
    case "$1" in
        -h|--help ) usage; exit 0
            ;;
        -i )        reads_dir=$(readlink -f "$2")
            shift
            ;;
        -o )        out_dir=$(readlink -f "$2")
            shift
            ;;
        -t )        threads="$2"
            shift
            ;;
        -k )        kmers="$2"
            shift
            ;;
        * )         echo "Option '$1' not recognized"; exit 1
            ;;
    esac
    shift
done

# SPAdes #

wget http://cab.spbu.ru/files/release3.12.0/SPAdes-3.12.0-Linux.tar.gz
tar -xzf SPAdes-3.12.0-Linux.tar.gz #extraer archivo
# cd SPAdes-3.12.0-Linux/bin/
# export PATH="directory_path:$PATH" #add SPAdes installation directory to the PATH variable


#Para verificar la instalación, remember add spades to PATH variable
#spades.py --test

# Run metaSPAdes #

for i in "$reads_dir"/*f-paired*.fq.gz;do

	echo " PE assembly"
	spades.py --meta \
    -o "$out_dir" `# Se asume que ya se agregó el directorio a la variable PATH` \
    -1 "$i" `# Forward files` \
    -2 $(echo "$i" | sed 's/.1.gz/.2.gz/') `# Reverse sequences` \
    -s "$unpaired_unaligned" `# unpaired. "$reads_dir"/*unpaired_unaligned.fq.gz` \
    -t "$2" `#16` \
    -k "$kmers" `# list of k-mer sizes to be use separeted comma(no more than 128)`

done

