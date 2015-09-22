#!/usr/bin/env bash

# Take variables
metaDir=$1 # Metagenomes directory
genomeDir=$2 # Genomes directory
L=$3 # Minimal length
S=$4 # Similarity
WL=$5 # Word length
EXT=$6 # Files extension

# Metagenome files array
metaFiles=()
x=0

# Genome files array
genoFiles=()
y=0

# Check bad call error
if [ $# != 6 ]; then
	echo "***ERROR*** Use: $0 metagenomesDir genomesDir L(200) S(40) K(8) fastaFilesExtension"
	exit -1
fi

# Take bin directory
BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo "Metagenome files:"

# Take all metagenome files
for elem in $(ls -d ${metaDir}/*.$EXT | awk -F "/" '{print $NF}' | awk -F ".$EXT" '{print $1}')
do
	metaFiles[$x]=$elem
	x=`expr $x + 1`
	echo "   $elem"
done

echo "Genome files:" 

# Take all genome files
for elem in $(ls -d ${genomeDir}/*.$EXT | awk -F "/" '{print $NF}' | awk -F ".$EXT" '{print $1}')
do
	genoFiles[$y]=$elem
	y=`expr $y + 1`
	echo "   $elem"
done

a