#!/usr/bin/env bash

# This script is used to compare a metagenome with a genome using MetaGecko
# Parameters to call it are:
#  @param metagenome relative or absolute path to metagenome fasta file
#  @param genome relative or absolute path to genome fasta file
#  @param outFileName name of output file. Extension ".frags" will be added.
#  @param S minimum similarity of fragments.
#  @param L minimum length of fragments.
# This script will generate the following files and folders:
#  @file {genome}-revercomp is the reverse complementary fasta file of genome file given.
#  @folder dictionaries are the fodler where genome and metagenome dictionaries will be stored.
#    @file {metagenome}.(d2hP/d2hW) are the files that conform the metagenome dictionary. 
#    @file {genome}.(d2hP/d2hW) are the files that conform the genome dictionary.
#    @file {genome}.words.sort it's a file that Gecko uses to generate the genome dictionary.
#    @file {genome}-revercomp.(d2hP/d2hW) are the files that conform the reverse complementary genome dictionary.
#    @file {genome}-revercomp.words.sort it's a file that Gecko uses to generate the reverse complementary genome dictionary.
#  @folder fragments is the folder where fragments will be stored.
#    @file {metagenome}-{genome}.frags is the fragments file generated.
# In some cases this last file could not appear. If the program ends without print errors in
# the screen look if exists one or any of the following files:
#  @file {metagenome}-{genome}-f.frags is the forward genome-metagenome fragments.
#  @file {metagenome}-{genome}-revercomp.frags is the reverse genome-metagenome fragments.
# When there are not any hit (any equal segment) between genome (f/r) and metagenome files, the
# FRAG file will not be generated. If any of this file doesn't exists, the combined fragment
# file ({metagenome}-{genome}.frags) will not be generated.
# NOTE: if any of this files exists check if the program didn't write any error message. If it
#       didn't, means that there're 0% of equal parts between metagenome and genome.
# Check arguments
if [ $# != 5 ]; then
        echo "***ERROR*** Use: $0 metagenome genome S L prefix "
        exit -1
fi


MGDIR=$(pwd)

metagenome=$(basename "$1")
genome=$(basename "$2")
metagExt="${genome##*.}"
genoExt="${metagenome##*.}"
metagenome="${metagenome%.*}"
genome="${genome%.*}"


#Copy Metagenome and Database to current folder
if [[ ! -e "$MGDIR/$(basename "$1")" ]]; then
    ln $1 $MGDIR/$(basename "$1")
fi 

if [[ ! -e "$MGDIR/$(basename "$2")" ]]; then
    ln $2 $MGDIR/$(basename "$2")
fi 


mkdir -p ${MGDIR}/dictionaries
mkdir -p ${MGDIR}/fragments


L=$4
S=$3



# Take binaries folder
BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Create database dictionarie
if [[ ! -f dictionaries/${genome}.d2hP || ! -f dictionaries/${genome}.d2hP ]];    then
        echo " Writing genome dict: ${genome}.${genoExt}"
        #New dict genome
        ${BINDIR}/dict ${genome}.${genoExt} dictionaries/${genome} 32
fi


# Create metagenome dictionarie
if [[ ! -f dictionaries/${metagenome}.d2hP || ! -f dictionaries/${metagenome}.d2hP ]];	then
	echo " Writting metagenome dict: ${metagenome}.${metagExt}"
	#New dict metagenome
	${BINDIR}/dict ${metagenome}.${metagExt} dictionaries/${metagenome} 32 f
fi

# Generate fragments
if [[ ! -f fragments/${metagenome}-${genome}-f.frags ]];	then
	echo "Writting fragments for: ${metagenome} - ${genome}"
	${BINDIR}/frag dictionaries/${metagenome} ${metagenome}.${metagExt} dictionaries/${genome} ${genome}.${genoExt} fragments/${metagenome}-${genome}-f $S $L $5
fi
