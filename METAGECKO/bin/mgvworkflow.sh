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

metagenome=$(basename "$1")
metagExt="${metagenome##*.}"
metagenome="${metagenome%.*}"
genome=$(basename "$2")
genoExt="${genome##*.}"
genome="${genome%.*}"
L=$4
S=$3

# Check arguments
if [ $# != 5 ]; then
	echo "***ERROR*** Use: $0 metagenome genome S L prefix"
	exit -1
fi

# Take binaries folder
BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Create necessary directories
mkdir -p dictionaries
mkdir -p fragments

# Check if dictionary exists
if [[ ! -f dictionaries/${genome}.d2hP || ! -f dictionaries/${genome}.d2hW ]];	then
	echo " Writting genome dict: ${genome}.${genoExt}"
	# Find words and order
	${BINDIR}/words ${genome}.${genoExt} dictionaries/${genome}.words.unsort
	${BINDIR}/sortWords 10000000 32 dictionaries/${genome}.words.unsort dictionaries/${genome}.words.sort

	# Create hash table in disk
	${BINDIR}/w2hd dictionaries/${genome}.words.sort dictionaries/${genome}
fi

# Now with reverse seq
if [[ ! -f dictionaries/${genome}-revercomp.d2hP || ! -f dictionaries/${genome}-revercomp.d2hW ]];	then
	# Generate reverse
	${BINDIR}/reverseComplement ${genome}.${genoExt} ${genome}-revercomp.${genoExt}
	echo " Writting genome dict: ${genome}-revercomp.${genoExt}"
	
	# Find words and order
	${BINDIR}/words ${genome}-revercomp.${genoExt} dictionaries/${genome}-revercomp.words.unsort
	${BINDIR}/sortWords 10000000 32 dictionaries/${genome}-revercomp.words.unsort dictionaries/${genome}-revercomp.words.sort

	# Create hash table in disk
	${BINDIR}/w2hd dictionaries/${genome}-revercomp.words.sort dictionaries/${genome}-revercomp
fi

# Create metagenome dictionarie
if [[ ! -f dictionaries/${metgenome}.d2hP || ! -f dictionaries/${metgenome}.d2hP ]];	then
	echo " Writting metagenome dict: ${metagenome}.${metagExt}"
	#New dict metagenome
	${BINDIR}/dict ${metagenome}.${metagExt} dictionaries/${metagenome} 32
fi

# Generate fragments
if [[ ! -f fragments/${metagenome}-${genome}-f.frags ]];	then
	echo "Writting fragments for: ${metagenome} - ${genome}"
	${BINDIR}/frag dictionaries/${metagenome} ${metagenome}.${metagExt} dictionaries/${genome} ${genome}.${genoExt} fragments/${metagenome}-${genome}-f $S $L f $5
fi

# Reverse fragments
if [[ ! -f fragments/${metagenome}-${genome}-revercomp.frags ]];	then
	echo "Writting new fragments for: ${metagenome} - ${genome}-revercomp"
	${BINDIR}/frag dictionaries/${metagenome} ${metagenome}.${metagExt} dictionaries/${genome}-revercomp ${genome}-revercomp.${genoExt} fragments/${metagenome}-${genome}-revercomp $S $L r $5
	# Fixe coordinates
	${BINDIR}/fixeReverseFrags ${genome}.${genoExt} ${genome}-revercomp.${genoExt} fragments/${metagenome}-${genome}-revercomp
fi

# Combine fragment
if [[ -f fragments/${metagenome}-${genome}-revercomp.frags && -f fragments/${metagenome}-${genome}-f.frags ]];	then
	echo "Combining forward and reverse: ${metagenome}-${genome}"
	${BINDIR}/combineFragments fragments/${metagenome}-${genome}-f.frags fragments/${metagenome}-${genome}-revercomp.frags fragments/${metagenome}-${genome}.frags 1
fi