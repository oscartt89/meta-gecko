#!/usr/bin/env bash

# This script is used to compare a metagenome set with a genome set using MetaGecko
# Parameters to call it are:
#  @param metagenomesDir is the relative or absolute path to metagenomes folder.
#  @param genomesDir is the relative or absolute path to genomes folder.
#  @param L is the minimum length of fragments.
#  @param S minimum similarity of fragments.
#  @param MetagExtension is the extension of metagenomes files.
#  @param GenoExtension is the extension of genomes files.
#  @param prefix is the length of the prefix (in bytes) of the word wanted. Ex: want use k-mers of length 20 you must use prefix
#         5 because 1 byte = 4 letters => 20 letters = 5 bytes
# This script will generate the following files and folders:
#  @folder {metagenomesDir}/dictionaries is the folder where metagenome dictionaries will be generated. 
#    @file {metagenome}.(d2hP/d2hW) are the files that conform the metagenome dictionary. 
#  @folder {genomeDir}/complementaries is the folder where reverse complementary genome fasta files will be stored.
#    @file {genome}-revercomp is the reverse complementary fasta file of genome file given.
#  @folder {genomesDir}/dictionaries is the folder where genome dictionaries will be generated.
#    @file {genome}.(d2hP/d2hW) are the files that conform the genome dictionary.
#    @file {genome}-revercomp.(d2hP/d2hW) are the files that conform the reverse complementary genome dictionary.
#    @file {genome}.words.sort it's a file that Gecko uses to generate the genome dictionary.
#    @file {genome}-revercomp.words.sort it's a file that Gecko uses to generate the reverse complementary genome dictionary.
#  @folder ./fragments is the folder where fragments generated will be stored.
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



metaDir=$1
genomeDir=$2
L=$3
S=$4
EXT=$5
GEXT=$6

# Metagenomes array
metagenomes=()
x=0

# Genomes array
genomes=()
y=0

if [ $# != 7 ]; then
	echo "***ERROR*** Use: $0 metagenomesDir genomesDir L S MetagExtension GenoExtension prefix"
	exit -1
fi

# Take binaries folder
BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Take metagenome files
for elem in $(ls -d ${metaDir}/*.$EXT | awk -F "/" '{print $NF}' | awk -F ".$EXT" '{print $1}')
do
	metagenomes[$x]=$elem
	x=`expr $x + 1`
done

# Take genome files
for elem in $(ls -d ${genomeDir}/*.$GEXT | awk -F "/" '{print $NF}' | awk -F ".$GEXT" '{print $1}')
do
	genomes[$y]=$elem
	y=`expr $y + 1`
done

# Create necessary directories
mkdir -p ${metaDir}/dictionaries
mkdir -p ${genomeDir}/dictionaries
mkdir -p ${genomeDir}/complementaries
mkdir -p fragments


# Write all genome dictionaries
for ((k=0 ; k < ${#genomes[@]} ; k++))
do
	genome=${genomes[$k]}
	if [[ ! -f $genomeDir/dictionaries/${genome}.d2hP || ! -f $genomeDir/dictionaries/${genome}.d2hW ]];	then
		echo " Writting genome dict: ${genome}.${GEXT}"
		# Find words and order
		${BINDIR}/words $genomeDir/${genome}.${GEXT} $genomeDir/dictionaries/${genome}.words.unsort
		${BINDIR}/sortWords 10000000 32 $genomeDir/dictionaries/${genome}.words.unsort $genomeDir/dictionaries/${genome}.words.sort

		# Create hash table in disk
		${BINDIR}/w2hd $genomeDir/dictionaries/${genome}.words.sort $genomeDir/dictionaries/${genome}
	fi

	# Now with reverse seq
	if [[ ! -f $genomeDir/dictionaries/${genome}-revercomp.d2hP || ! -f $genomeDir/dictionaries/${genome}-revercomp.d2hW ]];	then
		${BINDIR}/reverseComplement $genomeDir/${genome}.${GEXT} $genomeDir/complementaries/${genome}-revercomp.${GEXT}
		echo " Writting genome dict: ${genome}-revercomp.${GEXT}"
		# Find words and order
		${BINDIR}/words $genomeDir/complementaries/${genome}-revercomp.${GEXT} $genomeDir/dictionaries/${genome}-revercomp.words.unsort
		${BINDIR}/sortWords 10000000 32 $genomeDir/dictionaries/${genome}-revercomp.words.unsort $genomeDir/dictionaries/${genome}-revercomp.words.sort

		# Create hash table in disk
		${BINDIR}/w2hd $genomeDir/dictionaries/${genome}-revercomp.words.sort $genomeDir/dictionaries/${genome}-revercomp
	fi
done

# Generate metagenome dictionaries and fragments
for ((i=0 ; i < ${#metagenomes[@]} ; i++))
do
	seqX=${metagenomes[$i]}

	# Create metagenome dictionarie
	if [[ ! -f $metaDir/dictionaries/${seqX}.d2hP || ! -f $metaDir/dictionaries/${seqX}.d2hW ]];	then
		echo " Writting metagenome dict: ${seqX}.${EXT}"
		# Metagenome dict
		${BINDIR}/dict $metaDir/${seqX}.${EXT} $metaDir/dictionaries/${seqX} 32
	fi

	# Generate fragments
	for ((j=0 ; j < ${#genomes[@]} ; j++))
	do
		seqY=${genomes[$j]}
		# Forward
		if [[ ! -f fragments/${seqX}-${seqY}-f.frags ]];	then
			echo "Writting new fragments for: ${seqX} - ${seqY}"
			${BINDIR}/frag $metaDir/dictionaries/${seqX} $metaDir/${seqX}.${EXT} $genomeDir/dictionaries/${seqY} $genomeDir/${seqY}.${GEXT} fragments/${seqX}-${seqY}-f $S $L f $7
		fi
		# Reverse
		if [[ ! -f fragments/${seqX}-${seqY}-revercomp.frags ]];	then
			echo "Writting new fragments for: ${seqX} - ${seqY}-revercomp"
			${BINDIR}/frag $metaDir/dictionaries/${seqX} $metaDir/${seqX}.${EXT} $genomeDir/dictionaries/${seqY}-revercomp $genomeDir/complementaries/${seqY}-revercomp.${GEXT} fragments/${seqX}-${seqY}-revercomp $S $L r $7
			# Fixe fragments
			${BINDIR}/fixeReverseFrags $genomeDir/${seqY}.${GEXT} $genomeDir/complementaries/${seqY}-revercomp.${GEXT} fragments/${seqX}-${seqY}-revercomp.frags
		fi

		# Combine fragment
		if [[ -f fragments/${seqX}-${seqY}-revercomp.frags && -f fragments/${seqX}-${seqY}-f.frags ]];	then
			echo "Combining forward and reverse: ${seqX}-${seqY}"
			${BINDIR}/combineFragments fragments/${seqX}-${seqY}-f.frags fragments/${seqX}-${seqY}-revercomp.frags fragments/${seqX}-${seqY}.frags 1
		fi
	done
done

