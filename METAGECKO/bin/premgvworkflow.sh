#!/usr/bin/env bash

# This script compares a metagenome against a database of genomes
# Input parameters:
#  @param metagenome 	relative or absolute path to metagenome file.
#  @param database 		relative or absolute path to database file.
#  @param S 			minimum similarity to report a fragment.
#  @param L 			minimum length to report a fragment.
#  @param prefix		size of mers to look for.
#
# This script will generate the following files and folders:
#  @file {genome}-revercomp 				The reverse complementary of the genome.
#  @folder dictionaries 					Folder where the generated dictionaries will be stored.
#    @file {metagenome}.(d2hP/d2hW) 		Dictionary of the metagenome.
#    @file {genome}.(d2hP/d2hW) 			Dictionary of the database.
#    @file {genome}-revercomp.(d2hP/d2hW) 	Dictionary of the reverse complement of the database.
#    @file {genome}.words.sort 				File that contains all k-mers in the database.
#    @file {genome}-revercomp.words.sort 	Sorted file of the above.
#  @folder fragments 						The folder where fragments will be stored.
#    @file {metagenome}-{genome}.frags 		The generated fragments file.
#
# If no fragments file is generated:
# (1) There might be no hits between the metagenome and the database. Lower the similarity and length thresholds and try a lower prefix size.
# (2) Check the terminal output to verify that there were no errors reported.
#



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

# Run the preselection method

echo "Finding seeds for database reduction"
${BINDIR}/quickhits ${metagenome}.${metagExt} ${genome}.${genoExt} dictionaries/${genome}.hseq 3
${BINDIR}/dbReduce ${genome}.${genoExt} dictionaries/${genome}.hseq

# Set new extension
genoExt=${genoExt}.presec


# Create database dictionarie
if [[ ! -f dictionaries/${genome}.d2hP || ! -f dictionaries/${genome}.d2hP ]];    then
        echo " Writing genome dict: ${genome}.${genoExt}"
        #New dict genome
        ${BINDIR}/dict ${genome}.${genoExt} dictionaries/${genome} 32
fi


# Create metagenome dictionarie
if [[ ! -f dictionaries/${metagenome}.d2hP || ! -f dictionaries/${metagenome}.d2hP ]];	then
	echo " Writing metagenome dict: ${metagenome}.${metagExt}"
	#New dict metagenome
	${BINDIR}/dict ${metagenome}.${metagExt} dictionaries/${metagenome} 32 f
fi

# Generate fragments
if [[ ! -f fragments/${metagenome}-${genome}-f.frags ]];	then
	echo "Writing fragments for: ${metagenome} - ${genome}"
	${BINDIR}/frag dictionaries/${metagenome} ${metagenome}.${metagExt} dictionaries/${genome} ${genome}.${genoExt} fragments/${metagenome}-${genome}-f $S $L $5
fi

rm -rf ${metagenome}.${metagExt}
rm -rf ${genome}.${genoExt}

