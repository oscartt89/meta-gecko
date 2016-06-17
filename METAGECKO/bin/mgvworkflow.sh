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


MGDIR=$(dirname "$1")

metagenome=$(basename "$1")
genome=$(basename "$2")
metagExt="${genome##*.}"
genoExt="${metagenome##*.}"
metagenome="${metagenome%.*}"
genome="${genome%.*}"


#Copy Database to metagenome folder if they are not in the same folder
if [[ ! -e "$MGDIR/$(basename "$2")" ]]; then
    ln $2 $MGDIR/$(basename "$2")
fi 


mkdir -p ${MGDIR}/dictionaries
mkdir -p ${MGDIR}/fragments


L=$4
S=$3


# Check arguments
if [ $# != 5 ]; then
	echo "***ERROR*** Use: $0 metagenome genome S L prefix"
	rm -rf $MGDIR
	exit -1
fi



# Take binaries folder
BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


# Check if dictionary exists
if [[ ! -f ${MGDIR}/dictionaries/${genome}.d2hP || ! -f ${MGDIR}/dictionaries/${genome}.d2hW ]];	then
	echo " Writing genome dict: ${MGDIR}/${genome}.${genoExt}"
	# Find words and order

	${BINDIR}/words ${MGDIR}/${genome}.${genoExt} ${MGDIR}/dictionaries/${genome}.words.unsort

	${BINDIR}/sortWords 10000000 32 ${MGDIR}/dictionaries/${genome}.words.unsort ${MGDIR}/dictionaries/${genome}.words.sort
	
	# Create hash table in disk
	${BINDIR}/w2hd ${MGDIR}/dictionaries/${genome}.words.sort ${MGDIR}/dictionaries/${genome}
fi


#: <<'END'

# Now with reverse seq
if [[ ! -f ${MGDIR}/dictionaries/${genome}-revercomp.d2hP || ! -f ${MGDIR}/dictionaries/${genome}-revercomp.d2hW ]];	then
	# Generate reverse
	${BINDIR}/reverseComplement ${MGDIR}/${genome}.${genoExt} ${MGDIR}/${genome}-revercomp.${genoExt}
	echo " Writing genome dict: ${MGDIR}/${genome}-revercomp.${genoExt}"
	
	# Find words and order
	${BINDIR}/words ${MGDIR}/${genome}-revercomp.${genoExt} ${MGDIR}/dictionaries/${genome}-revercomp.words.unsort
	${BINDIR}/sortWords 10000000 32 ${MGDIR}/dictionaries/${genome}-revercomp.words.unsort ${MGDIR}/dictionaries/${genome}-revercomp.words.sort

	# Create hash table in disk
	${BINDIR}/w2hd ${MGDIR}/dictionaries/${genome}-revercomp.words.sort ${MGDIR}/dictionaries/${genome}-revercomp
fi

# Create metagenome dictionarie
if [[ ! -f ${MGDIR}/dictionaries/${metgenome}.d2hP || ! -f ${MGDIR}/dictionaries/${metgenome}.d2hP ]];	then
	echo " Writing metagenome dict: ${MGDIR}/${metagenome}.${metagExt}"
	#New dict metagenome
	${BINDIR}/dict ${MGDIR}/${metagenome}.${metagExt} ${MGDIR}/dictionaries/${metagenome} 32
fi

# Generate fragments
if [[ ! -f ${MGDIR}/fragments/${metagenome}-${genome}-f.frags ]];	then
	echo "Writing fragments for: ${MGDIR}/${metagenome} - ${genome}"
	${BINDIR}/frag ${MGDIR}/dictionaries/${metagenome} ${MGDIR}/${metagenome}.${metagExt} ${MGDIR}/dictionaries/${genome} ${MGDIR}/${genome}.${genoExt} ${MGDIR}/fragments/${metagenome}-${genome}-f $S $L f $5
fi

# Reverse fragments
if [[ ! -f ${MGDIR}/fragments/${metagenome}-${genome}-revercomp.frags ]];	then
	echo "Writing new fragments for: ${MGDIR}/${metagenome} - ${genome}-revercomp"
	${BINDIR}/frag ${MGDIR}/dictionaries/${metagenome} ${MGDIR}/${metagenome}.${metagExt} ${MGDIR}/dictionaries/${genome}-revercomp ${MGDIR}/${genome}-revercomp.${genoExt} ${MGDIR}/fragments/${metagenome}-${genome}-revercomp $S $L r $5
	# Fixe coordinates
	
	${BINDIR}/fixeReverseFrags ${MGDIR}/${genome}.${genoExt} ${MGDIR}/${genome}-revercomp.${genoExt} ${MGDIR}/fragments/${metagenome}-${genome}-revercomp.frags
	
fi

# Combine fragment
if [[ -f ${MGDIR}/fragments/${metagenome}-${genome}-revercomp.frags && -f ${MGDIR}/fragments/${metagenome}-${genome}-f.frags ]];	then
	echo "Combining forward and reverse: ${MGDIR}/${metagenome}-${genome}"
	${BINDIR}/combineFragments ${MGDIR}/fragments/${metagenome}-${genome}-f.frags ${MGDIR}/fragments/${metagenome}-${genome}-revercomp.frags ${MGDIR}/fragments/${metagenome}-${genome}.frags 1
fi
#END

#Remove hard link
rm $MGDIR/$(basename "$2")



