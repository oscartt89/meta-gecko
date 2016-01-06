#!/bin/bash 

# This script launch a set of experiments where, each time, a metagenome
# is artificially generated using a given genome. After that a dictionary
# of the metagenome is created using the old gecko version and the new 
# geckoMGV saving times spended on this actions in a file given.

BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


if [ $# != 9 ]; then
   echo " ==== ERROR ... you called this script inappropriately."
   echo ""
   echo "   usage: $0 genome.fasta startR startL maxR maxL incrR incrL d_dataFile f_dataFile"
   echo ""
   exit -1
fi

# Store info
startR=$2
startL=$3
maxR=$4
maxL=$5
incrR=$6
incrL=$7

# Create necessary directories
mkdir -p ExpOldNew/old
mkdir -p ExpOldNew/new
mkdir -p ExpOldNew/genoDic

# Create genome dictionary
${BINDIR}/words $1 ExpOldNew/geno.words.unsort
${BINDIR}/sortWords 10000000 32 ExpOldNew/geno.words.unsort ExpOldNew/geno.words.sort
${BINDIR}/w2hd ExpOldNew/geno.words.sort ExpOldNew/genoDic/gd
rm -f ExpOldNew/geno.words.sort

for (( R=$startR; R <= $maxR; R=R+$incrR ))
do
	for (( L=$startL; L <= $maxL; L=L+$incrL ))
	do
		java -jar ${BINDIR}/MetagGenerator.jar ExpOldNew/fakeM.fasta $1 $R $L

######################################### NEW #############################################

		dc_newInit=$(date -u +"%s")

		#New dict
		${BINDIR}/dic ExpOldNew/fakeM.fasta ExpOldNew/new/newDic

		dc_newEnd=$(date -u +"%s")

		dc_newTime=$(($dc_newEnd-$dc_newInit))

		echo "NEW DIC: $((dc_newTime / 60)) min $((dc_newTime % 60)) sec"


		fr_newInit=$(date -u +"%s")

		#New fragment
		${BINDIR}/frags ExpOldNew/genoDic/ ExpOldNew/new/ 100 32 32 ExpOldNew/FTest

		fr_newEnd=$(date -u +"%s")

		fr_newTime=$(($fr_newEnd-$fr_newInit))

		echo "NEW FRAG: $((fr_newTime / 60)) min $((fr_newTime % 60)) sec"

######################################### OLD #############################################

		# Take old time
		dc_oldInit=$(date -u +"%s")

		# find words and order
		#echo "${BINDIR}/words"
		${BINDIR}/words ExpOldNew/fakeM.fasta ExpOldNew/old/old.words.unsort
		#echo "${BINDIR}/sortWords 10000000 32"
		${BINDIR}/sortWords 10000000 32 ExpOldNew/old/old.words.unsort ExpOldNew/old/old.words.sort

		# Create hash table in disk
		#echo "${BINDIR}/w2hd Ds.words.sort Ds"
		${BINDIR}/w2hd ExpOldNew/old/old.words.sort ExpOldNew/old/oldDic

		dc_oldEnd=$(date -u +"%s")

		dc_oldTime=$(($dc_oldEnd-$dc_oldInit))

		echo "OLD DIC: $((dc_oldTime / 60)) min $((dc_oldTime % 60)) sec"

		fr_oldInit=$(date -u +"%s")

		#Old fragment
		#Hits
		${BINDIR}/hits ExpOldNew/old/oldDic ExpOldNew/genoDic/gd ExpOldNew/oldF.hits 1000 32
		#Sort hits
		${BINDIR}/sortHits 10000000 32 ExpOldNew/oldF.hits ExpOldNew/oldF.hits.sorted
		#FilterHits
		${BINDIR}/filterHits ExpOldNew/oldF.hits.sorted ExpOldNew/oldF.hits.sorted.filtered 32
		#FragHits
		${BINDIR}/FragHits ExpOldNew/fakeM.fasta $1 ExpOldNew/oldF.hits.sorted.filtered ExpOldNew/oldF.frags 32 100 32 1 f



		fr_oldEnd=$(date -u +"%s")

		fr_oldTime=$(($fr_oldEnd-$fr_oldInit))

		echo "OLD FRAG: $((fr_oldTime / 60)) min $((fr_oldTime % 60)) sec"

######################################## WRITE ############################################

		${BINDIR}/write $8 $R $L $dc_newTime $dc_oldTime
		${BINDIR}/write $9 $R $L $fr_newTime $fr_oldTime
	done
done

######################################## DELETE ###########################################
rm -f fakeM.fasta
rm -f ExpOldNew/new/newDic.metag.d2hP
rm -f ExpOldNew/new/newDic.metag.d2hR
rm -f ExpOldNew/new/newDic.metag.d2hW
rm -f ExpOldNew/old/old.words.sort
rm -f ExpOldNew/old/oldDic.d2hP
rm -f ExpOldNew/old/oldDic.d2hW
rm -rf ExpOldNew
