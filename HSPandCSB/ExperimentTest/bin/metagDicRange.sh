#!/bin/bash 

# This script launch a set of experiments where, each time, a metagenome
# is artificially generated using a given genome. After that a dictionary
# of the metagenome is created using the old gecko version and the new 
# geckoMGV saving times spended on this actions in a file given.

BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


if [ $# != 8 ]; then
   echo " ==== ERROR ... you called this script inappropriately."
   echo ""
   echo "   usage: $0 genome.fasta startR startL maxR maxL incrR incrL dataFile"
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

for (( R=$startR; R <= $maxR; R=R+$incrR ))
do
	for (( L=$startL; L <= $maxL; L=L+$incrL ))
	do
		java -jar ${BINDIR}/MetagGenerator.jar fakeM.fasta $1 $R $L

######################################### NEW #############################################

		newInit=$(date -u +"%s")

		#New dict
		${BINDIR}/dic fakeM.fasta newDic

		newEnd=$(date -u +"%s")

		newTime=$(($newEnd-$newInit))

		echo "NEW: $((newTime / 60)) min $((newTime % 60)) sec"

######################################### OLD #############################################

		# Take old time
		oldInit=$(date -u +"%s")

		# find words and order
		#echo "${BINDIR}/words"
		${BINDIR}/words fakeM.fasta old.words.unsort
		#echo "${BINDIR}/sortWords 10000000 32"
		${BINDIR}/sortWords 10000000 32 old.words.unsort old.words.sort

		# Create hash table in disk
		#echo "${BINDIR}/w2hd Ds.words.sort Ds"
		${BINDIR}/w2hd old.words.sort oldDic

		oldEnd=$(date -u +"%s")

		oldTime=$(($oldEnd-$oldInit))

		echo "OLD: $((oldTime / 60)) min $((oldTime % 60)) sec"

######################################## WRITE ############################################

		${BINDIR}/write $8 $R $L $newTime $oldTime
	done
done

######################################## DELETE ###########################################
rm -f fakeM.fasta
rm -f oldDic.d2hP
rm -f oldDic.d2hW
rm -f old.words.sort
rm -f newDic.metag.d2hP
rm -f newDic.metag.d2hR
rm -f newDic.metag.d2hW
