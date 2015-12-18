#!/bin/bash 

# This script launch a set of experiments where, each time, a metagenome
# is artificially generated using a given genome. After that a dictionary
# of the metagenome is created using the old gecko version and the new 
# geckoMGV saving times spended on this actions in a file given.
# WARNING! : actually the names of files and directories aren't calculeted. 
# 		Replace it directories for your's

BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

maxR=1000000
maxL=1000

for (( R=10000; R <= $maxR; R=R+10000 ))
do
	for (( L=100; L <= $maxL; L=L+100 ))
	do
		java -jar ${BINDIR}/MetagGenerator.jar ../ExperimentoSecuencial/Ds.fasta ../Genome/Ds-W501-2L.fasta $R $L


		newInit=$(date -u +"%s")

		#New dict
		${BINDIR}/dic ../ExperimentoSecuencial/Ds.fasta ../ExperimentoSecuencial/DsNew

		newEnd=$(date -u +"%s")

		newTime=$(($newEnd-$newInit))

		echo "NEW: $((newTime / 60)) min $((newTime % 60)) sec"



		# Take old time
		oldInit=$(date -u +"%s")

		# find words and order
		#echo "${BINDIR}/words"
		${BINDIR}/words ../ExperimentoSecuencial/Ds.fasta ../ExperimentoSecuencial/Ds.words.unsort
		#echo "${BINDIR}/sortWords 10000000 32"
		${BINDIR}/sortWords 10000000 32 ../ExperimentoSecuencial/Ds.words.unsort ../ExperimentoSecuencial/Ds.words.sort

		# Create hash table in disk
		#echo "${BINDIR}/w2hd Ds.words.sort Ds"
		${BINDIR}/w2hd ../ExperimentoSecuencial/Ds.words.sort ../ExperimentoSecuencial/Ds

		oldEnd=$(date -u +"%s")

		oldTime=$(($oldEnd-$oldInit))

		echo "OLD: $((oldTime / 60)) min $((oldTime % 60)) sec"
		${BINDIR}/write ../ExperimentoSecuencial/data.txt $R $L $newTime $oldTime
	done
done


