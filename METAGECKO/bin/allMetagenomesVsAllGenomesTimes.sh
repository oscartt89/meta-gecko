#!/usr/bin/env bash
metaDir=$1
genomeDir=$2
L=$3
S=$4
EXT=$5
GEXT=$6
GDFile="genomeDictTimes.txt"
MDFile="metagenomeDictTimes.txt"
FFile="fragmentsTimes.txt"

metagenomes=()
x=0

genomes=()
y=0

if [ $# != 7 ]; then
	echo "***ERROR*** Use: $0 metagenomesDir genomesDir L(200) S(40) MetagExtension GenoExtension prefix"
	exit -1
fi

BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

for elem in $(ls -d ${metaDir}/*.$EXT | awk -F "/" '{print $NF}' | awk -F ".$EXT" '{print $1}')
do
	metagenomes[$x]=$elem
	x=`expr $x + 1`
	#echo "X: $elem"
done

for elem in $(ls -d ${genomeDir}/*.$GEXT | awk -F "/" '{print $NF}' | awk -F ".$GEXT" '{print $1}')
do
	genomes[$y]=$elem
	y=`expr $y + 1`
	#echo "X: $elem"
done

# Create necessary directories
mkdir -p ${metaDir}/dictionaries
mkdir -p ${genomeDir}/dictionaries
mkdir -p ${genomeDir}/complementaries
mkdir -p fragments

# Write info headers if it's necessary
if [ ! -f  $GDFile ]
then
	${BINDIR}/write $GDFile Genome DictTime Day
fi

if [ ! -f  $MDFile ]
then
	${BINDIR}/write $MDFile Metagenome DictTime Day
fi

if [ ! -f  $FFile ]
then
	${BINDIR}/write $FFile Fragments FragsTime Day
fi

# Write all genome dictionaries
for ((k=0 ; k < ${#genomes[@]} ; k++))
do
	genome=${genomes[$k]}
	if [[ ! -f $genomeDir/dictionaries/${genome}.d2hP ]];	then
		echo " Writting genome dict: ${genome}.${GEXT}"
		
		# Take init time
		genoDTime_Init=$(date -u +"%s")

		# Find words and order
		${BINDIR}/words $genomeDir/${genome}.${GEXT} $genomeDir/dictionaries/${genome}.words.unsort
		${BINDIR}/sortWords 10000000 32 $genomeDir/dictionaries/${genome}.words.unsort $genomeDir/dictionaries/${genome}.words.sort

		# Create hash table in disk
		${BINDIR}/w2hd $genomeDir/dictionaries/${genome}.words.sort $genomeDir/dictionaries/${genome}

		genoDTime_End=$(date -u +"%s")

		genoDTime=$(($genoDTime_End-$genoDTime_Init))

		echo "Time: $((genoDTime / 60)) min $((genoDTime % 60)) sec"

		# Write times
		NOW=$(date +"%m-%d-%Y")
		${BINDIR}/write $GDFile ${genome} ${genoDTime} ${NOW} 
	fi

	# Now with reverse seq
	if [[ ! -f $genomeDir/dictionaries/${genome}-revercomp.d2hP ]];	then
		${BINDIR}/reverseComplement $genomeDir/${genome}.${GEXT} $genomeDir/complementaries/${genome}-revercomp.${GEXT}

		echo " Writting genome dict: ${genome}-revercomp.${GEXT}"
		
		# Take init time
		genoDTime_Init=$(date -u +"%s")

		# Find words and order
		${BINDIR}/words $genomeDir/complementaries/${genome}-revercomp.${GEXT} $genomeDir/dictionaries/${genome}-revercomp.words.unsort
		${BINDIR}/sortWords 10000000 32 $genomeDir/dictionaries/${genome}-revercomp.words.unsort $genomeDir/dictionaries/${genome}-revercomp.words.sort

		# Create hash table in disk
		${BINDIR}/w2hd $genomeDir/dictionaries/${genome}-revercomp.words.sort $genomeDir/dictionaries/${genome}-revercomp

		genoDTime_End=$(date -u +"%s")

		genoDTime=$(($genoDTime_End-$genoDTime_Init))

		echo "Time: $((genoDTime / 60)) min $((genoDTime % 60)) sec"

		# Write times
		NOW=$(date +"%m-%d-%Y")
		${BINDIR}/write $GDFile ${genome}-revercomp ${genoDTime} ${NOW}
	fi
done

for ((i=0 ; i < ${#metagenomes[@]} ; i++))
do
	seqX=${metagenomes[$i]}

	# Create metagenome dictionarie
	if [[ ! -f $metaDir/dictionaries/${seqX}.d2hP ]];	then
		echo " Writting metagenome dict: ${seqX}.${EXT}"

		metagDTime_Init=$(date -u +"%s")

		#New dict metagenome
		${BINDIR}/dict $metaDir/${seqX}.${EXT} $metaDir/dictionaries/${seqX} 32

		metagDTime_End=$(date -u +"%s")

		metagDTime=$(($metagDTime_End-$metagDTime_Init))

		echo "Time: $((metagDTime / 60)) min $((metagDTime % 60)) sec"

		# Write times
		NOW=$(date +"%m-%d-%Y")
		${BINDIR}/write $MDFile ${seqX} ${metagDTime} ${NOW}
	fi

	# Generate fragments
	for ((j=0 ; j < ${#genomes[@]} ; j++))
	do
		seqY=${genomes[$j]}
		if [[ ! -f fragments/${seqX}-${seqY}.frags ]];	then
			if [[ ! -f fragments/${seqX}-${seqY}-f.frags ]];	then
				echo "Writting new fragments for: ${seqX} - ${seqY}"

				fr_newInit=$(date -u +"%s")

				${BINDIR}/frag $metaDir/dictionaries/${seqX} $metaDir/${seqX}.${EXT} $genomeDir/dictionaries/${seqY} $genomeDir/${seqY}.${GEXT} fragments/${seqX}-${seqY}-f $S $L f $7

				fr_newEnd=$(date -u +"%s")

				fr_newTime=$(($fr_newEnd-$fr_newInit))

				echo "Time: $((fr_newTime / 60)) min $((fr_newTime % 60)) sec"
				# Write times
				NOW=$(date +"%m-%d-%Y")
				${BINDIR}/write $FFile ${seqX}-${seqY} ${fr_newTime} ${NOW}
			fi
			# Reverse
			if [[ ! -f fragments/${seqX}-${seqY}-revercomp.frags ]];	then
				echo "Writting new fragments for: ${seqX} - ${seqY}-revercomp"

				fr_newInit=$(date -u +"%s")

				${BINDIR}/frag $metaDir/dictionaries/${seqX} $metaDir/${seqX}.${EXT} $genomeDir/dictionaries/${seqY}-revercomp $genomeDir/complementaries/${seqY}-revercomp.${GEXT} fragments/${seqX}-${seqY}-revercomp $S $L r $7

				# Fixe fragments
				${BINDIR}/fixeReverseFrags $genomeDir/${seqY}.${GEXT} $genomeDir/complementaries/${seqY}-revercomp.${GEXT} fragments/${seqX}-${seqY}-revercomp.frags

				fr_newEnd=$(date -u +"%s")

				fr_newTime=$(($fr_newEnd-$fr_newInit))

				echo "Time: $((fr_newTime / 60)) min $((fr_newTime % 60)) sec"
				# Write times
				NOW=$(date +"%m-%d-%Y")
				${BINDIR}/write $FFile ${seqX}-${seqY}-revercomp ${fr_newTime} ${NOW}
			fi

			# Combine fragment
			if [[ -f fragments/${seqX}-${seqY}-revercomp.frags && -f fragments/${seqX}-${seqY}-f.frags ]];	then
				echo "Combining forward and reverse: ${seqX}-${seqY}"

				comb_Init=$(date -u +"%s")

				${BINDIR}/combineFragments fragments/${seqX}-${seqY}-f.frags fragments/${seqX}-${seqY}-revercomp.frags fragments/${seqX}-${seqY}.frags 1

				comb_End=$(date -u +"%s")
				comb_Time=$(($comb_End-$comb_Init))
				${BINDIR}/write $FFile Combine:${seqX}-${seqY} ${comb_Time}
			fi
		fi
	done
done

