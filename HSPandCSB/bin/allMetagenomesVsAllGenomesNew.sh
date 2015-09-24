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


############################## TAKE BINARIES PATH
# Take bin directory
BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo "----Taking files"
echo "Metagenome files:"


############################## TAKE GENOME AND METAGENOME FILES
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


############################## CALCULATE GENOME DICTIONARIES
echo ""
echo "----Calculating dictionaries"

# Prepare workspace
if [[ ! -d intermediateFiles ]]; then
	mkdir intermediateFiles
	mkdir intermediateFiles/dictionaries
elif [[ ! -d intermediateFiles/dictionaries ]];	then
	mkdir intermediateFiles/dictionaries
fi

	cd intermediateFiles/dictionaries # Move to dictionaries folder

anyDictionaryCalculated=0

for ((i=0 ; i < ${#genoFiles[@]} ; i++))
do
	genoF=${genoFiles[$i]}

	# Check if dictionary already exist the dictionary
	didntExists=0
	if [[ ! -f ${genoF}.d2hP ]]; then
		if [[ ! -f ${genoF}.d2hW ]]; then
			if [[ ! -f ${genoF}-*.d2hP ]]; then
				if [[ ! -f ${genoF}-*.d2hW ]]; then
					didntExists=1
					anyDictionaryCalculated=1
					${BINDIR}/dictionary.sh ../../${genomeDir}/${genoF}.$EXT $WL &
					#echo "Dictionary created: ${genoF}"
				fi
			fi
		fi
	fi
	if [ $didntExists -eq 0 ]; then # 
		echo "Dictionary already exists: ${genoF}"
	fi
done

if [ $anyDictionaryCalculated -eq 1 ];	then
	echo "Waiting for the calculation of the dictionaries"

	for job in `jobs -p`
	do
    	#echo $job
    	wait $job
	done
fi

cd ../../ # Go back to init directory


############################## STUDY EACH {METAGENOME Vs GENOME} CASE
for ((i=0 ; i < ${#metaFiles[@]} ; i++))
do
	for ((j=0 ; j < ${#genoFiles[@]} ; j++))
	do
		${BINDIR}/metagenomeVsGenome.sh ${metaDir}/${metaFiles[$i]}.$EXT ${genomeDir}/${genoFiles[$j]}.$EXT ${L} ${S} ${WL}
	done
done