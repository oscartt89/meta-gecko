#!/bin/bash 

FL=1000   # frequency limit
BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ $# != 3 ]; then
   echo " ==== ERROR ... you called this script inappropriately."
   echo ""
   echo "   usage:  $0 seqXName.fasta WL outDir"
   echo ""
   exit -1
fi

WL=$2     # wordSize
seqName=$(basename "$1")
extension="${seqName##*.}"
seqName="${seqName%.*}"

# find words and order
echo "${BINDIR}/words $1 $3/${seqName}.words.unsort"
${BINDIR}/words $1 $3/${seqName}.words.unsort
echo "${BINDIR}/sortWords 10000000 32 $3/${seqName}.words.unsort $3/${seqName}.words.sort"
${BINDIR}/sortWords 10000000 32 $3/${seqName}.words.unsort $3/${seqName}.words.sort

# Create hash table in disk
echo "${BINDIR}/w2hd $3/${seqName}.words.sort $3/${seqName} ${WL}"
${BINDIR}/w2hd $3/${seqName}.words.sort $3/${seqName}

