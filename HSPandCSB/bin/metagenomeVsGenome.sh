#!/usr/bin/env bash

# Take variables
metagenome=$1
genome=$2
length=$3
similarity=$4
wordLength=$5

# Check bad call error
if [ $# -lt 5 ]; then
   echo " ==== ERROR ... you called this script inappropriately."
   echo ""
   echo "   usage:  $0 metagenome genome lenght similarity WL"
   echo ""
   exit -1
fi