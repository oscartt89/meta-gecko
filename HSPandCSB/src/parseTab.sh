# parseTab

# This program takes gene features files in fasta format and parses it into a bin file.




folderIn=$1 # folder with geneFeatures
output=$2 # folder output files

path=../bin


if [ $# -ne 2 ]; then
	echo "***ERROR*** Use: $0 folderIn FolderOut"
	exit -1
fi

array=()
x=0
for elem in $(ls -d $1/*.geneFeature | awk -F "/" '{print $NF}'|awk -F ".geneFeature" '{print $1}' )
#for elem in $folder/*
do
	#if [[ -f $elem ]]; then
		array[$x]=$elem
		x=`expr $x + 1`
		echo "X: $elem"
	#fi
done

for ((i=0 ; i < ${#array[@]} ; i++))
	do
		grep ">" $1/${array[$i]}.geneFeature | sed -r 's/>lcl\|//g' | sed -r 's/\[gene=//g' | sed -r 's/\[locus_tag=//g' | sed -r 's/\[location=//g' | sed -r 's/<//g' | sed -r 's/>//g' | sed -r 's/\]//g' > $output/${array[$i]}.tmp
		$path/parseGeneFeatureTab $output/${array[$i]}.tmp $output/${array[$i]}.gene $output/${array[$i]}.gene.csv
		rm $output/${array[$i]}.tmp
	done





