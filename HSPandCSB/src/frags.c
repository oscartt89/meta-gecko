/* @file frags.c
 * @author Fernando Moreno Jabato <jabato@uma.es>
 * @description: this file contains de code to compare and create
 * 		common frags between a genome set and a metagenome using 
 * 		their dictionaries.
 */

#include "frags.h"

int main(int ac, char** av){
	if(ac!=4){
		fprintf(stderr,"USE: %s genomeSetFolder metagenomeFolder out\n",av[0]);
		return -1;
	}

	// Variables
	dictionaryG* dicGSet;
	dictionaryM* dicMSet;
	int numGenomes,numMetags;

	// Load genome set dictionaries
	if((numGenomes=readGenomeSet(av[1],dicGSet))<0) return -1;
	// Load metagenome set dictionaries
	if((numMetags=readMetagenomeSet(av[2],dicMSet))<0) return -1;

	fprintf(stdout, "%d\n", numMetags);

	return 0;
}