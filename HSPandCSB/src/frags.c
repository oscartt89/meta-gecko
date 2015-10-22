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
	wentry** dictionaries;
	int* lengths;

	// Load genome set dictionaries
	if(readGenomeSet(av[1],dictionaries,lengths)<0) return -1;
	

	return 0;
}