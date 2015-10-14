/* @file dictionary.c
 * @author Fernando Moreno Jabato <jabato@uma.es>
 * @description This program create a dictionary of a mono/multi sequence
 *    file of fasta format. 
 */

#include "dictionary.h"

int main(int ac, char** av){
 	if(ac!=3){
		fprintf(stderr,"USE: dictionary seqFile.IN dictionary.OUT\n");
		return -1;
	}

	// Variables
	wentry *words; // Array of words
	int numWords = -1;

	if((words = (wentry*) malloc(sizeof(wentry)*MAX_WORDS))==NULL){
		fprintf(stderr, "Error initializing words array\n");
		free(words);
		return -1;
	}

	// Read file and take kmers
	if((numWords = takeWords(words,av[1])) < 0) return -1; 
	// Sort words
	if(quickSort(words,0,numWords-1) < 0) return -1;

	return 0;
}