/* @author Fernando Moreno Jabato <jabato@uma.es>
 * @date 30-ept-2015
 * @description This file contains necessary implementation for create metagenome
 *     dictionaries from a file given.
 * @license all rights reserved to BitLAB (http://www.bitlab-es.com/bitlab/)
 *     and to author. 
 */
#include "metag.h"

int main(int ac, char **av){
	// Variables
	FILE *metagenome, *dictionary, *tempDic;
	int wordLength;
	char *read;

	// Check bad call error
	if(ac != 4){
		fprintf(stderr, "Error: bad call error. Use: %s metagenome outputFile wordLength\n", av[0]);
		return -1;
	}

	// Open streams
	if((metagenome = fopen(av[1],"rt")) == NULL){
		fprintf(stderr, "Error opening metagenome file. [%s]\n", av[1]);
		return -1;
	}

	if((dictionary = fopen(av[2],"wt")) == NULL){
		fprintf(stderr, "Error opening output file. [%s]\n", av[2]);
		return -1;
	}

	// Allocate read space on disk
	if((read = (char*) malloc(sizeof(char)*MAXREADLENGTH))==NULL){
		fprintf(stderr, "Error: problem reserving space from reads.\n");
		return -1;
	}

	// Take word length
	wordLength = atoi(av[3]);

	// Start to read sequences
	while(takeRead(metagenome,read)){
		tempDic = createDictionary(read,wordLength);
	}

} // END MAIN