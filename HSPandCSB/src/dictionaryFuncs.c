/* @author Fernando Moreno Jabato <jabato@uma.es>
 * @date 30-Sept-2015
 * @description This file contains necessary functions for 
 *     handle and create dictionaries.
 * @license all rights reserved to BitLAB (http://www.bitlab-es.com/bitlab/)
 *     and to author. 
 */
#include "dictionary.h"



void createDictionary(char *sequence,int wordLength){
	int seqLength = strlen(sequence);
	int numWords = 0, pos = 0;

//	if(seqLength < wordLength)
//		writeDictionary(NULL);
	
	struct WordList unsortedHead = ;
		// Init word list
		unsorted.nextWord = NULL;
		initWord(unsorted->word,wordLength);

	struct WordList *last = &unsorted;
	struct WordList temp;

	// First word
	unsorted->word->sequence[0] = sequence[0];
	unsorted->word->sequence[1] = sequence[1];
	unsorted->word->sequence[2] = sequence[2];
	unsorted->word.pos = pos;

	// Read words
	while(pos < (seqLength - 2)){
		// Init temporal instance
		temp = 

		// Link words

	}
}



int initWord(struct Word *word, int wordLength){
	if((word.sequence = (char*) malloc(sizeof(char)*wordLength)==NULL){
		fprintf(stderr, "Error reserving space from word on disk. [Length:%d]\n", wordLength);
		return -1;
	}
	return 1;
}

struct WordList createWordListInstance(int wordLength){
	struct WordList *wl = (struct *WordList) malloc(sizeof(struct WordList));

	initWord(wl->)
}