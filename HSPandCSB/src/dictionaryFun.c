#include "dictionary.h"

/*
 */
void shift_word(word* w){
	int i;
	for(i=0;i<BYTES_IN_WORD-1;i++){
		w->b[i]<<=2;
		w->b[i]|=(w->b[i+1]>>6);
	}
	w->b[BYTES_IN_WORD-1]<<=2;
}


/* This method is used to store a wentry instance into an array of wentry. 
 * @param wArr is the array of wentry** where the new wentry will be stored.
 * @param word is the new word to be stored.
 * @param length is the new length of the array.
 * @return 0 if the process ends correctly. If something got wrong return a
 * 		negative value.
 */
int storeWord(wentry** wArr,wentry *word,int length){
	// Take enough memory on array
	if((length % MAX_WORDS) == 0){ // Limit of words
		if(realloc(wArr,sizeof(wentry*)*(length + MAX_WORDS)) == NULL){
			fprintf(stderr, "Error reallocating the words array. New size: %d\n", length);
			free(wArr);
			return -1;
		}
		fprintf(stdout, "Memoria realocada\t%d\n",length);
	}

	// Take memory for new word
	if((wArr[length-1] = (wentry*) malloc(sizeof(wentry)))==NULL){
		fprintf(stderr, "Error allocating space for the new word\n");
		free(wArr);
		return -1;
	}

	// Copy the new word
	memcpy(wArr[length-1],word,sizeof(wentry));

	return 0;
}