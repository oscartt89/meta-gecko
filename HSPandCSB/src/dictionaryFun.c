#include "dictionary.h"

/*
 */
void shift_word(word * w){
	int i;
	for(i=0;i<BYTES_IN_WORD-1;i++){
		w->b[i]<<=2;
		w->b[i]|=(w->b[i+1]>>6);
	}
	w->b[BYTES_IN_WORD-1]<<=2;
}


/*
 */
int storeWord(wentry** wArr,wentry *word,int newLength){
	// Take enough memory on array
	if(newLength == 1){ // Initialize array
		if((wArr = malloc(sizeof(wentry*)*newLength))==NULL){
			fprintf(stderr, "Error initializing words array\n");
			free(wArr);
			return -1;
		}
	}else{
		if(realloc(wArr,sizeof(wentry*)*newLength) == NULL){
			fprintf(stderr, "Error reallocating the words array. New size: %d\n", newLength);
			free(wArr);
			return -1;
		}
	}

	// Take memory for new word
	if((wArr[newLength-1] = (wentry*) malloc(sizeof(wentry)))==NULL){
		fprintf(stderr, "Error allocating space for the new word\n");
		free(wArr);
		return -1;
	}

	fprintf(stdout, "Ready to copy memory.\n");


	// Copy the new word
	memcpy(wArr[newLength-1],word,sizeof(wentry));

	return 0;
}