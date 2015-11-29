#include "dictionary.h"


/* This function is used to shift bits in a unsigned char array
 *	@param w: word structure where char array to be shifted is stored.
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
int storeWord(wentry* wArr,wentry *word,int length){
	// Take enough memory on array
	if((length % MAX_WORDS) == 0){ // Limit of words
		if(realloc(wArr,sizeof(wentry*)*(length + MAX_WORDS)) == NULL){
			fprintf(stderr, "Error reallocating the words array. New size: %d\n", length);
			free(wArr);
			return -1;
		}
		fprintf(stdout, "Memoria realocada\t%d\n",length);
	}

	// Copy the new word
	memcpy(&wArr[length-1],word,sizeof(wentry));

	return 0;
}


/* This function is used to write an array of words in a new dictionary.
 *	@param words: array of wentry instances that will be written on the dictionary.
 *	@param numWords: number of instances on words array.
 *	@param wDic: words dictionary file pointer.
 *	@param pDic: positions dictionary file pointer.
 *	@param rDic: reads dictionary file pointer.
 */
void writeDic(wentry* words,int numWords,FILE *wDic,FILE *pDic,FILE *rDic){
	//Variables
	hashentry he;
	Read re;
	int i;

	// Read info
	re.readIndex = words[0].seq;
	re.num = 0;
	re.pos = ftell(wDic);
	// First word
	memcpy(&he.w.b[0],&words[0].w.b[0],BYTES_IN_WORD);
	he.pos=ftell(pDic);
	he.num=0;

	// Write all
	for(i=0 ; i<numWords; ++i){
		if(wordcmp(&he.w.b[0],&words[i].w.b[0],BYTES_IN_WORD)!=0){ // New sequence
			fwrite(&he,sizeof(hashentry),1,wDic); // Write kmer info
			memcpy(&he.w.b[0],&words[i].w.b[0],BYTES_IN_WORD); // Tke new "current" kmer
			he.pos=ftell(pDic); // Take position of kmer locations on pDic
			he.num=0;
			re.num++;
		}

		// Write new location
		fwrite(&words[i].pos,sizeof(uint64_t),1,pDic);
		he.num++;
	}

	// Write kmer stored on buffer
	fwrite(&he,sizeof(hashentry),1,wDic);
	re.num++;

	// Write Read information
	fwrite(&re,sizeof(Read),1,rDic);
}