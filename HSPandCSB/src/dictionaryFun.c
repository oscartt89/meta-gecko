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


/* 	
 */
int quickSort(wentry* words, int left, int right){
	int j;

	if(left < right){
		// divide and conquer
		if((j = partition(words, left, right))<0) return -1;
		quickSort(words, left, j-1);
		quickSort(words, j+1, right);
	}
	return 0;
}


/*
 */
int partition(wentry* words, int left, int right){
	int i = left;
	int j = right+1;
	wentry *t;

	if((t = (wentry*) malloc(sizeof(wentry)))==NULL){
		fprintf(stderr, "Error allocating memory for auxiliar variable.\n");
		return -1;
	}

	// left sera el pivote
	// y contendra la mediana de left, right y (left+right)/2
	int mid = (int) (left+right)/2;

	if(wordComparator(&words[mid],&words[right]))
		SWAP(&words[mid],&words[right],t);

	if(wordComparator(&words[mid],&words[left]))
		SWAP(&words[mid],&words[left],t);

	if(wordComparator(&words[left],&words[right]))
		SWAP(&words[left],&words[right],t);

	while(1){
		do{
			++i;
		}while(!wordComparator(&words[i],&words[left]) && i <= right);

		do{
			--j;
		}while(wordComparator(&words[j],&words[left]) && j >= left);

		if( i >= j ) break;

		SWAP(&words[i],&words[j],t);
	}

	SWAP(&words[left],&words[j],t);

	return j;
}


/* This function is used to compare two wentry instances. The criterion
 * used is:
 * 		1 - Compare sequences (alphabetically).
 * 		2 - Compare position on sequence.
 * @param w1 word to be compared.
 * @param w2 word to be compared.
 * @return a positive number if w1 is greater than w2, a negative number
 * 		if w2 is greater than w1 and zero if both are equal.
 */
int wordComparator(wentry* w1,wentry* w2){
	if(w1->w.b[0] > w2->w.b[0]) return 1;
	else if(w1->w.b[0] < w2->w.b[0]) return -1;

	if(w1->w.b[1] > w2->w.b[1]) return 1;
	else if(w1->w.b[1] < w2->w.b[1]) return -1;

	if(w1->w.b[2] > w2->w.b[2]) return 1;
	else if(w1->w.b[2] < w2->w.b[2]) return -1;

	if(w1->w.b[3] > w2->w.b[3]) return 1;
	else if(w1->w.b[3] < w2->w.b[3]) return -1;

	if(w1->w.b[4] > w2->w.b[4]) return 1;
	else if(w1->w.b[4] < w2->w.b[4]) return -1;

	if(w1->w.b[5] > w2->w.b[5]) return 1;
	else if(w1->w.b[5] < w2->w.b[5]) return -1;

	if(w1->w.b[6] > w2->w.b[6]) return 1;
	else if(w1->w.b[6] < w2->w.b[6]) return -1;

	if(w1->w.b[7] > w2->w.b[7]) return 1;
	else if(w1->w.b[7] < w2->w.b[7]) return -1;

	if(w1->pos > w2->pos) return 1;
	return 0;
}


/*
 */
void writeDic(wentry* words,int numWords,FILE *wDic,FILE *pDic,FILE *rDic){
	//Variables
	hashentry he;
	read re;
	int i;
	uint64_t loc;

	// Read info
	re.readIndex = words[0].seq;
	re.num = 0;
	re.pos = ftell(wDic);
	// First word
	memcpy(&he.w.b[0],&words[0].w.b[0],8);
	he.pos=ftell(pDic);
	he.num=0;

	// Write all
	for(i=0 ; i<numWords; ++i){
		loc = words[i].pos;
		if(wordcmp(&he.w.b[0],&words[i].w.b[0],WORD_SIZE)!=0){ // New sequence
			fwrite(&he,sizeof(hashentry),1,wDic); // Write kmer info
			memcpy(&he.w.b[0],&words[i].w.b[0],8); // Tke new "current" kmer
			he.pos=ftell(pDic); // Take position of kmer locations on pDic
			he.num=0;
			re.num++;
		}

		// Write new location
		fwrite(&loc,sizeof(uint64_t),1,pDic);
		he.num++;
	}

	// Write kmer stored on buffer
	fwrite(&he,sizeof(hashentry),1,wDic);
	re.num++;

	// Write Read information
	fwrite(&re,sizeof(read),1,rDic);
}


/*
 */
int wordcmp(unsigned char *w1, unsigned char*w2, int n) {
	int i;
	for (i=0;i<n;i++) {
		if (w1[i]<w2[i]) return -1;
		if (w1[i]>w2[i]) return +1;
	}
	return 0;
}


inline void SWAP(wentry *w1,wentry *w2,wentry *t){
	memcpy(t,w1,sizeof(wentry));
	memcpy(w1,w2,sizeof(wentry));
	memcpy(w2,t,sizeof(wentry));
}