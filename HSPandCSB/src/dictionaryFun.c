#include "dictionary.h"


/* This method implements the process of open the file, read the
 * sequences and take the words.
 * @param words is an array of words.
 * @param IN is the path of the sequence file.
 * @return the number of elements stored on words if everything 
 * 		was OK and a negative number in other cases.
 */
int takeWords(wentry *words, char *IN){
 	// Variables
	FILE *f;
	char c;

	// Open sequence file
	if((f = fopen(IN,"rt"))==NULL){
		fprintf(stderr, "Error openening sequence file\n");
		return -1;
	}

	// Fasta files starts with a comment. Avoid it.
	c = fgetc(f);
	while(feof(f) & c != '\n')
		c = fgetc(f);

	// Check
	if(feof(f)){
		fprintf(stderr, "Error reading sequence file. Not sequence found.\n");
		return -1;
	}

	// Start to read sequence
	// Prepare workspace
	unsigned long index = 0;
	unsigned long inEntry = 0; // Length of the well formed sequence stored on the buffer
	unsigned long NW = 0; // Number of well done sequences
		unsigned long totC = 0; // Total of characters read
		unsigned long NoACGT = 0; // Total of NoACGT character read on sequence
		unsigned long NoSeq = 0; // Number of no-sequence lines
	wentry temp; // Buffer
	temp.seq = 0;

	// Start to read
	c = fgetc(f);
	while(!feof(f)){
		if (!isupper(toupper(c))){ // Comment, empty or quality (+) line
			if(c=='>'){ // Comment line
				c = fgetc(f);
				while (c != '\n')
					c = fgetc(f);
				temp.seq++; // New sequence
				inEntry = 0; // Reset buffer length
				index++; //
			}
			NoSeq++;
			c=fgetc(f);
			continue;
		}
		shift_word(&temp.w); // Shift bits sequence
		// Add new nucleotid
		switch (c) {
			case 'A': // A = 00 
				inEntry++;
				break;
			case 'C': // C = 01
				temp.w.b[BYTES_IN_WORD-1]|=1;
				inEntry++;
				break;
			case 'G': // G = 10
				temp.w.b[BYTES_IN_WORD-1]|=2;
				inEntry++;
				break;
			case 'T': // T = 11
				temp.w.b[BYTES_IN_WORD-1]|=3;
				inEntry++;
				break;
			default : // Bad formed sequence
				inEntry=0; NoACGT++; break;
		}
		index++;
		totC++;
		if(inEntry >= (unsigned long)WORD_SIZE){ // Full well formed sequence 
			temp.pos=index-WORD_SIZE;
			NW++;
			if(storeWord(words,&temp,NW) < 0) return -1;
		}
		c=fgetc(f);
	}

	fclose(f); // Close input stream

	return NW;
}


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
	int  i, j;
	wentry* pivot;
	wentry* temp;

	if((pivot = (wentry*) malloc(sizeof(wentry)))==NULL){
		fprintf(stderr, "Error allocating memory for pivot.\n");
		return -1;
	}
	if((temp = (wentry*) malloc(sizeof(wentry)))==NULL){
		fprintf(stderr, "Error allocating memory for auxiliar variable.\n");
		return -1;
	}

	memcpy(pivot,&words[left],sizeof(wentry));
	i = left;
	j = right+1;
		
	while(i >= j){
		do ++i; while(wordComparator(&words[i],pivot)<=0 && i <= right);
		do --j; while(wordComparator(&words[j],pivot) > 0);
		if(i < j){
			memcpy(temp,&words[i],sizeof(wentry));
			memcpy(&words[i],&words[j],sizeof(wentry));
			memcpy(&words[j],temp,sizeof(wentry));
		}
	}

	memcpy(temp,&words[left],sizeof(wentry));
	memcpy(&words[left],&words[j],sizeof(wentry));
	memcpy(&words[j],temp,sizeof(wentry));

	free(temp);
	free(pivot);
	
	return j;
}


/* This function is used to compare two wentry instances. The criterion
 * used is:
 * 		1 - Compare sequence index.
 * 		2 - Compare sequences (alphabetically).
 * 		3 - Compare position on sequence.
 * @param w1 word to be compared.
 * @param w2 word to be compared.
 * @return a positive number if w1 is greater than w2, a negative number
 * 		if w2 is greater than w1 and zero if both are equal.
 */
int wordComparator(wentry* w1,wentry* w2){
	if(w1->seq > w2->seq) return 1;
	else if(w1->seq < w2->seq) return -1;

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