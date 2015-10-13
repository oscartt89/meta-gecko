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
	wentry **words; // Array of words

	if((words = malloc(sizeof(wentry*)*MAX_WORDS))==NULL){
		fprintf(stderr, "Error initializing words array\n");
		free(words);
		return -1;
	}

	// Read file and take kmers
	if(takeWords(words,av[1])<0) return -1; 

	return 0;
}


/* This method implements the process of open the file, read the
 * sequences and take the words.
 * @param words is an array of words.
 * @param IN is the path of the sequence file.
 * @return the number of elements stored on words if everything 
 * 		was OK and a negative number in other cases.
 */
int takeWords(wentry **words, char *IN){
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
}//


