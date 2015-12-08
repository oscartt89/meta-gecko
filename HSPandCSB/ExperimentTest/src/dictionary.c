/* @file dictionary.c
 * @author Fernando Moreno Jabato <jabato@uma.es>
 * @description This program create a dictionary of a mono/multi sequence
 *    file of metagenome fasta format. 
 */

#include "dictionary.h"

int main(int ac, char** av){
 	if(ac!=3){
		fprintf(stderr,"USE: dictionary seqFile.IN dictionary.OUT\n");
		return -1;
	}

	// Read file and take kmers
	if((createDictionary(av[1],av[2])) < 0) return -1; 

	return 0;
}


/*
*/
int createDictionary(char *fIN, char *fOUT){
	// Variables
	wentry *words; // Array of words
	FILE *metag, *wDic, *pDic, *rDic;
	char c;

	char *fname;
	if((fname = (char*) malloc(sizeof(char)*MAX_FILE_LENGTH))==NULL){
		fprintf(stderr, "createDictionary:: Error allocating memory for file name.\n");
		return -1;
	}

	strcpy(fname,fOUT);


	// Space for array of words
	if((words = (wentry*) malloc(sizeof(wentry)*MAX_WORDS))==NULL){
		fprintf(stderr, "createDictionary:: Error initializing words array\n");
		free(words);
		return -1;
	}

	// Open input file
	if((metag = fopen(fIN,"rt"))==NULL){
		fprintf(stderr, "createDictionary:: Error opening metagenome file.\n");
		return -1;
	}

	// Open word dicitonary file
	if((wDic = fopen(strcat(fname,".metag.d2hW"),"wb"))==NULL){
		fprintf(stderr, "createDictionary:: Error opening words dictionary file.\n");
		return -1;
	}
	strcpy(fname,fOUT);

	// Open position dicitonary file
	if((pDic = fopen(strcat(fname,".metag.d2hP"),"wb"))==NULL){
		fprintf(stderr, "createDictionary:: Error opening position dictionary file.\n");
		return -1;
	}
	strcpy(fname,fOUT);

	// Open read dicitonary file
	if((rDic = fopen(strcat(fname,".metag.d2hR"),"wb"))==NULL){
		fprintf(stderr, "createDictionary:: Error opening reads dictionary file.\n");
		return -1;
	}
	free(fname);

	// READ SEQUENCES
	// Fasta files starts with a comment or read info line. Avoid it.
//	c = fgetc(metag);
//	while(feof(metag) & c != '\n')
//		c = fgetc(metag);

	// Check
	if(feof(metag)){
		fprintf(stderr, "createDictionary:: Error reading metagenome file. No sequence found.\n");
		return -1;
	}

	// Necessary variables
	unsigned long index = 0; // Index on sequence
	unsigned long inEntry = 0; // Length of the well formed sequence stored on the buffer
	uint64_t NW = 0; // Number of well done sequences
	wentry temp; // Buffer
	temp.seq = -1;

	// Read sequences
	c = fgetc(metag);
	while(!feof(metag)){
		if (!isupper(toupper(c))){ // Comment, empty or quality (+) line
			if(c=='>'){ // Comment line
				c = fgetc(metag);
				while (c != '\n')
					c = fgetc(metag);
				temp.seq++; // New sequence
				inEntry = 0; // Reset buffer length
				index=0; // Reset index

				if(NW > 0){
					// Store read kmers taken
					quickSort_W(words,0,NW); // Sort kmers
					writeDic(words,NW,wDic,pDic,rDic);
					NW = 0;
				}
			}

			c=fgetc(metag); // First char of next sequence
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
				inEntry=0; break;
		}
		index++;
		if(inEntry >= (unsigned long)WORD_SIZE){ // Full well formed sequence 
			temp.pos=index-WORD_SIZE;
			NW++;
			if(storeWord(words,&temp,NW) < 0) return -1;
		}
		c=fgetc(metag);
	}
	// Store buffered kmers
	if(NW > 0){
		quickSort_W(words,0,NW); // Sort kmers
		writeDic(words,NW,wDic,pDic,rDic);
	}

	free(words);
	fclose(metag);
	fclose(wDic);
	fclose(pDic);
	fclose(rDic);

	return 0;
}