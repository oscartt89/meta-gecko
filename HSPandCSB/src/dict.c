/*
 * @author Fernando Moreno Jabato <jabato@uma.es>
 * @description This file encodes the workflow of GECKO for create
 *    metagenome dictionaries.
 * @licence all rights reserved to the author and BitLAB group (University
 *    of Malaga).
 */
#include "dict.h"

/**
 * This main function encodes the workflow of metagenome dictionary
 * creation. The workflow is:
 *   1 - Read a char of a read.
 *   2 - IF sequence is not long and good enough GO TO 1.
 *   3 - Store sequence (word) on words buffer.
 *   4 - IF bufer is full, write buffer sorted on intermediate files.
 *   5 - IF there are more Reads or chars unread on current Read GO 
 *       TO 1.
 *   6 - Write buffered words on intermediate files after sort it.
 *   7 - Read intermediate files word per word and write final dictionary
 *       files.
 */
int main(int ac, char** av){
	// Check arguments
	if(ac!=4){
		fprintf(stderr, "Bad call error.\nUSE: dict metag.IN dictName WL\n");
		return -1;
	}else if(atoi(av[3])%4 != 0){
		fprintf(stderr, "Error WL must be a 4 multiple and it's \"%s\".\n", av[3]);
		return -1;
	}

	// Variables
	int WL = atoi(av[3]); // Word length
	BYTES_IN_WORD = WL / 4;
	wentry *buffer; //Buffer of read words
	FILE *metag; // Input files
	FILE *wDic,*pDic;  // Output files
	FILE *bIndx, *wrds; // Intermediate files
	uint64_t /*readW = 0,*/ wordsInBuffer = 0, maxWordsStored = 0; // Absolute and buffer read words and control varaible
	uint32_t numBuffWritten = 0;
	char *fname;

	// Allocate necessary memory
	// Memory for buffer
	if((buffer = (wentry*) malloc(sizeof(wentry)*BUFFER_LENGTH))==NULL){
		fprintf(stderr, "Error allocating memory for words buffer.\n");
		return -1;
	}

	// Memory for file names handler
	if((fname = (char*) malloc(sizeof(char)*MAX_FILE_LENGTH))==NULL){
		fprintf(stderr, "Error allocating memory for file names handler.\n");
		return -1;
	}

	// Open current necessary files
	// Open metagenome file
	if((metag = fopen(av[1],"rt"))==NULL){
		fprintf(stderr, "Error opening metagenome file.\n");
		return -1;
	}

	//Open intermediate files
	strcpy(fname,av[2]); // Copy outDic name
	if((bIndx = fopen(strcat(fname,".bindx"),"wb"))==NULL){
		fprintf(stderr, "Error opening buffer index file.\n");
		return -1;
	}
	// Write first line of buffer index file (WordLength)
	fwrite(&WL,sizeof(int),1,bIndx);


	// Open words repo
	strcpy(fname,av[2]);
	if((wrds = fopen(strcat(fname,".wrds"),"wb"))==NULL){
		fprintf(stderr, "Error opening words repository.\n");
		return -1;
	}

	// Free unnecessary memory
	free(fname);

	// START WORKFLOW
	// Read metagenome file
	// Necessary variables
	uint64_t seqPos = 0, crrSeqL = 0;
	char c; // Auxiliar -> Read char per char
	wentry temp; // Auxiliar word variable
	if((temp.w.b = (unsigned char*) malloc(sizeof(unsigned char)*BYTES_IN_WORD))==NULL){
		fprintf(stderr, "Error allocating memory for temp.\n");
		return -1;
	}
	temp.seq = 0;

	// Start to read
	c = fgetc(metag);
	while(!feof(metag)){
		// Check if it's a special line
		if(!isupper(toupper(c))){ // Comment, empty or quality (+) line
			if(c=='>'){ // Comment line
				c = fgetc(metag);
				while(c != '\n') // Avoid comment line
					c = fgetc(metag);
				temp.seq++; // New sequence
				crrSeqL = 0; // Reset buffered sequence length
				seqPos = 0; // Reset index
			}
			c=fgetc(metag); // First char of next sequence
			continue;
		}
		shift_word(&temp.w); // Shift bits sequence
		// Add new nucleotid
		switch (c) {
			case 'A': // A = 00 
				crrSeqL++;
				break;
			case 'C': // C = 01
				temp.w.b[BYTES_IN_WORD-1]|=1;
				crrSeqL++;
				break;
			case 'G': // G = 10
				temp.w.b[BYTES_IN_WORD-1]|=2;
				crrSeqL++;
				break;
			case 'T': // T = 11
				temp.w.b[BYTES_IN_WORD-1]|=3;
				crrSeqL++;
				break;
			default : // Bad formed sequence
				crrSeqL = 0; break;
		}
		seqPos++;
		if(crrSeqL >= (uint64_t)WL){ // Full well formed sequence 
			temp.pos = seqPos - WL; // Take position on read
			// Store the new word
			storeWord(buffer,temp,++wordsInBuffer,&maxWordsStored);
			if(wordsInBuffer == BUFFER_LENGTH){ // Buffer is full
				if(writeBuffer(buffer,bIndx,wrds,wordsInBuffer) < 0){
					return -1;
				}else{
					// Update info
					//readW += wordsInBuffer;
					wordsInBuffer = 0;
				}
			}
		}
		c = fgetc(metag);
	}

	// Free&Close unnecessary varaibles
	free(temp.w.b);
	fclose(metag);

	// Write buffered words
	if(wordsInBuffer != 0)
		if(writeBuffer(buffer,bIndx,wrds,wordsInBuffer) < 0) return -1;

	// Free buffer space
	freeBuffer(buffer,maxWordsStored);












	fclose(bIndx);
	fclose(wrds);

	// Everything finished. All it's ok.
	return 0;
}