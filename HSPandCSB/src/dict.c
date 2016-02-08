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
	if(ac<4 || ac>5){
		fprintf(stderr, "Bad call error.\nUSE: dict seq.IN dictName WL\nOR: dict seq.IN dictName WL metag/geno\n");
		return -1;
	}else if(atoi(av[3])%4 != 0){
		fprintf(stderr, "Error WL must be a 4 multiple and it's \"%s\".\n", av[3]);
		return -1;
	}

	// Take extensions
	bool isMetag = true;
	if(ac==5 && strcmp(av[4],(const char*)&("geno"))==0) isMetag = false;

	// Variables
	WL = (uint16_t) atoi(av[3]); // Word length
	BYTES_IN_WORD = WL / 4;
	wentry *buffer; //Buffer of read words
	FILE *seqFile; // Input files
	FILE *wDic,*pDic;  // Output files
	FILE *bIndx, *wrds; // Intermediate files
	uint64_t readW = 0, wordsInBuffer = 0,i; // Absolute and buffer read words and auxiliar variable(i)
	uint32_t numBuffWritten = 0;
	char *fname;
	unsigned char *MemoryBlock;
	bool removeIntermediataFiles = true; // Config it if you want save or not the itnermediate files

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
	// Open sequence file
	if((seqFile = fopen(av[1],"rt"))==NULL){
		fprintf(stderr, "Error opening sequences file.\n");
		return -1;
	}

	//Open intermediate files
	strcpy(fname,av[2]); // Copy outDic name
	if((bIndx = fopen(strcat(fname,".bindx"),"wb"))==NULL){
		fprintf(stderr, "Error opening buffer index file (write).\n");
		return -1;
	}
	// Write first line of buffer index file (WordLength)
	fwrite(&WL,sizeof(uint16_t),1,bIndx);


	// Open words repo
	strcpy(fname,av[2]);
	if((wrds = fopen(strcat(fname,".wrds"),"wb"))==NULL){
		fprintf(stderr, "Error opening words repository (write).\n");
		return -1;
	}

	// START WORKFLOW
	// Read sequence file
	// Necessary variables
	uint64_t seqPos = 0, crrSeqL = 0;
	char c; // Auxiliar -> Read char per char
	wentry temp; // Auxiliar word variable
	if((temp.seq = (unsigned char*) malloc(sizeof(unsigned char)*BYTES_IN_WORD))==NULL){
		fprintf(stderr, "Error allocating memory for temp.\n");
		return -1;
	}
	temp.seqIndex = 0;

	// Memory for buffer words
	if((MemoryBlock = (unsigned char*)malloc(sizeof(unsigned char)*BYTES_IN_WORD*BUFFER_LENGTH))==NULL){
		fprintf(stderr, "Error allocating space for words.\n");
		return -1;
	}
	uint64_t mbIndex=0; // Memory block index
	for(i=0;i<BUFFER_LENGTH;++i, mbIndex+=BYTES_IN_WORD)
		buffer[i].seq = &MemoryBlock[mbIndex];

	// Start to read
	c = fgetc(seqFile);
	while(!feof(seqFile)){
		// Check if it's a special line
		if(!isupper(toupper(c))){ // Comment, empty or quality (+) line
			if(c=='>'){ // Comment line
				c = fgetc(seqFile);
				while(c != '\n') // Avoid comment line
					c = fgetc(seqFile);
				temp.seqIndex++; // New sequence
				crrSeqL = 0; // Reset buffered sequence length
				seqPos = 0; // Reset index
			}
			c=fgetc(seqFile); // First char of next sequence
			continue;
		}
		shift_word(temp.seq); // Shift bits sequence
		// Add new nucleotid
		switch (c) {
			case 'A': // A = 00 
				crrSeqL++;
				break;
			case 'C': // C = 01
				temp.seq[BYTES_IN_WORD-1]|=1;
				crrSeqL++;
				break;
			case 'G': // G = 10
				temp.seq[BYTES_IN_WORD-1]|=2;
				crrSeqL++;
				break;
			case 'T': // T = 11
				temp.seq[BYTES_IN_WORD-1]|=3;
				crrSeqL++;
				break;
			default : // Bad formed sequence
				crrSeqL = 0; break;
		}
		seqPos++;
		if(crrSeqL >= (uint64_t)WL){ // Full well formed sequence 
			temp.pos = seqPos - WL; // Take position on read
			// Store the new word
			storeWord(&buffer[wordsInBuffer],temp);
			readW++; wordsInBuffer++;
			if(wordsInBuffer == BUFFER_LENGTH){ // Buffer is full
				if(writeBuffer(buffer,bIndx,wrds,wordsInBuffer) < 0){
					return -1;
				}else{
					// Update info
					//readW += wordsInBuffer;
					wordsInBuffer = 0;
					numBuffWritten++;
				}
			}
		}
		c = fgetc(seqFile);
	}

	// Free&Close unnecessary varaibles
	fclose(seqFile);

	// Write buffered words
	if(wordsInBuffer != 0){
		if(writeBuffer(buffer,bIndx,wrds,wordsInBuffer) < 0) return -1;
		else numBuffWritten++;
	}

	// Free buffer space
	free(MemoryBlock);
	free(buffer);

	// Read Intermediate files and create final dictionary
		// Close write buffer and open read buffers
		fclose(bIndx);
		fclose(wrds);

		strcpy(fname,av[2]); // Copy outDic name
		if((bIndx = fopen(strcat(fname,".bindx"),"rb"))==NULL){
			fprintf(stderr, "Error opening buffer index file (read).\n");
			return -1;
		}
		fread(&WL,sizeof(uint16_t),1,bIndx); // Read header = word length

		// Open words repo
		strcpy(fname,av[2]);
		if((wrds = fopen(strcat(fname,".wrds"),"rb"))==NULL){
			fprintf(stderr, "Error opening words repository (read).\n");
			return -1;
		}

		// Open positions file
		strcpy(fname,av[2]); // Copy outDic name
		if((pDic = fopen(strcat(fname,isMetag?".md2hP":".gd2hP"),"wb"))==NULL){
			fprintf(stderr, "Error opening positions dictionary file.\n");
			return -1;
		}

		// Open dictionary words file 
		strcpy(fname,av[2]); // Copy outDic name
		if((wDic = fopen(strcat(fname,isMetag?".md2hW":".gd2hW"),"wb"))==NULL){
			fprintf(stderr, "Error opening words dictionary file.\n");
			return -1;
		}
		fwrite(&WL,sizeof(uint16_t),1,wDic); // Word length

////////////////////////////////////////////////////////////////////
fprintf(stderr, "TEST1\n");
////////////////////////////////////////////////////////////////////

	// Prepare necessary variables
	uint64_t arrPos[numBuffWritten]; // Where stat to read buffer i
	int64_t wordsUnread[numBuffWritten]; // Words unread of buffer i
	uint32_t index[numBuffWritten];
	uint64_t wrdsInBuff[numBuffWritten];
	wentry **words; // Buffers
	uint32_t reps = 0, activeBuffers = numBuffWritten; // Auxiliar variable
	wentry *WentryBlock;
////////////////////////////////////////////////////////////////////
fprintf(stderr, "TEST2\n");
////////////////////////////////////////////////////////////////////

	// Allocate space for buffers
	if((words = (wentry **) malloc(sizeof(wentry*)*numBuffWritten))==NULL){
		fprintf(stderr, "Error allocating matrix of words.\n");
		return -1;
	}
	if((WentryBlock = (wentry*) malloc(sizeof(wentry)*MERGE_BUFFER_LENGTH*numBuffWritten))==NULL){
		fprintf(stderr, "Error allocating wentry block.\n");
		return -1;
	}
////////////////////////////////////////////////////////////////////
fprintf(stderr, "TEST3\n");
////////////////////////////////////////////////////////////////////

	// Read info about buffers
	i = 0;
	uint64_t aux64;
////////////////////////////////////////////////////////////////////
fprintf(stderr, "TEST4\n");
////////////////////////////////////////////////////////////////////
	do{
		fread(&arrPos[i],sizeof(uint64_t),1,bIndx); // Position on words file
		fread(&aux64,sizeof(uint64_t),1,bIndx); // Number of words on set
		wordsUnread[i] = (int64_t) aux64;
		++i;
	}while(i<numBuffWritten);
////////////////////////////////////////////////////////////////////
fprintf(stderr, "TEST5\n");
////////////////////////////////////////////////////////////////////

	// Take memory for words & read first set of words
	if((MemoryBlock = (unsigned char*)malloc(sizeof(unsigned char)*BYTES_IN_WORD*activeBuffers))==NULL){
		fprintf(stderr, "Error allocating space for bufffer words.\n");
		return -1;
	}
////////////////////////////////////////////////////////////////////
fprintf(stderr, "TEST6\n");
////////////////////////////////////////////////////////////////////

	// Load matrix
	mbIndex = 0;
	uint64_t wbIndex = 0,j;
	for(i=0 ;i<activeBuffers ;++i, wbIndex+=MERGE_BUFFER_LENGTH){
		words[i] = &WentryBlock[wbIndex];
		wrdsInBuff[i]=0;

		fseek(wrds,arrPos[i],SEEK_SET);

		for(j=0;j<MERGE_BUFFER_LENGTH && wordsUnread[i]>0;++j, mbIndex+=BYTES_IN_WORD){
			words[i][j].seq = &MemoryBlock[mbIndex];
			loadWord(&words[i][j],wrds);
			wordsUnread[i]--;
			wrdsInBuff[i]++;
		}
		// Update info
		arrPos[i] = (uint64_t) ftell(wrds);
	}
////////////////////////////////////////////////////////////////////
fprintf(stderr, "TEST7\n");
////////////////////////////////////////////////////////////////////

	// Sort buffers
	sortBuffers(&words,activeBuffers,&index[0],&wrdsInBuff[0]);
////////////////////////////////////////////////////////////////////
fprintf(stderr, "TEST8\n");
////////////////////////////////////////////////////////////////////

	// First entrance
	fwrite(words[0][0].seq,sizeof(unsigned char),BYTES_IN_WORD,wDic);
	uint64_t pos = (uint64_t)ftell(pDic);
	fwrite(&pos,sizeof(uint64_t),1,wDic); // position on pDic
	fwrite(&words[0][0].seqIndex,sizeof(uint32_t),1,pDic); // Read index
	fwrite(&words[0][0].pos,sizeof(uint64_t),1,pDic); // Position on read
	reps++; // Increment number of repetitions
	index[0]++;
////////////////////////////////////////////////////////////////////
fprintf(stderr, "TEST9\n");
////////////////////////////////////////////////////////////////////
	storeWord(&temp,words[0][0]); // Update last word written
/*	if(wordsUnread[i] > 0){
		fseek(wrds,arrPos[i],SEEK_SET);
		loadWord(&words[i],wrds);
		lastLoaded = i;
		arrPos[i] = (uint64_t) ftell(wrds);
		wordsUnread[i]--;
	}
*/
////////////////////////////////////////////////////////////////////
fprintf(stderr, "TEST10\n");
////////////////////////////////////////////////////////////////////
	checkOrder(&words,&activeBuffers,&index[0],&wrdsInBuff[0]);
////////////////////////////////////////////////////////////////////
fprintf(stderr, "TEST11\n");
////////////////////////////////////////////////////////////////////

	// Write final dictionary file
	while(!finished(&wordsUnread[0],numBuffWritten)){
		// Reload matrix if it's necessary
		if(activeBuffers <= 0)
			loadMatrix(words,&wordsUnread[0],numBuffWritten,&index[0],&wrdsInBuff[0],&activeBuffers,&arrPos[0],wrds);	

		while(activeBuffers > 0){
			// Store word in buffer
			writeWord(&words[0][index[0]],wDic,pDic,wordcmp(words[0][index[0]].seq,temp.seq,BYTES_IN_WORD)!=0? false:true,&reps);
			storeWord(&temp,words[0][index[0]]); // Update last word written
			index[0]++;
			// Check order and shift if it's necessary
			checkOrder(&words,&activeBuffers,&index[0],&wrdsInBuff[0]);
		}
	}
////////////////////////////////////////////////////////////////////
fprintf(stderr, "TEST12\n");
////////////////////////////////////////////////////////////////////

	// Write last word index
	fwrite(&reps,sizeof(uint32_t),1,wDic); // Write num of repetitions

	// Deallocate words buffer
	free(MemoryBlock);

	if(removeIntermediataFiles){
		strcpy(fname,av[2]);
		remove(strcat(fname,".bindx"));
		strcpy(fname,av[2]);
		remove(strcat(fname,".wrds"));
	}

	// Free space
	free(fname);
	free(temp.seq);
	free(WentryBlock);
	free(words);

	// Everything finished. All it's ok.
	return 0;
}