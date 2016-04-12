/*
 * @author Fernando Moreno Jabato <jabato@uma.es>
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

	/////////////////////////////////////////////////////////////////
    fprintf(stdout, "\tDict: Starting dictionaries program.\n");
    /////////////////////////////////////////////////////////////////

	// Variables
	uint16_t WL = (uint16_t) atoi(av[3]); // Word length
	BYTES_IN_WORD = WL / 4;
	wentry *buffer; //Buffer of read words
	FILE *metag; // Input files
	FILE *wDic,*pDic;  // Output files
	FILE *bIndx, *wrds; // Intermediate files
	uint64_t readW = 0, wordsInBuffer = 0,i; // Absolute and buffer read words and auxiliar variable(i)
	uint32_t numBuffWritten = 0;
	char *fname;
	unsigned char *WordsBlock;
	bool removeIntermediataFiles = true; // Config it if you want save or not the itnermediate files

    /////////////////////////////////////////////////////////////////
    fprintf(stdout, "\tDict: Opening/creating necessary files.");
    fflush(stdout);
    /////////////////////////////////////////////////////////////////

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

	// Memory for sequences
	if((WordsBlock = (unsigned char *)malloc(sizeof(unsigned char)*BYTES_IN_WORD*BUFFER_LENGTH))==NULL){
		fprintf(stderr, "Error allocating memory for sequences.\n");
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

    /////////////////////////////////////////////////////////////////
    fprintf(stdout, " (Done)\n");
    fprintf(stdout, "\tDict: Reading k-mers.");
    fflush(stdout);
    /////////////////////////////////////////////////////////////////

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
	temp.w.WL = WL;

	// Memory for buffer words
	uint64_t blockIndex = 0;
	for(i=0;i<BUFFER_LENGTH;++i,blockIndex+=BYTES_IN_WORD)
		buffer[i].w.b = &WordsBlock[blockIndex];


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
		c = fgetc(metag);
	}

	// Free&Close unnecessary varaibles
	fclose(metag);
	// Write buffered words
	if(wordsInBuffer != 0){
		if(numBuffWritten == 0){ // Special case, only one buffer
            /////////////////////////////////////////////////////////////////
            fprintf(stdout, " (Done)\n");
            fprintf(stdout, "\tDict: Writting dictionary.");
            fflush(stdout);
            /////////////////////////////////////////////////////////////////

			// Sort buffer 
			quicksort_W(buffer,0,wordsInBuffer-1);

			// Close unnecesary files
			fclose(bIndx);
			fclose(wrds);

			// Open necessary files
			// Open positions file
			strcpy(fname,av[2]); // Copy outDic name
			if((pDic = fopen(strcat(fname,".d2hP"),"wb"))==NULL){
				fprintf(stderr, "Error opening positions dictionary file.\n");
				return -1;
			}

			// Open dictionary words file 
			strcpy(fname,av[2]); // Copy outDic name
			if((wDic = fopen(strcat(fname,".d2hW"),"wb"))==NULL){
				fprintf(stderr, "Error opening words dictionary file.\n");
				return -1;
			}
			fwrite(&WL,sizeof(uint16_t),1,wDic); // Word length

			// Declare necessary varibles
			uint32_t reps = 0;

			// First entrance
			int k;
			for(k=0;k<BYTES_IN_WORD;++k)
				fwrite(&buffer[0].w.b[k],sizeof(unsigned char),1,wDic); // write first word
			uint64_t pos = (uint64_t)ftell(pDic);
			fwrite(&pos,sizeof(uint64_t),1,wDic); // position on pDic
			fwrite(&buffer[0].seq,sizeof(uint32_t),1,pDic); // Read index
			fwrite(&buffer[0].pos,sizeof(uint64_t),1,pDic); // Position on read
			reps++; // Increment number of repetitions
			storeWord(&temp,buffer[0]); // Update last word written

			// Write final dictionary file
			uint64_t index;
			for(index=1; index<wordsInBuffer; ++index){	
				// Store word in buffer
				writeWord(&buffer[index],wDic,pDic,wordcmp(buffer[index].w,temp.w,BYTES_IN_WORD)!=0? false:true,&reps);
				storeWord(&temp,buffer[index]); // Update last word written
			}
			// Write last word index
			fwrite(&reps,sizeof(uint32_t),1,wDic); // Write num of repetitions

            /////////////////////////////////////////////////////////////////
            fprintf(stdout, " (Done)\n");
            fprintf(stdout, "\tDict: Closing the program.\n");
            fflush(stdout);
            /////////////////////////////////////////////////////////////////

			// Close files
			fclose(pDic);
			fclose(wDic);

			// Free memory
			free(WordsBlock);
			free(buffer);
			free(fname);
			free(temp.w.b);

			// Remove itnermediate files if it's necessary
			if(removeIntermediataFiles){
				strcpy(fname,av[2]);
				remove(strcat(fname,".bindx"));
				strcpy(fname,av[2]);
				remove(strcat(fname,".wrds"));
			}
			// End program
			return 0;
		}else if(writeBuffer(buffer,bIndx,wrds,wordsInBuffer) < 0){ return -1;
		}else numBuffWritten++;
	}else if(wordsInBuffer == 0 && numBuffWritten == 0){ // Special case
        /////////////////////////////////////////////////////////////////
        fprintf(stdout, "\tDict: Any k-mer found.\n");
        fprintf(stdout, "\tDict: Closing the program.\n");
        /////////////////////////////////////////////////////////////////

		free(WordsBlock);
		free(buffer);
		free(fname);
		free(temp.w.b);

		fclose(bIndx);
		fclose(wrds);

		// Remove itnermediate files if it's necessary
		if(removeIntermediataFiles){
			strcpy(fname,av[2]);
			remove(strcat(fname,".bindx"));
			strcpy(fname,av[2]);
			remove(strcat(fname,".wrds"));
		}

		return 0;
	}

	// Free buffer space
	free(WordsBlock);
	free(buffer);

    /////////////////////////////////////////////////////////////////
    fprintf(stdout, " (Done)\n");
    fprintf(stdout, "\tDict: Writting dictionary.");
    fflush(stdout);
    /////////////////////////////////////////////////////////////////

	// Read Intermediate files and create final dictionary
		// Close write buffer and open read buffers
		fclose(bIndx);
		fclose(wrds);

		strcpy(fname,av[2]); // Copy outDic name
		if((bIndx = fopen(strcat(fname,".bindx"),"rb"))==NULL){
			fprintf(stderr, "Error opening buffer index file (read).\n");
			return -1;
		}
		
		if(fread(&WL,sizeof(uint16_t),1,bIndx) != 1){ // Read header = word length
			fprintf(stderr, "Error reading word length at header. Value not found.\n");
			return -1;
		} 

		// Open words repo
		strcpy(fname,av[2]);
		if((wrds = fopen(strcat(fname,".wrds"),"rb"))==NULL){
			fprintf(stderr, "Error opening words repository (read).\n");
			return -1;
		}

		// Open positions file
		strcpy(fname,av[2]); // Copy outDic name
		if((pDic = fopen(strcat(fname,".d2hP"),"wb"))==NULL){
			fprintf(stderr, "Error opening positions dictionary file.\n");
			return -1;
		}

		// Open dictionary words file 
		strcpy(fname,av[2]); // Copy outDic name
		if((wDic = fopen(strcat(fname,".d2hW"),"wb"))==NULL){
			fprintf(stderr, "Error opening words dictionary file.\n");
			return -1;
		}
		fwrite(&WL,sizeof(uint16_t),1,wDic); // Word length

	// Prepare necessary variables
	uint64_t arrPos[numBuffWritten];
	int64_t wordsUnread[numBuffWritten];
	node_W *words = NULL;
	uint32_t reps = 0;
	uint64_t lastLoaded = -1, activeBuffers = numBuffWritten;
	unsigned char *BuffWordsBlock;
	wentry *WentryBlock;

	// Memory for sequences
	if((BuffWordsBlock = (unsigned char *)malloc(sizeof(unsigned char)*activeBuffers*BYTES_IN_WORD*READ_BUFF_LENGTH))==NULL){
		fprintf(stderr, "Error allocating memory for sequences (merge).\n");
		return -1;
	}

	if((WentryBlock = (wentry*) malloc(sizeof(wentry)*activeBuffers*READ_BUFF_LENGTH))==NULL){
		fprintf(stderr, "Error allocating memory for wentry block.\n");
		return -1;
	}

	// Read info about buffers
	i = 0;
	uint64_t aux64;
	do{
		if(fread(&arrPos[i],sizeof(uint64_t),1,bIndx)!=1){ // Position on words file
			fprintf(stderr, "Error reading position on P file.\n");
			return -1;
		}
		if(fread(&aux64,sizeof(uint64_t),1,bIndx)!=1){// Number of words on set
			fprintf(stderr, "Error reading number of repetitions at W file.\n");
			return -1;
		}
		wordsUnread[i] = (int64_t) aux64;
		++i;
	}while(i<activeBuffers);
	fclose(bIndx);

	// Take memory for words & read first set of words
	blockIndex = 0;
	uint64_t wblockIndex = 0, j,read;
	node_W *currNode = NULL; 
	for(i=0 ;i<activeBuffers ;++i, wblockIndex += READ_BUFF_LENGTH){
		currNode = (node_W*) malloc(sizeof(node_W)); // New node
		currNode->word = &WentryBlock[wblockIndex]; // Wentry array
		for(j=0;j<READ_BUFF_LENGTH;++j,blockIndex += BYTES_IN_WORD) // Words memory
			currNode->word[j].w.b = &BuffWordsBlock[blockIndex];
		currNode->next = words;
		fseek(wrds,arrPos[i],SEEK_SET);
		read = loadWord(&currNode->word,wrds,wordsUnread[i]);
		currNode->buff = i;
		currNode->index = 0;
		currNode->words_loaded = read;
		words = currNode;
		// Update info
		arrPos[i] = (uint64_t) ftell(wrds);
		wordsUnread[i]-=read;
	}
	words = currNode;

	// Sort array
	sortList(&words);

	// First entrance
	int k;
	for(k=0;k<BYTES_IN_WORD;++k)
		fwrite(&words->word[words->index].w.b[k],sizeof(unsigned char),1,wDic); // write first word
	uint64_t pos = (uint64_t)ftell(pDic);
	fwrite(&pos,sizeof(uint64_t),1,wDic); // position on pDic
	fwrite(&words->word[words->index].seq,sizeof(uint32_t),1,pDic); // Read index
	fwrite(&words->word[words->index].pos,sizeof(uint64_t),1,pDic); // Position on read
	reps++; // Increment number of repetitions
	storeWord(&temp,words->word[words->index]); // Update last word written
	words->index += 1;
	if(words->index >= words->words_loaded){
		if(wordsUnread[words->buff] > 0){
			fseek(wrds,arrPos[words->buff],SEEK_SET);
			read = loadWord(&words->word,wrds,wordsUnread[words->buff]);
			words->index = 0;
			words->words_loaded = read;
			lastLoaded = words->buff;
			arrPos[i] = (uint64_t) ftell(wrds);
			wordsUnread[i]-=read;
			checkOrder(&words,false);
		}else{
			checkOrder(&words,true);
			activeBuffers--;
			if(activeBuffers <= 0){
				free(BuffWordsBlock);
				free(WentryBlock);
				free(fname);
				free(temp.w.b);
				fclose(wDic);
				fclose(pDic);
				fclose(wrds);
				// Delete intermediate files
				if(removeIntermediataFiles){
					strcpy(fname,av[2]);
					remove(strcat(fname,".bindx"));
					strcpy(fname,av[2]);
					remove(strcat(fname,".wrds"));
				}

				return 0;
			}
		}
	}

	// Write final dictionary file
	while(activeBuffers > 0){	
		// Store word in buffer
		writeWord(&words->word[words->index],wDic,pDic,wordcmp(words->word[words->index].w,temp.w,BYTES_IN_WORD)!=0? false:true,&reps);
		storeWord(&temp,words->word[words->index]); // Update last word written
		words->index += 1;
		if(words->index >= words->words_loaded){
			// Load next word if it's possible
			if(wordsUnread[words->buff] > 0){
				if(words->buff != lastLoaded){
					fseek(wrds,arrPos[words->buff],SEEK_SET);
					lastLoaded = words->buff;
				}
				read = loadWord(&words->word,wrds,wordsUnread[words->buff]);
				words->index = 0;
				words->words_loaded = read;
				arrPos[words->buff] = (uint64_t) ftell(wrds);
				wordsUnread[words->buff]-= read;
				checkOrder(&words,false);
			}else{
				checkOrder(&words,true);
				activeBuffers--;
			}
		}else checkOrder(&words,false);
	}
	// Write last word index
	fwrite(&reps,sizeof(uint32_t),1,wDic); // Write num of repetitions

    /////////////////////////////////////////////////////////////////
    fprintf(stdout, " (Done)\n");
    fprintf(stdout, "\tDict: Closing the program.\n");
    fflush(stdout);
    /////////////////////////////////////////////////////////////////

	// Deallocate words buffer
	free(BuffWordsBlock);
	free(WentryBlock);

	fclose(wDic);
	fclose(pDic);
	fclose(wrds);

	if(removeIntermediataFiles){
		strcpy(fname,av[2]);
		remove(strcat(fname,".bindx"));
		strcpy(fname,av[2]);
		remove(strcat(fname,".wrds"));
	}

	// Free space
	free(fname);
	free(temp.w.b);
	// Everything finished. All it's ok.
	return 0;
}