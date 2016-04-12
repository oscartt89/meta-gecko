/*
 * @author Fernando Moreno Jabato <jabato@uma.es>
 * @description This file encodes the workflow of GECKO for create
 *    metagenome dictionaries.
 * @licence all rights reserved to the author and BitLAB group (University
 *    of Malaga).
 */
#include "frags.h"

int main(int ac, char** av){
	// Check arguments
	if(ac!=10){
		fprintf(stderr, "Bad call error.\nUSE: frags metagDic metagFile genoDic genoFile out minS minL f/r prefix\n");
		return -1;
	}

	/////////////////////////////////////////////////////////////////
	fprintf(stdout, "\tFrags: Starting fragments program.\n");
	/////////////////////////////////////////////////////////////////
	
	// Variables
	FILE *mW,*mP,*gW,*gP; // Dictionaries
	FILE *hIndx,*hts; // Intermediate files
	FILE *fr; // Fragments file
	Hit *buffer;
	uint64_t hitsInBuffer = 0, genomeLength, nStructs, metagenomeLength;
	uint16_t mWL;
	uint16_t BytesGenoWord = 8, BytesMetagWord;
	buffersWritten = 0; // Init global variable (frags.h)1
	S_Threshold = (uint64_t) atoi(av[6]); // Similarity threshold
	L_Threshold = (uint64_t) atoi(av[7]); // Length threshold
	prefixSize = atoi(av[9]); // Prefix array legnth
	Sequence *genome; // Sequence for genome
	Reads *metagenome; // Short sequence array for metagenome
	char *fname; // File names handler
	bool removeIntermediataFiles = true; // Internal variable to delete intermediate files	

	// Allocate necessary memory
	// Memory for buffer
	if((buffer = (Hit*) malloc(sizeof(Hit)*MAX_BUFF))==NULL){
		fprintf(stderr, "Error allocating memory for hits buffer.\n");
		return -1;
	}

	// Memory for file names handler
	if((fname = (char*) malloc(sizeof(char)*MAX_FILE_LENGTH))==NULL){
		fprintf(stderr, "Error allocating memory for file names handler.\n");
		return -1;
	}

	/////////////////////////////////////////////////////////////////
	fprintf(stdout, "\tFrags: Opening/creating necessary files.");
	fflush(stdout);
	/////////////////////////////////////////////////////////////////

	// Open current necessary files
	// Open metagenome positions file
	strcpy(fname,av[1]);
	if((mP = fopen(strcat(fname,".d2hP"),"rb"))==NULL){
		fprintf(stderr, "Error opening metagenome positions dictionaries.\n");
		return -1;
	}

	// Open metagenome words file
	strcpy(fname,av[1]);
	if((mW = fopen(strcat(fname,".d2hW"),"rb"))==NULL){
		fprintf(stderr, "Error opening metagenome words dictionaries.\n");
		return -1;
	}
	// Read words header = WordLength
		// Check
		if(fread(&mWL,sizeof(uint16_t),1,mW)!=1){ 
			fprintf(stderr, "Error, couldn't find word length.\n");
			return -1;
		}else if(mWL % 4 != 0){
			fprintf(stderr, "Error, word length of metagenome dictionary isn't a 4 multiple.\n");
			return -1;
		}else{
			BytesMetagWord = mWL/4;
			if(BytesMetagWord < prefixSize || BytesGenoWord < prefixSize){
				fprintf(stderr, "Error: prefix is too long.\n");
				return -1;
			}
		}

	// Open genome postions file
	strcpy(fname,av[3]);
	if((gP = fopen(strcat(fname,".d2hP"),"rb"))==NULL){
		fprintf(stderr, "Error opening genome positions dictionaries.\n");
		return -1;
	}

	// Open genome words file
	strcpy(fname,av[3]);
	if((gW = fopen(strcat(fname,".d2hW"),"rb"))==NULL){
		fprintf(stderr, "Error opening genome words dictionaries.\n");
		return -1;
	}

	//Open intermediate files
	strcpy(fname,av[5]); // Copy outDic name
	if((hIndx = fopen(strcat(fname,".hindx"),"wb"))==NULL){
		fprintf(stderr, "Error opening buffer index file.\n");
		return -1;
	}

	// Open hits repo
	strcpy(fname,av[5]);
	if((hts = fopen(strcat(fname,".hts"),"wb"))==NULL){
		fprintf(stderr, "Error opening hits repository.\n");
		return -1;
	}

	// Search hits
		// Prepare necessary variables
		WordEntry we[2]; // [0]-> Metagenome [1]-> Genome
		uint64_t lastFirstHit = (uint64_t)ftell(gW);
		bool firstmatch = true;
		int cmp;

		// Take memory
		if((we[0].seq = (unsigned char *) malloc(sizeof(unsigned char)*BytesMetagWord))==NULL){
			fprintf(stderr, "Error allocating memory for metagenome entrance.\n");
			return -1;
		}else we[0].WB = BytesMetagWord;

		if((we[1].seq = (unsigned char *) malloc(sizeof(unsigned char)*BytesGenoWord))==NULL){
			fprintf(stderr, "Error allocating memory for metagenome entrance.\n");
			return -1;
		}else we[1].WB = BytesGenoWord;

	// Read first entrances
	if(readWordEntrance(&we[0],mW,BytesMetagWord)<0) return -1;

	readHashEntry(&we[1],gW);

	/////////////////////////////////////////////////////////////////
	fprintf(stdout, " (Done)\n");
	fprintf(stdout, "\tFrags: Generating seeds.");
	fflush(stdout);
	/////////////////////////////////////////////////////////////////

	// Search
	if(prefixSize == BytesMetagWord && prefixSize == BytesGenoWord){
		while(!feof(mW) && !feof(gW)){
			if((cmp = wordcmp(we[0].seq,we[1].seq,BytesGenoWord))==0) // Hit
				generateHits(buffer,we[0],we[1],mP,gP,hIndx,hts,&hitsInBuffer,BytesGenoWord);
			// Load next word
			if(cmp >= 0) // New genome word is necessary
				readHashEntry(&we[1],gW);
			if(cmp <= 0) // New metagenome word is necessary
				if(readWordEntrance(&we[0],mW,BytesMetagWord)<0) return -1;
		}
	}else{
		while(!feof(mW)){
			// Check hit
			if((cmp = wordcmp(we[0].seq,we[1].seq,prefixSize))==0){ // Hit
				generateHits(buffer,we[0],we[1],mP,gP,hIndx,hts,&hitsInBuffer,prefixSize);
				if(firstmatch){
					lastFirstHit = (uint64_t)(ftell(gW) - sizeof(hashentry));
					firstmatch = false;
				}
			}
			// Check if could be more
			if(cmp >= 0){ // Could be more
				// Load next genome word
				readHashEntry(&we[1],gW);
				if(feof(gW)){ // End of genome file
					// Load next metagenome word
					if(readWordEntrance(&we[0],mW,BytesMetagWord)<0) return -1;
					// Reset values and come back at dict
					firstmatch = true;
					fseek(gW,lastFirstHit,SEEK_SET); // Reset geno dict
					readHashEntry(&we[1],gW);
				}
			}else if(cmp < 0){ // No more matches, take next metag word
				// Load next metagenome word
				if(readWordEntrance(&we[0],mW,BytesMetagWord)<0) return -1;
				// Reset values and come back at dict
				firstmatch = true;
				fseek(gW,lastFirstHit,SEEK_SET); // Reset geno dict
				readHashEntry(&we[1],gW);
			}
		}
	}

	/////////////////////////////////////////////////////////////////
	fprintf(stdout, " (Generated)\n");
	fprintf(stdout, "\tFrags: Loading sequences of genome.");
	fflush(stdout);
	/////////////////////////////////////////////////////////////////

	// Load sequences
	genome = LeeSeqDB(av[4], &genomeLength, &nStructs);

	/////////////////////////////////////////////////////////////////
	fprintf(stdout, " (Loaded)\n");
	fprintf(stdout, "\tFrags: Loading sequences of metagenome.");
	fflush(stdout);
	/////////////////////////////////////////////////////////////////

	metagenome = LoadMetagenome(av[2],&metagenomeLength);

	/////////////////////////////////////////////////////////////////
	fprintf(stdout, " (Loaded)\n");
	/////////////////////////////////////////////////////////////////

	// Write buffered hits
	if(hitsInBuffer > 0){
		if(buffersWritten > 0)
			writeHitsBuff(buffer,hIndx,hts,hitsInBuffer);
		else{ // Only one buffer
			// Sort buffer
			quicksort_H(buffer,0,hitsInBuffer-1);
			
			// Close unnecesary files
			fclose(mW); fclose(gW);
			fclose(mP); fclose(gP);
			fclose(hIndx);
			fclose(hts);

			// Free unnecesary variables 
			free(we[0].seq);
			free(we[1].seq);

			// Declare necessary variables
			FragFile frag;
			uint64_t index;
			int64_t distX,distY;
			
			// Init first fragment
			frag.block = 0;
			frag.strand = av[8][0];

			// Write first fragment
			Reads *currRead;
				// Search first read
				currRead = metagenome;
				while(currRead->seqIndex != buffer[0].seqX){
					if(currRead->next == NULL){
						fprintf(stderr, "Error searching first read.\n");
						return -1;
					}
					currRead = currRead->next;
				}
			// Open final files
			// Open final fragments file
			strcpy(fname,av[5]); // Copy outDic name
			if((fr = fopen(strcat(fname,".frags"),"wb"))==NULL){
				fprintf(stderr, "Error opening fragments final file.\n");
				return -1;
			}

			// Write headers
			writeSequenceLength(&metagenomeLength, fr);
			writeSequenceLength(&genomeLength, fr);

			/////////////////////////////////////////////////////////////////
			fprintf(stdout, "\tFrags: Extending seeds .");
			fflush(stdout);
			/////////////////////////////////////////////////////////////////

			// Generate first fragment
			FragFromHit(&frag,&buffer[0],currRead,genome,genomeLength,nStructs,fr);
			
			// Generate fragments
			for(index=1; index < hitsInBuffer; ++index){
				if(buffer[index].diag == frag.diag &&
						buffer[index].seqX == frag.seqX && 
						buffer[index].seqY == frag.seqY){ // Possible fragment
					// Check if are collapsable
					distX = (int64_t)(buffer[index].posX - frag.xEnd);
					distY = (int64_t)(buffer[index].posY - frag.yEnd);
					if(distX > 0 || distY > 0){ // Not collapsable by extension
						// Generate fragment 
						FragFromHit(&frag, &buffer[index],currRead,genome,genomeLength,nStructs,fr);
					}
				}else{ // New fragment
					// Check correct read index
					//currRead = metagenome;
					while(currRead->seqIndex != buffer[index].seqX){
						if(currRead->next == NULL){
							fprintf(stderr, "Error searching read index.\n");
							return -1;
						}
						currRead = currRead->next;
					}
					// Generate new fragment
					FragFromHit(&frag, &buffer[index],currRead,genome,genomeLength,nStructs,fr);
				}
			}

			/////////////////////////////////////////////////////////////////
			fprintf(stdout, " (Generated)\n");
			fprintf(stdout, "\tFrags: Closing the program.\n");
			/////////////////////////////////////////////////////////////////

			// Close output file
			fclose(fr);
//			freeReads(&metagenome);
			// Free unnecesary memory
			free(buffer);
			// Remove intermediate files
			if(removeIntermediataFiles){
				strcpy(fname,av[5]);
				remove(strcat(fname,".hts"));
				strcpy(fname,av[5]);
				remove(strcat(fname,".hindx"));
			}

			free(fname);

			// End program
			return 0;
		}
	}else if(hitsInBuffer == 0 && buffersWritten == 0){
		/////////////////////////////////////////////////////////////////
		fprintf(stdout, "\tFrags: Any match found.\n");
		fprintf(stdout, "\tFrags: Closing the program.\n");
		/////////////////////////////////////////////////////////////////

		// Free auxiliar buffers
		free(we[0].seq);
		free(we[1].seq);
		free(buffer);

		// Close files
		fclose(mW); fclose(gW);
		fclose(mP); fclose(gP);
		fclose(hIndx);
		fclose(hts);
		freeReads(&metagenome);

		// Remove intermediate files
		if(removeIntermediataFiles){
			strcpy(fname,av[5]);
			remove(strcat(fname,".hts"));
			strcpy(fname,av[5]);
			remove(strcat(fname,".hindx"));
		}

		free(fname);

		// End program
		return 0;
	}

	// Free auxiliar buffers
	free(we[0].seq);
	free(we[1].seq);
	free(buffer);

	// Close files
	fclose(mW); fclose(gW);
	fclose(mP); fclose(gP);
	fclose(hIndx);
	fclose(hts);

	// Open necessary files
	// Open intermediate files
	// Index file
	strcpy(fname,av[5]); // Copy outDic name
	if((hIndx = fopen(strcat(fname,".hindx"),"rb"))==NULL){
		fprintf(stderr, "Error opening buffer index file (read).\n");
		return -1;
	}

	// Open hits repo
	strcpy(fname,av[5]);
	if((hts = fopen(strcat(fname,".hts"),"rb"))==NULL){
		fprintf(stderr, "Error opening hits repository (read).\n");
		return -1;
	}

	// Open final fragments file
	strcpy(fname,av[5]); // Copy outDic name
	if((fr = fopen(strcat(fname,".frags"),"wb"))==NULL){
		fprintf(stderr, "Error opening fragments final file.\n");
		return -1;
	}

	// Write header
	writeSequenceLength(&metagenomeLength, fr);
	writeSequenceLength(&genomeLength, fr);

	// Prepare necessary variables
	node_H *hitsList = NULL;
	uint64_t hitsUnread[buffersWritten];
	uint64_t positions[buffersWritten];
	uint64_t lastLoaded = -1, activeBuffers = buffersWritten;
	Hit *HitsBlock;
	FragFile frag;

	// Take memory for hits
	if((HitsBlock = (Hit*) malloc(sizeof(Hit)*activeBuffers*READ_BUFF_LENGTH))==NULL){
		fprintf(stderr, "Error allocating memory for hits block.\n");
		return -1;
	}

	// Read buffers info
	uint64_t i = 0;
	do{
		if(fread(&positions[i],sizeof(uint64_t),1,hIndx)!=1){ // Take position on hits file
			fprintf(stderr, "Error reading position at hIndx file.\n");
			return -1;
		}
		if(fread(&hitsUnread[i],sizeof(uint64_t),1,hIndx)!=1){
			fprintf(stderr, "Error reading hits in buffer at hIndx file.\n");
			return -1;
		}
		++i;
	}while(i < activeBuffers);

	// Load first hits
	node_H *currNode = NULL;
	uint64_t read, blockIndex = 0;
	for(i=0 ;i<activeBuffers; ++i, blockIndex += READ_BUFF_LENGTH){
		currNode = (node_H*) malloc(sizeof(node_H));
		currNode->next = hitsList;
		currNode->hits = &HitsBlock[blockIndex];
		currNode->buff = i;
		fseek(hts,positions[i],SEEK_SET);
		read = loadHit(&currNode->hits,hts,hitsUnread[i]);
		currNode->index = 0;
		currNode->hits_loaded = read;
		// Update info
		positions[i] = (uint64_t) ftell(hts);
		hitsUnread[i]-=read;
		lastLoaded = i;
		hitsList = currNode;
	}

	// Assign head
	hitsList = currNode;

	// Sort hits
	sortList(&hitsList);	

	/////////////////////////////////////////////////////////////////
	fprintf(stdout, "\tFrags: Extending seeds.");
	fflush(stdout);
	/////////////////////////////////////////////////////////////////

	// Init fragment info
	frag.block = 0;
	frag.strand = av[8][0];

	// Write first fragment
	Reads *currRead;
		// Search first read
		currRead = metagenome;
		while(currRead->seqIndex != hitsList->hits[0].seqX){
			if(currRead->next == NULL){
				fprintf(stderr, "Error searching first read.\n");
				return -1;
			}
			currRead = currRead->next;
		}

	// Generate first fragment
	FragFromHit(&frag,&hitsList->hits[0],currRead,genome,genomeLength,nStructs,fr);

	// Move to next
	hitsList->index +=1;
	// Load new hit
	if(hitsList->index >= hitsList->hits_loaded){
		if(hitsUnread[hitsList->buff] > 0){
			if(hitsList->buff != lastLoaded){
				fseek(hts,positions[hitsList->buff],SEEK_SET);
				lastLoaded = hitsList->buff;
			}
			read = loadHit(&hitsList->hits,hts,hitsUnread[hitsList->buff]);
			hitsList->index = 0;
			hitsList->hits_loaded = read;
			positions[hitsList->buff] = (uint64_t) ftell(hts);
			hitsUnread[hitsList->buff]-=read;
			checkOrder(&hitsList,false);
		}else{
			checkOrder(&hitsList,true);
			activeBuffers--;
		}
	}else{
		checkOrder(&hitsList,false);
	}

	// Search hits and generate fragmetents
	int64_t distX,distY;

	// Read hits & generate fragments
	while(activeBuffers > 0 && hitsList != NULL){
		if(hitsList->hits[hitsList->index].diag == frag.diag &&
			hitsList->hits[hitsList->index].seqX == frag.seqX && 
			hitsList->hits[hitsList->index].seqY == frag.seqY){ // Possible fragment
			// Check if are collapsable
			distX = (int64_t)(hitsList->hits[hitsList->index].posX - frag.xEnd);
			distY = (int64_t)(hitsList->hits[hitsList->index].posY - frag.yEnd);
			if(distX > 0 || distY > 0){ // Not collapsable by xtension
				// Generate fragment 
				FragFromHit(&frag, &hitsList->hits[hitsList->index],currRead,genome,genomeLength,nStructs,fr);
			}
		}else{ // Different diag or seq
			// Check correct read index
			//currRead = metagenome;
			while(currRead->seqIndex != hitsList->hits[hitsList->index].seqX){
				if(currRead->next == NULL){
					fprintf(stderr, "Error searching read index.\n");
					return -1;
				}
				currRead = currRead->next;
			}
			// Generate new fragment
			FragFromHit(&frag, &hitsList->hits[hitsList->index],currRead,genome,genomeLength,nStructs,fr);
		}

		// Move to next
		hitsList->index +=1;
		// Load new hit
		if(hitsList->index >= hitsList->hits_loaded){
			if(hitsUnread[hitsList->buff] > 0){
				if(hitsList->buff != lastLoaded){
					fseek(hts,positions[hitsList->buff],SEEK_SET);
					lastLoaded = hitsList->buff;
				}
				read = loadHit(&hitsList->hits,hts,hitsUnread[hitsList->buff]);
				hitsList->index = 0;
				hitsList->hits_loaded = read;
				positions[hitsList->buff] = (uint64_t) ftell(hts);
				hitsUnread[hitsList->buff]-=read;
				checkOrder(&hitsList,false);
			}else{
				checkOrder(&hitsList,true);
				activeBuffers--;
			}
		}else{
			checkOrder(&hitsList,false);
		}
	}

	/////////////////////////////////////////////////////////////////
	fprintf(stdout, " (Generated)\n");
	/////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////
	fprintf(stdout, "\tFrags: Closing the program.\n");
	/////////////////////////////////////////////////////////////////

	// Close files
	fclose(hIndx);
	fclose(hts);
	fclose(fr);

	// Remove intermediate files
	if(removeIntermediataFiles){
		strcpy(fname,av[5]);
		remove(strcat(fname,".hts"));
		strcpy(fname,av[5]);
		remove(strcat(fname,".hindx"));
	}

	// Free malloc blocks
	free(HitsBlock);
	free(fname);

	// Free linked list
//	freeReads(&metagenome);

	// Everything finished OK
	return 0;
}