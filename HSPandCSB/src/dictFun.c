/*
 * @author Fernando Moreno Jabato <jabato@uma.es>
 * @licence all rights reserved to the author and BitLAB group (University
 *    of Malaga).
 */
#include "dict.h"


/* This function is used to shift bits in a unsigned char array
 *	@param w: word array to e shifted.
 */
void shift_word(unsigned char* w){
	int i;
	for(i=0;i<BYTES_IN_WORD-1;i++){
		w[i]<<=2;
		w[i]|=(w[i+1]>>6);
	}
	w[BYTES_IN_WORD-1]<<=2;
}


/* This function is used to store a given word in another word variable.
 *  @param container variable where word will be stored.
 *  @param word that will be stored.
 */
inline void storeWord(wentry* container,wentry word){
	container->pos = word.pos;
	container->seqIndex = word.seqIndex;
	int i;
	for(i=0;i<BYTES_IN_WORD;++i)
		container->seq[i] = word.seq[i];
}


/* This function is used to write a given set of words in two intermediate files.
 *  @param buff set of words to be written.
 *  @param index file of intermediate files.
 *  @param words file of intermediate files.
 *  @param numWords number of instances on buff set.
 *  @return a negative number if any error happens or a non negative number if
 *     the process finished correctly.
 */
int writeBuffer(wentry* buff,FILE* index,FILE* words,uint64_t numWords){
	// Sort buffer
	quicksort_W(buff,0,numWords-1);
	
	// Write buffer info on buffer index file
	uint64_t pos = (uint64_t) ftell(words); // Buffer start position
	fwrite(&pos,sizeof(uint64_t),1,index);
	fwrite(&numWords,sizeof(uint64_t),1,index); // Number of words
	// Write words on words file
	for(pos=0;pos<numWords;++pos){
		fwrite(&buff[pos].pos,sizeof(uint64_t),1,words); 
		fwrite(&buff[pos].seqIndex,sizeof(uint32_t),1,words);
		//fwrite(&buff[pos].w.WL,sizeof(uint16_t),1,words);
//		int i;
//		for(i=0;i<BYTES_IN_WORD;++i)
//			fwrite(&buff[pos].seq[i],sizeof(unsigned char),1,words);
fwrite(buff[pos].seq,sizeof(unsigned char),BYTES_IN_WORD,words);
	}
	// Buffer correctly written on intermediate files
	return 0;
}


/* This function compare two arrays of unsigned chars with the same length.
 *  @param w1: first array to be compared.
 *  @param w2: second array to be compared.
 *  @param n: length of BOTH arrays.
 *  @retun a positive number if w1>w2, a negative number if w1>w2 and zero if they are equal.
 */
int wordcmp(unsigned char *w1, unsigned char *w2, int n){
	int i;
	for(i=0;i<n;i++)
		if(w1[i] < w2[i]) return -1;
		else if(w1[i] > w2[i]) return +1;

	return 0;
}


/* This function is used to compare two wentry instances. The criterion
 * used is:
 * 		1 - Compare sequences (alphabetically).
 *      2 - Compare Read index.
 *		3 - Compare position on sequence.
 * @param w1 word to be compared.
 * @param w2 word to be compared.
 * @return a positive number if w1 is greater than w2, a negative number
 * 		if w2 is greater than w1 and zero if both are equal.
 */
int wordComparator(wentry* w1,wentry* w2){
	int wComp;
	if((wComp = wordcmp(w1->seq,w2->seq,BYTES_IN_WORD)) != 0) return wComp;

	if(w1->seqIndex > w2->seqIndex) return 1;
	else if(w1->seqIndex < w2->seqIndex) return -1;

	if(w1->pos > w2->pos) return 1;
	else if(w1->pos < w2->pos) return -1;
	
	return 0;
}


/* Function used to compare two wentry variables
 *  @param w1 word to be compared.
 *  @param w2 word to be compared
 *  @return zero if w2 are greater or equal and a positive number if
 *     w1 is greater.
 */
int GT(wentry w1, wentry w2){
	int i;
	for(i=0;i<BYTES_IN_WORD;i++)
		if(w1.seq[i] < w2.seq[i]) return 0;
		else if(w1.seq[i] > w2.seq[i]) return 1;

	if(w1.seqIndex > w2.seqIndex) return 1;
	else if(w1.seqIndex < w2.seqIndex) return 0;

	if(w1.pos > w2.pos) return 1;
	return 0;
}


/* This function is necessary for quicksort functionality.
 *  @param arr array to be sorted.
 *  @param left inde of the sub-array.
 *  @param right index of the sub-array.
 */
int partition(wentry* arr, int left, int right) {
  int i = left;
  int j = right+1;
  wentry t;

  // left sera el pivote
  // y contendra la mediana de left, r y (l+r)/2
  int mid = (left+right)/2;

  if(GT(arr[mid],arr[right])){
		SWAP_W(arr[mid],arr[right],t);
  }

  if(GT(arr[mid],arr[left])){
		SWAP_W(arr[mid],arr[left],t);
  }

  if(GT(arr[left],arr[right])){
		SWAP_W(arr[left],arr[right],t);
	}

	while(1){
		do{
			++i;
		}while(!GT(arr[i],arr[left]) && i <= right);

		do{
			--j;
		}while(GT(arr[j],arr[left]) && j >= left);

		if(i >= j) break;

		SWAP_W(arr[i],arr[j],t)
	}

	SWAP_W(arr[left],arr[j],t)

	return j;
}


/* This function is used to sort a wentry array.
 *  @param arr array to be sorted.
 *  @param left index where start to sort.
 *  @param right index where end sorting action.
 *
 */
int quicksort_W(wentry* arr, int left,int right) {
   int j;

	if( left < right ) {
 	// divide and conquer
       j = partition( arr, left, right);
       //  j=(left+r)/2;
       quicksort_W( arr, left, j-1);
       quicksort_W( arr, j+1, right);
   }
   return 0;
}


/* This function is used to check if all words has been read from intermediate file.
 *  @param unread array of unread words of each buffer segment.
 *  @param length of the unread array.
 *  @return true if there are not unread word and false in other cases.
 */
bool finished(int64_t *unread, uint64_t length){
	uint64_t i;
	for(i=0; i<length; ++i)
		if(unread[i] >= 0) return false;
	return true;
}


/* This function is used to load a word from a words intermediate file.
 *  @param word varaible where loaded wentry ill be stored.
 *  @param wFile pointer to words intermediate file.
 */
inline void loadWord(wentry *word,FILE* wFile){
	fread(&word->pos,sizeof(uint64_t),1,wFile); 
	fread(&word->seqIndex,sizeof(uint32_t),1,wFile);
	//fread(&word->w.WL,sizeof(uint16_t),1,wFile);
//	int i;
//	for(i=0;i<BYTES_IN_WORD;++i)
//		fread(&word->seq[i],sizeof(unsigned char),1,wFile);
fread(word->seq,sizeof(unsigned char),BYTES_IN_WORD,wFile);
}


/* This function is used to sort all buffers.
 *  @param 
 */
void sortBuffers(wentry*** buffers,uint64_t numBuff,uint32_t *indexes,uint64_t *wrdsInBuff){
	uint64_t i;
	bool placed = false;
	node *firstNode, *currNode;
////////////////////////////////////////////////////////////////////
fprintf(stderr, "\tTEST1\n");
////////////////////////////////////////////////////////////////////

	// First node
	node n1;
		n1.prev = NULL;
		n1.next = NULL;
		n1.buff = *buffers[0];
		n1.index = indexes[0];
		n1.numWords = wrdsInBuff[0];

	firstNode = &n1;
////////////////////////////////////////////////////////////////////
fprintf(stderr, "\tTEST2\n");
////////////////////////////////////////////////////////////////////

	// Sort using linked list
	for(i=1; i<numBuff; ++i){
		currNode = firstNode;
		placed = false;
////////////////////////////////////////////////////////////////////
fprintf(stderr, "\tTEST3\n");
////////////////////////////////////////////////////////////////////
		while(!placed){
			if(wordComparator(&(*buffers[i][indexes[i]]),&currNode->buff[currNode->index])>0){
				if(currNode->next == NULL){ // Place as tail
					node n;
					n.buff = *buffers[i];
					n.index = indexes[i];
					n.numWords = wrdsInBuff[i];
					n.prev = currNode;
					n.next = NULL;
					currNode->next = &n;
					placed = true;
				}else{ // Continue
					currNode = currNode->next;
				}
			}else{ // Place here (insert)
				node n;
				n.buff = *buffers[i];
				n.index = indexes[i];
				n.numWords = wrdsInBuff[i];
				n.prev = currNode->prev;
				if(n.prev != NULL)
					currNode->prev->next = &n;
				n.next = currNode;
				currNode->prev = &n;
				placed = true;
			}
		}
	}
////////////////////////////////////////////////////////////////////
fprintf(stderr, "\tTEST4\n");
////////////////////////////////////////////////////////////////////

	// Put right order in webtry matrix
	currNode = firstNode;
	*buffers[0] = currNode->buff;
	indexes[0] = currNode->index;
	wrdsInBuff[0] = currNode->numWords;
////////////////////////////////////////////////////////////////////
fprintf(stderr, "\tTEST5\n");
////////////////////////////////////////////////////////////////////
	for(i=1;i<numBuff;++i){
////////////////////////////////////////////////////////////////////
fprintf(stderr, "\tT\n");
////////////////////////////////////////////////////////////////////
		currNode = currNode->next;
////////////////////////////////////////////////////////////////////
fprintf(stderr, "\te\n");
////////////////////////////////////////////////////////////////////
		*buffers[i] = currNode->buff;
////////////////////////////////////////////////////////////////////
fprintf(stderr, "\ts\n");
////////////////////////////////////////////////////////////////////
		indexes[i] = currNode->index;
////////////////////////////////////////////////////////////////////
fprintf(stderr, "\tt\n");
////////////////////////////////////////////////////////////////////
		wrdsInBuff[i] = currNode->numWords;
	}
////////////////////////////////////////////////////////////////////
fprintf(stderr, "\tTEST6\n");
////////////////////////////////////////////////////////////////////

}


/* This function is used to write an entrance on dictionary files. The order of
 * each entrance on dictionaries are:
 *     - WDictionary : Word<unsigned char*BYTES_IN_WORD> + PosOnPDic<uint64_t> + NumReps<uint16_t>
 *     - PDictionary : ReadIndex<uint32_t> + PosOnRead<uint64_t>
 *  @param word to be written.
 *  @param wDic words dictioanry.
 *  @param pDic positions dictionary.
 *  @param sameThanLastWord boolean value that indicate if the current word is the same than the last written.
 *  @param words equal than last written.
 */
inline void writeWord(wentry *word, FILE* w, FILE* p, bool sameThanLastWord, uint32_t *words){
	if(!sameThanLastWord){ // Write new word
		uint64_t aux;
		fwrite(words,sizeof(uint32_t),1,w); // Write num of repetitions
//		for(aux=0;aux<BYTES_IN_WORD;++aux)
//			fwrite(&word->seq[aux],sizeof(unsigned char),1,w); // Write new word
fwrite(&word->seq,sizeof(unsigned char),BYTES_IN_WORD,w);
		aux = (uint64_t) ftell(p);
		fwrite(&aux,sizeof(uint64_t),1,w); // Write new postions on positions dictionary
		*words = 0; // Update value
	}
	fwrite(&word->seqIndex,sizeof(uint32_t),1,p); // Read index
	fwrite(&word->pos,sizeof(uint64_t),1,p); // Position on read
	*words+=1; // Increment number of repetitions
}


/* This function is used to check order and sort if it's necessary.
 *  
 */
void checkOrder(wentry*** buff,uint32_t *activeBuffers,uint32_t *indexes,uint64_t *wrdsInBuff){
	int newPos;
////////////////////////////////////////////////////////////////////
fprintf(stderr, "\tTEST1\n");
////////////////////////////////////////////////////////////////////

	// Exception -> Buffer ended
	if(indexes[0] >= wrdsInBuff[0]){  // ESTE IF GENERA UN ERROR "Violacion del segmento core"
////////////////////////////////////////////////////////////////////
fprintf(stderr, "\t\tT\n");
////////////////////////////////////////////////////////////////////

		wentry *temp = *buff[0];
		uint32_t tempIndx = indexes[0];
		uint64_t tempWrds = wrdsInBuff[0];
		int i;
		// Shift all
////////////////////////////////////////////////////////////////////
fprintf(stderr, "\t\tT\n");
////////////////////////////////////////////////////////////////////

		for(i=0; i<*activeBuffers; ++i){
			*buff[i] = *buff[i+1];
			indexes[i] = indexes[i+1];
			wrdsInBuff[i] = wrdsInBuff[i+1];
		}
////////////////////////////////////////////////////////////////////
fprintf(stderr, "\t\te\n");
////////////////////////////////////////////////////////////////////
		*buff[*activeBuffers-1] = temp;
		indexes[*activeBuffers-1] = tempIndx;
		wrdsInBuff[*activeBuffers-1] = tempWrds;
		activeBuffers--; // Update
////////////////////////////////////////////////////////////////////
fprintf(stderr, "\t\ts\n");
////////////////////////////////////////////////////////////////////

		return;
	}
////////////////////////////////////////////////////////////////////
fprintf(stderr, "\tTEST2\n");
////////////////////////////////////////////////////////////////////

	// Seek new position
	for(newPos=1; newPos < *activeBuffers-1;++newPos)
		if(wordComparator(&(*buff[0][indexes[0]]),&(*buff[newPos][indexes[newPos]])) <= 0)
			break;
////////////////////////////////////////////////////////////////////
fprintf(stderr, "\tTEST3\n");
////////////////////////////////////////////////////////////////////

	// Shift positions
	if(newPos > 0){
		wentry *temp;
		uint32_t tempIndx = indexes[0];
		uint64_t tempWrds = wrdsInBuff[0];
		int i;
		for(i=0;i<newPos-1;++i){
			*buff[i] = *buff[i+1];
			indexes[i] = indexes[i+1];
			wrdsInBuff[i] = wrdsInBuff[i+1];
		}
		*buff[newPos-1] = temp;
		indexes[newPos-1] = tempIndx;
		wrdsInBuff[newPos-1] = tempWrds;
	}	
////////////////////////////////////////////////////////////////////
fprintf(stderr, "\tTEST4\n");
////////////////////////////////////////////////////////////////////

}


/*
 */
void loadMatrix(wentry** buff,int64_t *unread,uint32_t numBuffs,uint32_t *indexes,uint64_t *wrdsInBuff,uint32_t *activeBuff,uint64_t *positions,FILE *words){
	uint32_t i,j=0,k;
	*activeBuff = 0;
	for(i=0;i<numBuffs;++i){
		if(unread[i]>0){
			wrdsInBuff[j] = 0;
			indexes[j] = 0;
			fseek(words,positions[i],SEEK_SET);

			// Load words
			for(k=0;k<MERGE_BUFFER_LENGTH && unread[i]>0;++k){
				loadWord(&buff[j][k],words);
				unread[i]--;
				wrdsInBuff[j]++;
			}

			// Update info
			positions[i] = (uint64_t) ftell(words);
			++j;
			*activeBuff+=1;
		}
	}
	sortBuffers(&buff,*activeBuff,indexes,wrdsInBuff);
}



void showWord(unsigned char* w, int wsize) {
	char Alf[] = { 'A', 'C', 'G', 'T' };
	char ws[wsize*4+4];
	int i;
	unsigned char c;
	fprintf(stdout, "Word:");
	for (i = 0; i < wsize; i++) {
		c = w[i];
		c = c >> 6;
		ws[4*i] = Alf[(int) c];
		fprintf(stdout,"%c",Alf[(int) c]);
		c = w[i];
		c = c << 2;
		c = c >> 6;
		ws[4*i+1] = Alf[(int) c];
		fprintf(stdout,"%c",Alf[(int) c]);
		c = w[i];
		c = c << 4;
		c = c >> 6;
		ws[4*i+2] = Alf[(int) c];
		fprintf(stdout,"%c",Alf[(int) c]);
		c = w[i];
		c = c << 6;
		c = c >> 6;
		ws[4*i+3] = Alf[(int) c];
		fprintf(stdout,"%c",Alf[(int) c]);
	}
	ws[wsize*4+4] = '\0';
	fprintf(stdout, "\n");
}
