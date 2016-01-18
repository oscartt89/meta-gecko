/*
 * @author Fernando Moreno Jabato <jabato@uma.es>
 * @licence all rights reserved to the author and BitLAB group (University
 *    of Malaga).
 */
#include "dict.h"

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


/* This function is used to store a given word in another word variable.
 *  @param container variable where word will be stored.
 *  @param word taht will be stored.
 */
inline void storeWord(wentry* container,wentry word){
	container->pos = word.pos;
	container->seq = word.seq;
	container->w.WL = word.w.WL;
	int i;
	for(i=0;i<BYTES_IN_WORD;++i)
		container->w.b[i] = word.w.b[i];
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
///////////////////////////////////////////////////////////////////////////
int i;
for(i=0;i<numWords;++i){
	fprintf(stdout, "S:%"PRIu32"  P:%"PRIu64" ",buff[i].seq,buff[i].pos);
	showWord(&buff[i].w,BYTES_IN_WORD);
}
//////////////////////////////////////////////////////////////////////////
	// Sort buffer
	quickSort_W(buff,0,numWords);
///////////////////////////////////////////////////////////////////////////
fprintf(stderr, "||\n");
for(i=0;i<numWords;++i){
	fprintf(stdout, "S:%"PRIu32"  P:%"PRIu64" ",buff[i].seq,buff[i].pos);
	showWord(&buff[i].w,BYTES_IN_WORD);
}
//////////////////////////////////////////////////////////////////////////
	
	// Write buffer info on buffer index file
	uint64_t pos = (uint64_t) ftell(words); // Buffer start position
	fwrite(&pos,sizeof(uint64_t),1,index);
	fwrite(&numWords,sizeof(uint64_t),1,index); // Number of words
	// Write words on words file
	for(pos=0;pos<numWords;++pos){
		fwrite(&buff[pos].pos,sizeof(uint64_t),1,words); 
		fwrite(&buff[pos].seq,sizeof(uint32_t),1,words);
		//fwrite(&buff[pos].w.WL,sizeof(uint16_t),1,words);
		int i;
		for(i=0;i<BYTES_IN_WORD;++i)
			fwrite(&buff[pos].w.b[i],sizeof(unsigned char),1,words);
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
int wordcmp(word w1, word w2, int n){
	int i;
	for(i=0;i<n;i++)
		if(w1.b[i] < w2.b[i]) return -1;
		else if(w1.b[i] > w2.b[i]) return +1;

	return 0;
}


/* This function is used to swap/interchange two wentry instances.
 *  @param w1: wentry that will be swapped.
 *  @param w2: wentry that will be swapped.
 */
inline void SWAP_W(wentry *w1,wentry *w2){
	wentry t;
	if((t.w.b = (unsigned char*)malloc(sizeof(unsigned char)*BYTES_IN_WORD))==NULL){
		fprintf(stderr, "SWAP_W:: Error allocating memory.\n");
	}else{
		storeWord(&t,*w1);
		storeWord(w1,*w2);
		storeWord(w2,t);
		free(t.w.b);
	}
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
	if((wComp = wordcmp(w1->w,w2->w,BYTES_IN_WORD)) != 0) return wComp;

	if(w1->seq > w2->seq) return 1;
	else if(w1->seq < w2->seq) return -1;

	if(w1->pos > w2->pos) return 1;
	else if(w1->pos < w2->pos) return -1;
	
	return 0;
}


/* Function used to sort an array of wentrys.
 *  @param words: array of wentry that will be sorted.
 *  @param start: index where start to sort.
 *  @param length: length of the array to sort (starting on "start" index).
 */
void quickSort_W(wentry *words, uint64_t start, uint64_t length){
	//Check exceptions
  if (length < 2) return;
  else if(length == 2){
  	if(wordComparator(&words[start],&words[start+1])>0)
  		SWAP_W(&words[start],&words[start+1]);
  	return;
  }
  if(start < 0) return;
	uint64_t right, left, pivot, end = start + length - 1;
	int changed = 0;

  pivot = rand() % (end + 1 - start) + start;
  right = start;
  left = end;

  while(right != pivot || left != pivot) {
  	// While right is lower than pivot, continue
      while(wordComparator(&words[right],&words[pivot])<0 && right <= pivot) right++;
      // While left y greater than pivot, continue
      while(wordComparator(&words[pivot],&words[left])<0 && left >= pivot) left--;
      // All array checked -> end
      if(right>=left)
          break;
      
      if(right > pivot){ // Right sub-array is sorted
      	quickSort_W(words,pivot,length-pivot);
      	changed = 0;
      }else if(left < pivot){ // Left sub-array is sorted
      	quickSort_W(words,start,pivot);
      	changed = 0;
      }else{
	    // Swap selected values
	    SWAP_W(&words[right],&words[left]);
	    changed = 1;

	    if(pivot == right && pivot != end) pivot++;
	  	else if(pivot == left && pivot > 0) pivot--;
	  }
      
      if(right < pivot) right++;
      if(left > pivot) left--;
  }

  if(changed && right<left){ // If anything changed, don't continue, it's sorted
  	quickSort_W(words,start,right);
  	quickSort_W(words,right,length-right);
  }else if(changed) quickSort_W(words,start,end);


  return;
}


/* This function is used to check if all words has been read from intermediate file.
 *  @param unread array of unread words of each buffer segment.
 *  @param length of the unread array.
 *  @return true if there are not unread word and false in other cases.
 */
bool finished(uint64_t *unread, uint64_t length){
	uint64_t i;
	for(i=0; i<length; ++i)
		if(unread[i] > 0) return false;
	return true;
}


/* This function is used to load a word from a words intermediate file.
 *  @param word varaible where loaded wentry ill be stored.
 *  @param wFile pointer to words intermediate file.
 */
inline void loadWord(wentry *word,FILE* wFile){
	fread(&word->pos,sizeof(uint64_t),1,wFile); 
	fread(&word->seq,sizeof(uint32_t),1,wFile);
	//fread(&word->w.WL,sizeof(uint16_t),1,wFile);
	int i;
	for(i=0;i<BYTES_IN_WORD;++i)
		fread(&word->w.b[i],sizeof(unsigned char),1,wFile);
}


/* This function is used to return the index of the lower word on a wentry array.
 *  @param words array where search.
 *  @param length of words array.
 *  @return the index of the lowest wentry on words array.
 */
uint64_t lowestWord(wentry *words,uint64_t length){
	uint64_t i,j=0;
	for(i=1;i<length;++i)
		if(wordComparator(&words[j],&words[i]) > 0) j = i;
	return j;
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
inline void writeWord(wentry *word, FILE* w, FILE* p, bool sameThanLastWord, uint16_t *words){
	uint64_t aux;
	if(!sameThanLastWord){ // Write new word
		fwrite(words,sizeof(uint16_t),1,w); // Write num of repetitions
		fwrite(word->w.b,sizeof(unsigned char),BYTES_IN_WORD,w); // Write new word
		aux = (uint64_t) ftell(p);
		fwrite(&aux,sizeof(uint64_t),1,w); // Write new postions on positions dictionary
		*words = 0; // Update value
	}
	fwrite(&word->seq,sizeof(uint32_t),1,p); // Read index
	fwrite(&word->pos,sizeof(uint64_t),1,p); // Position on read
	*words+=1; // Increment number of repetitions
}


/* This funtion is used to deallocate a wentry array in disk.
 *  @param arr wentry array to be deallocated.
 *  @param length of wentry array.
 */
inline void freeWArray(wentry *arr,uint64_t length){
	uint64_t i;
	for(i=0;i<length;++i)
		free(arr[i].w.b);
	free(arr);
}


void showWord(word* w, int wsize) {
	char Alf[] = { 'A', 'C', 'G', 'T' };
	char ws[wsize*4+4];
	int i;
	unsigned char c;
	fprintf(stdout, "Word:");
	for (i = 0; i < wsize; i++) {
		c = w->b[i];
		c = c >> 6;
		ws[4*i] = Alf[(int) c];
		fprintf(stdout,"%c",Alf[(int) c]);
		c = w->b[i];
		c = c << 2;
		c = c >> 6;
		ws[4*i+1] = Alf[(int) c];
		fprintf(stdout,"%c",Alf[(int) c]);
		c = w->b[i];
		c = c << 4;
		c = c >> 6;
		ws[4*i+2] = Alf[(int) c];
		fprintf(stdout,"%c",Alf[(int) c]);
		c = w->b[i];
		c = c << 6;
		c = c >> 6;
		ws[4*i+3] = Alf[(int) c];
		fprintf(stdout,"%c",Alf[(int) c]);
	}
	ws[wsize*4+4] = '\0';
	fprintf(stdout, "\n");
}