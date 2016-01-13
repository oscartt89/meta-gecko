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


/* This method is used to store a wentry instance into an array of wentry. 
 * @param wArr is the array of wentry** where the new wentry will be stored.
 * @param word is the new word to be stored.
 * @param length is the new length of the array.
 * @param maxStored is a variable used to know if memory allocation is necessary or not.
 * @return 0 if the process ends correctly. If something got wrong return a
 * 		negative value.
 */
int storeWord(wentry *wArr,wentry word,uint64_t length,uint64_t *maxStored){
	// Allocate memory if it's necessary
	if((length > *maxStored) == 0){ // Limit of words
		if((wArr[length-1].w.b = (unsigned char*)malloc(sizeof(unsigned char)*BYTES_IN_WORD))==NULL){
			fprintf(stderr, "storeWord:: Error allocating space for word\n");
			return -1;
		}
		*maxStored++;
	}
	// Copy the new word
	memcpy(&wArr[length-1],&word,sizeof(wentry));
	return 0;
}


/* This function is used to write a given set of words in two intermediate files.
 *  @param buff set of words to be written.
 *  @param index file of intermediate files.
 *  @param positions file of intermediate files.
 *  @param numWords number of instances on buff set.
 *  @return a negative number if any error happens or a non negative number if
 *     the process finished correctly.
 */
int writeBuffer(wentry* buff,FILE* index,FILE* positions,uint64_t numWords){
	// Sort buffer
	quickSort_W(buff,0,numWords);
	
	// Write buffer info on buffer index file
	uint64_t pos = (uint64_t) ftell(positions); // Buffer start position
	fwrite(&pos,sizeof(uint64_t),1,index);
	fwrite(&numWords,sizeof(uint64_t),1,index); // Number of words
	// Write words on words file
	fwrite(buff,sizeof(buff[0]),numWords,positions);

	// Buffer correctly written on intermediate files
	return 0;
}


/* This function compare two arrays of unsigned chars with the same length.
 *  @param w1: first array to be compared.
 *  @param w2: second array to be compared.
 *  @param n: length of BOTH arrays.
 *  @retun a positive number if w1>w2, a negative number if w1>w2 and zero if they are equal.
 */
int wordcmp(unsigned char *w1, unsigned char*w2, int n){
	int i;
	for(i=0;i<n;i++)
		if(w1[i] < w2[i]) return -1;
		else if(w1[i] > w2[i]) return +1;

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
		memcpy(&t,w1,sizeof(wentry));
		memcpy(w1,w2,sizeof(wentry));
		memcpy(w2,&t,sizeof(wentry));
		free(t.w.b);
	}
}


/* This function is used to compare two wentry instances. The criterion
 * used is:
 * 		1 - Compare sequences (alphabetically).
 * 		2 - Compare position on sequence.
 *      3 - Compare Read index.
 * @param w1 word to be compared.
 * @param w2 word to be compared.
 * @return a positive number if w1 is greater than w2, a negative number
 * 		if w2 is greater than w1 and zero if both are equal.
 */
int wordComparator(wentry* w1,wentry* w2){
	int wComp;
	if((wComp = wordcmp(&w1->w.b[0],&w2->w.b[0],BYTES_IN_WORD)) != 0) return wComp;

	if(w1->pos > w2->pos) return 1;
	else if(w1->pos < w2->pos) return -1;

	if(w1->seq > w2->seq) return 1;
	else if(w1->seq < w2->seq) return -1;

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
      // Swap selected values
      SWAP_W(&words[right],&words[left]);
      changed = 1;
      if(right < pivot) right++;
      if(left > pivot) left--;
  }

  if(changed && right<left){ // If anything changed, don't continue, it's sorted
  	quickSort_W(words,start,right);
  	quickSort_W(words,right,length-right);
  }

  return;
}


/* This function is used to deallocate a buffer from disk
 *  @param buff buffer to be deallocated.
 *  @param maxStored length of the buffer.
 */
void freeBuffer(wentry* buff, uint64_t maxStored){
	int k;
	for(k=0; k<maxStored; ++k)
		free(buff[k].w.b);
	free(buff);
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


/* This function is used to write an entrance on dictionary files
 *  @param word to be written.
 *  @param wDic words dictioanry.
 *  @param pDic positions dictionary.
 *  @param sameThanLastWord boolean value that indicate if the current word is the same than the last written.
 *  @param words equal than last written.
 */
inline void writeWord(wentry *word, FILE* w, FILE* p, bool sameThanLastWord, uint16_t *words){
	//if()
}