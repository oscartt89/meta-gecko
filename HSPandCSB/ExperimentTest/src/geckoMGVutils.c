#include "geckoMGVutils.h"

/* This function is used to swap/interchange two wentry instances.
 *  @param w1: wentry that will be swapped.
 *  @param w2: wentry that will be swapped.
 */
inline void SWAP_W(wentry *w1,wentry *w2){
	wentry t;
	memcpy(&t,w1,sizeof(wentry));
	memcpy(w1,w2,sizeof(wentry));
	memcpy(w2,&t,sizeof(wentry));
}


/* This function is used to swap/interchange two hits instances.
 *  @param h1: hit that will be swapped.
 *  @param h2: hit that will be swapped.
 */
inline void SWAP_H(hit *h1,hit *h2){
  hit t;
  memcpy(&t,h1,sizeof(hit));
  memcpy(h1,h2,sizeof(hit));
  memcpy(h2,&t,sizeof(hit));
}


/* This function compare two arrays of unsigned chars with the same length.
 *  @param w1: first array to be compared.
 *  @param w2: second array to be compared.
 *  @param n: length of BOTH arrays.
 *  @retun a positive number if w1>w2, a negative number if w1>w2 and zero if they are equal.
 */
int wordcmp(unsigned char *w1, unsigned char*w2, int n) {
	int i;
	for (i=0;i<n;i++) {
		if (w1[i]<w2[i]) return -1;
		if (w1[i]>w2[i]) return +1;
	}
	return 0;
}


/* This function is used to compare two wentry instances. The criterion
 * used is:
 * 		1 - Compare sequences (alphabetically).
 * 		2 - Compare position on sequence.
 * @param w1 word to be compared.
 * @param w2 word to be compared.
 * @return a positive number if w1 is greater than w2, a negative number
 * 		if w2 is greater than w1 and zero if both are equal.
 */
int wordComparator(wentry* w1,wentry* w2){
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
	else if(w1->pos < w2->pos) return -1;

	return 0;
}


/* This function is used to compare two hit instances. The criterion
 * used is:
 *    1 - Compare sequence 1 index.
 *    2 - Compare sequence 2 index.
 *    3 - Compare sequence 1 start.
 *    4 - Compare sequence 2 start.
 *    2 - Compare length.
 * @param h1 hit to be compared.
 * @param h2 hit to be compared.
 * @return a positive number if h1 is greater than h2, a negative number
 *    if h2 is greater than h1 and zero if both are equal.
 */
int hitComparator(hit* h1,hit* h2){
  if(h1->seq1 > h2->seq1) return 1;
  else if(h1->seq1 < h2->seq1) return -1;

  if(h1->seq2 > h2->seq2) return 1;
  else if(h1->seq2 < h2->seq2) return -1;

  if(h1->start1 > h2->start1) return 1;
  else if(h1->start1 < h2->start1) return -1;

  if(h1->start2 > h2->start2) return 1;
  else if(h1->start2 < h2->start2) return -1;

  if(h1->length > h2->length) return 1;
  else if(h1->length < h2->length) return -1;

  return 0;
}


/* Function used to sort an array of wentrys.
 *  @param words: array of wentry that will be sorted.
 *  @param start: index where start to sort.
 *  @param length: length of the array to sort (starting on "start" index).
 */
void quickSort_W(wentry* words, uint64_t start, uint64_t length){
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


/* Function used to sort an array of hits.
 *  @param hits: array of hit that will be sorted.
 *  @param start: index where start to sort.
 *  @param length: length of the array to sort (starting on "start" index).
 */
void quickSort_H(hit* hits, uint64_t start, uint64_t length){
  //Check exceptions
  if (length < 2) return;
  else if(length == 2){
    if(hitComparator(&hits[start],&hits[start+1])>0)
      SWAP_H(&hits[start],&hits[start+1]);
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
      while(hitComparator(&hits[right],&hits[pivot])<0 && right <= pivot) right++;
      // While left y greater than pivot, continue
      while(hitComparator(&hits[pivot],&hits[left])<0 && left >= pivot) left--;
      // All array checked -> end
      if(right>=left)
          break;
      // Swap selected values
      SWAP_H(&hits[right],&hits[left]);
      changed = 1;
      if(right < pivot) right++;
      if(left > pivot) left--;
  }

  if(changed && right<left){ // If anything changed, don't continue, it's sorted
    quickSort_H(hits,start,right);
    quickSort_H(hits,right,length-right);
  }

  return;
}