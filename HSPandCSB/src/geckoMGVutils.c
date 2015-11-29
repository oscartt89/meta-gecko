#include "geckoMGVcommon.h"

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
 *  @param left: index where start to sort on left side.
 *  @param right: index where start to sort on right side.
 *  @return: a negative number if something get wrong. Zero in other cases.
 */
int quickSort_W(wentry* words, int left, int right){
	int j;

	if(left < right){
		// divide and conquer
		if((j = partition_W(words, left, right))<0) return -1;
		quickSort_W(words, left, j-1);
		quickSort_W(words, j+1, right);
	}
	return 0;
}


/* Auxiliar function for quicsort_W function.
 *  @param words: array of wentrys that will be sorted.
 *  @param left: index where start to sort on left side.
 *  @param right: index where start to sort on right side.
 *  @return: a negative number if something get wrong. Zero in other cases.
 */
int partition_W(wentry* words, int left, int right){
	int i = left;
	int j = right+1;

	// Mid will be the pivot
  	// and will store the median = (left+right)/2
	int mid = (int) (left+right)/2;

	if(wordComparator(&words[mid],&words[right]))
		SWAP_W(&words[mid],&words[right]);

	if(wordComparator(&words[mid],&words[left]))
		SWAP_W(&words[mid],&words[left]);

	if(wordComparator(&words[left],&words[right]))
		SWAP_W(&words[left],&words[right]);

	while(1){
		do{
			++i;
		}while(!wordComparator(&words[i],&words[left]) && i <= right);

		do{
			--j;
		}while(wordComparator(&words[j],&words[left]) && j >= left);

		if( i >= j ) break;

		SWAP_W(&words[i],&words[j]);
	}

	SWAP_W(&words[left],&words[j]);

	return j;
}


/* Function used to sort an array of hits.
 *  @param hits: array of hits that will be sorted.
 *  @param left: index where start to sort on left side.
 *  @param right: index where start to sort on right side.
 *  @return: a negative number if something get wrong. Zero in other cases.
 */
int quickSort_H(hit* hits, int left, int right){
  int j;

  if(left < right){
    // divide and conquer
    if((j = partition_H(hits, left, right))<0) return -1;
    quickSort_H(hits, left, j-1);
    quickSort_H(hits, j+1, right);
  }
  return 0;
}


/* Auxiliar function for quicsort function.
 *  @param hits: array of hits that will be sorted.
 *  @param left: index where start to sort on left side.
 *  @param right: index where start to sort on right side.
 *  @return: a negative number if something get wrong. Zero in other cases.
 */
int partition_H(hit* hits, int left, int right){
  int i = left;
  int j = right+1;

  // Mid will be the pivot
  // and will store the median = (left+right)/2
  int mid = (int) (left+right)/2;

  if(hitComparator(&hits[mid],&hits[right]))
    SWAP_H(&hits[mid],&hits[right]);

  if(hitComparator(&hits[mid],&hits[left]))
    SWAP_H(&hits[mid],&hits[left]);

  if(hitComparator(&hits[left],&hits[right]))
    SWAP_H(&hits[left],&hits[right]);

  while(1){
    do{
      ++i;
    }while(!hitComparator(&hits[i],&hits[left]) && i <= right);

    do{
      --j;
    }while(hitComparator(&hits[j],&hits[left]) && j >= left);

    if( i >= j ) break;

    SWAP_H(&hits[i],&hits[j]);
  }

  SWAP_H(&hits[left],&hits[j]);

  return j;
}