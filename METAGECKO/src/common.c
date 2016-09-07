#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>


/* 	This function compares two arrays of unsigned chars with the same length.
 *  @param w1: first array to be compared.
 *  @param w2: second array to be compared.
 *  @param bytes_to_check: length of BOTH arrays in bytes
 *  @retun a positive number if w1>w2, a negative number if w1<w2 and zero if they are equal.
 */
 
int wordcmpbytearray(const unsigned char *w1, const unsigned char *w2, int bytes_to_check) {
	int i;
	for (i=0;i<bytes_to_check;i++) {
		if (w1[i]<w2[i]) return -1;
		if (w1[i]>w2[i]) return +1;
	}
	return 0;
}


/* 	This function is used to shift left bits in a unsigned char array
 *	@b: char array representing the word using compressed 2-bit characters
 */
 
void shift_byte_array_left(unsigned char * b, unsigned int bytes_to_check) {
    unsigned int i;
    for (i = 0; i < bytes_to_check - 1; i++) {
        b[i] <<= 2;
        b[i] |= (b[i + 1] >> 6);
    }
    b[bytes_to_check - 1] <<= 2;
}

/* 	This function is used to shift right bits in a unsigned char array
 *	@b: char array representing the word using compressed 2-bit characters
 */
void shift_byte_array_right(unsigned char * b, unsigned int bytes_to_check) {
    unsigned int i;
    for (i = bytes_to_check - 1; i > 0; i--) {
        b[i] >>= 2;
        b[i] |= (b[i - 1] << 6);
    }
    b[i] >>= 2;
}

/*  Sorts a vector of unsigned char bytes where each individual value takes up "bytes" bytes
 *  @array: The unsigned char vector to sort
    @x:     Starting position to sort the vector (use 0)
    @y:     Ending position to sort
    @bytes: How many bytes will be used for each individual value in the vector
 */

void QuickSortByteArray(unsigned char * array, uint64_t x, uint64_t y, unsigned int bytes) {

    unsigned char pivote[bytes], aux[bytes];
    uint64_t x1, y1;

    memcpy(&pivote[0], array+(((x+y)/2)*bytes), bytes);
    x1 = x;
    y1 = y;

    do{ 
        while (wordcmpbytearray(pivote, array+x1*bytes, bytes) > 0) x1++;
        while (wordcmpbytearray(pivote, array+y1*bytes, bytes) < 0) y1--;
        if (x1 < y1) { 

            memcpy(&aux[0], array+x1*bytes, bytes);
            memcpy(array+x1*bytes, array+y1*bytes, bytes);
            memcpy(array+y1*bytes, &aux[0], bytes);
            x1++;
            y1--;
        }else if (x1 == y1) x1++;
    } while (x1 <=y1);

    if (x<y1) QuickSortByteArray(array,x,y1,bytes);
    if (x1<y) QuickSortByteArray(array,x1,y,bytes);
}

/* Translates an unsigned char word of ACTG letters compressed in 2 bits to a readable char
 *  @b: Unsigned array of chars compressed
 *  @kmer: Char array to receive readable translation
 *  @WORD_LENGTH: Length in bits of the word to translate
 */

void showByteArray(const unsigned char * b, char * kmer, uint16_t WORD_LENGTH) {
    char Alf[] = { 'A', 'C', 'G', 'T' };
    int i;
    int wsize = WORD_LENGTH/4;
    unsigned char c;
    for (i = 0; i < wsize; i++) {
        c = b[i];
        c = c >> 6;
        kmer[4*i] = Alf[(int) c];
        c = b[i];
        c = c << 2;
        c = c >> 6;
        kmer[4*i+1] = Alf[(int) c];
        c = b[i];
        c = c << 4;
        c = c >> 6;
        kmer[4*i+2] = Alf[(int) c];
        c = b[i];
        c = c << 6;
        c = c >> 6;
        kmer[4*i+3] = Alf[(int) c];
    }
    kmer[32]='\0';
}

/* Performs a binary search on an unsigned char array comparing "bytes" bytes at a time
 *  @byteword:  A pointer to the word that will be compared
 *  @bytearray: The vector holding the words
 *  @bytes:     How many bytes per value
    @t_first:     index (multiplied by the number of bytes) of the first kmer
    @t_last:      index (multiplied by the number of bytes) of the last kmer

    Returns
        The position of the byteword in the array if it was found, if not, -1 is returned
 */

int64_t binarySearchByteArray(unsigned char * byteword, unsigned char * bytearray, unsigned int bytes, uint64_t t_first, uint64_t t_last) {

    //char kmer1[32], kmer2[32];
    //showByteArray(byteword, kmer1, 32);


    int64_t first = t_first;
    int64_t last = t_last - 1;
    int64_t middle = (first+last)/2;
    int compare;

    while (first <= last) {
        
        //showByteArray(bytearray+(middle*bytes), kmer2, 32);

        //printf("Comparing:\n%s\n%s  ", kmer1, kmer2);
        //printf("accessing %"PRIu64" first:%"PRIu64", last: %"PRIu64"\n", middle*bytes, first, last);
        compare = wordcmpbytearray(byteword, bytearray+(middle*bytes), bytes);
        //printf("yields: %d\n", compare);
        if (compare == 0) return middle;
        if (compare < 0) last=middle - 1; 
        else first = middle + 1;

        middle = (first + last)/2;
    }
    return -1;
}




/* Prints the kmer table values converted to readable format (ascii)
 *  @array: The vector holding the words
    @nSequences:    The number of words
    @bytes: Bytes per word
    @out:   File handler to where printing will be redirected
 *  
 */

void printTable(unsigned char * array, uint64_t nSequences, unsigned int bytes, FILE * out){
    uint64_t i;
    char kmer[32];
    for(i=0;i<nSequences;i++){
        showByteArray((unsigned char *)array+(i*bytes), kmer, 32);
        fprintf(out, "%s->%"PRIu64"\n", kmer, i);
    }
}

/*  Converts ascii to uint64_t
    @text:  The char vector containing the number in ASCII

    Returns the value converted to uint64_t
*/

uint64_t asciiToUint64(const char *text){
    uint64_t number=0;

    for(;*text;text++){
        char digit=*text-'0';           
        number=(number*10)+digit;
    }
    return number;
}

