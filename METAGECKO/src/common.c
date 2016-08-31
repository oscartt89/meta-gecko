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
    kmer[WORD_LENGTH]='\0';
}


int64_t binarySearchByteArray(unsigned char * byteword, unsigned char * bytearray, unsigned int bytes, uint64_t nMers) {

    //char kmer1[32], kmer2[32];

    int64_t first = 0;
    int64_t last = nMers - 1;
    int64_t middle = (first+last)/2;
    int compare;

    while (first <= last) {
        //showByteArray(byteword, kmer1, 32);
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


void printTable(unsigned char * array, uint64_t nSequences, unsigned int bytes, FILE * out){
    uint64_t i;
    char kmer[32];
    for(i=0;i<nSequences/bytes;i++){
        showByteArray((unsigned char *)array+(i*bytes), kmer, 32);
        fprintf(out, "%s\n", kmer);
    }
}

