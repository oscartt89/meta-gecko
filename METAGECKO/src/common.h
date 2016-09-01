#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>

int wordcmpbytearray(const unsigned char *w1, const unsigned char *w2, int bytes_to_check);
void shift_byte_array_left(unsigned char * b, unsigned int bytes_to_check);
void shift_byte_array_right(unsigned char * b, unsigned int bytes_to_check);
void QuickSortByteArray(unsigned char * array, uint64_t x, uint64_t y, unsigned int bytes);
void showByteArray(const unsigned char * b, char * kmer, uint16_t WORD_LENGTH);
int64_t binarySearchByteArray(unsigned char * byteword, unsigned char * bytearray, unsigned int bytes, uint64_t nMers);
void printTable(unsigned char * array, uint64_t nSequences, unsigned int bytes, FILE * out);
uint64_t asciiToUint64(const char *text);