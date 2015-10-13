// LIBRARIES
#include <stdio.h> //IO stream
#include <stdlib.h> //memory functions
#include <errno.h>
#include <string.h> //string tools
#include <ctype.h>
#include <stdint.h> //unit64_t

// VARIABLES
static const int  WORD_SIZE = 32;
static const int BYTES_IN_WORD = 8;
static const uint64_t MAX_WORDS = 1000000;
static const int MAX_NUM_PROCS = 2;

//STRUCTS
typedef struct {
	//Each letter is stored using 2 bits
	//We have 4 letters per byte and a
	//maximum of 32 in 'b'
    unsigned char b[8];
} word;

typedef struct {
	//Word compressed in binary format
    word w;
    //Ocurrence position in the sequence
    uint64_t pos;
    //For multiple sequence files this var
    //reflects in what sequence occurs the
    //word
    uint64_t seq;
} wentry;


// FUNCTIONS
int takewords(wentry**,char*);
void shift_word(word*);
int storeWord(wentry**,wentry*,int);
int partition(wentry**,int,int);
int quickSort(wentry**,int,int);
int wordComparator(wentry*,wentry*);