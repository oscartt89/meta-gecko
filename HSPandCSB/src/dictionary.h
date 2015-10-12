// LIBRARIES
#include <stdio.h> //IO stream
#include <stdlib.h> //memory functions
#include <errno.h>
#include <string.h> //string tools
#include <ctype.h>
#include <stdint.h> //unit64_t

// VARIABLES
static int WORD_SIZE = 32;
static int BYTES_IN_WORD = 8;

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
// dictionary.c 
void takewords(wentry**,char*);

//DictionaryFun.c
void shift_word(word*);
int storeWord(wentry**,wentry*,int);