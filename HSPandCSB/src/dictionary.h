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
static const int MAX_FILE_LENGTH = 1024;

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
    uint32_t seq;
} wentry;

typedef struct {
    //Word compressed in binary format
    word w;
    //Position of first location on Positions Dictionary
    uint64_t pos;
    //Number of ocurrences inside the
    //sequence. This is used to know the
    //number of locations stored in the
    //positions file
    uint16_t num;
} hashentry;

typedef struct {
    // Index of read
    uint32_t readIndex;
    // Position on word dictionary
    uint64_t pos;
    // Number of different kmers
    uint16_t num;
} read;

// FUNCTIONS
int createDictionary(char*,char*);
//
void shift_word(word*);
int storeWord(wentry*,wentry*,int);
int partition(wentry*,int,int);
int quickSort(wentry*,int,int);
int wordComparator(wentry*,wentry*);
void writeDic(wentry*,int,FILE*,FILE*,FILE*);
int wordcmp(unsigned char*,unsigned char*,int); 
inline void SWAP(wentry*,wentry*,wentry*);