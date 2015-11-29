#include <stdio.h> //IO stream
#include <stdlib.h> //memory functions
#include <errno.h> // errors
#include <string.h> //string tools
#include <stdint.h> //unit64_t

// VARAIBLES
#define BYTES_IN_WORD 8
#define MAX_FILE_LENGTH 1024
#define MAX_WORDS 10000
#define BITS_NUCLEOTIDE 2

//STRUCTS
typedef struct {
	//Each letter is stored using 2 bits
	//We have 4 letters per byte and a
	//maximum of 32 in 'b'
    unsigned char b[BYTES_IN_WORD];
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
    //Ocurrence position in the position dictionary
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
} Read;

typedef struct{
    // Start position of hit on seq1
    uint64_t start1;
    // Start position of hit on seq2
    uint64_t start2;
    // Length of hit
    uint64_t length;
    // Sequence 1
    uint32_t seq1;
    // Sequence 2;
    uint32_t seq2;
} hit;

// FUNCTIONS
// Swap functions
inline void SWAP_W(wentry*,wentry*);
inline void SWAP_H(hit*,hit*);
// Comparation functions
int wordcmp(unsigned char*,unsigned char*,int);
int wordComparator(wentry*,wentry*);
int hitComparator(hit*,hit*);
// Sort functions
int quickSort_W(wentry*,int,int);
int partition_W(wentry*,int,int);
int quickSort_H(hit*,int,int);
int partition_H(hit*,int,int);