/*
 * @author Fernando Moreno Jabato <jabato@uma.es>
 * @licence all rights reserved to the author and BitLAB group (University
 *    of Malaga).
 */

#include <stdio.h> // IO functions
#include <stdlib.h> // Memory functions
#include <errno.h> // errors
#include <stdint.h> // unit64_t ...
#include <inttypes.h> 
#include <stdbool.h> // Boolean varaibles

// VARAIBLES
#define MAX_BUFF 10000000
#define MAX_FILE_LENGTH 1024


// STRUCTS
typedef struct{
    bool metag; // True -> Metagenome entrance; False -> Genome entrance
    uint16_t WB; // Word bytes
    unsigned char *seq; // Must be an unsigned char [WB]
    uint64_t pos; // Position on locations file
    uint32_t reps; // Number of instances of this word
} WordEntry;

typedef struct {
	// Diagonal where the hit is located
	// This value is calculated as:
	// posX - posY
	int64_t diag;
    // Ocurrence position in sequence X
    uint64_t posX;
    // Ocurrence position in sequence Y
    uint64_t posY;
    // For multiple sequence files this var
    // reflects in what sequence of X file
    // occurs the word
    uint32_t seqX;
    // For multiple sequence files this var
    // reflects in what sequence of Y file
    // occurs the word
    uint32_t seqY;
    // Length of the hit
    uint64_t length;
} Hit;

// OLD Struct for GECKO genome dictionaries
typedef struct {
    //Each letter is stored using 2 bits
    //We have 4 letters per byte and a
    //maximum of 32 in 'b'
    unsigned char b[8];
} word;

// OLD Struct for GECKO genome dictionaries
typedef struct {
    //Word compressed in binary format
    word w;
    //Ocurrence position in the locations dictionary
    uint64_t pos;
    //Number of locations stored in the
    //positions file
    uint64_t num;
} hashentry;

// OLD Struct for GECKO genome dictionaries
typedef struct {
    //Ocurrence position in the sequence
    uint64_t pos;
    //For multiple sequence files this var
    //reflects in what sequence occurs the
    //word
    uint64_t seq;
} location;

// FUNCTIONS
int wordcmp(unsigned char*,unsigned char*,int);
int WEComparer(WordEntry w1, WordEntry w2);
int readHashEntry(WordEntry*,FILE*);
int readWordEntrance(WordEntry*,FILE*,uint16_t);
int generateHits(Hit*,WordEntry,WordEntry,FILE*,FILE*,FILE*,FILE*,uint64_t*,int,int);