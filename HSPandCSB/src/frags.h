// LIBRARIES
#include <stdio.h> //IO stream
#include <stdlib.h> //memory functions
#include <errno.h> //errors
#include <stdint.h> //unit64_t
#include <dirent.h> //directory functions
#include <string.h> //string tools


// VARIABLES
#define MAX_NAME_L 1024
#define MAX_GENOME_SET 20
static const uint64_t MAX_WORDS = 1000000; // Common with dictionary.h variable


// STRUCTS
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

typedef struct {
    //Word compressed in binary format
    word w;
    //Ocurrence position in the sequence
    uint64_t pos;
    //Number of ocurrences inside the
    //sequence. This is used to know the
    //number of locations stored in the
    //positions file
    uint64_t num;
} hashentry;

typedef struct {
    // Index of read
    uint64_t readIndex;
    // Position on word dictionary
    uint64_t pos;
    // Number of different kmers
    uint64_t num;
} READ; // "read" is used on dirent.h

typedef struct {
	// Name of dictionaries
	char name[MAX_NAME_L];
	// Words dictionary
	char W[MAX_NAME_L];
	// Positions dictionary
	char P[MAX_NAME_L];
} dictionaryG; // Genome dictionary

typedef struct {
	// Name of dictionaries
	char name[MAX_NAME_L];
	// Words dictionary
	char W[MAX_NAME_L];
	// Positions dictionary
	char P[MAX_NAME_L];
	// Reads dictionary
	char R[MAX_NAME_L];
} dictionaryM; // Metagenome dictionary

// FUNCTIONS
int readGenomeSet(char*,wentry**,int*);