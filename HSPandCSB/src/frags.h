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
#include <string.h> // String functions
#include <stdbool.h> // Boolean varaibles

// VARAIBLES
#define MAX_BUFF 10000000
#define MAX_FILE_LENGTH 1024

// LINE FUNCTIONS
#define SWAP(a,b,t) t=a; a=b; b=t;

// GLOBAL VARAIBLES
uint64_t buffersWritten;
uint64_t S_Threshold; // Similarity threshold

// STRUCTS
typedef struct{
    bool metag; // True -> Metagenome entrance; False -> Genome entrance
    uint16_t WB; // Word bytes
    unsigned char *seq; // Must be an unsigned char [WB]
    uint64_t pos; // Position on locations file
    uint32_t reps; // Number of instances of this word
} WordEntry;

typedef struct{
    uint32_t seq; // Sequence index
    uint64_t pos; // Position on sequence
} LocationEntry;

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

typedef struct{
    //Diagonal where the frag is located
    //This value is calculated as:
    //posX - posY
    int64_t diag;
    //Start position in sequence X
    uint64_t xStart;
    //Start position in Sequence Y
    uint64_t yStart;
    //End position in Sequence X
    uint64_t xEnd;
    //End position in Sequence Y
    uint64_t yEnd;
    //Fragment Length
    //For ungaped aligment is:
    //xEnd-xStart+1
    uint64_t length;
    //Number of identities in the
    //fragment
    uint64_t ident;
    //Score of the fragment. This
    //depends on the score matrix
    //used
    uint64_t score;
    //Percentage of similarity. This
    //is calculated as score/scoreMax
    //Where score max is the maximum
    //score possible
    float similarity;
    //sequence number in the 'X' file
    /*uint64_t*/uint32_t seqX;
    //sequence number in the 'Y' file
    /*uint64_t*/uint32_t seqY;
    //synteny block id
    int64_t block;
    //'f' for the forward strain and 'r' for the reverse
    char strand;
} FragFile;

//
struct item{
    Hit hits;
    uint64_t buff;
    struct item *next;
};

typedef struct item node; 

// FUNCTIONS
int wordcmp(unsigned char*,unsigned char*,int);
int HComparer(Hit w1, Hit w2);
void readHashEntry(WordEntry*,FILE*);
int readWordEntrance(WordEntry*,FILE*,uint16_t);
int generateHits(Hit*,WordEntry,WordEntry,FILE*,FILE*,FILE*,FILE*,uint64_t*);
void loadLocationEntrance(LocationEntry*,FILE*,uint32_t,bool); 
inline void storeHit(Hit*,LocationEntry,LocationEntry,uint64_t);
void writeHitsBuff(Hit*,FILE*,FILE*,uint64_t);
int GT(Hit,Hit);
int partition(Hit*,int,int);
void quicksort_H(Hit*,int,int);
inline void loadHit(Hit*,FILE*);
uint64_t lowestHit(Hit*,uint64_t,int64_t*);
bool finished(int64_t*,uint64_t);
inline void writeFragment(FragFile,FILE*);
void SWAP_H(Hit*,Hit*,Hit);
void copyHit(Hit*,Hit);
void push(node**,node**);
void move(node**,node**);
void sortList(node**);
