// LIBRARIES
#include <stdio.h> //IO stream
#include <stdlib.h> //memory functions
#include <errno.h> //errors
#include <stdint.h> //unit64_t
#include <dirent.h> //directory functions
#include <string.h> //string tools


// VARIABLES
#define MAX_NAME_L 1024
#define MAX_GENOME_SET 40
#define MAX_METAGENOME_SET 10
static const uint64_t MAX_WORDS = 1000000; // Common with dictionary.h variable
#define BYTES_WORD 8
#define BITS_NUCLEOTIDE 2

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
    uint32_t seq;
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
    //Word compressed in binary format
    word w;
    //Ocurrence position in the sequence
    uint64_t pos;
    //Number of ocurrences inside the
    //sequence. This is used to know the
    //number of locations stored in the
    //positions file
    uint16_t num;
} hashentryNew;

typedef struct {
    // Index of read
    uint32_t readIndex;
    // Position on word dictionary
    uint64_t pos;
    // Number of different kmers
    uint16_t num;
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

typedef struct {
    //Ocurrence position in the sequence
    uint64_t pos;
    //For multiple sequence files this var
    //reflects in what sequence occurs the
    //word
    uint64_t seq;
} location; // Struct for genome d2hP files

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
    uint64_t seqX;
    //sequence number in the 'Y' file
    uint64_t seqY;
    //synteny block id
    int64_t block;
    //'f' for the forward strain and 'r' for the reverse
    char strand;
}FragFile;

// FUNCTIONS
int readGenomeSet(char*,dictionaryG**);
int readMetagenomeSet(char*,dictionaryM**);
uint64_t loadRead(FILE*,FILE*,FILE*,wentry**,int);
uint64_t loadGenome(dictionaryG,wentry**,int);
uint64_t hits(wentry*,wentry*,hit**,uint64_t,uint64_t,int);
int wordcmp(unsigned char *,unsigned char*,int); // Copied from dictionaryFun.c
int quickSort(hit*,int,int); // Copied from ""
int partition(hit*,int,int);// "" "" ""
int hitComparator(hit*,hit*); // "" "" ""
inline void SWAP(hit*,hit*,hit*);
uint64_t groupHits(hit*,uint64_t);
int calculateFragments(hit*,uint64_t,int,int,FILE*);