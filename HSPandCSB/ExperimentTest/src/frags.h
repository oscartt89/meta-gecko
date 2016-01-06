// LIBRARIES
#include "geckoMGVutils.h"
#include <dirent.h> //directory functions


// VARIABLES
#define MAX_GENOME_SET 40
#define MAX_METAGENOME_SET 10
#define MAX_HITS 1000
#define SCORE 4
#define MAX_FILES 30

// STRUCTS
typedef struct {
    //Word compressed in binary format
    word w;
    //Ocurrence position in the position dictionary
    uint64_t pos;
    //Number of ocurrences inside the
    //sequence. This is used to know the
    //number of locations stored in the
    //positions file
    uint64_t num;
} hashentryOld;

typedef struct{
    // Word compressed in binary format
    word w;
    // Sequence
    uint32_t seq;
    // Num of instances on location array
    uint16_t num;
    // Locations
    uint64_t* locations;
}HE;

typedef struct {
	// Name of dictionaries
	char name[MAX_FILE_LENGTH];
	// Words dictionary
	char W[MAX_FILE_LENGTH];
	// Positions dictionary
	char P[MAX_FILE_LENGTH];
} dictionaryG; // Genome dictionary

typedef struct {
	// Name of dictionaries
	char name[MAX_FILE_LENGTH];
	// Words dictionary
	char W[MAX_FILE_LENGTH];
	// Positions dictionary
	char P[MAX_FILE_LENGTH];
	// Reads dictionary
	char R[MAX_FILE_LENGTH];
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
// Files and directories hsndler functions
int readGenomeSet(char*,dictionaryG**);
int readMetagenomeSet(char*,dictionaryM**);
int startsWithDot(char*);
int endsWith(char*,char*);
int listFiles(char*, char***);
// Loaders
uint64_t loadRead(FILE*,FILE*,FILE*,HE**,int);
int64_t loadGenome(dictionaryG,HE**,int);
// Hits functions
int64_t hits(HE*,HE*,hit**,uint64_t,uint64_t,int);
int64_t groupHits(hit*,uint64_t);
// Fragment functions
int calculateFragments(hit*,uint64_t,int,int,FILE*);
inline void storeFragFile(FragFile*,hit*,float);
// Space handlers
inline void free_HE(HE*,int);
inline void free_Files(char**,int);
// Sort functions
inline void SWAP_S(char*,char*);
void quickSort_S(char**,uint64_t,uint64_t);