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
#include <ctype.h> // Char functions: isupper,...
#include <stdbool.h> // Boolean varaibles

// VARIABLES
#define BUFFER_LENGTH 10000000
#define MAX_FILE_LENGTH 1024
#define FLAG 0
#define MERGE_BUFFER_LENGTH 1000
// FUNCTIONS
#define SWAP_W(a,b,t) t=a; a=b; b=t;

int BYTES_IN_WORD;
uint16_t WL;

// STRUCTS
typedef struct {
	// Word compressed in binary format
    unsigned char *seq;
    // Ocurrence position in the sequence (Read)
    uint64_t pos;
    // Read index
    uint32_t seqIndex;
} wentry;

typedef struct node{
    // Wentry
    wentry *buff;
    // Index where start to read
    uint32_t index;
    // Words in buffer
    uint64_t numWords;
    // Previous
    struct node *prev;
    // Next
    struct node *next; 
}node;

// FUNCTIONS
inline void shift_word(unsigned char*);
inline void storeWord(wentry*,wentry);
int writeBuffer(wentry*,FILE*,FILE*,uint64_t);
int wordcmp(unsigned char*,unsigned char*,int);
int wordComparator(wentry*,wentry*);
int partition(wentry*,int,int);
int quicksort_W(wentry*,int,int);
bool finished(int64_t*,uint64_t);
inline void loadWord(wentry*,FILE*);
void sortBuffers(wentry***,uint64_t,uint32_t*,uint64_t*);
inline void writeWord(wentry*,FILE*,FILE*,bool,uint32_t*);
void checkOrder(wentry***,uint32_t*,uint32_t*,uint64_t*);
void loadMatrix(wentry**,int64_t*,uint32_t,uint32_t*,uint64_t*,uint32_t*,uint64_t*,FILE*);
void showWord(unsigned char*,int);