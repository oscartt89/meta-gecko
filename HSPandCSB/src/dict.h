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
#define READ_BUFF_LENGTH 10000
// FUNCTIONS
#define SWAP_W(a,b,t) t=a; a=b; b=t;

int BYTES_IN_WORD;

// STRUCTS
typedef struct {
    uint16_t WL;
	// Each letter is stored using 2 bits
	// We have 4 letters per byte and 
	// each unsigned char size is 1 byte.
    unsigned char *b;
} word;

typedef struct {
	// Word compressed in binary format
    word w;
    // Ocurrence position in the sequence (Read)
    uint64_t pos;
    // Read index
    uint32_t seq;
} wentry;

struct item{
    wentry *word;
    uint64_t buff;
    uint64_t index;
    uint64_t words_loaded;
    struct item *next;
};

typedef struct item node; 


// FUNCTIONS
inline void shift_word(word*);
inline void storeWord(wentry*,wentry);
int writeBuffer(wentry*,FILE*,FILE*,uint64_t);
int wordcmp(word,word,int);
int wordComparator(wentry*,wentry*);
int partition(wentry*,int,int);
int quicksort_W(wentry*,int,int);
long int searchInIndex(word,FILE*,long int,bool*);
bool finished(int64_t*,uint64_t);
uint64_t loadWord(wentry**,FILE*,int64_t);
inline void writeWord(wentry*,FILE*,FILE*,bool,uint32_t*);
void checkOrder(node**,bool);
int GT(wentry,wentry);
void push(node**,node**);
void move(node**,node**);
void sortList(node**);
void showWord(word*,int);
