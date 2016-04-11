/*
 * @author Fernando Moreno Jabato <jabato@uma.es>
 * @licence all rights reserved to the author and BitLAB group (University
 *    of Malaga).
 */

#include <stdio.h> // IO functions
#include <stdlib.h> // Memory functions
#include <errno.h> // errors

#include <string.h> // String functions
#include <ctype.h> // Char functions: isupper,...
#include <stdbool.h> // Boolean varaibles

#include "metag_structs.h"

// VARIABLES
#define BUFFER_LENGTH 10000000
#define MAX_FILE_LENGTH 1024
#define FLAG 0
#define READ_BUFF_LENGTH 10000
// FUNCTIONS
#define SWAP_W(a,b,t) t=a; a=b; b=t;

int BYTES_IN_WORD;

 


// FUNCTIONS
inline void shift_word(Word*);
inline void storeWord(wentry*,wentry);
int writeBuffer(wentry*,FILE*,FILE*,uint64_t);
int wordcmp(Word,Word,int);
int wordComparator(wentry*,wentry*);
int partition(wentry*,int,int);
int quicksort_W(wentry*,int,int);
uint64_t loadWord(wentry**,FILE*,int64_t);
inline void writeWord(wentry*,FILE*,FILE*,bool,uint32_t*);
void checkOrder(node_W**,bool);
int GT(wentry,wentry);
void push(node_W**,node_W**);
void move(node_W**,node_W**);
void sortList(node_W**);
