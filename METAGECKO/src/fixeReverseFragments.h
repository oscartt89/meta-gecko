/*
 * @author Fernando Moreno Jabato <jabato@uma.es>
 * @licence all rights reserved to the author and BitLAB group (University
 *    of Malaga).
 */

#include <stdio.h> // IO functions
#include <stdlib.h> // Memory functions
#include <string.h> // String functions
#include <ctype.h>
#include <arpa/inet.h>
#include <stdbool.h> // Boolean varaibles
#include "metag_common.h"

#define FRAG_BUFFER_SIZE 100000
#define MAX_FILE_LENGTH 1024 // Already exists
#define READ_BUFF_LENGTH 10000 // Already exists


// Functions
void readFragment(FragFile*,FILE*);
void writeFragment(FragFile,FILE*);
int fragmentComparator(FragFile,FragFile);
void readSequenceLength(uint64_t*,FILE*);
void writeSequenceLength(uint64_t*,FILE*);
void endianessConversion(char*,char*,int);
uint64_t loadStats(Read**,char*);
void changeOrder(Read**);
inline void seekStats(uint64_t,Read*,Read*,Read**,Read**);
void writeFragBuffer(FragFile*,FILE*,FILE*,uint64_t);
int partition(FragFile*,int,int);
void quicksort_F(FragFile*,int,int);
void SWAP_F(FragFile*,FragFile*,FragFile);
void copyFFile(FragFile*,FragFile);
void freeReadList(Read**);
void checkOrder(fnode**,bool);
void push(fnode**,fnode**);
void move(fnode**,fnode**);
void sortList(fnode**);
uint64_t loadFrag(FragFile**,FILE*,int64_t);