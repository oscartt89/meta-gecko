/*
 * @author Fernando Moreno Jabato <jabato@uma.es>
 * @licence all rights reserved to the author and BitLAB group (University
 *    of Malaga).
 */

#include <stdio.h> // IO functions
#include <stdlib.h> // Memory functions
#include <errno.h> // errors
#include <string.h> // String functions
#include <stdbool.h> // Boolean varaibles
#include <ctype.h>
#include <arpa/inet.h>
#include "metag_structs.h"

// VARAIBLES
#define MAX_BUFF 10000000
#define MAX_FILE_LENGTH 1024
#define READ_BUFF_LENGTH 10000
#define Eq_Value 4
#define Dif_Value -4
#define Score_Threshold 0


// LINE FUNCTIONS
#define SWAP(a,b,t) t=a; a=b; b=t;

// GLOBAL VARAIBLES
uint64_t buffersWritten;
uint64_t S_Threshold; // Similarity threshold
uint64_t L_Threshold; // Length threshold
uint64_t hitLength;
int prefixSize; // Word size


// FUNCTIONS
int wordcmp(unsigned char*,unsigned char*,int);
int HComparer(Hit w1, Hit w2);
void readHashEntry(WordEntry*,FILE*);
int readWordEntrance(WordEntry*,FILE*,uint16_t);
int generateHits(Hit*,WordEntry,WordEntry,FILE*,FILE*,FILE*,FILE*,uint64_t*,int);
void loadLocationEntrance(LocationEntry*,FILE*,uint32_t,bool); 
inline void storeHit(Hit*,LocationEntry,LocationEntry,uint64_t);
void writeHitsBuff(Hit*,FILE*,FILE*,uint64_t);
int GT(Hit,Hit);
int partition(Hit*,int,int);
void quicksort_H(Hit*,int,int);
uint64_t loadHit(Hit**,FILE*,int64_t);
uint64_t lowestHit(Hit*,uint64_t,int64_t*);
bool finished(int64_t*,uint64_t);
void writeFragment(FragFile,FILE*);
void SWAP_H(Hit*,Hit*,Hit);
void copyHit(Hit*,Hit);
void push(node_H**,node_H**);
void move(node_H**,node_H**);
void sortList(node_H**);
void checkOrder(node_H**,bool);
void FragFromHit(FragFile*,Hit*,Reads*,Sequence*,uint64_t,uint64_t,FILE*);
char getValue(Sequence*,uint64_t,int);
Sequence* LeeSeqDB(char*,uint64_t*,uint64_t*);
Reads* LoadMetagenome(char*,uint64_t*);
inline void freeReads(Reads**);
void writeSequenceLength(uint64_t*,FILE*);
void endianessConversion(char*,char*,int);