/*
 * @author Fernando Moreno Jabato <jabato@uma.es>
 * @licence all rights reserved to the author and BitLAB group (University
 *    of Malaga).
 */

#include "metag_common.h"

// VARAIBLES
#define MAX_BUFF 10000000
#define MAX_FILE_LENGTH 1024
#define READ_BUFF_LENGTH 10000
#define Eq_Value 4
#define Dif_Value -4
#define Score_Threshold 0


// LINE FUNCTIONS
#define SWAP(a, b, t) t=a; a=b; b=t;

// FUNCTIONS

int HComparer(Hit w1, Hit w2);


uint64_t readHashEntrance(HashEntry *we, FILE *wD, uint16_t SeqBytes, char * byteBufferHits, int nextMWorGW, uint64_t * currHit, uint64_t loadFrom, uint64_t * tRead);

int generateHits(Hit *, HashEntry, HashEntry, FILE *, FILE *, FILE *, FILE *, uint64_t *, int, uint64_t *, uint64_t, uint64_t);

inline void loadLocationEntrance(LocationEntry *, FILE *, uint32_t);

inline void storeHit(Hit *, LocationEntry, LocationEntry, uint64_t, uint64_t);

void writeHitsBuff(Hit *, FILE *, FILE *, uint64_t, int, uint64_t *);

int GT(Hit, Hit);

int partition(Hit *, int, int);

void quicksort_H(Hit *, int, int);

inline uint64_t loadHit(Hit *, FILE *, int64_t);

uint64_t lowestHit(Hit *, uint64_t, int64_t *);

bool finished(int64_t *, uint64_t);

void writeFragment(FragFile, FILE *);

void SWAP_H(Hit *, Hit *, Hit *);

void push(node_H **, node_H **);

void move(node_H **, node_H **);

void sortList(node_H **);

void checkOrder(node_H **, bool);

int FragFromHit(FragFile *, Hit *, Reads *, Sequence *, uint64_t, uint64_t, FILE *, unsigned int, unsigned int, float);

inline char getValueOnRead(Reads *r, uint64_t pos);

char getValue(Sequence *, uint64_t, int);

Sequence *LeeSeqDB(char *, uint64_t *, uint64_t *);

Reads *LoadMetagenome(char *, uint64_t *);

inline void freeReads(Reads *);

void freeGenomes(Sequence * genomes);

void writeSequenceLength(uint64_t *, FILE *);

void endianessConversion(char *, char *, int);

inline int filteredHit(Hit h1, Hit h2, int prefixSize);

inline uint32_t getSizeOfIndexEntry(uint16_t WB);