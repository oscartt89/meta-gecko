// LIBRARIES
#include "geckoMGVutils.h"
#include <ctype.h> // Char functions: isupper,...

// VARIABLES
static const int  WORD_SIZE = 32;

// FUNCTIONS
int createDictionary(char*,char*);
//
void shift_word(word*);
int storeWord(wentry*,wentry*,int);
void writeDic(wentry*,int,FILE*,FILE*,FILE*); 