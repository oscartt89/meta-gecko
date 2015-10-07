/* @author Fernando Moreno Jabato <jabato@uma.es>
 * @date 30-Sept-2015
 * @description This file is a header that includes all code 
 *     to handle and create dictionaries.
 * @license all rights reserved to BitLAB (http://www.bitlab-es.com/bitlab/)
 *     and to author. 
 */

#include "stdio.h"
#include "stdlib.h"
#include "string.h"

struct Word{
	char *sequence;
	int pos;
};

struct WordList{
	Word *nextWord;
	Word *word;
};

// Dictionary prototipes
void createDictionary(char*,int);
void initWord(struct Word*, int wordLength);