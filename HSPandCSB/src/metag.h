/* @author Fernando Moreno Jabato <jabato@uma.es>
 * @date 24-Sept-2015
 * @description This file is a header that includes all code 
 *     to handle metagenomes.
 * @license all rights reserved to BitLAB (http://www.bitlab-es.com/bitlab/)
 *     and to author. 
 */

#include "stdio.h"
#include "stdlib.h"
#include "string.h"

#define MAXREADLENGTH 100000


struct Read{
	char *ID;    // Read ID.
	char *seq;   // Sequence (read).
};

// Prototipes
int countReads(FILE*);