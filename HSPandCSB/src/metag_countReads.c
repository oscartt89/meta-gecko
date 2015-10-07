/* @author Fernando Moreno Jabato <jabato@uma.es>
 * @date 28-Sept-2015
 * @description This file contains necessary implementation for return the
 *     number of reads on a metagenome file given.
 * @license all rights reserved to BitLAB (http://www.bitlab-es.com/bitlab/)
 *     and to author. 
 */
#include "metag.h"

int main(int ac, char **av){
	if(ac!=2){
		fprintf(stderr, "Bad call error. Use: %s metagenomFile\n", av[0]);
		return -1;
	}

	FILE *meta;

	// Open input stream
	if((meta = fopen(av[1],"r"))==NULL) {
		fprintf(stderr,"Error opening metagenome file. [%s]\n",av[1]);
		return -1;
	}

	return countReads(meta);
}