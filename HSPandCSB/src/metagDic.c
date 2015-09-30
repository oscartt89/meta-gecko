/* @author Fernando Moreno Jabato <jabato@uma.es>
 * @date 30-ept-2015
 * @description This file contains necessary implementation for create metagenome
 *     dictionaries from a file given.
 * @license all rights reserved to BitLAB (http://www.bitlab-es.com/bitlab/)
 *     and to author. 
 */
#include "metag.h"

int main(int ac, char **av){
	// Variables
	FILE *metagenome, *dictionary;
	int wordLength;

	// Check bad call error
	if(ac != 3){
		fprintf(stderr, "Error: bad call error. Use: %s metagenome wordLength\n", av[0]);
		return -1;
	}

	
}