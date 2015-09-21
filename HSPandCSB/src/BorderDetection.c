/* BorderDetection  

Given tow CSB, it searches end and start areas and look for hits in these areas.
The output are tow histogram of hits in these areas. 

 
 ------------------arjona@uma.es----Oct/2014-------------*/

#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include "commonFunctions.h"
#include "structs.h"
#include "commonFunctions.h"
#include "fragmentv3.h"
#include "JAMfunctions.h"

#define POINT 4



int main(int ac, char **av) {

	if (ac != 7){
		printf("ac:%d\n",ac);
		terror("usage: BorderDetection Fragments.csb.frag csbA csbF csbG csbB output.csv");
		
	}
	// Get parameters
	int bA,bF,bG,bB;
	bA = atoi(av[2]);
	bF = atoi(av[3]);
	bG = atoi(av[4]);
	bB = atoi(av[5]);
	
	
	/* Read Fragments */
	Fragmentv3* f;
	int nf;
	int xtotal,ytotal;
		
	f=readFragmentsv3(av[1],&nf,&xtotal,&ytotal);
	/******************/
	int i=0;
	
	
	//printf("%d\n",(int)pow(2,8*sizeof(int)));exit(-1);
	/* Debug
	printf("f.xIni\tf.yIni\tf.xFin\tf.yFin\tf.length\tf.strand\tf.score\tf.block\n");
	for(i=0;i<nf;i++){
		printf("%d\t%d\t%d\t%d\t%d\t%c\t%d\t%d\n",(int)f[i].xIni,(int)f[i].yIni,(int)f[i].xFin,(int)f[i].yFin,(int)f[i].length,f[i].strand,(int)f[i].score,(int)f[i].block);
	}
	*/
	/*******************/
	
	/* Calculating A,F,G,B as two fragment */
	Fragmentv3 A,F,G,B;
	resetFragmentv3(&A);
	resetFragmentv3(&F);
	resetFragmentv3(&G);
	resetFragmentv3(&B);
	for(i=0;i<nf;i++){
		if ((int)f[i].block==bA){
			updateFragmentv3(&A,f[i]);
		}else if ((int)f[i].block==bF){
			updateFragmentv3(&F,f[i]);
		}else if ((int)f[i].block==bG){
			updateFragmentv3(&G,f[i]);
		}else if ((int)f[i].block==bB){
			updateFragmentv3(&B,f[i]);
		}
	}
	
	// Debug
	
	printFragv3(A);
	printFragv3(F);
	printFragv3(G);
	printFragv3(B);
	
	/****************************/
	
	
	

	return 0;
}
