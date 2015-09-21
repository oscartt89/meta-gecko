
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <math.h>

#include "fragmentv3.h"
#include "fragmentv2.h"



void updateFragmentv3(Fragmentv3* f, Fragmentv3 g);
void resetFragmentv3(Fragmentv3* f);
void printFragv3(Fragmentv3 f);
int overlapFragFile ( struct FragFile f, struct FragFile g,char c,int sol);
int overlapFragmentv3 ( Fragmentv3 f, Fragmentv3 g,char c,int sol);


int MIN (long int a, long int b);
int MAX (long int a, long int b);
