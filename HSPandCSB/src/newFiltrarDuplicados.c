/*
 * filtrarDuplicados.c
 *
 *  Created on: 17/03/2014
 *      Author: jarjonamedina
 * 
 
 We identify fragments that are overlapped in X or Y axis. 
 If two fragment are overlapped we compare its score value. The lesser we mark as a duplicated. 

 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <math.h>

#include "fragmentv3.h"



  #define MAX_LEVELS 1000000


int solape ( Fragmentv3 f, Fragmentv3 g,char c); // Overlapping function
void quickSortX(Fragmentv3 *f, int elements);
void quickSortY(Fragmentv3 *f, int elements,int ytotal);

int max(int a, int b);
int min(int a, int b);

int main(int ac,char** av){


if(ac<3){
	printf("Use: filtrarDuplicados fragsv3 fragsFiltr.fil frags.dup m[op] \n");
	exit(-1);
}





/* Open file */
FILE* fout;
	if((fout=fopen(av[2],"wb"))==NULL){
		printf("***ERROR Opening input file %s\n",av[2]);
		exit(-1);
	}
/* Open file */
FILE* foutDUP;
	if((foutDUP=fopen(av[3],"wb"))==NULL){
		printf("***ERROR Opening input file %s\n",av[3]);
		exit(-1);
	}



/* Read Fragments */
	Fragmentv3* f;
	int nf,xtotal,ytotal;
	f=readFragmentsv3(av[1],&nf,&xtotal,&ytotal);

int i,j;


// Reset frags
for(i=0;i<nf;i++){
	f[i].gap=0;
	f[i].events=0;
	f[i].duplications=0;
}
/*  Empezamos */
for(i=0;i<nf;i++){
	if(f[i].block>=0){
		if(!f[i].gap){
			for(j=i+1;j<nf;j++){
				if(!f[j].gap && !f[i].gap && f[j].block>=0){
					if( solape(f[i],f[j],'x') ){
						if(f[i].score>f[j].score){
							f[j].gap=1;
						}else{
							f[i].gap=1;
						}
					}
				}
			
			}
		
		}
	}
}

/***********************************/

//**********************/

/* Write fragments */ ////score>1
fwrite(&xtotal,sizeof(int),1,fout);
fwrite(&ytotal,sizeof(int),1,fout);

fwrite(&xtotal,sizeof(int),1,foutDUP);
fwrite(&ytotal,sizeof(int),1,foutDUP);

int ndup;
int nfrags;
ndup=nfrags=0;

//

for(i=0;i<nf;i++){

	
	if(f[i].gap==0){
		f[i].duplications=0;
	}else{
		f[i].block=-2;
	}
	
	fwrite(&f[i],sizeof(Fragmentv3),1,fout);
	printf("%d\t%ld\t%ld\t%ld\t%ld\t%c\t%ld\t%ld\t%ld\t%d\t%Le\t%d\t%d\n",f[i].block,f[i].xIni,f[i].yIni,f[i].xFin,f[i].yFin,f[i].strand,f[i].score,f[i].ident,f[i].length,f[i].gap,f[i].pvalue,f[i].duplications,f[i].block);

}


/*********************** END OF FILTER **************/
fclose(fout);
fclose(foutDUP);
return 0;
}

/*******************************************
*****************************************/


int solape ( Fragmentv3 f, Fragmentv3 g,char c){ 

// it returns  1 if overlapped

int sol;
sol=(int)c;

int fyIni,fyFin,gyIni,gyFin;
int fxIni,fxFin,gxIni,gxFin;
sol=0;

// Set data
fxIni=f.xIni;
fxFin=f.xFin;
gxIni=g.xIni;
gxFin=g.xFin;

fyIni=f.yIni;
fyFin=f.yFin;
gyIni=g.yIni;
gyFin=g.yFin;

if(f.strand=='r'){
	fyIni=f.yFin;
	fyFin=f.yIni;
}
if(g.strand=='r'){
	gyIni=g.yIni;
	gyFin=g.yFin;
}
// Program

if ( (gxIni<=fxIni) && (gxFin>=fxIni) )return 1;
if ( (gxFin>=fxFin) && (gxIni<fxFin) )return 1;
if ( (gxIni>=fxIni) && (gxFin<=fxFin) )return 1;

if ( (gyIni<=fyIni) && (gyFin>=fyIni) )return 1;
if ( (gyFin>=fyFin) && (gyIni<fyFin) )return 1;
if ( (gyIni>=fyIni) && (gyFin<=fyFin) )return 1;


//
	return sol;
	
}
//    quickSort(&myArray[3],5); // sorts elements 3, 4, 5, 6, and 7


void quickSortX(Fragmentv3 *f, int elements) {



  int  beg[MAX_LEVELS], end[MAX_LEVELS], i=0, L, R, swap ;
	Fragmentv3 piv;
  beg[0]=0; end[0]=elements;
  while (i>=0) {
    L=beg[i]; R=end[i]-1;
    if (L<R) {
      piv=f[L];
      while (L<R) {
        while (f[R].xIni>=piv.xIni && L<R) R--; if (L<R) f[L++]=f[R];
        while (f[L].xIni<=piv.xIni && L<R) L++; if (L<R) f[R--]=f[L]; }
      f[L]=piv; beg[i+1]=L+1; end[i+1]=end[i]; end[i++]=L;
      if (end[i]-beg[i]>end[i-1]-beg[i-1]) {
        swap=beg[i]; beg[i]=beg[i-1]; beg[i-1]=swap;
        swap=end[i]; end[i]=end[i-1]; end[i-1]=swap; }}
    else {
      i--; }}}



/*****************/
/*****************/

void quickSortY(Fragmentv3 *f, int elements,int ytotal) {

  int  beg[MAX_LEVELS], end[MAX_LEVELS], i=0, L, R, swap ;
	Fragmentv3 piv;
	long ytemp;
	ytemp=ytotal;
// Change Y component

	for(i=0;i<elements;i++){
		if(f[i].strand=='r'){
			// we have already done this in fragv2tov3.c
			//f[i].yIni=ytotal-f[i].yIni;
			//f[i].yFin=ytotal-f[i].yFin;

			ytemp=f[i].yIni;
			f[i].yIni=f[i].yFin;
			f[i].yFin=ytemp;	
		}
	}
//************/
  beg[0]=0; end[0]=elements;
  while (i>=0) {
    L=beg[i]; R=end[i]-1;
    if (L<R) {
      piv=f[L];
      while (L<R) {
        while (f[R].yIni>=piv.yIni && L<R) R--; if (L<R) f[L++]=f[R];
        while (f[L].yIni<=piv.yIni && L<R) L++; if (L<R) f[R--]=f[L]; }
      f[L]=piv; beg[i+1]=L+1; end[i+1]=end[i]; end[i++]=L;
      if (end[i]-beg[i]>end[i-1]-beg[i-1]) {
        swap=beg[i]; beg[i]=beg[i-1]; beg[i-1]=swap;
        swap=end[i]; end[i]=end[i-1]; end[i-1]=swap; }}
    else {
      i--; }}

// Change Y component

	for(i=0;i<elements;i++){
		if(f[i].strand=='r'){
			ytemp=f[i].yIni;
			f[i].yIni=f[i].yFin;
			f[i].yFin=ytemp;	

			//f[i].yIni=ytotal-f[i].yIni;
			//f[i].yFin=ytotal-f[i].yFin;


		}
	}
//*************/

}
/****************************************/
/*************************************/



int max(int a, int b){if (a>=b)return a;if(b>a)return b;return 0;}
int min(int a, int b){if (a>=b)return b;if(b>a)return a;return 0;}

