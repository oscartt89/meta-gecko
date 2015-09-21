/* mgFoIview2TXT-R2G: produces a text files to be displayed with the 
 *                best matches for each read (genomnes over reads)
 *                complements the program mgFoIview2TXT who matches reads over genomes
 *
 *     syntax: ./mgFoIview2TXT-R2G LISTgenomes reads.File reads.IND frags.file  matchingFILE  out.TXT 
 * LISTgenomes : contains the list of genome files usewd in the matching process (in this version is hard-coded)
 * reads*.IND : mapped reads (and associted index) file (obtain the index with mgReadsIndex prog)
 * frags.file  : fragments from the comparisson process (comes from FragHits program)
 * matchingFILE: produced by bestGenomeInRead program
 * out.TXT     : name of the output file (cisualized with the Java program) 
 *
 * Note: ONLY the matching part of the genome represented 
 * ortrelles at uma.es
 * date: November the 1st, 2014
 * ------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <ctype.h>
#include <time.h>

#include "structs.h" //For the FragFile struct
#include "fragmentv2.h"
#include "comparisonFunctions.h"
#include "commonFunctions.h"

#define max(a,b)    (((a)>(b)) ? (a):(b))
#define min(x,y)    (((x) < (y)) ? (x) : (y))
#define MAXnFRG   1000
#define RuleSTEP   100
#define SEP " "
#define MAXCACHE    40

/*
#define NEGRO  0
#define AZUL   1
#define ROJO   2
#define NARANJ 3
#define VERDE  4
#define AMARI  5
#define MORAO  6
#define BLANC  7
#define MAXFRAME 800
#define MaxWIDTH 100  // numer of nucleotids (from +100)
*/

//Struct for reads index tuple
struct rIndex {
	    uint64_t pos;      //Start position of read
	    uint64_t Lac;      // accumulated length
};

/*struct genomeChunk{
    char ident[MAXLID+1];
    int  from;
    int  to;
    char *datos;
};*/

struct Files {
   char fname[1024];
   char shortName[40];
};

struct SequenceCache {
	uint64_t timesUsed;
	uint64_t lastUse;
	int64_t sequenceNumber;
	uint64_t GL;
	struct Sequence *seq;
};

int howManyFiles(char *fname);
struct Files* LoadListOfFiles(char*fListName, int *nF);
struct Sequence *leeRead(FILE *, int *, int, uint64_t *);
void printReadGuide(FILE *, struct Sequence *, int, uint64_t, int, int);

//Sequence cache functions
struct SequenceCache *search(struct SequenceCache *sc, int64_t sequenceNumber);
void insert(struct SequenceCache *sc, struct Sequence *s, uint64_t sequenceNumber, uint64_t GL);

struct SequenceCache *search(struct SequenceCache *sc, int64_t sequenceNumber){
	uint64_t i;
	int64_t pos = -1;

	for(i=0; i < MAXCACHE; i++){
		if(sc[i].sequenceNumber==sequenceNumber){
			sc[i].timesUsed++;
			sc[i].lastUse = 0;
			pos = i;
		} else {
			sc[i].lastUse++;
		}
	}

	return (pos != -1)? &sc[pos] : NULL;
}

uint64_t lessUsed(struct SequenceCache *sc){
	uint64_t i, pos = 0, minTimesUsed=UINTMAX_MAX;
	for(i=0; i<MAXCACHE; i++){
		if(sc[i].timesUsed<minTimesUsed){
			minTimesUsed = sc[i].timesUsed;
			pos = i;
		}
	}
	return pos;
}

uint64_t mostUsed(struct SequenceCache *sc){
	uint64_t i, pos = 0, maxTimesUsed=0;
	for(i=0; i<MAXCACHE; i++){
		if(sc[i].timesUsed>maxTimesUsed){
			maxTimesUsed = sc[i].timesUsed;
			pos = i;
		}
	}
	return pos;
}

uint64_t lru(struct SequenceCache *sc){
	uint64_t i, pos = 0, maxLastUse=0;
	for(i=0; i<MAXCACHE; i++){
		if(sc[i].lastUse>maxLastUse){
			maxLastUse = sc[i].lastUse;
			pos = i;
		}
	}
	return pos;
}

uint64_t mru(struct SequenceCache *sc){
	uint64_t i, pos = 0, minLastUse=UINTMAX_MAX;
	for(i=0; i<MAXCACHE; i++){
		if(sc[i].lastUse<minLastUse){
			minLastUse = sc[i].lastUse;
			pos = i;
		}
	}
	return pos;
}

uint64_t randomC(struct SequenceCache *sc){
	int64_t tmp = sc[0].lastUse;
	tmp++;
	return rand() % MAXCACHE;
}

//This function assumes the sequence is not present
void insert(struct SequenceCache *sc, struct Sequence *s, uint64_t sequenceNumber, uint64_t GL){
	uint64_t pos = 0;
	//Search for the place to be stored
	while(pos<MAXCACHE && sc[pos].sequenceNumber>=0){
		pos++;
	}

	//if there is no place, use STRATEGY to find the entry to be replaced
	if(pos>=MAXCACHE){
		pos = STRATEGY(sc);
	}

	sc[pos].timesUsed=0;
	sc[pos].lastUse = 0;
	free(sc[pos].seq);
	//We have the position
	sc[pos].sequenceNumber=sequenceNumber;
	sc[pos].seq=s;
	sc[pos].GL=GL;
}

/* / get a genome chunk from a genome seq
struct genomeChunk *loadGenomeChunk(struct Sequence *G, int From, int To) {
        char c;
        int i,len=0,j, k=0;
        struct genomeChunk *Gch;
        int maxLen;

        if (( Gch= (struct genomeChunk*)calloc(1,sizeof(struct genomeChunk)) )==NULL)
           terror("memory for sctruct genomeChunk (1)");
        Gch->from = From;
        Gch->to   = To;
        maxLen  = To - From + 1;
        if ( (Gch->datos=(char *)calloc(1,(maxLen + 1)*sizeof(char)))==NULL)
           terror("memory for datos in sctruct genomeChunk (2)");

        strcpy(Gch->ident, G->ident);

	for (j=0, i=From; i<=To; i++, j++)
	   Gch->datos[j]= G->datos[i];
	Gch->datos[maxLen]=0x00;

        return Gch;
}
*/

int howManyFiles(char *fname){
        FILE *f;
        char line[1024];
        int nL=0;

        if ((f=fopen(fname,"rt"))==NULL){
            fprintf(stdout,"Opening frags file LIST");
	    exit(-1);
	}

        fgets(line,1024,f);
        while(!feof(f)) {
                if (line[0]!='#' && (int)strlen(line)!=0) nL++;
                fgets(line,1024,f);
        }
        fclose(f);
        return nL;
}

// Load to memory a list of files --------------------------------------------------
// datafile format: fileName[tab]nSeq[tab][format][newLINE]

struct Files* LoadListOfFiles(char*fListName, int *nF) {

        FILE *f;
        struct Files*L=NULL;
        char line[1024];
        int N=0,nFiles;

        nFiles = howManyFiles(fListName);
        if ((L=(struct Files*) calloc(nFiles,sizeof(struct Files)))==NULL)
           terror("memory for list of files");

        if ((f=fopen(fListName,"rt"))==NULL)
                terror("opening List of Files");

        fgets(line,1024,f);
        while(!feof(f)) {
                if (line[(int)strlen(line)-1]=='\n') line[(int)strlen(line)-1]=0x00;

                if (line[0]!='#' && (int)strlen(line)>2) {
			sscanf(line, "%s %s",L[N].fname,L[N].shortName);  
                }
		N++;
                fgets(line,1024,f);
        }
        fclose(f);
        (*nF) = nFiles;
        return L;
}


struct Sequence *leeRead(FILE *f, int *rL, int from, uint64_t *nChunks) {
	char c = '\0';
	uint64_t lon = 0, k = 0;

	fseek(f, from, SEEK_SET);

	struct Sequence *sX; //sX will be the first elem

	//Initialize
	*nChunks = 0;

	//Memory
	if ((sX = (struct Sequence*) malloc(sizeof(struct Sequence))) == NULL)
		terror("Memory...");

	while ((c = getc(f)) != '>' && !feof(f))
		; //start seq

	c = getc(f);
	while (k < MAXLID && c != '\n' && c != ' ') {
		if (feof(f))
			return 0;

		sX->ident[k++] = c;
		c = getc(f);
	}

	sX->ident[k] = 0; //end of data.

	while (c != '\n')
		c = getc(f);
	c = getc(f);

	//start list with sX2
	while (c != '>' && !feof(f)) {
		c = toupper(c);
		if (isupper(c))
			sX->datos[lon++] = c;
		c = getc(f);
	}

	if (lon < MAXLS)
		sX->datos[lon] = 0x00;

	//fprintf(stdout, "sx-datos: %s \n", sX->datos);
	*nChunks = 1;
	*rL = (int) lon;
	
    return sX;

}

// draw Read and Rule header (for each read & Read info)
void printReadGuide(FILE *fOut, struct Sequence *R, int rL, uint64_t nChunks, int rN, int nFG){
    int i, j, k=RuleSTEP;
    j= rL /RuleSTEP;
    fprintf(fOut,">Read/Genome [rNumber=%-8d](L:%-9d nF=%-3d): ",rN,rL,nFG);
    for (i=0; i< rL; i+=RuleSTEP){
	fprintf(fOut,"....-....1....-....2....-....3....-....4....-....5....-....6....-....7....-....8....-....9....-.");
	fprintf(fOut,"%4d",k);
        k+=RuleSTEP;
    }
    fprintf(fOut,"\n");
    fprintf(fOut,"%50s : ", R->ident);
    for(i=0;i<rL;i++){
        fprintf(fOut, "%c", getValue(R, i, nChunks));
    }
    fprintf(fOut,"\n");
}

int main(int ac, char **av) {
    int nF;
    FILE *fR, *fRI,*fFoI, *fG, *fOut;
    //FILE *fF;
    uint64_t nx,ny;
    //long t;
    struct FragFile *F, *f;

    int nFp;

    int i,j;
    int xI,xE, yI,yE,xL;
    int rL, tLen, nReads;
    struct rIndex *rI;
    int r, rpos, rOffset; //FrOffset, LastFrag = -1;
    struct Sequence *R=NULL, *G=NULL;
    uint64_t nChunksR, nChunksG, tmp;
    int yIni, yEnd,numR;
    int nG=0, GL=0;
    //int sizeF= sizeofFragment();
    // to manage the new mapping file
    int readNant=-1;
    int readN,genomeN, maxX,score, totIDENT,totLEN, nFG;
    float Totcov, porID;
    int frags[MAXnFRG];
    struct SequenceCache sc[MAXCACHE];
    struct SequenceCache *ce;
  //======================================================


    if (ac!=7) 
  	  terror("syntax: mgFoIview2TXT-R2G LISTgenomes reads.File reads.IND fragsFILE matchingFILE  out.TXT");

    //Initialize random seed
    srand(time(NULL));

    // open OUT  file
    if ((fOut=fopen(av[6],"wt"))==NULL)
	    terror("Opening OUTPUT file");

    // Load list of Genomes
    struct Files *GTable;
    GTable = LoadListOfFiles(av[1], &nG);

    for (i=0; i<nG; i++)
	fprintf(fOut, "Genome [%-4d] [%10s] %s\n",i,GTable[i].shortName,GTable[i].fname);
    fprintf(fOut,"Number of Genomes=%d\n",nG);


    // load reads index (only those containing reads in the FoI)
    // if the number of reads is too large it would be possible
    // to manage the index in DISK
    // open reads Index

    //fprintf(stdout, "Reads index Start\n");
    if ((fRI=fopen(av[3],"rb"))==NULL)
	terror("Opening rweads Index file");

    fseek(fRI,0,SEEK_END);
    tLen=ftell(fRI);
    nReads = tLen/sizeof(struct rIndex);
    if ((rI = (struct rIndex*)calloc(nReads, sizeof(struct rIndex)))==NULL)
	terror("memory for READS index");
    fseek(fRI,0,SEEK_SET);
    fread(rI, nReads, sizeof(struct rIndex), fRI);
    fclose(fRI);
    //fprintf(stdout, "Reads index End\n");

    // open reads file------------------------
    if ((fR=fopen(av[2],"rt"))==NULL)
	    terror("Opening reads file");

    // open fragments file-----
    //if ((fF=fopen(av[4],"rt"))==NULL)
    //    terror("Opening fragments file");

    //store fragments in memory
    //fprintf(stdout, "Fragments read start\n");
    f = readFragments(av[4], &nF, &nx, &ny); 

    //fprintf(stdout, "Fragments read End\n");

    // number of frags (Check OSCAR .... if valid
    //fseek(fF, 0L, SEEK_END);
    //t= ftell(fF);
    //nF = (t - 2*sizeof(uint64_t))/sizeF;
//    fprintf(stdout, "Number of frags: %d\n",nF);
    //fseek(fF, 0L, SEEK_SET);
    //readSequenceLength(&nx, fF);
    //readSequenceLength(&ny, fF);
//{int xnx=nx, xny=ny; printf("nX=%d ny=%d sizeF=%d===========================================\n",xnx, xny,sizeF);} 

    // open MATCHING file
    if ((fFoI=fopen(av[5],"rt"))==NULL)
        terror("Opening MATCHING Frags file");

    for(i=0;i<MAXCACHE;i++){
	sc[i].lastUse = 0;
        sc[i].timesUsed=0;
	sc[i].sequenceNumber=-1;
	sc[i].GL = 0;
	sc[i].seq = NULL;
    }
	    
    // kick-off================================================

    fscanf(fFoI, "%d\t%d\t%f\t%d\t%d\t%d\t%d\t%f\t%d\t",&readN,&genomeN,&Totcov, &maxX,&score, &totIDENT,&totLEN,&porID, &nFG);

    while(!feof(fFoI)){

      for (i=0;i<nFG; i++)  // read the frag numbers
         fscanf(fFoI,"\t%d",&frags[i]);

//printf("r=%d G=%d nFr=%d----------------------\n",readN, genomeN, nFG);

      if (readNant!= readN) { // new read
	//if (readNant!=-1) ;
	fprintf(fOut, "\n");

        // Load the new read 

        r      = (F->strand == 'f')? F->seqY : (nReads-1)-F->seqY;
        rpos   = rI[readN].pos;
        rOffset= rI[readN].Lac; // not used ...now relative positions inside the read
	if(R!=NULL){
		free(R);
		R=NULL;
	}
	R     = leeRead(fR, &rL, rpos, &nChunksR);
	readNant=readN;

	printReadGuide(fOut, R, rL, nChunksR, readN, nFG);
      }

      
      if ((ce = search(sc, (int64_t)genomeN)) == NULL) {
	//if (genAnt!= -1) ;//fprintf(fOut,"\n");

	  //load genome sequence-------
	    if ((fG=fopen(GTable[genomeN].fname,"rt"))==NULL){
	        fprintf(stderr, "Genome file: %s\n", GTable[genomeN].fname);
		terror("Opening genome file");
	}

	G = LeeSeqDB(fG, &tmp, &nChunksG, 0);
	insert(sc, G, genomeN, tmp);
	if(G==NULL)terror("Error loading genome");
	ce = search(sc, (int64_t)genomeN);
	GL=(int)tmp;
//printf("G:%s:",GTable[genomeN].shortName);for(j=0;j<30;j++) printf("%c",G.datos[j]); printf(" (L=%d)\n",GL);

	    fclose(fG);
	  // end load Genome (can/MUST be optimised)
	  // to optimise, store the genomes (ony data) as a cotinuos string
	  // to allow direct access from Disk

      }

      G = ce->seq;
      GL = ce->GL;

      // fragments LOOP---------------------
      for (i=0;i<nFG; i++) { // read the frag numbers
	
	// read the fragment-------- 
        // NOT SURE IF THIS WILL WORK--- otherwise we need a FRAGs-INDEX

	//FrOffset = sizeF * (frags[i] - LastFrag - 1);
	//fseek(fF, FrOffset, SEEK_CUR);
	//readFragment(&F, fF);
        //
	//LastFrag = frags[i];

        F = &f[frags[i]];

        yIni=F->yStart; yEnd=F->yEnd; numR=F->seqY;
        xI=F->xStart; xE=F->xEnd; yI=F->yStart; yE=F->yEnd; xL=F->length;
	// chek tge coordinates X ir Y+

/*        colIni = (int)F->xStart; // - From;
        colEnd = (int)F->xEnd;
        Lprint = (int)F->length;
        rIni   = (F->strand=='f')?(int)F->yStart:(int)(ny-F->yStart-1); // start point of the genome
*/

//printf("Frg(%d):xIni=%d yIni=%d yEnd=%d xL=%d %c  ..rIni=%d\n",frags[i],xI,yI, yE, xL, F->strand, rIni);
/* working over the fragment coordinates---relative, therefore, not needed
        if (colIni < 0) colIni=0;
        if (colEnd > To) colEnd=To;
        colEnd = (colEnd - From);
        Lprint = colEnd - colIni +1;
*/

//        if((colEnd-colIni)<=1) colEnd+=2;  // for SHORT frags

        fprintf(fOut,"r=%4d g=%7d L=%4d [%c] %%I=%3.2f %14s :",xI,yI, xL, F->strand, ((float)F->ident)/F->length, GTable[genomeN].shortName); 

//	fprintf(stdout, "g=%s rP=%4d gP=%7d gE=%d L=%4d [%c]",GTable[genomeN].shortName, xI,yI, yE, xL, F->strand); fflush(stdout);

	fprintf(fOut, " ");
        for (j=0; j< xI /*colIni*/; j++) fprintf(fOut, SEP);


        // HERE PRINT GENOMES chunks----
	// G   contains the genome sequence; GL the length
	// Gch contain the genome chunk specified in the Fragment coordinates (DELETED

// Load genomeChunk into memory (DELETED-------------working over the genome sequence)
//        Gch = loadGenomeChunk(&G, From, To);
//printf("Gch:%s:[d=%d  h=%d] ",Gch->ident, From, To);for(j=0;j<30;j++) printf("%c",Gch->datos[j]); printf("\n");

        if(F->strand == 'f'){
			if (yI+xL >GL){  fprintf(fOut, "GLen=%d <----bad coordinates----->\n",GL); continue;}
				else for (j=0; j<xL; j++) fprintf(fOut, "%c", getValue(G, yI+j, nChunksG)); // [rIni-rOffset+j]);
	} else {
			char tmp;
			for (j=0; j<xL; j++){
				if(((GL-1)-yI-xL)<0){
					fprintf(fOut, "GLen=%d Position=%d <----bad coordinates----->\n",GL,((GL-1)-yI-xL)); continue;
				} else {
					tmp = getValue(G,(GL-1)-yI-j,nChunksG);
					switch(tmp){
						case 'A': tmp = 'T';
						break;
						case 'C': tmp = 'G';
						break;
						case 'G': tmp = 'C';
						break;
						case 'T': tmp = 'A';
						break;
					}
					fprintf(fOut, "%c",tmp);
				}
			}	

	}

        fprintf(fOut,"\n");
     }

/*
        if (Gch->datos) free(Gch->datos);
        if (Gch) free(Gch);
*/
	// next matching read-genome
	fscanf(fFoI, "%d\t%d\t%f\t%d\t%d\t%d\t%d\t%f\t%d\t",&readN,&genomeN,&Totcov, &maxX,&score, &totIDENT,&totLEN,&porID, &nFG);
	nFp++;
    }
//    fprintf(stdout,"-----------done\n");
    fclose(fFoI);
    fclose(fR);
    //fclose(fF);
    fclose(fOut);
    return 0;
}
