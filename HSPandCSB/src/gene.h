


#ifndef GENE_H_
#define GENE_H_
/*****/
// Gene Struct

typedef struct {
	char id[50];
	char locus[50];
	char gene[50];
	char strand;
	int start;
	int end;
	// NewField
	int length;	
	int* cov;
	char product[1000];
}GeneFeatures;
/****/

GeneFeatures* readGenes(char* s,int* nf);
int getLocusLong(GeneFeatures* g, int ini, int end,int max,char c);
int sameStrand(char g,char s);
void setCOV(GeneFeatures* G, int ini,int end,int n);
#endif /* GENE_H_ */
