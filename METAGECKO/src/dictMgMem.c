#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"
#include "metag_common.h"

static uint32_t ksizeTables[3] = {2,4,8};


int main(int argc, char ** av){
    
    if(argc != 3) terror("USE: dictMgMem <metagenome> <ksize[1=8,2=16,3=32]>");
    
    
    
    //Database to read kmers from
    FILE * database;
    
        
    //Open database
    database = fopen64(av[1], "rt");
    if(database == NULL) terror("Could not open database");
    
    //Ksize for seeds
    uint32_t ksizeidx = (uint32_t) atoi(av[2]);
    
    
    //Get size of db
    fseeko64(database, 0L, SEEK_END);
    uint64_t totalSize = ftello64(database);
    fseeko64(database, 0L, SEEK_SET);
    
    //Allocate memory, as much as the size of the file for a table that will be sorted
    unsigned char * tableForward = (unsigned char *) malloc(ksizeTables[ksizeidx]*totalSize);
    unsigned char * tableReverse = (unsigned char *) malloc(ksizeTables[ksizeidx]*totalSize);
    if(tableForward == NULL) terror("Could not allocate memory for forward table");
    if(tableReverse == NULL) terror("Could not allocate memory for reverse table");
    
    
    //Variables to read kmers
    char c = 'N'; //Char to read characters

    unsigned char b[8], br[8];
    memset(b, 0, ksizeTables[ksizeidx]);
    memset(br, 0, ksizeTables[ksizeidx]);
    
    //Variables to account for positions
    int strandF = 1, strandR = 1;
    uint64_t pos = 0, crrSeqL = 0; //total length, current sequence length
    uint64_t idx = 0; //index for table access
    
    //Print info
    fprintf(stdout, "[INFO] Computing tree of mers\n");
    
    c = fgetc(database);
    while(!feof(database)){
        // Check if it's a special line
        if (!isupper(toupper(c))) { // Comment, empty or quality (+) line
            if (c == '>') { // Comment line
                c = fgetc(database);
                while (c != '\n') c = fgetc(database); //Avoid comment line

                crrSeqL = 0; // Reset buffered sequence length

                pos++; // Absolute coordinates: add one for the "*"
            }
            c = fgetc(database); // First char of next sequence
            continue;
        }
        
        if (strandF) shift_byte_array_left(b, ksizeTables[ksizeidx]); // Shift bits sequence
        if (strandR) shift_byte_array_right(br, ksizeTables[ksizeidx]); // Shift bits sequence

        // Add new nucleotide
        switch (c) {
            case 'A': // A = 00
                crrSeqL++;
                if (strandR) br[0] |= 192;
                break;
            case 'C': // C = 01
                if (strandF) b[ksizeTables[ksizeidx] - 1] |= 1;
                if (strandR) br[0] |= 128;
                crrSeqL++;
                break;
            case 'G': // G = 10
                if (strandF) b[ksizeTables[ksizeidx] - 1] |= 2;
                if (strandR) br[0] |= 64;
                crrSeqL++;
                break;
            case 'T': // T = 11
                if (strandF) b[ksizeTables[ksizeidx] - 1] |= 3;
                crrSeqL++;
                break;
            default : // Bad formed sequence
                crrSeqL = 0;
                break;
        }
        pos++;
        if (crrSeqL >= (uint64_t) ksizeTables[ksizeidx]*4) { // Full well formed sequence
            if(strandF){
                memcpy(tableForward+idx, b, ksizeTables[ksizeidx]);
            }
            if(strandR){
                memcpy(tableReverse+idx, b, ksizeTables[ksizeidx]);
            }
            idx += ksizeTables[ksizeidx];
            
        }
        c = fgetc(database);

    }
    fprintf(stdout, "[INFO] Sequence of length %"PRIu64" has %"PRIu64" mers of size k=%d\n", pos, pos-ksizeTables[ksizeidx]*4, ksizeTables[ksizeidx]*4);
    fprintf(stdout, "[INFO] Used %"PRIu64" bytes, which are %"PRIu64" MegaBytes\n", idx, idx/(1024*1024));
    



    fclose(database);
    

    free(tableForward);
    free(tableReverse);
    
    
    
    
    return 0;
}