/*
 * @author Fernando Moreno Jabato <jabato@uma.es>
 * @licence all rights reserved to the author and BitLAB group (University
 *    of Malaga).
 */

#include "fixeReverseFragments.h"


// Functions
void readFragment(FragFile*,FILE*);
void writeFragment(FragFile,FILE*);
int fragmentComparator(FragFile,FragFile);
void readSequenceLength(uint64_t*,FILE*);
void writeSequenceLength(uint64_t*,FILE*);
void endianessConversion(char*,char*,int);
uint64_t loadStats(Read**,char*);
void changeOrder(Read**);
inline void seekStats(uint64_t,Read*,Read*,Read**,Read**);
void writeFragBuffer(FragFile*,FILE*,FILE*,uint64_t);
int partition(FragFile*,int,int);
void quicksort_F(FragFile*,int,int);
void SWAP_F(FragFile*,FragFile*,FragFile);
void copyFFile(FragFile*,FragFile);
void freeReadList(Read**);
void checkOrder(fnode**,bool);
void push(fnode**,fnode**);
void move(fnode**,fnode**);
void sortList(fnode**);
uint64_t loadFrag(FragFile**,FILE*,int64_t);


// MAIN
int main(int ac, char** av){
    // Check arguments
    if(ac!=4){
        fprintf(stderr, "Bad call error.\nUSE: fixeReverseFragments fastaFile fastaReverseFile reverseFrags\n");
        return -1;
    }

    /////////////////////////////////////////////////////////////////
    fprintf(stdout, "\tReverse: Starting reverse fragments program.\n");
    /////////////////////////////////////////////////////////////////

    // Variables
    FILE *frags, *fIndx, *fBuff;
    Read *f_reads, *r_reads, *curr_f, *curr_r;
    uint64_t num_f_reads, num_r_reads, metagLength, genoLength;
    uint64_t fragsInBuffer, buffersWritten = 0;
    FragFile *buffer;
    char *fname;
    bool removeIntermediataFiles = true;

    // Take stats
    num_f_reads = loadStats(&f_reads,av[1]);
    num_r_reads = loadStats(&r_reads,av[2]);

    if(num_f_reads != num_r_reads){
        fprintf(stderr, "Error:: diferent number of reads on fasta files:\n\tForward:%"PRIu64"\n\tReverse:%"PRIu64"\n", num_f_reads,num_r_reads);
        return -1;
    }

    // Change forward indexes
    curr_f = f_reads;
    do{
        curr_f->index = num_r_reads - curr_f->index - 1;
        curr_f = curr_f->next;
    }while(curr_f != NULL);

    // Re-order reverse-reads
    if(num_r_reads > 1)
        changeOrder(&f_reads);
 

    /////////////////////////////////////////////////////////////////
    fprintf(stdout, "\tReverse: Opening/creating necessary files.");
    fflush(stdout);
    /////////////////////////////////////////////////////////////////


    // Re-calculate coordinates
    // Open fragments file
    if((frags = fopen(av[3],"rb"))==NULL){
        fprintf(stderr, "Error opening fragments file.\n");
        return -1;    
    }
        // Read headers
        readSequenceLength(&metagLength,frags);
        readSequenceLength(&genoLength,frags);

    // Take memory for frags buffer
    if((buffer = (FragFile*) malloc(sizeof(FragFile)*FRAG_BUFFER_SIZE))==NULL){
        fprintf(stderr, "Error allocating memory for fragments buffer.\n");
        return -1;
    }
    // Memory for file names handler
    if((fname = (char*) malloc(sizeof(char)*MAX_FILE_LENGTH))==NULL){
        fprintf(stderr, "Error allocating memory for file names handler.\n");
        return -1;
    } 

    // Open intermediate files
    strcpy(fname,av[3]);
    if((fIndx = fopen(strcat(fname,".fixe.findx"),"wb"))==NULL){
        fprintf(stderr, "Error opening intermediate buffer index file.\n");
        return -1;
    }
    strcpy(fname,av[3]);
    if((fBuff = fopen(strcat(fname,".fixe.fbuff"),"wb"))==NULL){
        fprintf(stderr, "Error opening intermediate buffer file.\n");
        return -1;
    }


    /////////////////////////////////////////////////////////////////
    fprintf(stdout, " (Done)\n");
    fprintf(stdout, "\tReverse: Loading original fragments.");
    fflush(stdout);
    /////////////////////////////////////////////////////////////////

    // Start to read fragments and re-calculate values
    fragsInBuffer = 0;

    // Read first fragment
    readFragment(&buffer[fragsInBuffer],frags);
    
    // Seek stats
    seekStats(buffer[fragsInBuffer].seqY,f_reads,r_reads,&curr_f,&curr_r);

    // Change index
    buffer[fragsInBuffer].seqY = num_r_reads - buffer[fragsInBuffer].seqY - 1;

    /* Change coordinates process:
     *   1 - Re-calculate start => YStart = reverse_offset - Y + forward_offset
     *   2 - Start in forward is the end => YEnd = YStart
     *   3 - Calculate the real YStart => YStart = YEnd - length
     */
    buffer[fragsInBuffer].yStart = (curr_r->Lac + curr_r->length) - buffer[fragsInBuffer].yStart + curr_f->Lac;
    buffer[fragsInBuffer].yEnd = buffer[fragsInBuffer].yStart;
    buffer[fragsInBuffer].yStart = buffer[fragsInBuffer].yEnd - buffer[fragsInBuffer].length;

    // Update number of fragments read
    fragsInBuffer++;

    // Read fragment
    readFragment(&buffer[fragsInBuffer],frags);

    while(!feof(frags)){
        // Seek stats
        if(curr_f->index > buffer[fragsInBuffer].seqY)
            seekStats(buffer[fragsInBuffer].seqY,f_reads,r_reads,&curr_f,&curr_r);
        else if(curr_f->index < buffer[fragsInBuffer].seqY)
            seekStats(buffer[fragsInBuffer].seqY,curr_f,curr_r,&curr_f,&curr_r);

        // Change index
        buffer[fragsInBuffer].seqY = num_r_reads - buffer[fragsInBuffer].seqY - 1;

        // Change coordinates
        buffer[fragsInBuffer].yStart = (curr_r->Lac + curr_r->length) - buffer[fragsInBuffer].yStart + curr_f->Lac;
        buffer[fragsInBuffer].yEnd = buffer[fragsInBuffer].yStart;
        buffer[fragsInBuffer].yStart = buffer[fragsInBuffer].yEnd - buffer[fragsInBuffer].length;

        // Update number of fragments
        fragsInBuffer++;

        // Check buffer
        if(fragsInBuffer >= FRAG_BUFFER_SIZE){
            writeFragBuffer(buffer,fIndx,fBuff,fragsInBuffer); // Write buffer
            // Update info
            buffersWritten++;
            fragsInBuffer = 0;
        }

        // Read new fragment
        readFragment(&buffer[fragsInBuffer],frags);
    }

    // Check buffered fragments
    if(fragsInBuffer > 0){
        if(buffersWritten > 0){
           writeFragBuffer(buffer,fIndx,fBuff,fragsInBuffer);
        }else{ // Only one buffer
            /////////////////////////////////////////////////////////////////
            fprintf(stdout, " (Done)\n");
            fprintf(stdout, "\tReverse: Writting re-calculated fragments.");
            fflush(stdout);
            /////////////////////////////////////////////////////////////////

            // Sort buffer
            quicksort_F(buffer,0,fragsInBuffer-1);

            // Close intermediate files
            fclose(fIndx);
            fclose(fBuff);
            fclose(frags);

            // Open fragments
            if((frags = fopen(av[3],"wb"))==NULL){
                fprintf(stderr, "Error opening fragments file (w).\n");
                return -1;    
            }

            // Write headers
            writeSequenceLength(&metagLength,frags);
            writeSequenceLength(&genoLength,frags);

            // Start to write fragments
            uint64_t pos;
            for(pos = 0; pos < fragsInBuffer; ++pos)
                writeFragment(buffer[pos],frags);

            // Free and close
            fclose(frags);
            free(buffer);

            // Remove intermediate files
            if(removeIntermediataFiles){
                strcpy(fname,av[3]);
                remove(strcat(fname,".fixe.findx"));
                strcpy(fname,av[3]);
                remove(strcat(fname,".fixe.fbuff"));
            }

            /////////////////////////////////////////////////////////////////
            fprintf(stdout, " (Done)\n");
            fprintf(stdout, "\tReverse: Closing the program.\n");
            fflush(stdout);
            /////////////////////////////////////////////////////////////////

            free(fname);
            freeReadList(&f_reads);
            freeReadList(&r_reads);

            return 0;
        }
    }else if(fragsInBuffer == 0 && buffersWritten == 0){
        /////////////////////////////////////////////////////////////////
        fprintf(stdout, "\tReverse: Any fragment found.\n");
        fprintf(stdout, "\tReverse: Closing the program.\n");
        /////////////////////////////////////////////////////////////////
        // Close all
        free(buffer);

        freeReadList(&f_reads);
        freeReadList(&r_reads);        

        fclose(fIndx);
        fclose(fBuff);
        fclose(frags);

        // Remove intermediate files
        if(removeIntermediataFiles){
            strcpy(fname,av[3]);
            remove(strcat(fname,".fixe.findx"));
            strcpy(fname,av[3]);
            remove(strcat(fname,".fixe.fbuff"));
        }

        free(fname);

        return 0;
    }

    /////////////////////////////////////////////////////////////////
    fprintf(stdout, " (Done)\n");
    fprintf(stdout, "\tReverse: Writting re-calculated fragments.");
    fflush(stdout);
    /////////////////////////////////////////////////////////////////


    // Close unnecessary variables
    fclose(fIndx);
    fclose(fBuff);
    fclose(frags);

    // Free unnecessaru variables
    free(buffer);
    freeReadList(&f_reads);
    freeReadList(&r_reads);

    // Open necessary files
    if((frags = fopen(av[3],"wb"))==NULL){
        fprintf(stderr, "Error opening fragments file (w).\n");
        return -1;    
    }
        // Read headers
        writeSequenceLength(&metagLength,frags);
        writeSequenceLength(&genoLength,frags);

    // Open intermediate files
    strcpy(fname,av[3]);
    if((fIndx = fopen(strcat(fname,".fixe.findx"),"rb"))==NULL){
        fprintf(stderr, "Error opening intermediate buffer index file (r).\n");
        return -1;
    }
    strcpy(fname,av[3]);
    if((fBuff = fopen(strcat(fname,".fixe.fbuff"),"rb"))==NULL){
        fprintf(stderr, "Error opening intermediate buffer file (r).\n");
        return -1;
    }

    // Prepare necessary variables
    fnode *fragsList = NULL;
    uint64_t fragsUnread[buffersWritten];
    uint64_t positions[buffersWritten];
    uint64_t lastLoaded, activeBuffers = buffersWritten;
    FragFile *FragsBlock;

    // Take memory for hits
    if((FragsBlock = (FragFile*) malloc(sizeof(FragFile)*activeBuffers*READ_BUFF_LENGTH))==NULL){
        fprintf(stderr, "Error allocating memory for fragments block.\n");
        return -1;
    }

    // Read buffers info
    uint64_t i = 0;
    do{
        if(fread(&positions[i],sizeof(uint64_t),1,fIndx)!=1){
            fprintf(stderr, "Error reading position at fIndx.\n");
            return -1;
        }
        if(fread(&fragsUnread[i],sizeof(uint64_t),1,fIndx)!=1){
            fprintf(stderr, "Error reading unread fragments at fIndx.\n");
            return -1;
        }
        ++i;
    }while(i < activeBuffers);

    // Load first hits
    fnode *currNode = NULL;
    uint64_t read, blockIndex = 0;
    for(i=0 ;i<activeBuffers; ++i, blockIndex += READ_BUFF_LENGTH){
        currNode = (fnode*) malloc(sizeof(fnode));
        currNode->next = fragsList;
        currNode->frags = &FragsBlock[blockIndex];
        currNode->buff = i;
        fseek(fBuff,positions[i],SEEK_SET);
        read = loadFrag(&currNode->frags,fBuff,fragsUnread[i]); //
        currNode->index = 0;
        currNode->frags_loaded = read;
        // Update info
        positions[i] = (uint64_t) ftell(fBuff);
        fragsUnread[i]-=read;
        lastLoaded = i;
        fragsList = currNode;
    }

    // Assign head
    fragsList = currNode;

    // Sort hits
    sortList(&fragsList);    

    // Read hits & generate fragments
    while(activeBuffers > 0){
        // Write fragments
        writeFragment(fragsList->frags[fragsList->index],frags);

        // Move to next
        fragsList->index +=1;
        // Load new hit
        if(fragsList->index >= fragsList->frags_loaded){
            if(fragsUnread[fragsList->buff] > 0){
                if(fragsList->buff != lastLoaded){
                    fseek(fBuff,positions[fragsList->buff],SEEK_SET);
                    lastLoaded = fragsList->buff;
                }
                read = loadFrag(&fragsList->frags,fBuff,fragsUnread[fragsList->buff]); //
                fragsList->index = 0;
                fragsList->frags_loaded = read;
                positions[fragsList->buff] = (uint64_t) ftell(fBuff);
                fragsUnread[fragsList->buff] -= read;
                checkOrder(&fragsList,false);
            }else{
                checkOrder(&fragsList,true);
                activeBuffers--;
            }
        }else{
            checkOrder(&fragsList,false);
        }
    }

    /////////////////////////////////////////////////////////////////
    fprintf(stdout, " (Done)\n");
    fprintf(stdout, "\tReverse: Closing the program.\n");
    fflush(stdout);
    /////////////////////////////////////////////////////////////////

    // Close files
    free(FragsBlock);
    fclose(fIndx);
    fclose(fBuff);
    fclose(frags);

    // Remove intermediate files
    if(removeIntermediataFiles){
        strcpy(fname,av[3]);
        remove(strcat(fname,".fixe.findx"));
        strcpy(fname,av[3]);
        remove(strcat(fname,".fixe.fbuff"));
    }

    free(fname);

    // End
    return 0;
}



/*
 */
void writeFragment(FragFile frag, FILE *f){
    char tmpArray[8];
    if(htons(1)==1){
        //Big endian
        fwrite(&frag.diag, sizeof(int64_t), 1, f);
        fwrite(&frag.xStart, sizeof(uint64_t), 1, f);
        fwrite(&frag.yStart, sizeof(uint64_t), 1, f);
        fwrite(&frag.xEnd, sizeof(uint64_t), 1, f);
        fwrite(&frag.yEnd, sizeof(uint64_t), 1, f);
        fwrite(&frag.length, sizeof(uint64_t), 1, f);
        fwrite(&frag.ident, sizeof(uint64_t), 1, f);
        fwrite(&frag.score, sizeof(uint64_t), 1, f);
        fwrite(&frag.similarity, sizeof(float), 1, f);
        fwrite(&frag.seqX, sizeof(uint64_t), 1, f);
        fwrite(&frag.seqY, sizeof(uint64_t), 1, f);
        fwrite(&frag.block, sizeof(int64_t), 1, f);
        fputc(frag.strand, f);
    } else {
        //Little endian
        endianessConversion((char *)(&frag.diag), tmpArray, sizeof(int64_t));
        fwrite(tmpArray, sizeof(int64_t), 1, f);
        endianessConversion((char *)(&frag.xStart), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *)(&frag.yStart), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *)(&frag.xEnd), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *)(&frag.yEnd), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *)(&frag.length), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *)(&frag.ident), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *)(&frag.score), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *)(&frag.similarity), tmpArray, sizeof(float));
        fwrite(tmpArray, sizeof(float), 1, f);
        endianessConversion((char *)(&frag.seqX), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *)(&frag.seqY), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *)(&frag.block), tmpArray, sizeof(int64_t));
        fwrite(tmpArray, sizeof(int64_t), 1, f);
        fputc(frag.strand, f);
    }
}


/**
 * Function to read a fragment from the specified file
 */
void readFragment(FragFile *frag, FILE *f){
    char tmpArray[8];

    if(htons(1)==1){
        //big endian
        if(fread(&frag->diag, sizeof(int64_t), 1, f)!=1){
            if(feof(f))return;
            fprintf(stderr,"Error reading the HSP diagonal");
        }
        if(fread(&frag->xStart, sizeof(uint64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP xStart");
        }
        if(fread(&frag->yStart, sizeof(uint64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP yStart");
        }
        if(fread(&frag->xEnd, sizeof(uint64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP xEnd");
        }
        if(fread(&frag->yEnd, sizeof(uint64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP yEnd");
        }
        if(fread(&frag->length, sizeof(uint64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP length");
        }
        if(fread(&frag->ident, sizeof(uint64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP identities");
        }
        if(fread(&frag->score, sizeof(uint64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP score");
        }
        if(fread(&frag->similarity, sizeof(float), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP similarity");
        }
        if(fread(&frag->seqX, sizeof(uint64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP seqX");
        }
        if(fread(&frag->seqY, sizeof(uint64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP seqY");
        }
        if(fread(&frag->block, sizeof(int64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP block");
        }
        frag->strand = fgetc(f);
    } else {
        //little endian
        if(fread(tmpArray, sizeof(int64_t), 1, f)!=1){
            if(feof(f))return;
            fprintf(stderr,"Error reading the HSP diagonal");
        }
        endianessConversion(tmpArray, (char *)(&frag->diag), sizeof(int64_t)); 
        if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP xStart");
        }
        endianessConversion(tmpArray, (char *)(&frag->xStart), sizeof(uint64_t)); 
        if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP yStart");
        }
        endianessConversion(tmpArray, (char *)(&frag->yStart), sizeof(uint64_t)); 
        if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP xEnd");
        }
        endianessConversion(tmpArray, (char *)(&frag->xEnd), sizeof(uint64_t)); 
        if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP yEnd");
        }
        endianessConversion(tmpArray, (char *)(&frag->yEnd), sizeof(uint64_t)); 
        if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP length");
        }
        endianessConversion(tmpArray, (char *)(&frag->length), sizeof(uint64_t)); 
        if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP identity");
        }
        endianessConversion(tmpArray, (char *)(&frag->ident), sizeof(uint64_t)); 
        if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP score");
        }
        endianessConversion(tmpArray, (char *)(&frag->score), sizeof(uint64_t)); 
        if(fread(tmpArray, sizeof(float), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP float");
        }
        endianessConversion(tmpArray, (char *)(&frag->similarity), sizeof(float)); 
        if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP seqX");
        }
        endianessConversion(tmpArray, (char *)(&frag->seqX), sizeof(uint64_t)); 
        if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP seqY");
        }
        endianessConversion(tmpArray, (char *)(&frag->seqY), sizeof(uint64_t)); 
        if(fread(tmpArray, sizeof(int64_t), 1, f)!=1){
            fprintf(stderr,"Error reading the HSP block");
        }
        endianessConversion(tmpArray, (char *)(&frag->block), sizeof(int64_t)); 
        frag->strand = fgetc(f);
    }
}


/* This function is used to comapre two fragments. The order of comparation is:
 *     1 - SeqX
 *     2 - SeqY
 *     3 - diag
 *     4 - xStart
 *     5 - length
 *     6 - yStart
 *     7 - block
 *     8 - strand
 *     9 - score
 *  @param f1 fragment to be compared.
 *  @param f2 fragment to be compared.
 *  @return negative number if f1<f2, positive number if f1>f2 or zero if f1=f2.
 */
int fragmentComparator(FragFile f1, FragFile f2){
    if(f1.seqX > f2.seqX) return 1;
    else if(f1.seqX < f2.seqX) return -1;

    if(f1.seqY > f2.seqY) return 1;
    else if(f1.seqY < f2.seqY) return -1;

    if(f1.diag > f2.diag) return 1;
    else if(f1.diag < f2.diag) return -1;

    if(f1.xStart > f2.xStart) return 1;
    else if(f1.xStart < f2.xStart) return -1;

    if(f1.length > f2.length) return 1;
    else if(f1.length < f2.length) return -1;

    if(f1.yStart > f2.yStart) return 1;
    else if(f1.yStart < f2.yStart) return -1;

    if(f1.block > f2.block) return 1;
    else if(f1.block < f2.block) return -1;

    if(f1.strand == 'f' && f2.strand == 'r') return 1;
    else if(f1.strand == 'r' && f2.strand == 'f') return -1;

    if(f1.score > f2.score) return 1;
    else if(f1.score < f2.score) return -1;

    return 0;
}

/*
 */
int GT(FragFile f1, FragFile f2){
    if(f1.seqX > f2.seqX) return 1;
    else if(f1.seqX < f2.seqX) return 0;

    if(f1.seqY > f2.seqY) return 1;
    else if(f1.seqY < f2.seqY) return 0;

    if(f1.diag > f2.diag) return 1;
    else if(f1.diag < f2.diag) return 0;

    if(f1.xStart > f2.xStart) return 1;
    else if(f1.xStart < f2.xStart) return 0;

    if(f1.length > f2.length) return 1;
    else if(f1.length < f2.length) return 0;

    if(f1.yStart > f2.yStart) return 1;
    else if(f1.yStart < f2.yStart) return 0;

    if(f1.block > f2.block) return 1;
    else if(f1.block < f2.block) return 0;

    if(f1.strand == 'f' && f2.strand == 'r') return 1;
    else if(f1.strand == 'r' && f2.strand == 'f') return 0;

    if(f1.score > f2.score) return 1;
    else return 0;
}


/**
 * Function to read the sequence length
 */
void readSequenceLength(uint64_t *length, FILE *f){
    char tmpArray[8];
    if(htons(1)==1){
        //big endian
        if(fread(length, sizeof(uint64_t), 1, f)!=1){
            fprintf(stderr,"Error reading sequence length");
        }
    } else {
        //little endian
        if(fread(tmpArray, sizeof(uint64_t), 1, f)!=1){
            fprintf(stderr,"Error reading sequence length");
        }
        endianessConversion(tmpArray, (char *)length, sizeof(uint64_t));
    }
}

/**
 * Function to write the sequence length
 */
void writeSequenceLength(uint64_t *length, FILE *f){
    char tmpArray[8];
    if(htons(1)==1){
        //big endian
        fwrite(length, sizeof(uint64_t), 1, f);
    } else {
        //little endian
        endianessConversion((char *)length, tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
    }
}


/*
 */
void endianessConversion(char *source, char *target, int numberOfBytes){
    int i,j;
    for(i=numberOfBytes-1;i>=0;i--){
        j=numberOfBytes-1-i;
        target[j]=source[i];
    }
}


/* This function is used to load the stats of a fasta multi-sequence file.
 *  @note based on mgReadsIndex.c program. 
 *  @param head is the Read linked list first element.
 *  @param filePath is the absolute/relative path to fasta file.
 *  @return the total of reads loaded.
 */
uint64_t loadStats(Read **head,char *filePath){
    // Variables
    FILE *file;
    char c;
    uint64_t readLen, totalLen=0, nReads=0;
    *head = NULL;
    Read *currRead = NULL,*lastRead = NULL;

    // Open file
    if((file = fopen(filePath,"rt"))==NULL){
        fprintf(stderr, "Error opening stats file: %s\n", filePath);
        return '\0';
    }    

    // initialize first read
    currRead = (Read*) malloc(sizeof(Read));
    currRead->Lac = 0;
    currRead->index = nReads;

    // Read
    c = fgetc(file);
    while(!feof(file)){
        readLen=0; //Reset read length
        // Seek next header
        while(c!='>') c = (char)fgetc(file);
        // Save read init position on file
        currRead->pos = ftell(file) - 1; 
        // Seek end of header
        while(c!='\n') c = (char)fgetc(file);
        // Read sequence
        while(c != '>' && !feof(file)){
            c = toupper(c);
            if(c > 64 && c < 91) readLen++;
            c = (char)fgetc(file);
        }
        // Store length
        currRead->length = readLen;

        // Store current read
        if(*head == NULL){ // First read
            *head = currRead; 
        }else{
            // Link with last node
            lastRead->next = currRead;
        }
        // Update last node
        lastRead = currRead;

        // Take memory for new read
        currRead = (Read*) malloc(sizeof(Read));
        // Update info
        currRead->Lac = lastRead->Lac + readLen + 1;
        currRead->next = NULL;
        currRead->index = ++nReads;

        totalLen += readLen;
    }

    // Free last node (it's not necessary)
    free(currRead);
    // Close file & end
    fclose(file);

    return nReads;
}


/* This function is used to change the linked list order to the reverse order.
 *  @param head is the Read linked list first element.
 */
void changeOrder(Read** head){
    // Variables
    Read *currRead,*last = NULL,*next;

    // Init
    currRead = *head;

    // Change order
    next = currRead->next;
    currRead->next = last;

    // Update
    last = currRead;
    currRead = next;

    while(currRead->next != NULL){
        // Take next
        next = currRead->next;

        // Change order
        currRead->next = last;

        // Update last & current
        last = currRead;    
        currRead = next;
    }

    // Link last & convert in head
    currRead->next = last;

    *head = currRead;
}


/* This function is used to seek a stats of a specific read using the read_index
 *  @warning both read linked list must be coordinated!!
 *  @param index is the index of the wanted read.
 *  @param start1 first element of one of the linked lists.
 *  @param start2 first element of one of the linked lists.
 *  @param stat1 stats wanted of the first linked list. If couldn't found, NULL will be stored.
 *  @param stat2 stats wanted of the second linked list. If couldn't found, NULL will be stored.
 */
inline void seekStats(uint64_t index,Read *start1,Read *start2,Read **stat1,Read **stat2){
    // Init
    *stat1 = start1;
    *stat2 = start2;

    // Check
    if((*stat1)->index == index && (*stat2)->index == index) return; // Found!

    // Start to search
    while((*stat1)->next != NULL || (*stat2)->next != NULL){
        // Move next
        *stat1 = (*stat1)->next;
        *stat2 = (*stat2)->next;
        // Check
        if((*stat1)->index == index && (*stat2)->index == index) return; // Found!
    }

    // Not found
    *stat1 = NULL;
    *stat2 = NULL;    
}


/* This function is used to write a buffer in the intermediate files. The order
 * in each dictioanry is:
 *   - Index: each entrance: Pos<uint64_t> FragsInBuff<uint64_t>
 *   - Frags: each entrance is a FragFile structure.
 *  @param buff buffer to be written.
 *  @param index intermediate file.
 *  @param frags intermediate file.
 *  @param fragsInBuff number of words stored on buffer.
 */
void writeFragBuffer(FragFile* buff,FILE* index,FILE* frags,uint64_t fragsInBuff){
    // Sort buffer
    quicksort_F(buff,0,fragsInBuff-1);

    // Write info on index file
    uint64_t pos = (uint64_t) ftell(frags);
    
    fwrite(&pos,sizeof(uint64_t),1,index);
    fwrite(&fragsInBuff,sizeof(uint64_t),1,index);
    
    // Write hits in hits file
    for(pos=0; pos<fragsInBuff; ++pos)
        writeFragment(buff[pos],frags);
}


/* This function is necessary for quicksort functionality.
 *  @param arr array to be sorted.
 *  @param left inde of the sub-array.
 *  @param right index of the sub-array.
 */
int partition(FragFile* arr, int left, int right){
   int i = left;
   int j = right + 1;
   FragFile t;

   // Pivot variable
   int pivot = (left+right)/2;

   if(GT(arr[pivot],arr[right]))
         SWAP_F(&arr[pivot],&arr[right],t);

   if(GT(arr[pivot],arr[left]))
         SWAP_F(&arr[pivot],&arr[left],t);

   if(GT(arr[left],arr[right]))
         SWAP_F(&arr[left],&arr[right],t);

    while(1){
        do{
            ++i;
        }while(!GT(arr[i],arr[left]) && i <= right);

        do{
            --j;
        }while(GT(arr[j],arr[left]) && j >= left);

        if(i >= j) break;

        SWAP_F(&arr[i],&arr[j],t);
    }

    SWAP_F(&arr[left],&arr[j],t);

    return j;
}


/* This function is used to sort a FragFile array.
 *  @param arr array to be sorted.
 *  @param left index where start to sort.
 *  @param right index where end sorting action.
 */
void quicksort_F(FragFile* arr, int left,int right){
    int j;

    if(left < right){
        // divide and conquer
        j = partition(arr,left,right);
        quicksort_F(arr,left,j-1);
        quicksort_F(arr,j+1,right);
   }
}


/* This function is used to swap two fragmentss variables
 *  @param h1 frag to be swapped.
 *  @param h2 frag to be swapped.
 *  @param t auxiliar frag.
 */
void SWAP_F(FragFile* h1, FragFile* h2, FragFile t){
    copyFFile(&t,*h1);
    copyFFile(h1,*h2);
    copyFFile(h2,t);
}

/* This function is used to copy a FragFile variable in another FragFile variable.
 *  @param toCopy where frag will be copied.
 *  @param copy frag to be copied.
 */
void copyFFile(FragFile* toCopy, FragFile copy){
    toCopy->diag = copy.diag;
    toCopy->xStart = copy.xStart;
    toCopy->xEnd = copy.xEnd;
    toCopy->yStart = copy.yStart;
    toCopy->yEnd = copy.yEnd;
    toCopy->length = copy.length;
    toCopy->ident = copy.ident;
    toCopy->score = copy.score;
    toCopy->similarity = copy.similarity;
    toCopy->seqX = copy.seqX;
    toCopy->seqY = copy.seqY;
    toCopy->block = copy.block;
    toCopy->strand = copy.strand;
}


/*
 */
void freeReadList(Read** head){
    Read *aux;
    while((*head)->next != NULL){
        aux = *head;
        *head = (*head)->next;
        free(aux);
    }
    free(*head);
}


/* This method push node B after A (A->C ==PUSH==> A->B->C)
 *  @param A node after B will be pushed.
 *  @param B node to be pushed.
 */
void push(fnode **A,fnode **B){
    (*B)->next = (*A)->next;
    (*A)->next = *B;
}


/* Move node after B to after A position and make linked list consistent.
 *  @param A reference node.
 *  @param B node after it will be moved.
 */
void move(fnode **A,fnode **B){
    fnode *temp = (*B)->next->next;
    push(A,&(*B)->next);
    (*B)->next = temp;
}


/* This emthod sort a linked list
 *  @param first node of the linked list.
 */
void sortList(fnode **first){
    if((*first)->next == NULL) return; // Linked list with only one element

    fnode *current = *first;
    fnode *aux;
    bool sorted = false;
    // Do until end
    while(!sorted){
        if(current->next == NULL) sorted = true;
        else if(GT(current->next->frags[current->index],current->frags[current->next->index])==0){ // Next is smaller
            // Search position
            if(GT(current->next->frags[current->next->index],(*first)->frags[(*first)->index])==0){ // New first node
                aux = current->next->next;
                current->next->next = *first;
                *first = current->next;
                current->next = aux;
            }else{ // Search position
                aux = *first;           
                while(GT(aux->next->frags[aux->next->index],current->next->frags[current->next->index])==1)
                    aux = aux->next;
                move(&aux,&current);
                // Chekc if it's the last node
                if(current->next == NULL) sorted = true;
            }
        }else{ // Go next
            current = current->next;
            if(current->next == NULL){ // End of the list
                // List sorted
                sorted = true;
            }
        }
    }
}


/* This function is used to check the correc order of the first node of a linked list.
 * If it's incorrect, this function sort it.
 *  @param list linked list to be checked.
 *  @param discardFirst a boolean value that indicate if first node should be deleted.
 */
void checkOrder(fnode** list,bool discardFirst){
    fnode *aux;
    if(discardFirst){
        aux = *list;
        *list = (*list)->next;
        free(aux);
    }else if((*list)->next != NULL){ // Check new position
        // Search new position
        if(GT((*list)->frags[(*list)->index],(*list)->next->frags[(*list)->next->index])==1){
            fnode *curr = (*list)->next;
            while(1){
                if(curr->next == NULL) break; // End of list
                else if(GT((*list)->frags[(*list)->index],curr->next->frags[curr->next->index])==0) break; // position found
                else curr = curr->next;
            }
            aux = (*list)->next;
            (*list)->next = curr->next;
            curr->next = *list;
            *list = aux;
        }
    }
}


/* This function is used to load a fragfile from a frags intermediate file.
 *  @param frag varaible where loaded frag will be stored.
 *  @param hbuff pointer to frags intermediate file.
 *  @param unread rest of frags on intermediate file.
 *  @return Number of frags read from intermediate file.
 */
uint64_t loadFrag(FragFile **frag,FILE* hbuff, int64_t unread){
    uint64_t j;
    for(j=0; j<READ_BUFF_LENGTH && unread > 0; ++j){
        readFragment(&(*frag)[j],hbuff);
        unread--;
    }
    return j;
}