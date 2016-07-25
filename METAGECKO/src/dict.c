/*
 * @author Fernando Moreno Jabato <jabato@uma.es>
 * @licence all rights reserved to the author and BitLAB group (University
 *    of Malaga).
 */
#include "dict.h"

/* This main contains the workflow to read a FASTA format file and take the sequences
 * words to generate a hashtable that links all found words (without repetitions) with 
 * the positions where it appears.
 * The way to invoke the program is:
 *   @use dict FASTA_File out WL
 * Where the parameters used are:
 *   @param FASTA_File is the FASTA format file where sequence that will be readed is 
 *          stored.
 *   @param out is the basename of the file where fragmetns will be stored.
 *   @param WL is the word length to use to create the dictionary. Only well formed
 *          words of length = WL will be stored on dictionary.
 *   @param strand is the direction of the sequence used to calculate the dictionary.
 *          can be 'f' (forward) or 'r' (reverse). NOTE: this parameter is optional. If
 *          it's not specified, both strands will be used to calculate the dictionary.
 * The result of the program will be the following:
 *   @file <out>.d2hW is a file that contains the words and a linkage to the second 
 *         final file. This linkage is formed by the number of positions that are linked
 *         to the current word and where start to read this positions in the second file.
 *   @file <out>.d2hP is a file that contains the pairs of coordinates. This pairs are
 *         formed by the sequence index and the position on the sequence where the
 *         linked word appears.
 * The program will write some messages in default output stream that shows the action that
 * is being done by the program each moment.
 * Warning: if any error is launched, the program automatically stops, if it happens is
 * normal that the following files appear in your output folder:
 *   @file <out>.bindx is an auxiliary file with information about some buffer used in the 
 *         program process.
 *   @file <out>.wrds is an auxilary file with information about the words read during 
 *         the read process.
 */
int main(int ac, char **av) {
    // Check arguments
    if (ac != 4 && ac != 5) {
        fprintf(stderr, "Bad call error.\nUSE: dict metag.IN dictName WL\nOR:  dict metag.IN dictName WL strand\n");
        return -1;
    }

    // Check arguments
    if (!exists(av[1])) { // Check if metagenome file exists
        fprintf(stderr, "Error:: Metagenome file specified doesn't exists.\n");
        return -1;
    }
    if (!is_int(av[3])) {
        fprintf(stderr, "Error:: Word length specified isn't a number.\n");
        return -1;
    } else if (atoi(av[3]) % 4 != 0 && atoi(av[3]) > 0) {
        fprintf(stderr, "Error:: Word length must be a 4 multiple and positive.\n");
        return -1;
    }
    if (ac == 5) {
        if (av[4][0] != 'f' && av[4][0] != 'r') { // Check forward/reverse argument
            fprintf(stderr, "Error:: Strand argument must be <f> or <r>\n");
            return -1;
        }
    }

    /////////////////////////// CHECKPOINT ///////////////////////////
    fprintf(stdout, "\tDict: Starting dictionaries program.\n");
    /////////////////////////// CHECKPOINT ///////////////////////////

    // Variables
    uint16_t WL = (uint16_t) atoi(av[3]); // Word length
    BYTES_IN_WORD = WL / 4;
    wentry *buffer; //Buffer of read words
    FILE *metag; // Input files
    FILE *wDic, *pDic;  // Output files
    FILE *bIndx, *wrds; // Intermediate files
    uint64_t wordsInBuffer = 0, i; // Absolute and buffer read words and auxiliar variable(i)
    uint32_t numBuffWritten = 0;
    char *fname;
    unsigned char *WordsBlock;
    bool removeIntermediataFiles = true; // Config it if you want save or not the itnermediate files
    bool strandF = false, strandR = false;

    if (ac == 5) {
        if (av[4][0] == 'f') strandF = true;
        else if (av[4][0] == 'r') strandR = true;
    } else {
        strandF = true;
        strandR = true;
    }

    // Check
    if (!strandF && !strandR) {
        fprintf(stderr, "Error, no strand specified.\n");
        return -1;
    }

    /////////////////////////// CHECKPOINT ///////////////////////////
    fprintf(stdout, "\tDict: Opening/creating necessary files.");
    fflush(stdout);
    /////////////////////////// CHECKPOINT ///////////////////////////

    // Allocate necessary memory
    // Memory for buffer
    if ((buffer = (wentry *) malloc(sizeof(wentry) * BUFFER_LENGTH)) == NULL) {
        fprintf(stderr, "Error allocating memory for words buffer.\n");
        return -1;
    }

    // Memory for file names handler
    if ((fname = (char *) malloc(sizeof(char) * MAX_FILE_LENGTH)) == NULL) {
        fprintf(stderr, "Error allocating memory for file names handler.\n");
        return -1;
    }

    // Memory for sequences
    if ((WordsBlock = (unsigned char *) malloc(sizeof(unsigned char) * BUFFER_LENGTH * BYTES_IN_WORD)) == NULL) {
        fprintf(stderr, "Error allocating memory for sequences.\n");
        return -1;
    }

    // Open current necessary files
    // Open metagenome file
    if ((metag = fopen(av[1], "rt")) == NULL) {
        fprintf(stderr, "Error opening metagenome file.\n");
        return -1;
    }

    //Open intermediate files
    strcpy(fname, av[2]); // Copy outDic name
    if ((bIndx = fopen(strcat(fname, ".bindx"), "wb")) == NULL) {
        fprintf(stderr, "Error opening buffer index file (write).\n");
        return -1;
    }

    // Write first line of buffer index file (WordLength)
    fwrite(&WL, sizeof(uint16_t), 1, bIndx);

    // Open words repo
    strcpy(fname, av[2]);
    if ((wrds = fopen(strcat(fname, ".wrds"), "wb")) == NULL) {
        fprintf(stderr, "Error opening words repository (write).\n");
        return -1;
    }

    /////////////////////////// CHECKPOINT ///////////////////////////
    fprintf(stdout, " (Done)\n");
    fprintf(stdout, "\tDict: Reading k-mers.");
    if (strandF && strandR) fprintf(stdout, "[f,r]");
    else if (strandF) fprintf(stdout, "[f]");
    else if (strandR) fprintf(stdout, "[r]");
    fflush(stdout);
    /////////////////////////// CHECKPOINT ///////////////////////////

    // START WORKFLOW
    // Read metagenome file
    // Necessary variables
    uint64_t seqPos = 0, crrSeqL = 0;
    char c; // Auxiliar -> Read char per char
    wentry temp, rev_temp; // Auxiliar word variable
    if ((temp.w.b = (unsigned char *) malloc(sizeof(unsigned char) * BYTES_IN_WORD)) == NULL) {
        fprintf(stderr, "Error allocating memory for temp.\n");
        return -1;
    }
    if ((rev_temp.w.b = (unsigned char *) malloc(sizeof(unsigned char) * BYTES_IN_WORD)) == NULL) {
        fprintf(stderr, "Error allocating memory for rev_temp.\n");
        return -1;
    }

    memset(temp.w.b, 0, BYTES_IN_WORD);
    memset(rev_temp.w.b, 0, BYTES_IN_WORD);

    temp.loc.seq = -1;
    temp.loc.strand = 'f';
    rev_temp.loc.strand = 'r';

    // Memory for buffer words
    uint64_t blockIndex = 0;
    for (i = 0; i < BUFFER_LENGTH; i++, blockIndex += BYTES_IN_WORD)
        buffer[i].w.b = &WordsBlock[blockIndex];


    // Memory and variables for reading buffer
    uint64_t posBuffer = READBUF + 1, tReadBuffer = 0;
    char *readBuffer = (char *) malloc(READBUF * sizeof(char));
    if (readBuffer == NULL) terror("Could not allocate memory for reading buffer");


    // Start to read
    c = buffered_fgetc(readBuffer, &posBuffer, &tReadBuffer, metag);
    while (!feof(metag) || (feof(metag) && posBuffer < tReadBuffer)) {
        // Check if it's a special line
        if (!isupper(toupper(c))) { // Comment, empty or quality (+) line
            if (c == '>') { // Comment line
                c = buffered_fgetc(readBuffer, &posBuffer, &tReadBuffer, metag);
                while (c != '\n') // Avoid comment line
                    c = buffered_fgetc(readBuffer, &posBuffer, &tReadBuffer, metag);

                temp.loc.seq++; // New sequence
                crrSeqL = 0; // Reset buffered sequence length

                seqPos++; // Absolute coordinates: add one for the "*"
            }
            c = buffered_fgetc(readBuffer, &posBuffer, &tReadBuffer, metag); // First char of next sequence
            continue;
        }
        if (strandF) shift_word_left(&temp.w); // Shift bits sequence
        if (strandR) shift_word_right(&rev_temp.w); // Shift bits sequence



        // Add new nucleotid
        switch (c) {
            case 'A': // A = 00
                crrSeqL++;
                if (strandR) rev_temp.w.b[0] |= 192;
                break;
            case 'C': // C = 01
                if (strandF) temp.w.b[BYTES_IN_WORD - 1] |= 1;
                if (strandR) rev_temp.w.b[0] |= 128;
                crrSeqL++;
                break;
            case 'G': // G = 10
                if (strandF) temp.w.b[BYTES_IN_WORD - 1] |= 2;
                if (strandR) rev_temp.w.b[0] |= 64;
                crrSeqL++;
                break;
            case 'T': // T = 11
                if (strandF) temp.w.b[BYTES_IN_WORD - 1] |= 3;
                crrSeqL++;
                break;
            default : // Bad formed sequence
                crrSeqL = 0;
                break;
        }
        seqPos++;
        if (crrSeqL >= (uint64_t) WL) { // Full well formed sequence
            if (strandF) {
                temp.loc.pos = seqPos - WL - 1; // Take position on read
                // Store the new word
                storeWord(&buffer[wordsInBuffer], temp);

                wordsInBuffer++;

                if (wordsInBuffer == BUFFER_LENGTH) { // Buffer is full
                    if (writeBuffer(buffer, bIndx, wrds, wordsInBuffer) < 0) {
                        return -1;
                    } else {
                        // Update info
                        //readW += wordsInBuffer;
                        wordsInBuffer = 0;
                        numBuffWritten++;
                    }
                }
            }
            if (strandR) {
                rev_temp.loc.pos = seqPos - 1; // Take position on read
                rev_temp.loc.seq = temp.loc.seq;

                // Store the new word
                storeWord(&buffer[wordsInBuffer], rev_temp);

                wordsInBuffer++;

                if (wordsInBuffer == BUFFER_LENGTH) { // Buffer is full
                    if (writeBuffer(buffer, bIndx, wrds, wordsInBuffer) < 0) {
                        return -1;
                    } else {
                        // Update info
                        //readW += wordsInBuffer;
                        wordsInBuffer = 0;
                        numBuffWritten++;
                    }
                }
            }
        }
        c = buffered_fgetc(readBuffer, &posBuffer, &tReadBuffer, metag);
    }


    // Free&Close unnecessary variables
    free(readBuffer);
    fclose(metag);
    // Write buffered words
    if (wordsInBuffer != 0) {
        if (numBuffWritten == 0) { // Special case, only one buffer
            /////////////////////////// CHECKPOINT ///////////////////////////
            fprintf(stdout, " (Done)\n");
            fprintf(stdout, "\tDict: Writting dictionary.");
            fflush(stdout);
            /////////////////////////// CHECKPOINT ///////////////////////////

            quicksort_W(buffer, 0, wordsInBuffer - 1);

            // Close unnecessary files
            fclose(bIndx);
            fclose(wrds);

            // Open necessary files
            // Open positions file
            strcpy(fname, av[2]); // Copy outDic name
            if ((pDic = fopen(strcat(fname, ".d2hP"), "wb")) == NULL) {
                fprintf(stderr, "Error opening positions dictionary file.\n");
                return -1;
            }

            // Open dictionary words file
            strcpy(fname, av[2]); // Copy outDic name
            if ((wDic = fopen(strcat(fname, ".d2hW"), "wb")) == NULL) {
                fprintf(stderr, "Error opening words dictionary file.\n");
                return -1;
            }
            fwrite(&WL, sizeof(uint16_t), 1, wDic); // Word length

            // Declare necessary variables
            uint32_t reps = 0;

            // First entrance
            fwrite(buffer[0].w.b, sizeof(unsigned char), BYTES_IN_WORD, wDic); // write first word
            uint64_t pos = (uint64_t) ftello(pDic);
            fwrite(&pos, sizeof(uint64_t), 1, wDic); // position on pDic
            fwrite(&buffer[0].loc, sizeof(LocationEntry), 1, pDic); // Location

            reps++; // Increment number of repetitions
            storeWord(&temp, buffer[0]); // Update last word written

            // Write final dictionary file
            uint64_t index;
            for (index = 1; index < wordsInBuffer; ++index) {
                int cmp = wordcmp(buffer[index].w.b, temp.w.b, BYTES_IN_WORD * 4);
                // Store word in buffer
                writeWord(&buffer[index], wDic, pDic,
                          cmp == 0 ? true : false, &reps);
                if (cmp != 0)
                    storeWord(&temp, buffer[index]); // Update last word written
            }

            // Write last word repetitions
            fwrite(&reps, sizeof(uint32_t), 1, wDic); // Write num of repetitions

            /////////////////////////// CHECKPOINT ///////////////////////////
            fprintf(stdout, " (Done)\n");
            fprintf(stdout, "\tDict: Closing the program.\n");
            fflush(stdout);
            /////////////////////////// CHECKPOINT ///////////////////////////

            // Close files
            fclose(pDic);
            fclose(wDic);

            // Free memory
            free(WordsBlock);
            free(buffer);
            free(temp.w.b);
            free(rev_temp.w.b);

            // Remove itnermediate files if it's necessary
            if (removeIntermediataFiles) {
                strcpy(fname, av[2]);
                remove(strcat(fname, ".bindx"));
                strcpy(fname, av[2]);
                remove(strcat(fname, ".wrds"));
            }

            free(fname);

            // End program
            return 0;
        } else if (writeBuffer(buffer, bIndx, wrds, wordsInBuffer) < 0) {
            return -1;
        } else numBuffWritten++;
    } else if (wordsInBuffer == 0 && numBuffWritten == 0) { // Special case
        /////////////////////////// CHECKPOINT ///////////////////////////
        fprintf(stdout, "\tDict: Any k-mer found.\n");
        fprintf(stdout, "\tDict: Closing the program.\n");
        /////////////////////////// CHECKPOINT ///////////////////////////

        free(WordsBlock);
        free(buffer);
        free(temp.w.b);
        free(rev_temp.w.b);

        fclose(bIndx);
        fclose(wrds);

        // Remove itnermediate files if it's necessary
        if (removeIntermediataFiles) {
            strcpy(fname, av[2]);
            remove(strcat(fname, ".bindx"));
            strcpy(fname, av[2]);
            remove(strcat(fname, ".wrds"));
        }

        free(fname);

        return 0;
    }

    // Free buffer space
    free(WordsBlock);
    free(buffer);

    /////////////////////////// CHECKPOINT ///////////////////////////
    fprintf(stdout, " (Done)\n");
    fprintf(stdout, "\tDict: Writing dictionary.");
    fflush(stdout);
    /////////////////////////// CHECKPOINT ///////////////////////////

    // Read Intermediate files and create final dictionary
    // Close write buffer and open read buffers
    fclose(bIndx);
    fclose(wrds);

    strcpy(fname, av[2]); // Copy outDic name
    if ((bIndx = fopen(strcat(fname, ".bindx"), "rb")) == NULL) {
        fprintf(stderr, "Error opening buffer index file (read).\n");
        return -1;
    }

    if (fread(&WL, sizeof(uint16_t), 1, bIndx) != 1) { // Read header = word length
        fprintf(stderr, "Error reading word length at header. Value not found.\n");
        return -1;
    }

    // Open words repo
    strcpy(fname, av[2]);
    if ((wrds = fopen(strcat(fname, ".wrds"), "rb")) == NULL) {
        fprintf(stderr, "Error opening words repository (read).\n");
        return -1;
    }

    // Open positions file
    strcpy(fname, av[2]); // Copy outDic name
    if ((pDic = fopen(strcat(fname, ".d2hP"), "wb")) == NULL) {
        fprintf(stderr, "Error opening positions dictionary file.\n");
        return -1;
    }

    // Open dictionary words file
    strcpy(fname, av[2]); // Copy outDic name
    if ((wDic = fopen(strcat(fname, ".d2hW"), "wb")) == NULL) {
        fprintf(stderr, "Error opening words dictionary file.\n");
        return -1;
    }

    fwrite(&WL, sizeof(uint16_t), 1, wDic); // Word length

    // Prepare necessary variables
    uint64_t arrPos[numBuffWritten];
    int64_t wordsUnread[numBuffWritten];
    node_W *words = NULL;
    uint32_t reps = 0;
    uint64_t lastLoaded = -1, activeBuffers = numBuffWritten;
    unsigned char *BuffWordsBlock;
    wentry *WentryBlock;

    // Memory for sequences
    if ((BuffWordsBlock = (unsigned char *) malloc(
            sizeof(unsigned char) * activeBuffers * BYTES_IN_WORD * READ_BUFF_LENGTH)) == NULL) {
        fprintf(stderr, "Error allocating memory for sequences (merge).\n");
        return -1;
    }

    if ((WentryBlock = (wentry *) malloc(sizeof(wentry) * activeBuffers * READ_BUFF_LENGTH)) == NULL) {
        fprintf(stderr, "Error allocating memory for wentry block.\n");
        return -1;
    }

    // Read info about buffers
    i = 0;
    uint64_t aux64;
    do {
        if (fread(&arrPos[i], sizeof(uint64_t), 1, bIndx) != 1) { // Position on words file
            fprintf(stderr, "Error reading position on P file.\n");
            return -1;
        }
        if (fread(&aux64, sizeof(uint64_t), 1, bIndx) != 1) {// Number of words on set
            fprintf(stderr, "Error reading number of repetitions at W file.\n");
            return -1;
        }
        wordsUnread[i] = (int64_t) aux64;
        ++i;
    } while (i < activeBuffers);
    fclose(bIndx);

    // Take memory for words & read first set of words
    blockIndex = 0;
    uint64_t wblockIndex = 0, j, read;
    node_W *currNode = NULL;
    for (i = 0; i < activeBuffers; ++i, wblockIndex += READ_BUFF_LENGTH) {
        currNode = (node_W *) malloc(sizeof(node_W)); // New node
        currNode->word = &WentryBlock[wblockIndex]; // Wentry array
        for (j = 0; j < READ_BUFF_LENGTH; ++j, blockIndex += BYTES_IN_WORD) // Words memory
            currNode->word[j].w.b = &BuffWordsBlock[blockIndex];
        currNode->next = words;
        fseeko(wrds, arrPos[i], SEEK_SET);
        read = loadWord(currNode->word, wrds, wordsUnread[i]);
        currNode->buff = i;
        currNode->index = 0;
        currNode->words_loaded = read;
        words = currNode;
        // Update info
        arrPos[i] = (uint64_t) ftello(wrds);
        wordsUnread[i] -= read;
    }
    words = currNode;

    // Sort array
    sortList(&words);

    // First entrance
    fwrite(words->word[words->index].w.b, sizeof(unsigned char), BYTES_IN_WORD, wDic); // write first word
    uint64_t pos = (uint64_t) ftello(pDic);
    fwrite(&pos, sizeof(uint64_t), 1, wDic); // position on pDic
    fwrite(&words->word[words->index].loc, sizeof(LocationEntry), 1, pDic); // Location

    reps++; // Increment number of repetitions
    storeWord(&temp, words->word[words->index]); // Update last word written
    words->index += 1;

    //REVIEW BUT I THINK IS DEAD CODE
//    if (words->index >= words->words_loaded) {
//        fprintf(stderr, "Entro\n");
//        if (wordsUnread[words->buff] > 0) {
//            fseeko(wrds, arrPos[words->buff], SEEK_SET);
//            read = loadWord(words->word, wrds, wordsUnread[words->buff]);
//            words->index = 0;
//            words->words_loaded = read;
//            lastLoaded = words->buff;
//            arrPos[words->buff] = (uint64_t) ftell(wrds);
//            wordsUnread[words->buff] -= read;
//            checkOrder(&words, false);
//        } else {
//            checkOrder(&words, true);
//            activeBuffers--;
//            if (activeBuffers <= 0) {
//                free(BuffWordsBlock);
//                free(WentryBlock);
//                free(fname);
//                free(temp.w.b);
//                free(rev_temp.w.b);
//                fclose(wDic);
//                fclose(pDic);
//                fclose(wrds);
//                // Delete intermediate files
//                if (removeIntermediataFiles) {
//                    strcpy(fname, av[2]);
//                    remove(strcat(fname, ".bindx"));
//                    strcpy(fname, av[2]);
//                    remove(strcat(fname, ".wrds"));
//                }
//
//                return 0;
//            }
//        }
//    }
    //END OF DEAD CODE

    // Write final dictionary file
    while (activeBuffers > 0) {
        int cmp = wordcmp(words->word[words->index].w.b, temp.w.b, BYTES_IN_WORD * 4);
        // Store word in buffer
        writeWord(&words->word[words->index], wDic, pDic,
                  cmp == 0 ? true : false, &reps);
        if (cmp != 0)
            storeWord(&temp, words->word[words->index]); // Update last word written
        words->index += 1;
        if (words->index >= words->words_loaded) {
            // Load next word if it's possible
            if (wordsUnread[words->buff] > 0) {
                if (words->buff != lastLoaded) {
                    fseeko(wrds, arrPos[words->buff], SEEK_SET);
                    lastLoaded = words->buff;
                }
                read = loadWord(words->word, wrds, wordsUnread[words->buff]);
                words->index = 0;
                words->words_loaded = read;
                arrPos[words->buff] = (uint64_t) ftello(wrds);
                wordsUnread[words->buff] -= read;
                checkOrder(&words, false);
            } else {
                checkOrder(&words, true);
                activeBuffers--;
            }
        } else checkOrder(&words, false);
    }

    // Write last word repetitions
    fwrite(&reps, sizeof(uint32_t), 1, wDic); // Write num of repetitions

    /////////////////////////// CHECKPOINT ///////////////////////////
    fprintf(stdout, " (Done)\n");
    fprintf(stdout, "\tDict: Closing the program.\n");
    fflush(stdout);
    /////////////////////////// CHECKPOINT ///////////////////////////

    // Deallocate words buffer
    free(BuffWordsBlock);
    free(WentryBlock);

    fclose(wDic);
    fclose(pDic);
    fclose(wrds);

    if (removeIntermediataFiles) {
        strcpy(fname, av[2]);
        remove(strcat(fname, ".bindx"));
        strcpy(fname, av[2]);
        remove(strcat(fname, ".wrds"));
    }

    // Free space
    free(fname);
    free(temp.w.b);
    free(rev_temp.w.b);
    // Everything finished. All it's ok.
    return 0;
}
