/*
 * @author Fernando Moreno Jabato <jabato@uma.es>
 * @licence all rights reserved to the author and BitLAB group (University
 *    of Malaga).
 */

// Necessary packages
#include <stdint.h> // write uint64_t ...
#include <inttypes.h>  // uint64_t ...
#include <unistd.h>  // System checks
#include <ctype.h> // Char and string functions
#include <arpa/inet.h> // some IO functions
#include <stdio.h> // IO functions
#include <stdlib.h> // Memory functions
#include <errno.h> // errors
#include <string.h> // String functions
#include <stdbool.h> // Boolean varaibles

//To disable padding and therefore reducing the size of the binary files
//Check if this affects execution time, probably a bit
#pragma pack(push, 1)


// Constants
#define MAXLID 200
#define MAX_READ_LENGTH 5000
#define MAXLS 100000000
#define READBUF 100000
#define PAR_THRESHOLD   10000

// STRUCTS
/* Structure used to store a sequence (k-mer) and its length.
 *  @program META-GECKO dictionaries.
 */
typedef struct {
	/* b is the sequence codified.
	 * Each aminoacid is stored using 2 bits
	 * We have 4 letters per byte and 
	 * each unsigned char size is 1 byte.
     */
    unsigned char *b;
} Word;

/* This structure is used to handle metagenome locations dictionaries.
 *  @program META-GECKO fragments. Loading metagenome dictionaries.
 */
typedef struct{
    uint32_t seq; // Sequence index
    uint64_t pos; // Position on sequence
    char strand; // Strand on the sequence
} LocationEntry;

/* This structure is used to handle k-mers in a sequence.
 *  @program META-GECKO dictionaries.
 */
typedef struct {
	// Word compressed in binary format
    Word w;
    //Location
    LocationEntry loc;
} wentry;


/* This structure is the base of a linked list of wentry structures and
 * information about what buffer stores the wentry instance.
 *  @program META-GECKO dictioanries. Sorting multiple buffers of wentry.
 */
struct item_W{
    wentry *word; // Words array
    uint64_t buff; // Buffer index
    uint64_t index; // Current word index in array 
    uint64_t words_loaded; // Number of words loaded on array
    struct item_W *next; // Next link [linked list]
};


/* This structure is used as node of a wentry linked list.
 *  @program META-GECKO dictioanries. Sorting multiple buffers of wentry.
 */
typedef struct item_W node_W;


/* This structure is used to handle metagenomes and genomes dictioanries.
 *  @program META-GECKO fragments. Loading dictionaries.
 */
typedef struct{
    uint16_t WB; // Word bytes
    unsigned char *seq; // Must be an unsigned char [WB]
    uint64_t pos; // Position on locations file
    uint32_t reps; // Number of instances of this word
} HashEntry;

/* This structure is used to handle hits (seeds).
 *  @program META-GECKO fragments. Generating seeds.
 */
typedef struct {
    // Diagonal where the hit is located
    // This value is calculated as:
    // posX - posY
    int64_t diag;
    // Ocurrence position in sequence X
    uint64_t posX;
    // Ocurrence position in sequence Y
    uint64_t posY;
    // For multiple sequence files this var
    // reflects in what sequence of X file
    // occurs the word
    uint32_t seqX;
    // For multiple sequence files this var
    // reflects in what sequence of Y file
    // occurs the word
    uint32_t seqY;
    // Length of the hit
    uint64_t length;
    // Strand on sequence X
    char strandX;
    // Strand on sequence Y
    char strandY;
} Hit;


/* Structure used to store a sequence (k-mer).
 *  @program META-GECKO fragments. Load genome dictinaries.
 *  @source structs.h (Gecko).
 */
typedef struct {
    //Each letter is stored using 2 bits
    //We have 4 letters per byte and a
    //maximum of 32 in 'b'
    unsigned char b[8];
} word;


/* This structure is used to handle k-mers in a sequence.
 *  @program META-GECKO fragments. Load genome dictionaries.
 *  @source structs.h (Gecko).
 */
typedef struct {
    //Word compressed in binary format
    word w;
    //Ocurrence position in the locations dictionary
    uint64_t pos;
    //Number of locations stored in the
    //positions file
    uint64_t num;
} hashentry;


/* This structure is used to handle genome locations dictionaries.
 *  @program META-GECKO fragments. Loading genome dictionaries.
 *  @source structs.h (Gecko).
 */
typedef struct {
    //Ocurrence position in the sequence
    uint64_t pos;
    //For multiple sequence files this var
    //reflects in what sequence occurs the
    //word
    uint64_t seq;
} location;


/* Structure used to handle fragments.
 *  @program META-GECKO fragments.
 *  @source structs.h
 */
typedef struct{
    //Diagonal where the frag is located
    //This value is calculated as:
    //posX - posY
    int64_t diag;
    //Start position in sequence X
    uint64_t xStart;
    //Start position in Sequence Y
    uint64_t yStart;
    //End position in Sequence X
    uint64_t xEnd;
    //End position in Sequence Y
    uint64_t yEnd;
    //Fragment Length
    //For ungaped aligment is:
    //xEnd-xStart+1
    uint64_t length;
    //Number of identities in the
    //fragment
    uint64_t ident;
    //Score of the fragment. This
    //depends on the score matrix
    //used
    uint64_t score;
    //Percentage of similarity. This
    //is calculated as score/scoreMax
    //Where score max is the maximum
    //score possible
    float similarity;
    //sequence number in the 'X' file
    uint64_t seqX;
    //sequence number in the 'Y' file
    uint64_t seqY;
    //synteny block id
    int64_t block;
    //Strand of the seuqnece Y.'f' for the forward strain and 'r' for the reverse
    char strand;
} FragFile;


/* This structure is the base of a linked list of hit structures and
 * information about what buffer stores the hit instance.
 *  @program META-GECKO fragments. Sorting multiple buffers of hits.
 */
struct item_H{
    Hit *hits; // Array of hits
    uint64_t buff; // Buffer index
    uint64_t index; // Current hit on hits array
    uint64_t hits_loaded; // Number of hits loaded on hits array
    struct item_H *next; // Link to next node [Linked list]
};


/* This structure is used as node of a hit linked list.
 *  @program META-GECKO fragments. Sorting multiple buffers of hit.
 */
typedef struct item_H node_H; 


/* This structure is used to store long (genome) sequences.
 *   @program META-GECKO fragments. Load genome.
 *   @source structs.h (Gecko).
 */
typedef struct{
    char ident[MAXLID+1];
    char datos[MAXLS];
} Sequence;


/* This structure is used to store short sequences (Reads).
 *  @program META-GECKO fragments. Load metagenome.
 */
struct ShortSequence{
    char * sequence; // Sequence
    uint32_t seqIndex; // Index of the sequence
    uint32_t length; // Sequence length
    uint64_t Lac;   //Accumulated length
    struct ShortSequence *next; // Used to create a linked list
};


/* This structure is used as node of a ShortSequence linked list.
 *  @program META-GECKO fragments.
 */
typedef struct ShortSequence Reads; // Linked list of short sequences


/* This structure is the base of a linked list of read statistics.
 *  @program META-GECKO fixeReverseFragments.
 */
struct rStats {
    uint64_t pos;         // Start position of read
    uint64_t Lac;         // Accumulated length == offset
    uint64_t length;      // Read length  
    uint64_t index;       // Read index   
    struct rStats *next;  // Next link (linked list)
};


/* This structure is used as node of a rStats linked list.
 *  @program META-GECKO fixeReverseFragments.
 */
typedef struct rStats Read;


/* This structure is the base of a linked list of FragFiles.
 *  @program META-GECKO fixeReverseFragments.
 */
struct fitem{
    FragFile *frags;
    uint64_t buff;
    uint64_t index;
    uint64_t frags_loaded;
    struct fitem *next;
};


/* This structure is used as node of a FragFile linked list.
 *  @program META-GECKO fixeReverseFragments.
 */
typedef struct fitem fnode;


/// Functions


/*
	Print the error message 's' and exit(-1)
*/
void terror(char *s);

/*
	Manually-controlled buffering to avoid getc's
*/

char buffered_fgetc(char *buffer, uint64_t *pos, uint64_t *read, FILE *f);

/*
	Print an n-bits word encoded using 2 bits per letter
*/

void showWord(Word *w, char *ws, uint16_t WORD_LENGTH);

/* 
	This function compares two arrays of unsigned chars with the same length.
 */
int wordcmp(unsigned char *w1, unsigned char *w2, int n);

/* 
	This function is used to check if a file exists or not.
 */
int exists(char *file);


/* 
	This function is used to check if a string given is an integer.
 */
int is_int(char const *str);


/* 
	This function is used to check if a string given is a float.
 */
int is_float(char const *str);

/* 
	Compute size of HashEntry
 */
uint16_t size_of_HashEntry(const uint16_t BYTES_IN_WORD);



