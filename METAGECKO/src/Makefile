CC=gcc
CFLAGS=-O3 -D_FILE_OFFSET_BITS=64 -Wall
BIN=../bin

all: gecko meta_gecko

gecko: words sortWords w2hd reverseComplement

meta_gecko: dict frags fixeReverse combineFragments reverseComplement

words: words.c
	$(CC) $(CFLAGS) words.c dictionaryFunctions.c commonFunctions.c -o $(BIN)/words

sortWords: sortWords.c
	$(CC) $(CFLAGS) -DBaseType=wentry sortWords.c dictionaryFunctions.c commonFunctions.c quicksort.c -lpthread -o $(BIN)/sortWords

w2hd: w2hd.c
	$(CC) $(CFLAGS) w2hd.c dictionaryFunctions.c commonFunctions.c -o $(BIN)/w2hd

reverseComplement: reverseComplement.c
	$(CC) $(CFLAGS) reverseComplement.c commonFunctions.c -o $(BIN)/reverseComplement

dict: dict.c
	$(CC) $(CFLAGS) dict.c dict.h dictFun.c -o $(BIN)/dict

frags: frags.c
	$(CC) $(CFLAGS) frags.c frags.h fragsFun.c -o $(BIN)/frag

fixeReverse: fixeReverseFragments.c
	$(CC) $(CFLAGS) fixeReverseFragments.c -o $(BIN)/fixeReverseFrags

combineFragments: combineFragments.c
	$(CC) $(CFLAGS) combineFragments.c -o $(BIN)/combineFragments