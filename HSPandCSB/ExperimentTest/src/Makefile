CC=gcc
BIN=../bin
FLAGS=-Wall -pg

all: dictionary fragment write-words write-frags

develope: dictionary-develope fragment-develope write-words write-frags

dictionary: dictionary.c
	$(CC) dictionary.c dictionary.h geckoMGVutils.h dictionaryFun.c geckoMGVutils.c -o $(BIN)/dic

fragment: frags.c
	$(CC) frags.c frags.h geckoMGVutils.h fragsFun.c geckoMGVutils.c -o $(BIN)/frags

dictionary-develope: dictionary.c
	$(CC) $(FLAGS) dictionary.c dictionary.h geckoMGVutils.h dictionaryFun.c geckoMGVutils.c -o $(BIN)/dic

fragment-develope: frags.c
	$(CC) $(FLAGS) frags.c frags.h geckoMGVutils.h fragsFun.c geckoMGVutils.c -o $(BIN)/frags	

write-words: writeWords.c
	$(CC) writeWords.c -o $(BIN)/WW

write-frags: fragB_T.c
	$(CC) fragB_T.c -o $(BIN)/WF