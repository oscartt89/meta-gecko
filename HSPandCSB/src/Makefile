CC=gcc
CXX=g++
CFLAGS=-O3 -D_FILE_OFFSET_BITS=64 -Wall -Wextra
BIN=../bin

all:  fuzzyGenomeAbundance-mapping bestGenomeInRead-mapping parseGFF parseGeneBankTab fuzzyGenomeAbundance genomeAbundance extractSeq BorderDetection getDuplicatedCSB fromTxt2MAT getMalla getCoverage4HSP improveCSBwithDuplicates reverseComplement words sortWords w2hd wordsStat hits sortHits filterHits FragHits combineFrags hitsStat fragStat sortFragsBySeqX fragv2tov3 fragv3tov2 filtroDup getCSB getFreqFasta kar2test getInfoCSB parseGeneFeature getCoverage4LocusTaginHSP mergeMetagenomeFrags bestGenomeInRead duplicatedMetaFrags getHistogramByLength mgFoIview2TXT-R2G mgReadsIndex
#borderDetectionOTS
png: af2pngrev aMf2pngrev

hits: hits.c
	$(CC) $(CFLAGS) hits.c comparisonFunctions.c dictionaryFunctions.c commonFunctions.c -o $(BIN)/hits

sortHits: sortHits.c
	$(CC) $(CFLAGS) -DBaseType=hit sortHits.c comparisonFunctions.c commonFunctions.c quicksort.c -lpthread -o $(BIN)/sortHits

filterHits: filterHits.c
	$(CC) $(CFLAGS) filterHits.c comparisonFunctions.c commonFunctions.c -o $(BIN)/filterHits

FragHits: FragHits.c
	$(CC) $(CFLAGS) FragHits.c comparisonFunctions.c commonFunctions.c -o $(BIN)/FragHits
	
combineFrags: combineFrags.c
	$(CC) $(CFLAGS) combineFrags.c comparisonFunctions.c commonFunctions.c -o $(BIN)/combineFrags
	
hitsStat: hitsStat.c
	$(CC) $(CFLAGS) hitsStat.c commonFunctions.c -o $(BIN)/hitsStat
	
fragStat: fragStat.c
	$(CC) $(CFLAGS) fragStat.c comparisonFunctions.c commonFunctions.c -o $(BIN)/fragStat
	
af2pngrev: af2pngrev.cpp
	$(CXX) $(CFLAGS) -Wno-write-strings af2pngrev.cpp commonFunctions.c comparisonFunctions.c -I/usr/include/freetype2  -L/usr/local/lib -lpngwriter -lz -lm -lfreetype -o $(BIN)/af2pngrev -lpng 

aMf2pngrev: aMf2pngrev.cpp
	$(CXX) $(CFLAGS) -Wno-write-strings aMf2pngrev.cpp commonFunctions.c comparisonFunctions.c -I/usr/include/freetype2  -L/usr/local/lib -lpngwriter -lz -lm -lfreetype -o $(BIN)/aMf2pngrev -lpng 

reverseComplement: reverseComplement.c
	$(CC) $(CFLAGS) reverseComplement.c commonFunctions.c -o $(BIN)/reverseComplement

words: words.c
	$(CC) $(CFLAGS) words.c dictionaryFunctions.c commonFunctions.c -o $(BIN)/words

sortWords: sortWords.c
	$(CC) $(CFLAGS) -DBaseType=wentry sortWords.c dictionaryFunctions.c commonFunctions.c quicksort.c -lpthread -o $(BIN)/sortWords

w2hd: w2hd.c
	$(CC) $(CFLAGS) w2hd.c dictionaryFunctions.c commonFunctions.c -o $(BIN)/w2hd

wordsStat: wordsStat.c
	$(CC) $(CFLAGS) wordsStat.c dictionaryFunctions.c commonFunctions.c -o $(BIN)/wordsStat

sortFragsBySeqX: sortFragsBySeqX.c
	$(CC) $(CFLAGS) -D'BaseType=struct FragFile' sortFragsBySeqX.c commonFunctions.c quicksort.c comparisonFunctions.c -lpthread -o $(BIN)/sortFragsBySeqX

newFragToBalazsVersion:
	$(CC) $(CFLAGS) commonFunctions.c comparisonFunctions.c newFragToBalazsVersion.c -o $(BIN)/newFragToBalazsVersion

kar2test:
	$(CC) $(CFLAGS) kar2test.c -lm -o $(BIN)/kar2test

fragv2tov3:
	$(CC) $(CFLAGS) fragv2tov3.c fragmentv3.c fragmentv2.c  -lm commonFunctions.c comparisonFunctions.c -o $(BIN)/fragv2tov3

fragv3tov2:
	$(CC) $(CFLAGS) fragv3tov2.c fragmentv3.c fragmentv2.c -lm commonFunctions.c comparisonFunctions.c -o $(BIN)/fragv3tov2

filtroDup:
	$(CC) $(CFLAGS) newFiltrarDuplicados.c fragmentv3.c -g -lm -o $(BIN)/filtroDup
getMalla:
	$(CC) $(CFLAGS) getMalla.c fragmentv3.c -g -lm -o $(BIN)/getMalla

getCSB:
	$(CC) $(CFLAGS) getCSB.c -lm postprocessingFunctions.c lista.c fragmentv3.c -o $(BIN)/getCSB
getDuplicatedCSB:
	$(CC) $(CFLAGS) getDuplicatedCSB.c -lm postprocessingFunctions.c lista.c fragmentv3.c -o $(BIN)/getDuplicatedCSB
getFreqFasta:
	$(CC) $(CFLAGS) getFreqFasta.c -lm -o $(BIN)/getFreqFasta

getInfoCSB:
	$(CC) $(CFLAGS) getInfoCSB.c fragmentv2.c commonFunctions.c comparisonFunctions.c -lm -o $(BIN)/getInfo

parseGeneFeature:
	$(CC) $(CFLAGS) parseGeneFeatureTab.c gene.c -lm -o $(BIN)/parseGeneFeatureTab
parseGFF:
	$(CC) $(CFLAGS) parseGFF.c gene.c -lm -o $(BIN)/parseGFF
parseGeneBankTab:
	$(CC) $(CFLAGS) parseGeneBankTab.c gene.c -lm -o $(BIN)/parseGeneBankTab
	
getCoverage4LocusTaginHSP:
	$(CC) $(CFLAGS) getCoverage4LocusTaginHSP.c gene.c fragmentv2.c commonFunctions.c comparisonFunctions.c -o $(BIN)/getCoverage4LocusTaginHSP

mergeMetagenomeFrags:
	$(CC) $(CFLAGS) mergeMetagenomeFrags.c commonFunctions.c comparisonFunctions.c -lm -o $(BIN)/mergeMetagenomeFrags

bestGenomeInRead:
	$(CC) $(CFLAGS) bestGenomeInRead.c commonFunctions.c comparisonFunctions.c fragmentv2.c -lm -o $(BIN)/bestGenomeInRead
bestGenomeInRead-mapping:
	$(CC) $(CFLAGS) bestGenomeInRead-mapping.c commonFunctions.c comparisonFunctions.c fragmentv2.c -lm -o $(BIN)/bestGenomeInRead-mapping
genomeAbundance:
	$(CC) $(CFLAGS) genomeAbundance.c commonFunctions.c comparisonFunctions.c fragmentv2.c -lm -o $(BIN)/genomeAbundance
fuzzyGenomeAbundance:
	$(CC) $(CFLAGS) fuzzyGenomeAbundance.c commonFunctions.c comparisonFunctions.c fragmentv2.c -lm -o $(BIN)/fuzzyGenomeAbundance
fuzzyGenomeAbundance-mapping:
	$(CC) $(CFLAGS) fuzzyGenomeAbundance-mapping.c commonFunctions.c comparisonFunctions.c fragmentv2.c -lm -o $(BIN)/fuzzyGenomeAbundance-mapping
duplicatedMetaFrags:
	$(CC) $(CFLAGS) duplicatedMetaFrags.c commonFunctions.c comparisonFunctions.c fragmentv2.c -lm -o $(BIN)/duplicatedMetaFrags

getHistogramByLength:
	$(CC) $(CFLAGS) getHistogramByLength.c fragmentv2.c  -lm commonFunctions.c comparisonFunctions.c -o $(BIN)/getHistogramByLength
getCoverage4HSP:
	$(CC) $(CFLAGS) getCoverage4HSP.c gene.c fragmentv2.c commonFunctions.c comparisonFunctions.c -o $(BIN)/getCoverage4HSP
BorderDetection:
	$(CC) $(CFLAGS) BorderDetection.c fragmentv3.c commonFunctions.c JAMfunctions.c -lm -o $(BIN)/BorderDetection
borderDetectionOTS:
	$(CC) $(CFLAGS) BorderDetection_OTS.c commonFunctions.c comparisonFunctions.c JAMfunctions.c fragmentv3.c -lm -o $(BIN)/BorderDetectionOTS
extractSeq:
	$(CC) $(CFLAGS) extractSeq.c -lm -o $(BIN)/extractSeq
fromTxt2MAT: 
	$(CC) $(CFLAGS) fromTxt2MAT.c -lm -o $(BIN)/fromTxt2MAT
improveCSBwithDuplicates: 
	$(CC) $(CFLAGS) improveCSBwithDuplicates.c JAMfunctions.c -lm postprocessingFunctions.c lista.c fragmentv3.c -o $(BIN)/improveCSBwithDuplicates

eevolutionary:
	$(CC) $(CFLAGS)	evolutionaryEvents.c lista.c fragmentv3.c -g -lm -o eevolutionary

mgFoIview2TXT-R2G:
	$(CC) $(CFLAGS) -DSTRATEGY=lessUsed mgFoIview2TXT-R2G.c comparisonFunctions.c commonFunctions.c fragmentv2.c -o $(BIN)/mgFoIview2TXT-R2G

mgReadsIndex:
	$(CC) $(CFLAGS) mgReadsIndex.c -o $(BIN)/mgReadsIndex

	#$(BIN)/borderDetectionOTS
clean:
	rm -rf  $(BIN)/fuzzyGenomeAbundance-mapping $(BIN)/bestGenomeInRead-mapping $(BIN)/parseGFF $(BIN)/parseGeneBankTab $(BIN)/fuzzyGenomeAbundance $(BIN)/genomeAbundance $(BIN)/extractSeq $(BIN)/BorderDetection $(BIN)/getDuplicatedCSB $(BIN)/getMalla $(BIN)/fromTxt2MAT $(BIN)/getCoverage4HSP $(BIN)/improveCSBwithDuplicates $(BIN)/fragv2tov3 $(BIN)/fragv3tov2 $(BIN)/filtroDup $(BIN)/getCSB $(BIN)/getFreqFasta $(BIN)/kar2test $(BIN)/getInfo $(BIN)/parseGeneFeatureTab $(BIN)/getCoverage4LocusTaginHSP $(BIN)/sortFragsBySeqX $(BIN)/mergeMetagenomeFrags $(BIN)/bestGenomeInRead $(BIN)/duplicatedMetaFrags $(BIN)/getHistogramByLength $(BIN)/mgFoIview2TXT-R2G $(BIN)/mgReadsIndex
	rm -rf $(BIN)/reverseComplement $(BIN)/words $(BIN)/sortWords $(BIN)/w2hd $(BIN)/wordsStat
	rm -rf $(BIN)/hits $(BIN)/sortHits $(BIN)/filterHits $(BIN)/FragHits $(BIN)/combineFrags $(BIN)/hitsStat $(BIN)/fragStat $(BIN)/af2pngrev $(BIN)/aMf2pngrev
