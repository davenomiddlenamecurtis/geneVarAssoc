#ifndef GETSEQUENCEHPP
#define GETSEQUENCEHPP
#include <stdio.h>


class faSequenceFile {
	FILE *fp;
	long startPos;
	int lineLength,gapLength;
public:
	faSequenceFile();
	~faSequenceFile() { fp && fclose(fp); }
	int inited() { return fp!=0; }
	int init(char *fn);
	int getSequence(char *s,long ePos,long len);
};

class aminoAcid {
public:
	char *name;
	char letter;
	char *code;
};

class codonReader {
	static char table[4][4][4];
	static int baseLookup[256];
	static char *nameTable[256];
public:
	char translate(char codon[3],char strand);
	codonReader();
};

char getCompBase(char b);
void getComplementarySequence(char *t,char *s); // space must be already allocated
char *RNA2DNA(char *seq); // conversion occurs in situ
// general utility functions, probably OK not to be in a class
#endif