#if 0
Copyright 2018 David Curtis

This file is part of the geneVarAssoc package.

geneVarAssoc is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

geneVarAssoc is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with geneVarAssoc.If not, see <http://www.gnu.org/licenses/>.
#endif

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
	const char *name;
	char letter;
	const char *code;
};

class codonReader {
	static char table[4][4][4];
	static int baseLookup[256];
	static const char *nameTable[256];
public:
	char translate(char codon[3],char strand);
	codonReader();
};

char getCompBase(char b);
void getComplementarySequence(char *t,char *s); // space must be already allocated
char *RNA2DNA(char *seq); // conversion occurs in situ
// general utility functions, probably OK not to be in a class
#endif