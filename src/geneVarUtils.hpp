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

#ifndef GENEVARUTILSHPP
#define GENEVARUTILSHPP 1
#include "getGene.hpp"
#include "masterLocusFile.hpp"

#define MAXVCFPERCC MAXVCFFILES
#define MAXDEPTH 5
#define	MAXFILENAMELENGTH 100

class gvaParams {
public:
	char geneListFn[MAXFILENAMELENGTH],baitFn[MAXFILENAMELENGTH],ccFn[2][MAXVCFPERCC][MAXFILENAMELENGTH],referencePath[MAXFILENAMELENGTH],
		geneName[100],sequencePath[100],posName[100], refAll[100], altAll[100];
	char bedFileFn[MAXFILENAMELENGTH], famFileFn[MAXFILENAMELENGTH], bimFileFn[MAXFILENAMELENGTH];
	char intervalListFn[MAXFILENAMELENGTH], testName[100]; // for intVarAssoc
	int useFreqs[2],nSubs[2],nCc[2],writeComments,writeScoreFile;
	int dontExtractVariants,keepTempFiles,doNotRun,onlyExtractVariants;
	int input(FILE *fp,analysisSpecs &spec);
	int readParms(int argc,char *argv[],analysisSpecs &spec);
	int getNextArg(char *nextArg,int argc,char *argv[],FILE *fp[MAXDEPTH],int *depth,int *argNum);
	int upstream,downstream,margin;
	int firstGeneNum,lastGeneNum;
};

int tbiExtractIntervals(char* tbiFilename, char* outFn, int appendToOld, int addChrInVCF, int removeSpaces, intervalList& iList);
int plinkExtractIntervals(char* bedFilename, char* famFilename, char* bimFilename, char* outFn, intervalList& iList, char* geneName);

#define MISSINGPHENOTYPE -999 // value used internally to indicate subject without a defined phenotype
#ifndef hereOK
#define hereOK() fprintf(stderr,"Got to line %d in %s OK\n",__LINE__,__FILE__)
#endif
#endif