#ifndef GENEVARUTILSHPP
#define GENEVARUTILSHPP 1
#include "getGene.hpp"
#include "masterLocusFile.hpp"

#define MAXVCFPERCC 10
#define MAXDEPTH 5

class gvaParams {
public:
	char geneListFn[100],baitFn[100],ccFn[2][MAXVCFPERCC][100],referencePath[100],geneName[100],sequencePath[100],posName[100];
	char intervalListFn[100], testName[100]; // for intVarAssoc
	int useFreqs[2],nSubs[2],nCc[2],writeComments,writeScoreFile;
	int dontExtractGene,keepTempFiles,doNotRun;
	int input(FILE *fp,analysisSpecs &spec);
	int readParms(int argc,char *argv[],analysisSpecs &spec);
	int getNextArg(char *nextArg,int argc,char *argv[],FILE *fp[MAXDEPTH],int *depth,int *argNum);
	int upstream,downstream,margin;
	int firstGeneNum,lastGeneNum;
	int wFunc;
	float wf;
};

#define hereOK() fprintf(stderr,"Got to line %d in %s OK\n",__LINE__,__FILE__)
#endif