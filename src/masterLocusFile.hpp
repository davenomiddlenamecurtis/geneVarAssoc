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

#ifndef masterLocusFileHPP
#define masterLocusFileHPP

#define MAXALL 10 // maximum number of alleles occurring at any locus
#define MAXALLLENGTH 200
// #define MAXALLLENGTH 300 // maximum length of the character string describing each REF allele 
// or all ALT alleles with commas separating them
// 100 was too short
// in theory one can have very long alleles in VCF files to describe chromosomal rearrangements
// the code for scanWord means that only a truncated version will be input
// but at present truncation causes an error
// 30/12/15 I am increasing this to 1000 to see if this fixes a problem I have
// when I do this I get a segmentation fault

#ifndef MAXVCFFILES
#define MAXVCFFILES 10
#endif
#define MAXFILENAMELENGTH 100
#define BUFFSIZE 1000000 // maximum length of e.g. line read in from VCF file was 900000
#define VCFFIELDLENGTH 2000 // maximum length for e.g. quality, format fields
// even 3000 was too short for a UK10K gene with many PolyPhen entries
#define MAXSCOREASSOCARGS 50
#define MAXEXPRESSIONLENGTH 2000

#include "dcindex.hpp"
#include "dcerror.hpp"
#include <stdio.h>
#include <string.h>
#include "getGene.hpp"
#include <string>
#include <map>
#include <list>
typedef std::pair<std::string,int> TStrIntPair;
typedef std::map<std::string,int> TStrIntMap;


#ifndef FILEPOSITION
#ifdef WIN32
#define FILEPOSITION __int64
#define FTELL(f) _ftelli64(f)
#define FSEEK(f,p,o) _fseeki64(f,p,o)
#else
#define FILEPOSITION long
#define FTELL(f) ftell(f)
#define FSEEK(f,p,o) fseek(f,p,o)
// for now
#endif
#endif

#if 0
Set things up so we can have different kinds of datafile.
masterLocus is general information about variant.
localLocus relates to one particular locusFile.
masterLocusFile has information about a number of locusFiles, possibly different types.
#endif

typedef char filenamestring[MAXFILENAMELENGTH];
typedef int allelePair[2];
typedef float probTriple[3];
#define MAXSTR 100
typedef char strEntry[MAXSTR+1];
#define DEFAULTNUMVCFFIELDSTOSKIP 9 // maybe one day allow this to vary

enum locusSNP { SNP_NO=0,SNP_MAYBE,SNP_YES } ;

class masterLocusFile;
class localLocus;

typedef int alleleMap[MAXALL];

class analysisSpecs {
public:
	analysisSpecs() 
	{ 
		unknownIfUntyped=1; // if there are no calls for a variant in the VCF file assume it has not been covered rather than all wildtype
		unknownIfNoPass=0; altIsCommon=0; 
		wildIfUnknown=0;
		onlycc01=1; // do not use subject unless cc status is 0 or 1
		sc=0; sp=0L; ec=25; ep=0L;
		skipIfNoPass=1;
		proportionCalledToPass=0.95;
		useConsequenceWeights=0;
		consequenceThreshold=NULL_CONSEQUENCE;
		useEnsembl=useEnsembl=willNeedEnsemblConsequence=willNeedInbuiltConsequence=0;
		GQThreshold=-1;
		hetDevThreshold=hetDevThresholdSq=-1;
		depthThreshold = -1;
		ABThreshold = -1;
		onlyUseSNPs=0;
		doRecessiveTest=0;
		weightThreshold=0;
		LDThreshold=1;
		phenotypes=NULL;
		showHapLocusNames=useHaplotypes=0;
		useProbs=0;
		count_hom_as_het=0;
		useTrios=0;
		dontMergeAlleles=ignoreAlleles=0;
		mergeAltAlleles = 1;
		omitIntrons = 0;
		spliceRegionSize = 100;
		debug=0;
		*alleleFreqStr=*alleleNumberStr=*alleleCountStr='\0';
		*commentExpression=*weightExpression='\0';
		wf=10;
		numVcfFieldsToSkip = DEFAULTNUMVCFFIELDSTOSKIP;
		removeVcfSpaces = 0;
	} 
int onlycc01,unknownIfUntyped,unknownIfNoPass,altIsCommon,sc,ec,skipIfNoPass,useConsequenceWeights,onlyUseSNPs,nExc,doRecessiveTest,addChrInVCF[MAXVCFFILES],useHaplotypes, showHapLocusNames,count_hom_as_het,useTrios,
ignoreAlleles,dontMergeAlleles,useProbs,wildIfUnknown,debug,omitIntrons,spliceRegionSize;
int useEnsembl,willNeedEnsemblConsequence,willNeedInbuiltConsequence,mergeAltAlleles,numVcfFieldsToSkip,removeVcfSpaces;
int *phenotypes;
TStrIntMap subPhenos;
long sp,ep;
float GQThreshold,proportionCalledToPass,hetDevThreshold,hetDevThresholdSq,depthThreshold,ABThreshold;
float weightThreshold,LDThreshold,wf;
float consequenceThreshold;
char exclusionStr[20][200];
char triosFn[200];
char alleleFreqStr[100],alleleNumberStr[100],alleleCountStr[100],
scoreassocArgs[MAXSCOREASSOCARGS][2][100];
int nScoreassocArgs;
char weightExpression[MAXEXPRESSIONLENGTH];
char commentExpression[MAXEXPRESSIONLENGTH];
char vepCommand[MAXEXPRESSIONLENGTH];
std::list<std::string> excludeExpressions;
};

#if 0
locusFile will be base class and derived classes will be unknown sizes
masterLocus will need to keep pointers to localLocus and use these sizes when reading and writing them
may not need to allocate memory each time localLocus is read and written
may be that read() and write() can assume that myLocalLocus[i] has constant size
#endif

typedef localLocus *LOCALLOCUSPTR;

enum locusFileType { NOFILETYPE,VCFFILE, SHAPEITHAPSFILE, IMPUTEDOSAGEFILE, PHASEDHAPSFILE };
#if 0
A masterLocus contains information about a particular variant. 
A localLocus contains information about that variant as it is contained in a particular locusDataFile.
There may be different kinds of locusDataFile containing information about the same variant.
However, different subsets of variants can be contained in different locusDataFiles.
If a locusDataFile does not contain information about a variant then the localLocus::filter field will contain "UNTYPED"
The masterLocusFile keeps track of which files are used.
#endif

class masterLocus {
protected:
	friend class masterLocusFile;
	int chr;
	long pos;
	locusSNP SNP;
	int nLocusFiles;
	consequenceType worstConsequenceType[MAXALL];
	char masterID[VCFFIELDLENGTH];
	char ref[MAXALLLENGTH+1];
	char alt[MAXALLLENGTH*MAXALL];
	int nAlls;
	char alls[MAXALL][MAXALLLENGTH+1];
	alleleMap *alleleMapping; // maps how allele in local file maps to global allele list for that locus
	char ensemblConsequence[MAXALL][50],quickConsequence[MAXALL][50];
	int genoCount[3],genoCcCount[2][3]; // may know this sometimes
	FILEPOSITION *locusPosInFile;
	LOCALLOCUSPTR *myLocalLocus; // keep a record of each one read in separately
	int read(FILE *fp);
	int write(FILE *fp);
	int writePredictorQuery(FILE *fp, analysisSpecs const &spec);
	int readQueryOutput(FILE *fp,analysisSpecs const &spec);
public:
	masterLocus(int nLF);
	~masterLocus();
	locusSNP isSNP() { return SNP; }
	int outputProbs(probTriple *prob,FILE *f,int whichFile,int nSubs,analysisSpecs const &spec);
	int outputAlleles(allelePair *all,FILE *f,int whichFile,int nSubs,analysisSpecs const &spec);
	int outputCalls(strEntry *call,FILE *f,int whichFile,int nSubs,analysisSpecs const &spec);
	int outputVcfGenotypes(FILE *fo,FILE *f,int whichFile,int nSubs,analysisSpecs const &spec);
	int print(FILE *fp);
	int printFeatures(FILE* fp,int showFreq=-1);
	const char *getID(); // get a human readable ID if one of myLocalLocus has one
//	int writeRiskVarInfo(char *s,int withFreqs=0);
//	bool readRiskVarInfo(char *s,int withFreqs=0);
	consequenceType getWorstConsequenceType(int a) { return worstConsequenceType[a]; }
	int getQuickFeature(refseqGeneInfo &r,int a);
	const char *reportQuickConsequence(int a) { return quickConsequence[a]; }
	// this is text of worst consequence after call to getQuickFeature()
	const char *reportEnsemblConsequence(int a) { return ensemblConsequence[a]; }
	int getGenoCount(int g) { return genoCount[g]; }
	void setLocusPosInFile(int f,FILEPOSITION l) { locusPosInFile[f] = l; }
	FILEPOSITION getLocusPosInFile(int f) { return locusPosInFile[f]; }
	int getChr() { return chr; }
	long getPos() { return pos; }
	const char *getAll(int a) { return alls[a]; }
};

class localLocus {
	friend class masterLocus;
	friend class masterLocusFile;
protected:
	int chr;
	long pos;
	char id[VCFFIELDLENGTH];
	char ref[MAXALLLENGTH];
	char alt[MAXALLLENGTH*MAXALLLENGTH];
	char alls[MAXALL][MAXALLLENGTH];
	float alleleFreq[MAXALL];
	int nAltAlls; // number of alt alleles for which frequency is given
	locusSNP SNP;
	char filter[VCFFIELDLENGTH];
	float AF,AC,AN;
	virtual int input(FILE *f, FILEPOSITION *locusPosInFile, analysisSpecs const &spec)=0; // not implemented, only for derived classes
	// this function is to read in information for the next valid locus, possibly with criterion that it must PASS
	virtual void clear();
	locusSNP isSNP();
	virtual int outputAlleles(allelePair *all,FILE *f,FILEPOSITION filePos,int nSubs,int *alleleMap,analysisSpecs const &spec)=0;
	virtual int outputProbs(probTriple *prob,FILE *f,FILEPOSITION filePos,int nSubs,int *alleleMap,analysisSpecs const &spec)=0;
public:
	virtual int outputCalls(strEntry *call, FILE *f, FILEPOSITION filePos,int nSubs, int *alleleMap, analysisSpecs const &spec) { return 0; }
	virtual int read(FILE *fp);
	virtual int write(FILE *fp);
	virtual locusFileType myType()=0;
	virtual int typeSpecificCopy(localLocus *src);
	localLocus();
	virtual ~localLocus() { ; }
	virtual const char *getInfo() { return ""; }
};

class locusFile {
protected:
	int nSubs;
	friend class masterLocusFile;
	FILE *fp;
	virtual locusFileType fileType()=0;
	virtual int outputSubNames(strEntry *subName, analysisSpecs &spec) { return 0; } // may not always be implemented
	virtual int readHeaderInfo()=0;
	// this function is guaranteed to do two things:
	// (1) count the number of subjects in this locusFile
	// (2) position fp so that it is pointing at the start of the information for the first locus
public:
	int getNSUbs() { return nSubs; }
	static char *buff; // anybody can use this
	locusFile() { if (buff==0) buff=new char[BUFFSIZE]; fp = 0; } // never gets freed
	virtual ~locusFile() { if (fp) fclose(fp);  }
};


typedef locusFile *LOCUSFILEPTR;
#define MAXLOCIINSCOREASSOCFILE 25000

class masterLocusFile  {
	dcIndex index;
	FILE *recordFile;
	masterLocus tempRecord;
	localLocus *tempLocus;
	int nLocusFiles;
	int currentLocusFile;
	filenamestring *lfFileNames;
	int *nSubs, *cc, *nAlls;
	LOCUSFILEPTR *locusFiles;
	locusFileType *fileTypes;
	int *holdsFreqs; // whether file provides information about allele frequencies
	int addLocus(FILE *f,analysisSpecs const &spec);
	FILEPOSITION currentRecPos;
public:
	void setFileType(int i, locusFileType t) { fileTypes[i] = t; }
	int providesFreqs(int i) { return holdsFreqs[i]; }
	locusFileType getFileType(int i) { return fileTypes[i]; }
	masterLocusFile(int nLF);
	~masterLocusFile();
	int openFiles(char *rfn,char *ifn);
	void closeFiles() { recordFile && fclose(recordFile); recordFile=0; index.close(); };
	int openLocusFiles();
	int closeLocusFiles();
	int setCurrentFile(char *fn);
	int addLocusFile(char *fn,locusFileType t);
	int readLocusFileEntries(char *fn,analysisSpecs const &spec,int aff); // must specify cc status
	int load(masterLocus &rec,FILEPOSITION pos);
	int loadFirst(analysisSpecs &spec);
	int loadNext(analysisSpecs &spec);
	int save(masterLocus &rec,FILEPOSITION pos);
	int fill(masterLocus &rec,localLocus *loc,FILEPOSITION locusPosInFile);
	int merge(masterLocus &rec,localLocus *loc,FILEPOSITION locusPosInFile);
	FILEPOSITION findFirstInRange(analysisSpecs const &spec);
	FILEPOSITION findFirstInRange(int chr,FILEPOSITION pos) { analysisSpecs newSpec; newSpec.sc=chr; newSpec.sp=pos; return findFirstInRange(newSpec); }
	int print(FILE* fp,analysisSpecs const &spec);
	int printFeatures(FILE* fp,analysisSpecs const &spec,int showFreq=-1);
	int outputCurrentAlleles(allelePair *all,analysisSpecs &spec);
	int outputCurrentProbs(probTriple *prob,analysisSpecs &spec);
	int outputAlleles(allelePair **all,analysisSpecs &spec);
	int outputProbs(probTriple **all,analysisSpecs &spec);
	int outputCalls(strEntry **call,analysisSpecs &spec);
	int outputSubNames(strEntry *subName,analysisSpecs &spec);
	int outputAffectionStatus(int *cc,analysisSpecs &spec);
	int outputMergedVCFHeader(FILE *fo);
	int outputMergedVcfGenotypes(FILE *fo,analysisSpecs const &spec);
	int outputAltFrequencies(float *freqs,int cc,analysisSpecs const &spec);
	int outputEurAltFrequencies(float *freqs,int cc,analysisSpecs const &spec);
	int outputSAInfo(int *useLocus,float *locusWeight,analysisSpecs const &spec);
	int getEnsemblConsequences(analysisSpecs const &spec);
	int getQuickConsequences(refseqGeneInfo &r,analysisSpecs const &spec,int redo=0);
	int writeScoreAssocFiles(char *root,float wf,int *useFreqs,int *suppliedNSubs,int writeNames,int writeComments,int writeScorefile,analysisSpecs &spec);
	int writeScoreAssocFiles(masterLocusFile &subFile,char *root,float wf,int *useFreqs,int *suppliedNSubs,int writeNames,int writeComments,int writeScorefile,analysisSpecs &spec);
//	int writeVars(char *fn,int *useFreqs,analysisSpecs &spec);
	int writeGenos(char *fn,int *useFreqs,analysisSpecs &spec);
	int writeGenoCounts(FILE *fo[2],char *geneName,long *varNum,analysisSpecs &spec,allelePair **a);
	int writeAltSubs(char *fn,analysisSpecs &spec);
	int gotoFirstInRange(analysisSpecs &spec);
	int gotoNextInRange(analysisSpecs &spec);
	int countNumberInRange(analysisSpecs &spec);
	int currentChr() { return tempRecord.chr; }
	FILEPOSITION currentPos() { return tempRecord.pos; }
	int currentIsSNP() { return tempRecord.isSNP(); }
	int currentNAlls() { return tempRecord.nAlls; }
	int getNSubs(int i) { return nSubs[i]; }
	int getTotalSubs();
	const char *currentAll(int a) { return tempRecord.alls[a]; }
	const char *currentID() { return tempRecord.masterID; }
	const char *currentQuickConsequence(int a) { return tempRecord.quickConsequence[a]; }
	consequenceType currentWorstConsequenceType(int a) { return tempRecord.worstConsequenceType[a]; }
//	int writeRiskVarInfo(char *s, int withFreqs = 0) { return tempRecord.writeRiskVarInfo(s, withFreqs); }
//	bool readRiskVarInfo(char *s,int withFreqs=0){ return tempRecord.readRiskVarInfo(s, withFreqs); }

};

bool scanWord(char **line,char *word,int maxLength,char token='\0');


#ifdef MSDOS
#define hereOK()
#else
#define hereOK() fprintf(stderr,"Got to line %d in %s OK\n",__LINE__,__FILE__)
#endif
#endif

