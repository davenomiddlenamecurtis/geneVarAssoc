#ifndef hapsLocusFileHPP
#define hapsLocusFileHPP

#include "masterLocusFile.hpp"

class hapsLocalLocus : localLocus {
	char rsName[100];
public:
	virtual int outputProbs(probTriple *all, FILE *f, FILEPOSITION filePos, int nSubs, int *alleleMap, analysisSpecs const &spec) { return 0; } // not implemented, though could be
	virtual int outputAlleles(allelePair *all,FILE *f,FILEPOSITION filePos,int nSubs,int *alleleMap,analysisSpecs const &spec);
	virtual int input(FILE *f, FILEPOSITION *locusPosInFile, analysisSpecs const &spec);
	virtual void clear() { rsName[0] = '\0'; }
	virtual int read(FILE *fp);
	virtual int write(FILE *fp);
	virtual int typeSpecificCopy(localLocus *src);
	hapsLocalLocus() { clear(); }
	~hapsLocalLocus() { ; }
	virtual locusFileType myType() { return SHAPEITHAPSFILE; }
};

class hapsLocusFile : locusFile {
	virtual locusFileType fileType() { return SHAPEITHAPSFILE; }
public:
	virtual int readHeaderInfo();
	virtual int outputSubNames(strEntry *subName, analysisSpecs &spec);
	hapsLocusFile() { ; }
};

#endif