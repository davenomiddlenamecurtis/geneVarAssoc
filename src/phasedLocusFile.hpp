#ifndef phasedLocusFileHPP
#define phasedLocusFileHPP

#include "masterLocusFile.hpp"

// read output of phase program

class phasedLocalLocus : localLocus {
public:
	virtual int outputAlleles(allelePair *all,FILE *f,long filePos,int nSubs,int *alleleMap,analysisSpecs const &spec);
	virtual int input(FILE *f, long *locusPosInFile, analysisSpecs const &spec);
	virtual void clear() { ; }
	virtual int read(FILE *fp);
	virtual int write(FILE *fp);
	virtual int typeSpecificCopy(localLocus *src);
	phasedLocalLocus() { clear(); }
	~phasedLocalLocus() { ; }
	virtual locusFileType myType() { return PHASEDHAPSFILE; }
};

class phasedLocusFile : locusFile {
	virtual locusFileType fileType() { return PHASEDHAPSFILE; }
public:
	virtual int readHeaderInfo();
	virtual int outputSubNames(strEntry *subName, analysisSpecs &spec);
	phasedLocusFile() { ; }
};

#endif