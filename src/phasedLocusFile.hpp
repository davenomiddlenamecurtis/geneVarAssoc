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