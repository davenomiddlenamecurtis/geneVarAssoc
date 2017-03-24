#include <stdlib.h>
#include <string.h>

#ifndef MSDOS 
#include <unistd.h>
#endif
#include "geneVarUtils.hpp"
#include "riskLocus.hpp"
#include <ctype.h>


	int outputGenotypes(FILE *fo, masterLocusFile *vf, analysisSpecs spec,int *cc,allelePair *a)
{
	int nSub,s,c,g,i;
	int gCount[5],gCcCount[2][5];
	char pathVarInfo[1000];
	nSub=vf->getTotalSubs();
	if (!vf->loadFirst(spec))
		return 0;
	vf->openLocusFiles();
	do {
		vf->outputCurrentAlleles(a,spec);
#if 0
		consequenceType t=vf->currentWorstConsequenceType();
		strcpy(comment,vf->currentQuickConsequence());
		for (ptr=comment;*ptr;++ptr)
				if (isspace(*ptr))
					*ptr='_';
		// fprintf(fo,":%d:%ld:%s:%d:%s:%s:\t",vf->currentChr(),vf->currentPos(),vf->currentID(),t,comment,vf->currentPolyPhen());
		sprintf(pathVarInfo,":%d:%ld:%s:%d:%s:%s:\t",vf->currentChr(),vf->currentPos(),vf->currentID(),t,comment,vf->currentPolyPhen());
#else
		vf->writeRiskVarInfo(pathVarInfo,0);
#endif
		for (c = 0; c < 5; ++c)
		{
			gCount[c] = 0;
			for (i=0;i<2;++i)
				gCcCount[i][c]=0;
		}
		for (s = 0; s < nSub; ++s)
		{
			int aa[2];
			for (i = 0; i < 2; ++i)
			{
				aa[i] = a[s][i];
				if (aa[i]>2)
					aa[i] = 2;
			}
			g=aa[0]+aa[1];
			++gCount[g];
			++gCcCount[cc[s]][g];
		}
		// fprintf(fo,":%d:%d:%d:\t",gCount[2],gCount[3],gCount[4]);
		sprintf(strchr(pathVarInfo,'\0'),":%d:%d:%d:\t",gCount[2],gCount[3],gCount[4]);
		sprintf(strchr(pathVarInfo,'\0'),":%d:%d:%d:%d:%d:%d:\t",
			gCcCount[0][2],gCcCount[0][3],gCcCount[0][4],gCcCount[1][2],gCcCount[1][3],gCcCount[1][4]);
		fprintf(fo,"%-200.200s",pathVarInfo);
		for (s=0;s<nSub;++s)
			fprintf(fo,"%d %d\t",a[s][0],a[s][1]);
		fprintf(fo,"\n");
	}
	while (vf->loadNext(spec));
	vf->closeLocusFiles();
	return 1;
}
