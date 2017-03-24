#include "masterLocusFile.hpp"
#include "geneVarUtils.hpp"
#ifndef MSDOS 
#include <unistd.h>
#endif
#include <ctype.h>
#include <string.h>

int masterLocusFile::writeAltSubs(char *fn, analysisSpecs &spec)
{
	int i,totalSub;
	FILE *fo;
	totalSub=0;
	for (i=0,totalSub=0;i<nLocusFiles;++i)
		totalSub+=nSubs[i];
	allelePair *all;
	all=(allelePair *)calloc(totalSub,sizeof(allelePair));
	strEntry *subName;
	subName=(strEntry *)calloc(totalSub,sizeof(strEntry));
	openLocusFiles();
	outputSubNames(subName,spec);
	if (gotoFirstInRange(spec))
		outputCurrentAlleles(all,spec);
	if (fn==0)
		fo=stdout;
	else
		fo=fopen(fn,"w");
	for (i=0;i<totalSub;++i)
		if (all[i][0]!=0 && (all[i][0]!=01|| all[i][1]!=1))
			fprintf(fo,"%s\t%d %d\n",subName[i],all[i][0],all[i][1]);
	if (fn!=0)
		fclose(fo);
	free(all);
	free (subName);
	return 1;
}

int main(int argc,char *argv[])
{
	char fn[100],fn2[100],line[1000],posSpec[100],vcfFn[100],vcfFnBuff[100],chrStr[20],*ptr;
	int i,totalSub,cc;
	FILE *fp;
	gvaParams gp;
	analysisSpecs spec;
	if (!gp.readParms(argc,argv,spec))
		exit(1);
	spec.useConsequenceWeights=0; // I am not going to annotate these variants
	masterLocusFile vf(gp.nCc[0]+gp.nCc[1]);
	if ((ptr=strchr(gp.posName,':'))==0)
		dcerror(1,"Usage: showAltSubs --arg-file something.arg --position 7:12139555");
	*ptr='\0';
	sprintf(posSpec,"%s%s:%s-%s",spec.addChrInVCF[0]?"chr":"",gp.posName,ptr+1,ptr+1);
	sprintf(fn,"gva.altSubs.db");
	sprintf(fn2,"gva.altSubs.vdx");
	unlink(fn);
	unlink(fn2);
	vf.openFiles(fn,fn2);
	for (cc=0;cc<2;++cc)
		for (i=0;i<gp.nCc[cc];++i)
		{
			strcpy(vcfFn,gp.ccFn[cc][i]);
			if (strchr(vcfFn,'*'))
			{
				strcpy(vcfFnBuff,vcfFn);
				*strchr(vcfFnBuff,'*')='\0';
				strcat(vcfFnBuff,gp.posName);
				strcat(vcfFnBuff,strchr(vcfFn,'*')+1);
				strcpy(vcfFn,vcfFnBuff);
			}
			sprintf(fn,"gva.altSubs.%s.%d.vcf",cc?"case":"cont",i+1);
			sprintf(line,"tabix -h %s %s >%s",vcfFn,posSpec,fn);
			system(line);
			vf.addLocusFile(fn,VCFFILE);
			vf.readLocusFileEntries(fn,spec,1); // this is just an easy way to make sure that files are known about
		}
	// sprintf(fn,"%s_%s.aso",argv[2],ptr+1);
	vf.writeAltSubs(0,spec);
	return 0;
}
