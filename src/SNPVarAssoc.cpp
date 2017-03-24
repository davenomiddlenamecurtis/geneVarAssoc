#include <stdlib.h>
#ifndef MSDOS 
#include <unistd.h>
#endif
#include "geneVarUtils.hpp"

int main(int argc,char *argv[])
{
	char fn[100],fn2[100],line[1000],SNPName[100],SNPNames[100],posSpec[100],*ptr;
	int i,first,cc;
	FILE *fp,*fg;
	gvaParams gp;
	analysisSpecs spec;
	fp=fopen(argv[1],"r");
	gp.input(fp,spec);
	fclose(fp);
	masterLocusFile vf(gp.nCc[0]+gp.nCc[1]);
	strcpy(SNPNames,argv[2]);
	
	sprintf(fn,"gva.SNPs.db");
	sprintf(fn2,"gva.SNPs.vdx");
	unlink(fn);
	unlink(fn2);
	vf.openFiles(fn,fn2);
	first=1;
	fg=fopen(SNPNames,"r");
	while (fgets(line,999,fg) && sscanf(line,"%s",SNPName)==1)
	{
		ptr=strchr(SNPName,':');
		*ptr='\0';
	int ff=0;
	for (cc=0;cc<2;++cc)
		for (i=0;i<gp.nCc[cc];++i)
		{
			sprintf(posSpec,"%s%s:%s-%s",spec.addChrInVCF[ff++]?"chr":"",SNPName,ptr+1,ptr+1);
			sprintf(fn,"gva.altSubs.%s.%d.vcf",cc?"case":"cont",i+1);
			sprintf(line,"tabix %s %s %s %s%s",first==1?"-h ":"",gp.ccFn[cc][i],posSpec,first==1?">":">>",fn);
			printf("%s\n",line);
			system(line);
		}
	first=0;
	}
	fclose(fg);
	for (cc=0;cc<2;++cc)
	for (i = 0; i < gp.nCc[cc]; ++i)
	{
		sprintf(fn, "gva.altSubs.%s.%d.vcf", cc ? "case" : "cont", i + 1);
		vf.addLocusFile(fn,VCFFILE);
		vf.readLocusFileEntries(fn, spec, cc);
	}
	sprintf(fn,"gva.%s","SNPs");
	vf.writeOldScoreAssocFiles(fn,gp.wf,gp.wFunc,gp.useFreqs,gp.nSubs,1,gp.writeComments,spec);

	sprintf(line,"scoreassoc gva.SNPs.par gva.SNPs.dat %s",argv[3]);
	if (gp.writeScoreFile==1)
		sprintf(strchr(line,'\0')," gva.%s.sco",fn);
	system(line);

	return 0;
}