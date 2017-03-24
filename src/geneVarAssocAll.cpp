#include <stdlib.h>
#ifndef MSDOS 
#include <unistd.h>
#endif
#include "geneVarUtils.hpp"

void getMLP(char *fn,float *mlp,float *recMlp,float *HWEMlp,float *HWEContMlp,float *recHOMMlp,float *HOMMlp,float *HOMContMlp,float *numYY)

{
	FILE *fp;
	char buff[1000];
	float thisYY;
	*mlp=*recMlp=*HWEMlp=*HWEContMlp=*HOMMlp=*HOMContMlp=0;
	*numYY=0;
	fp=fopen(fn,"r");
	while (*buff=0,fgets(buff,1999,fp) && strncmp("-log(p)",buff,strlen("-log(p)"))) ;
	sscanf(buff,"-log(p) =%f",mlp);
	while (*buff=0,fgets(buff,1999,fp) && strncmp("-log(p)",buff,strlen("-log(p)"))) 
	{
		if (strstr(buff,"always"))
			sscanf(strchr(buff,'(')+1,"%f",&thisYY);
		if (thisYY>*numYY)
			*numYY=thisYY;
	}
	sscanf(buff,"-log(p) =%f",recMlp);
	while (*buff=0,fgets(buff,1999,fp) && strncmp("-log(p)",buff,strlen("-log(p)"))) ;
	sscanf(buff,"-log(p) =%f",HWEMlp);
	while (*buff=0,fgets(buff,1999,fp) && strncmp("-log(p)",buff,strlen("-log(p)"))) ;
	sscanf(buff,"-log(p) =%f",HWEContMlp);
	while (*buff=0,fgets(buff,1999,fp) && strncmp("-log(p)",buff,strlen("-log(p)"))) ;
	sscanf(buff,"-log(p) =%f",recHOMMlp);
	while (*buff=0,fgets(buff,1999,fp) && strncmp("-log(p)",buff,strlen("-log(p)"))) ;
	sscanf(buff,"-log(p) =%f",HOMMlp);
	while (*buff=0,fgets(buff,1999,fp) && strncmp("-log(p)",buff,strlen("-log(p)"))) ;
	sscanf(buff,"-log(p) =%f",HOMContMlp);
	fclose(fp);
	return;
}

int main(int argc,char *argv[])
{
	geneExtractor gcase,gcont;
	refseqGeneInfo r;
	char fn[100],fn2[100],line[1000],geneName[100];
	FILE *fp,*fo;
	gvaParams gp;
	analysisSpecs spec;
	float mlp,recMlp,HWEMlp,HWEContMlp,recHOMMlp,HOMMlp,HOMContMlp,numYY;
	int result,problemExtractingGene;
	int i,k; k=0;
//	_CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );
	fp=fopen(argv[1],"r");
	if (fp == NULL)
		{ dcerror(1, "Could not open file %s", argv[1]); exit(1); }
	gp.firstGeneNum=atoi(argv[3]);
	gp.lastGeneNum=atoi(argv[4]);
	gp.input(fp,spec);
	fclose(fp);	
	masterLocusFile vf(gp.nCc[0]+gp.nCc[1]);
	r.setListFile(gp.geneListFn);
	r.setBaitsFile(gp.baitFn);
	if (gp.referencePath[0]!='\0')
		r.setReferencePath(gp.referencePath);
	r.setUpstream(gp.upstream);
	r.setUpstream(gp.downstream);
	r.setBaitMargin(gp.margin);
	r.goToStart();
	fo=fopen(argv[2],"w");
	while (r.getNextGene())
	{
		// printf("%s\n",r.getGene());
		++k;
		if (k>gp.lastGeneNum)
			break;
		if (k<gp.firstGeneNum)
			continue;
		unlink("gva.all.db");
		unlink("gva.all.vdx");
		vf.openFiles("gva.all.db","gva.all.vdx");
		problemExtractingGene=0;
		int ff;
		for (i=0;i<gp.nCc[0];++i)
		{
			gcont.setVariantFileName(gp.ccFn[0][i]);
			sprintf(fn,"gva.all.cont.%d.vcf",i+1);
			unlink(fn);
			if (gcont.extractGene(r,fn,0,spec.addChrInVCF[ff++])==0)
			{
				problemExtractingGene=1;
				break;
			}
			vf.addLocusFile(fn,VCFFILE);
			vf.readLocusFileEntries(fn,spec,0);
		}
		for (i=0;i<gp.nCc[1];++i)
		{
			gcase.setVariantFileName(gp.ccFn[1][i]);
			sprintf(fn,"gva.all.case.%d.vcf",i+1);
			unlink(fn);
			if (gcase.extractGene(r,fn,0,spec.addChrInVCF[ff++])==0)
			{
				problemExtractingGene=1;
				break;
			}
			vf.addLocusFile(fn,VCFFILE);
			vf.readLocusFileEntries(fn,spec,1);
		}
		if (problemExtractingGene==1)
		{
			vf.closeLocusFiles();
			vf.closeFiles();
			continue;
		}

		if (spec.consequenceThreshold!=0 || spec.useConsequenceWeights!=0)
		{
			if (spec.useEnsembl)
				vf.getEnsemblConsequences(spec);
			else
				vf.getQuickConsequences(r,spec);
		}
		unlink("gva.all.par");
		unlink("gva.all.dat");
		unlink("gva.all.sao");
		vf.writeOldScoreAssocFiles("gva.all",gp.wf,gp.wFunc,gp.useFreqs,gp.nSubs,1,gp.writeComments,spec);
		vf.closeLocusFiles();
		vf.closeFiles();
		system("scoreassoc gva.all.par gva.all.dat gva.all.sao");
		getMLP("gva.all.sao",&mlp,&recMlp,&HWEMlp,&HWEContMlp,&recHOMMlp,&HOMMlp,&HOMContMlp,&numYY);
		fprintf(fo," %6d %-5s %12d %-20s %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.0f %6.4f\n",
			k,r.getChr(),r.getStart(),r.getGene(),mlp,recMlp,HWEMlp,HWEContMlp,recHOMMlp,HOMMlp,HOMContMlp,numYY);
		printf(" %6d %-20s %6.4f\n",k,r.getGene(),mlp);
		if (mlp>=2.0 || (recMlp>=1 && HWEMlp>=2.0))
		{
			sprintf(fn,"gva.all.%s.%05.2f.sao",r.getGene(),mlp);
			unlink(fn);
			result=rename("gva.all.sao",fn);
			if (result!=0) // WTF?
			{
				sprintf(fn2,"type gva.all.sao > %s",fn);
				system(fn2);
			}

		}
	}
	fclose(fo);

	return 0;
}
