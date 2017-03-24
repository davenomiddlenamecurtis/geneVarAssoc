#include <stdlib.h>
#ifndef MSDOS 
#include <unistd.h>
#endif
#include "geneVarUtils.hpp"

// do gva analysis on list of genes
// C:\sequence\useScores\veryRare\SLPs>\msvc\vcf\geneVarAssocSome \sequence\sharedseq\gva.SSS.ct08.veryRare.score.par allGenes.txt
void getMLP(char *fn,float *slp,float *wilcSlp,float *recMlp,float *HWEMlp,float *HWEContMlp,float *HOMMlp,float *HOMContMlp,float *numYY)
{
	FILE *fp;
	char buff[1000];
	float thisYY;
	*slp=*wilcSlp=*recMlp=*HWEMlp=*HWEContMlp=*HOMMlp=*HOMContMlp=0;
	*numYY=0;
	fp=fopen(fn,"r");
#ifdef USEMLP
	while (*buff=0,fgets(buff,1999,fp) && strncmp("-log(p)",buff,strlen("-log(p)"))) ;
	sscanf(buff,"-log(p) =%f",mlp);
#else
	while (*buff=0,fgets(buff,1999,fp) && strncmp("SLP",buff,strlen("SLP"))) ;
	sscanf(buff,"SLP =%f",slp);
	while (*buff=0,fgets(buff,1999,fp) && strncmp("Wilcoxon SLP=",buff,strlen("Wilcoxon SLP"))) ;
	sscanf(buff,"Wilcoxon SLP=%f",wilcSlp);	
#endif
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
	FILE *fp,*fo,*fi;
	gvaParams gp;
	analysisSpecs spec;
	float slp,wilcSlp,recMlp,HWEMlp,HWEContMlp,HOMMlp,HOMContMlp,numYY;
	int result,problemExtractingGene;
	int i,k; k=0;
	fp=fopen(argv[1],"r");
	if (fp == NULL)
		{ dcerror(1, "Could not open file %s", argv[1]); exit(1); }
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
	fi=fopen(argv[2],"r");
	fo=fopen(argv[3],"w");
	dcerror.warn();
	while (fgets(line,100,fi) && sscanf(line,"%s",geneName)==1)
	{
		++k;
		if (!r.findGene(geneName) || !r.getNextGene())
			continue;
		unlink("gva.all.db");
		unlink("gva.all.vdx");
		vf.openFiles("gva.all.db","gva.all.vdx");
		problemExtractingGene=0;
		int ff=0;
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
	if (spec.consequenceThreshold!=0)
		sprintf(fn,"gva.%s.ct%02d",geneName,spec.consequenceThreshold);
	else if (spec.useConsequenceWeights!=0)
		sprintf(fn,"gva.%s.ucw",geneName);
	else
		sprintf(fn,"gva.%s",geneName);
		vf.writeOldScoreAssocFiles(fn,gp.wf,gp.wFunc,gp.useFreqs,gp.nSubs,1,gp.writeComments,spec);
		vf.closeLocusFiles();
		vf.closeFiles();
		sprintf(line,"scoreassoc %s.par %s.dat %s.sao",fn,fn,fn);
		if (gp.writeScoreFile == 1)
			sprintf(strchr(line, '\0'), " %s.sco", fn);
		system(line);
		sprintf(line,"%s.sao",fn);
		getMLP(line,&slp,&wilcSlp,&recMlp,&HWEMlp,&HWEContMlp,&HOMMlp,&HOMContMlp,&numYY);
		fprintf(fo," %6d %-5s %12d %-20s %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.0f\n",
			k,r.getChr(),r.getStart(),r.getGene(),slp,wilcSlp,recMlp,HWEMlp,HWEContMlp,HOMMlp,HOMContMlp,numYY);
		printf(" %6d %-20s %6.4f\n",k,r.getGene(),slp);
	}
	fclose(fo);
	fclose(fi);

	return 0;
}
