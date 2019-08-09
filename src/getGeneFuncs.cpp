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

#include "getGene.hpp"
#include "getSequence.hpp"
#include <ctype.h>

char refseqGeneInfo::geneLine[GENELINELENGTH];

#ifndef MSDOS
#define stricmp strcasecmp
#define strnicmp strncasecmp
#define strupr mystrupr
#include <ctype.h>
static char *mystrupr(char *s)
{
	char *t;
	t=s;
	while (*t)
	{
		if (isalpha(*t))
			*t=toupper(*t);
		++t;
	}
	return s;
}
#endif

// #bin	name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	score	name2	cdsStartStat	cdsEndStat	exonFrames
// use zero-based starts and one-based ends
// this is what we have in the refseq file downloaded from UCSC
// the coordinates we receive for pos will be one-based, as used by VCF and ensembl

int compInt ( const void * elem1, const void * elem2 ) 
{ return *(int *)elem1<*(int *)elem2?-1:*(int *)elem1==*(int *)elem2?0:1; }

void refseqTranscript::print3PrimeSequence(faSequenceFile &f,FILE *fp)
{
	char *threePrimeSeq;
	threePrimeSeq=(char *)malloc(txEnd-cdsEnd+1);
	f.getSequence(threePrimeSeq,cdsEnd+1,txEnd-cdsEnd);
	// strupr(threePrimeSeq);
	fprintf(fp,"%s",threePrimeSeq);
	free(threePrimeSeq);
}

int refseqTranscript::get3PrimeMatches(faSequenceFile &f,intervalList &matches,char *toMatch)
{
	char *threePrimeSeq,*site;
	int rv;
	site=(char*)malloc(strlen(toMatch)+1);
	if (strand=='+')
	{
		threePrimeSeq=(char *)malloc(txEnd-cdsEnd+1);
		f.getSequence(threePrimeSeq,cdsEnd+1,txEnd-cdsEnd);
		strupr(threePrimeSeq);
		strcpy(site,toMatch);
		rv=matches.getAllMatches(threePrimeSeq,site,cdsEnd+1);
	}
	else
	{
		threePrimeSeq=(char*)malloc(cdsStart-txStart+1);
		f.getSequence(threePrimeSeq,txStart,cdsStart-txStart);
		strupr(threePrimeSeq);
		getComplementarySequence(site,toMatch);
		rv=matches.getAllMatches(threePrimeSeq,site,txStart);
	}
	free(site); // one day make it so I can read in a  number of strings to match without having to get 3prime sequence again each time
	free(threePrimeSeq);
	return rv;
}

void refseqGeneInfo::getFivePrime(int *st,int *en,int t)
{
	if (strand=='+')
	{
		*st=transcript[t].txStart;
		*en=transcript[t].cdsStart; // I think this is OK because ends are 1-based, starts are 0-based
	}
	else
	{
		*st=transcript[t].cdsEnd;
		*en=transcript[t].txEnd;
	}
}

void refseqGeneInfo::getThreePrime(int *st,int *en,int t)
{
	if (strand=='+')
	{
		*st=transcript[t].cdsEnd;
		*en=transcript[t].txEnd;
	}
	else
	{
		*st=transcript[t].txStart;
		*en=transcript[t].cdsStart;
	}
}

int refseqTranscript::getEffect(faSequenceFile &f,int pos,char *all0,char *all1,int promoterLength,int downstreamLength,int getKozak)
{
	// if in coding region, check all alternative alleles for most deleterious
	int i;
	effect.clear();
	if (strand=='+')
	{
		if (pos<=exonStarts[0]-promoterLength || pos>exonEnds[exonCount-1]+downstreamLength)
		{
			effect.consequence[effect.nConsequence++]=INTERGENIC;
			return 1;
		}
	}
	else // if (strand=='-')
	{
		if (pos<=exonStarts[0]-downstreamLength||pos>exonEnds[exonCount-1]+promoterLength)
		{
			effect.consequence[effect.nConsequence++]=INTERGENIC;
			return 1;
		}
	}
// Kozak sequence
// (gcc)gccRccAUGG - A is at s, 0-coded
	if (getKozak && strand=='+' && cdsStart-pos<9 && pos-cdsStart<=5) // pos>cdsStart-9
		effect.consequence[effect.nConsequence++]=KOZAK;
	else if (getKozak && strand=='-' && pos-cdsEnd<=9 && cdsEnd-pos<5) // pos<=cdsEnd+9
		effect.consequence[effect.nConsequence++]=KOZAK;
	else if (strand=='+' && pos>exonStarts[0]-promoterLength && pos<=exonStarts[0])
		effect.consequence[effect.nConsequence++]=UPSTREAM;
	else if (strand=='-' && pos<=exonEnds[exonCount-1]+promoterLength && pos>exonEnds[exonCount-1])
		effect.consequence[effect.nConsequence++]=UPSTREAM;
	else if (strand=='-' && pos>exonStarts[0]-downstreamLength && pos<=exonStarts[0])
		effect.consequence[effect.nConsequence++]=DOWNSTREAM;
	else if (strand=='+' && pos<=exonEnds[exonCount-1]+downstreamLength && pos>exonEnds[exonCount-1])
		effect.consequence[effect.nConsequence++]=DOWNSTREAM;
	else for (i=0;i<exonCount-1;++i)
	{
		// put the (starts+1) in parentheses to make up for them being zero-based
		spliceSiteType sT;
		bool essential;
		int sSStart,sSEnd;
		if (strand == '+' && pos>exonEnds[i] - NSBDONOREXON && pos <= exonEnds[i] + NSBDONORINTRON)
		{
			essential=pos>exonEnds[i] && pos <= exonEnds[i] + 2;
			getSpliceSiteSequence(f, exonEnds[i] - NSBDONOREXON + 1, exonEnds[i] + NSBDONORINTRON, DONOR, essential, pos, all0, all1);
		}
		else if (strand == '+' && pos >= (exonStarts[i + 1] + 1) - NSBACCEPTORINTRON && pos<(exonStarts[i + 1] + 1) + NSBACCEPTOREXON)
		{
			essential=pos>=(exonStarts[i + 1] + 1) - 2 && pos<(exonStarts[i + 1] + 1);
			getSpliceSiteSequence(f, (exonStarts[i + 1] + 1) - NSBACCEPTORINTRON, (exonStarts[i + 1] + 1) + NSBACCEPTOREXON - 1, ACCEPTOR, essential, pos, all0, all1);
		}
		else if (strand == '-' && pos>exonEnds[i] - NSBACCEPTOREXON && pos <= exonEnds[i] + NSBACCEPTORINTRON)
		{
			essential=pos>exonEnds[i] && pos <= exonEnds[i] + 2;
			getSpliceSiteSequence(f, exonEnds[i] - NSBACCEPTOREXON + 1, exonEnds[i] + NSBACCEPTORINTRON, ACCEPTOR, essential, pos, all0, all1);
		}
		else if (strand == '-' && pos >= (exonStarts[i + 1] + 1) - NSBDONORINTRON && pos < (exonStarts[i + 1] + 1) + NSBDONOREXON)
		{
			essential=pos >= (exonStarts[i + 1] + 1) - 2 && pos < (exonStarts[i + 1] + 1);
			getSpliceSiteSequence(f, (exonStarts[i + 1] + 1) - NSBDONORINTRON, (exonStarts[i + 1] + 1) + NSBDONOREXON - 1, DONOR, essential, pos, all0, all1);
		}
		else if (pos>exonEnds[i]&&pos<=exonStarts[i+1])
			effect.consequence[effect.nConsequence++]=INTRONIC;
		// calling getSpliceSiteSequence sets effect.consequence[effect.nConsequence++]
	}

	// this will miss situation where insertion hits start of exon
	if (effect.nConsequence!=0)
		return 1;

	// by getting to here we have established position is in an exon
	if ((strand=='+' && pos<=cdsStart-(getKozak?9:0)) || (strand=='-' && pos>cdsEnd+(getKozak?9:0)))
		effect.consequence[effect.nConsequence++]=FIVEPRIME_UTR;
	else if ((strand=='-' && pos<=cdsStart) || (strand=='+' && pos>cdsEnd))
		effect.consequence[effect.nConsequence++]=THREEPRIME_UTR;
	if (effect.nConsequence!=0)
		return 1;

	if (strlen(all0)!=1 || strlen(all1)!=1)
	{
		// must cast to int here because strlen returns size_t, which is unsigned
		if ((int)(strlen(all1)-strlen(all0))%3!=0)
			effect.consequence[effect.nConsequence++]=FRAMESHIFT_CODING;
		else
			effect.consequence[effect.nConsequence++]=CODINGINDEL; 
	// may be same length but then something tricky
	}
	else
		getCodingEffect(f,pos,all0,all1);
	return 1;
}

void refseqGeneInfo::setReferencePath(char *s)
{
	char slashChar,*ptr;
#ifdef MSDOS
	slashChar='\\';
#else
	slashChar='/';
#endif
	strcpy(referencePath,s);
	ptr=strchr(referencePath,'\0');
	--ptr;
	if (*ptr!=slashChar)
	{
		ptr[1]=slashChar;
		ptr[2]='\0';
	}
}

codonReader cr;
consequenceType refseqGeneInfo::getEffect(int pos,char *all0,char *all1,int promoterLength,int downstreamLength,int getKozak)
{
	char fn[100];
	int t,c;
	faSequenceFile f;

	sprintf(fn,"%s%s.fa",referencePath,chr);
	if (!f.init(fn))
		return NULL_CONSEQUENCE;
	worstConsequence=NULL_CONSEQUENCE;
	for (t=0;t<nTranscript;++t)
	{
		transcript[t].getEffect(f,pos,all0,all1,promoterLength,downstreamLength,getKozak);
		for (c=0;c<transcript[t].effect.nConsequence;++c)
		if (worstConsequence<transcript[t].effect.consequence[c])
		{
			worstEffect=transcript[t].effect;
			worstConsequence=transcript[t].effect.consequence[c];
		}
	}
return worstConsequence;
}

const char *refseqGeneInfo::tellEffect()
{
	char featureBuff[10000];
	int c;
	for (c=0;c<NCONSEQUENCETYPES;++c)
		if (consequence[c].t==worstConsequence)
		{
			sprintf(featureBuff,"%s %s %s %s %s",consequence[c].str,geneName,worstEffect.aaStr,worstEffect.codonStr,worstEffect.spliceStr);
			break;
		}
	if (c==NCONSEQUENCETYPES)
	{
		dcerror(1,"No match for consequence type %d\n",worstConsequence);
		return NULL;
	}
	if (worstEffect.posInCDS != 0)
		sprintf(strchr(featureBuff, '\0'), " %d", worstEffect.posInCDS);
	if (worstEffect.aaPos != 0)
		sprintf(strchr(featureBuff, '\0'), " %d", worstEffect.aaPos);
	return featureBuff;
}

int refseqTranscript::getSpliceSiteSequence(faSequenceFile &f,int sSStart,int sSEnd,spliceSiteType sST,bool essential,int pos,char *a0,char *a1)
{
	char all0[2],all1[2];
	all0[0]=a0[0];
	all1[0]=a1[0];
	all0[1]=all1[1]='\0';
	//make local copies in case need to swap later
	char sSSeq[2][(NSBACCEPTOREXON+NSBACCEPTORINTRON)*2+1],buff[(NSBACCEPTOREXON+NSBACCEPTORINTRON)*2+1];
	int seqPos,i,j;
	f.getSequence(buff,sSStart,sSEnd-sSStart+1);
	for (i=0;i<=sSEnd-sSStart;++i)
		buff[i]=tolower(buff[i]);
	seqPos=(strand=='+')?0:sSEnd-sSStart;
	for (i=0;i<=sSEnd-sSStart;++i)
	{
		sSSeq[0][seqPos]=buff[i];
		if (sSStart+i==pos)
		{
			if (strlen(a0)==1 && strlen(a1)==1)
			{
				if (toupper(buff[i])!=all0[0])
					{
						dcerror(1,"base in site does not match reference allele %c at position %d - ",
							all0[0],pos);
						if (toupper(buff[i])!=all1[0])
							dcerror(1,"it does not match alt allele %c either\n",all1[0]);
						else
						{
							char a;
							dcerror(0,"Will swap and use alt allele %c as reference\n",all1[0]);
							a=all0[0]; all0[0]=all1[0]; all1[0]=a;
						}
					}
				sSSeq[1][seqPos]=all1[0];
			}
			else
				sSSeq[1][seqPos]='*'; // not a single base substitution
		}
		else
			sSSeq[1][seqPos]=buff[i];
		seqPos+=(strand=='+')?1:-1;
	}
	sSSeq[0][i]=sSSeq[1][i]='\0';
	if (strand=='-')
		for (i=0;i<=sSEnd-sSStart;++i)
			for (j=0;j<2;++j)
				sSSeq[j][i]=getCompBase(sSSeq[j][i]);
// M-A-G-[cut]-G-U-R-A-G-U (donor site)
// The splice donor site includes an almost invariant sequence GU at the 5' end of the intron
// C-U-R-[A]-Y (branch sequence 20-50 nucleotides upstream of acceptor site)
// the 3' end of the intron terminates the intron with an almost invariant AG sequence
// Upstream (5'-ward) from the AG there is a region high in pyrimidines (C and U), or polypyrimidine tract. 
// Y-rich-N-C-A-G-[cut]-G (acceptor site)
	if (sST==DONOR)
		for (i=NSBDONOREXON;i<NSBDONOREXON+2;++i)
			sSSeq[0][i]=toupper(sSSeq[0][i]); // aim to upper-case the invariant GT
	else
		for (i=NSBACCEPTORINTRON-2;i<NSBACCEPTORINTRON;++i)
			sSSeq[0][i]=toupper(sSSeq[0][i]); // aim to upper-case the invariant AG 
	sprintf(effect.spliceStr,"%s/%s",sSSeq[0],sSSeq[1]);
	effect.consequence[effect.nConsequence++]=essential?ESSENTIAL_SPLICE_SITE:SPLICE_SITE;
}

int refseqTranscript::getCodingEffect(faSequenceFile &f,int pos,char *a0,char *a1)
{
	consequenceType nsType;
	int posInGene,cdsStartExon, cdsEndExon,posExon,e,b,framePos;
	char oldAA,newAA;
	char codon[4],newCodon[4],buff[4],all0[2],all1[2];
	if (!strchr("GACT", *a0) || !strchr("GACT", *a1))
	{
		dcerror(1,"Bad allele pair %c %c at position %d\n", *a0, *a1, pos);
		return 0;
	}
	all0[0]=a0[0];
	all1[0]=a1[0];
	all0[1]=all1[1]='\0';
	//make local copies in case need to swap later
	posInGene=pos-cdsStart; // cdStart is base 0, pos and posInGene are base 1
	cdsStartExon=-1;
	for (e=0;e<exonCount;++e)
		if (cdsStart>=exonStarts[e] && cdsStart<exonEnds[e])
		{
			cdsStartExon=e;
			break;
		}
	if (cdsStartExon==-1)
	{
		dcerror(1,"Could not find cdsStartExon in %s",geneName);
		return 0;
	}
	posExon=-1;
	for (e=0;e<exonCount;++e)
		if (pos>exonStarts[e] && pos<=exonEnds[e])
		{
			posExon=e;
			break;
		}
	if (posExon==-1)
	{
		dcerror(1,"Could not find posExon in %s",geneName);
		return 0;
	}
	// for gene on the - strand we will find the position from cdsStart which is actually the "end" of the code
	// that is, posInGene is the position along the positive strand
	// this can give us framePos, also along the positive strand
	// remove lengths of intervening introns
	for (e=cdsStartExon;e<posExon;++e)
		posInGene-=exonStarts[e+1]-exonEnds[e]; 
	framePos=(posInGene-1)%3; // 0, 1 or 2 - framePos is zero-based
	f.getSequence(buff,pos-framePos,3);
	for (b=0;b<3;++b)
	{
		codon[b]=buff[b];
		if (pos-framePos+b==exonEnds[posExon]) // can never happen for last exon because would be out of frame
			f.getSequence(buff,exonStarts[posExon+1]-b,3);

	}
	codon[3]='\0';
	if (toupper(codon[framePos])!=all0[0])
	{
		dcerror(1,"base %d of codon %s does not match reference allele %c at position %d - ",
			framePos,codon,all0[0],pos);
		if (toupper(codon[framePos])!=all1[0])
			dcerror(1,"it does not match alt allele %c either\n",all1[0]);
		else
		{
			char a;
			dcerror(0,"Will swap and use alt allele %c as reference\n",all1[0]);
			a=all0[0]; all0[0]=all1[0]; all1[0]=a;
		}
	}
	strcpy(newCodon,codon);
	newCodon[framePos]=all1[0];
	for (b=0;b<3;++b)
	{
		if (b==framePos)
		{
		codon[b]=toupper(codon[b]);
		newCodon[b]=toupper(newCodon[b]);
		}
		else
		{
		codon[b]=tolower(codon[b]);
		newCodon[b]=tolower(newCodon[b]);
		}
	}
	oldAA=cr.translate(codon,strand); // calling translate() actually changes the codon if on - strand
	newAA=cr.translate(newCodon,strand);
	if (oldAA==newAA)
	{
		sprintf(effect.aaStr,"%c",oldAA);
		sprintf(effect.codonStr,"%s/%s",codon,newCodon);
		nsType=SYNONYMOUS_CODING;
	}
	else
	{
		if (newAA=='*')
			nsType=STOP_GAINED;
		else if (oldAA=='*')
			nsType=STOP_LOST;
		else
			nsType=NON_SYNONYMOUS_CODING;
		sprintf(effect.aaStr,"%c/%c",oldAA,newAA);
		sprintf(effect.codonStr,"%s/%s",codon,newCodon);
	}
	effect.consequence[effect.nConsequence++]=nsType;
	if (strand!='+')
	{
		posInGene = cdsEnd-pos +1; // cdsEnd is base 1, pos and posInGene are base 1
		cdsEndExon = -1;
		for (e = 0; e < exonCount; ++e)
			if (cdsEnd >= exonStarts[e] && cdsStart < exonEnds[e])
			{
				cdsEndExon = e;
				break;
			}
		if (cdsEndExon == -1)
		{
			dcerror(1, "Could not find cdsEndExon in %s", geneName);
			return 0;
		}
		for (e = posExon; e<cdsEndExon; ++e)
			posInGene -= exonStarts[e + 1] - exonEnds[e];
	}
	effect.posInCDS = posInGene; // allow comparison with ensembl - 1-based? - this is not matching VEP
	effect.aaPos = (effect.posInCDS - 1) / 3 + 1;
	return 1;
}

int refseqGeneInfo::getNextGene(int transcriptionStartCanVary)
{
	int firstLine,t;
	char lineName[100],lineChr[40],lineStrand;
	int lineTxStart,txStart,txEnd,lineTxEnd,lineCdsStart,cdsStart;
	long startPos;
	gotAllExons=allExonCount=0;
	if (geneListFile==0)
	{
		geneListFile=fopen(geneListFileName,"rb"); // binary so MSVC fseek() will work
		if (geneListFile==0)
		{
			dcerror(2,"Could not open file %s",geneListFileName);
			return 0;
		}
	fgets(geneLine,GENELINELENGTH-1,geneListFile); // miss header line
	}
	firstLine=1; // first line/transcript for this gene
	nTranscript=0;
	while (1)
	{
		startPos=ftell(geneListFile);
		if (!fgets(geneLine, GENELINELENGTH - 1, geneListFile))
		{
			if (firstLine)
				return 0; // ran out of file
			else
				break;
		}
		// check that both geneName and txStart have not changed - occasional genes are present in multiple locations (!)
		if (sscanf(geneLine,"%*d %*s %s %c %d %*d %d %d %*d %*s %*s %*d %s",
			lineChr,&lineStrand,&lineTxStart,&lineTxEnd,&lineCdsStart,lineName)!=6)
		{ dcerror(1,"Not enough parameters in line: %s",geneLine); return 0; }
		if (strchr(lineChr,'_'))
			continue; // ignore genes not assigned to chromosomes
		if (firstLine)
		{
			firstLine=0;
			strcpy(geneName,lineName);
			strcpy(chr,lineChr);
			txStart=lineTxStart;
			txEnd=lineTxEnd;
			cdsStart=lineCdsStart;
			strand=lineStrand;
		}
		if (strcmp(geneName, lineName) || strcmp(chr, lineChr) || (!transcriptionStartCanVary && txStart != lineTxStart && cdsStart != lineCdsStart)
			|| txStart>lineTxEnd || txEnd<lineTxStart)
			// DDX11L1 is on two different chromosomes!
			// SEC23B has different transcription starts
			// A2ML1 has different transcription starts and different coding starts - second transcript not translated?
			// MIR4315-1 is on two different bits of chr 17
			// do not allow transcripts that do not overlap at all
		{
			fseek(geneListFile,startPos,SEEK_SET); // go back to start of line for next gene
			break;
		}
		if (txStart>lineTxStart)
			txStart=lineTxStart;
		if (txEnd<lineTxEnd)
			txEnd=lineTxEnd;
		transcript[nTranscript++].read(geneLine);
	}
	firstExonStart=transcript[0].exonStarts[0];
	lastExonEnd=transcript[0].exonEnds[transcript[0].exonCount-1];
	for (t=1;t<nTranscript;++t)
	{
	if (firstExonStart>transcript[t].exonStarts[0])
		firstExonStart=transcript[t].exonStarts[0];
	if (lastExonEnd<transcript[t].exonEnds[transcript[t].exonCount-1])
		lastExonEnd=transcript[t].exonEnds[transcript[t].exonCount-1];
	}
	return 1;
}

char exonStartsStr[10000],exonEndsStr[10000];
int refseqTranscript::read(char *line)
{
	char *sPtr,*ePtr;
	int i;
	if (sscanf(line,"%*d %*s %s %c %d %d %d %d %d %s %s %*d %s",
			chr,&strand,&txStart,&txEnd,&cdsStart,&cdsEnd,&exonCount,exonStartsStr,exonEndsStr,geneName)!=10)
		{ dcerror(1,"Not enough parameters in line: %s",line); return 0; }
	for (i=0,sPtr=exonStartsStr,ePtr=exonEndsStr;i<exonCount;++i)
		{
			if (*sPtr=='"')
				++sPtr;
			if (*ePtr=='"')
				++ePtr;
			// above necessary if tab-delimited file saved from Excel
			// printf("%s\n",line);
			if (sscanf(sPtr,"%d",&exonStarts[i])!=1)
		{ dcerror(1,"Could not read enough exon starts from this line: %s",line); return 0; }
			sPtr=strchr(sPtr,',')+1;
			if (sscanf(ePtr,"%d",&exonEnds[i])!=1)
		{ dcerror(1,"Could not read enough exon ends from this line: %s",line); return 0; }
			ePtr=strchr(ePtr,',')+1;
			
		}
	return 1;
}
// #bin	name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	score	name2	cdsStartStat	cdsEndStat	exonFrames

int refseqGeneInfo::checkExonLengths()
{
	int t;
	for (t=0;t<nTranscript;++t)
		transcript[t].checkExonLengths();
	return 1;
}

int refseqTranscript::checkExonLengths()
{
	int i,len;
	for (i=0;i<exonCount;++i)
		printf("exon length mod 3 = %d\n",(exonEnds[i]+1-exonStarts[i])%3);
	len=cdsEnd-cdsStart+1;
	printf("coding region %d - %d, length = %d, mod 3 = %d\n",
		cdsStart,cdsEnd,len,len%3);
	len=0;
	for (i=0;i<exonCount;++i)
		if (exonEnds[i]>=cdsStart)
			break;
	len=exonEnds[i]-cdsStart+1;
	for (;i<exonCount;++i)
		if (exonEnds[i]>=cdsEnd)
		{
			len+=cdsEnd-exonStarts[i]+1;
			break;
		}
		else
			len+=exonEnds[i]-exonStarts[i]+1;

	printf("net length = %d, mod 3 = %d\n",
		len,len%3);

	return 1;
}

int refseqGeneInfo::tbiExtractGene(char *tbiFilename,char *outFn,int appendToOld,int addChrInVCF,int removeSpaces)
{
	int i,startPos,endPos,foundOne,systemStatus;
	char buff[1000],*tbiFn,*ptr,tbiFnBuff[1000];
	long lineStart;
	if ((ptr=strchr(tbiFilename,'*'))==0)
		tbiFn=tbiFilename;
	else
	{
		strcpy(tbiFnBuff,tbiFilename);
		ptr=strchr(tbiFnBuff,'*');
		*ptr='\0';
		strcat(tbiFnBuff,chr+3); // because first three chars of chr are "chr"
		ptr=strchr(tbiFilename,'*');
		strcat(tbiFnBuff,ptr+1);
		tbiFn=tbiFnBuff;
	}

#if 0
	// code was written originally to extract individual exons
	sprintf(line,"tabix -h %s ",tbiFn);
	for (i=0;i<exonCount;++i)
		sprintf(strchr(line,'\0'),"%s:%d-%d ",
		chr+3,
		exonStarts[i] - ((i==0&&strand=='+')?upstream:0) - ((i==0&&strand=='-')?downstream:0),
		exonEnds[i] + ((i==exonCount-1&&strand=='-')?upstream:0) + ((i==exonCount-1&&strand=='+')?downstream:0));
#endif
	// code to extract according to baits
	if (baitsFileName[0]!='\0')
	{
		if (baitsFile==0)
		{
			if ((baitsFile=fopen(baitsFileName,"rb"))==0)
			{
				dcerror(7,"Could not open baits file: %s",baitsFileName);
				return 0;
			}
		}
			// We will start slow and stupid
		fseek(baitsFile,0L,SEEK_SET);
			do
			{
				lineStart=ftell(baitsFile);
				if (!fgets(buff,999,baitsFile))
					return 0; // no bait found for this transcript
			} while (strncmp(buff,chr+3,strlen(chr+3)) || (sscanf(buff,"%*s %*d %d",&endPos),endPos-baitMargin<firstExonStart));

		fseek(baitsFile,lineStart,SEEK_SET);
		sprintf(geneLine,"b 692tabix %s%s ",tbiFn,appendToOld?"":" -h");

		foundOne=0; // what can happen is small transcript in refseq file may be missed completely
		while (fgets(buff,999,baitsFile)&&!strncmp(buff,chr+3,strlen(chr+3))&&(sscanf(buff,"%*s %d %d",&startPos,&endPos),startPos+baitMargin<=lastExonEnd))
			{
				foundOne=1;
#ifdef OLDTABIX
				sprintf(strchr(geneLine,'\0'),"%s:%d-%d ",chr+3,startPos+baitMargin,endPos+baitMargin);
#else
				sprintf(strchr(geneLine,'\0'),"%s:%d-%d ",chr+(addChrInVCF?0:3),startPos+baitMargin,endPos+baitMargin);
#endif
				if (strlen(geneLine)>3900) // line getting long, extract baits so far then go back for more - maximum is 4000?
				{
				sprintf(strchr(geneLine,'\0'),"%s %s %s",removeSpaces?"| sed s/' '/'_'/g ":"",appendToOld?">>":">",outFn);
				printf("Running command: %s\n",geneLine);
				checkSystem();
				systemStatus=system(geneLine);
				// printf("system returned %d\n",systemStatus);
				appendToOld=1;
				sprintf(geneLine,"tabix %s ",tbiFn);
				}
			}
		if (foundOne==0)
		{
			printf("Did not find any baits for %s\n",geneName);
			return 0;
		}
		if (strchr(geneLine+4,':')) // means we found another bait since last calling tabix, skipping start of line which may be C:
		{
		sprintf(strchr(geneLine,'\0'),"%s%s %s", removeSpaces ? "| sed s/' '/'_'/g " : "", appendToOld?">>":">",outFn);
		printf("Running command: %s\n",geneLine);
		checkSystem();
		systemStatus=system(geneLine);
		// printf("system returned %d\n",systemStatus);
		}
	}
	else
	{
	sprintf(geneLine,"tabix %s%s %s:%d-%d",tbiFn,appendToOld?"":" -h",
		chr+(addChrInVCF?0:3),
		firstExonStart - ((strand=='+')?upstream:downstream),
		lastExonEnd + ((strand=='-')?upstream:downstream));
	sprintf(strchr(geneLine,'\0'),"%s%s %s", removeSpaces ? "| sed s/' '/'_'/g " : "", appendToOld?">>":">",outFn);
	printf("Running command: %s\n",geneLine);
	systemStatus=system(geneLine);
	// printf("system returned %d\n",systemStatus);
	}
	return 1;
}

int refseqGeneInfo::goToStart()
{
	if (geneListFile==0)
	{
		geneListFile=fopen(geneListFileName,"rb"); // binary so MSVC fseek() will work
		if (geneListFile==0)
		{
			dcerror(2,"Could not open file %s",geneListFileName);
			return 0;
		}
	}
	else
		fseek(geneListFile,0L,SEEK_SET);
	fgets(geneLine,GENELINELENGTH-1,geneListFile); // miss header line
	return 1;
}


int refseqGeneInfo::findFirstGene(char *chrToFind,int posToFind)
{
	char lineChr[20],chrToFindStr[20];
	int lineStart;
	long startPos;
	sprintf(chrToFindStr,"chr%s",chrToFind);
	if (geneListFile==0)
	{
		geneListFile=fopen(geneListFileName,"rb"); // binary so MSVC fseek() will work
		if (geneListFile==0)
		{
			dcerror(2,"Could not open file %s",geneListFileName);
			return 0;
		}
	}
	else
		fseek(geneListFile,0L,SEEK_SET);
	fgets(geneLine,GENELINELENGTH-1,geneListFile); // miss header line
	do {
		startPos=ftell(geneListFile);
		if (!fgets(geneLine,GENELINELENGTH-1,geneListFile))
		{
			dcerror(3,"Could not find gene at position %s:%d in file",chrToFind,posToFind,geneListFileName);
			return 0;
		}
		if (sscanf(geneLine,"%*d %*s %s %*c %d",lineChr,&lineStart)!=2)
		{
			dcerror(4,"Could not read enough parameters from this line: %s",geneLine);
			return 0;
		}
	} while (stricmp(chrToFindStr,lineChr));
	// } while (stricmp(chrToFindStr,lineChr) || lineStart<posToFind); because transcripts are not in order along the chromosome
	fseek(geneListFile,startPos,SEEK_SET); // back to start of first line for this gene
	return 1;
}

int refseqGeneInfo::findGene(char *name)
{
	char lineName[100];
	if (!goToStart())
		return 0;
	long startPos;
	do {
		startPos=ftell(geneListFile);
		if (!fgets(geneLine,GENELINELENGTH-1,geneListFile) || sscanf(geneLine,"%s",lineName)!=1)
		{
			dcerror(3,"Could not find gene called %s in file %s\n",name,geneListFileName);
			return 0;
		}
		if (sscanf(geneLine,"%*d %*s %*s %*c %*d %*d %*d %*d %*d %*s %*s %*d %s",lineName)!=1)
		{
			dcerror(4,"Could not read enough parameters from this line: %s\n",geneLine);
			return 0;
		}
	} while (stricmp(name,lineName));
	fseek(geneListFile,startPos,SEEK_SET); // back to start of first line for this gene
	return 1;
}



int geneExtractor::downloadGene(refseqGeneInfo &g,char *outFn)
// download vcf file from 1000 genomes so need to specify chromosome number
{
	char tbiFn[200],firstHalf[200],secondHalf[200],*ptr;
	if (variantFileName[0]=='\0')
		setVariantFileName(DEFAULTVARIANTFILENAME);
	ptr=strstr(variantFileName,"chr");
	if (ptr==0)
	{
		dcerror(1,"Could not find string \"chr\" in filename %s",variantFileName);
		return 0;
	}
	strncpy(firstHalf,variantFileName,ptr-variantFileName);
	firstHalf[ptr-variantFileName]='\0';
	strcpy(secondHalf,ptr+3);
	sprintf(tbiFn,"%s%s%s",firstHalf,g.getChr(),secondHalf);
	return g.tbiExtractGene(tbiFn,outFn,0,0,0);
}

int geneExtractor::extractVariants(refseqGeneInfo &g,char *outFn,int appendToOld,int addChrInVCF,int removeSpaces)
// extract vcf file from local tbi file
{
	if (variantFileName[0]=='\0')
	{
		dcerror(1,"Need to give name for local vcf file in geneExtractor::extractVariants()");
		return 0;
	}
	else 
		return g.tbiExtractGene(variantFileName,outFn,appendToOld,addChrInVCF,removeSpaces);
}


