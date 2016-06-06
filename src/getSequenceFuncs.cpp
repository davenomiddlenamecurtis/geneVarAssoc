#include <ctype.h>
#include <string.h>
#include "getSequence.hpp"
#include "dcerror.hpp"

aminoAcid allAminoAcids[] =
{
{"Isoleucine",'I',"ATT"},
{"Isoleucine",'I',"ATC"},
{"Isoleucine",'I',"ATA"},
{"Leucine",'L',	"CTT"},
{"Leucine",'L',"CTC"},
{"Leucine",'L',"CTA"},
{"Leucine",'L',"CTG"},
{"Leucine",'L',"TTA"},
{"Leucine",'L',"TTG"},
{"Valine",'V',"GTT"},
{"Valine",'V',"GTC"},
{"Valine",'V',"GTA"},
{"Valine",'V',"GTG"},
{"Phenylalanine",'F',"TTT"},
{"Phenylalanine",'F',"TTC"},
{"Methionine",'M',"ATG"},
{"Cysteine",'C',"TGT"},
{"Cysteine",'C',"TGC"},
{"Alanine",'A',"GCT"},
{"Alanine",'A',"GCC"},
{"Alanine",'A',"GCA"},
{"Alanine",'A',"GCG"},
{"Glycine",'G',"GGT"},
{"Glycine",'G',"GGC"},
{"Glycine",'G',"GGA"},
{"Glycine",'G',"GGG"},
{"Proline",'P',"CCT"},
{"Proline",'P',"CCC"},
{"Proline",'P',"CCA"},
{"Proline",'P',"CCG"},
{"Threonine",'T',"ACT"},
{"Threonine",'T',"ACC"},
{"Threonine",'T',"ACA"},
{"Threonine",'T',"ACG"},
{"Serine",'S',"TCT"},
{"Serine",'S',"TCC"},
{"Serine",'S',"TCA"},
{"Serine",'S',"TCG"},
{"Serine",'S',"AGT"},
{"Serine",'S',"AGC"},
{"Tyrosine",'Y',"TAT"},
{"Tyrosine",'Y',"TAC"},
{"Tryptophan",'W',"TGG"},
{"Glutamine",'Q',"CAA"},
{"Glutamine",'Q',"CAG"},
{"Asparagine",'N',"AAT"},
{"Asparagine",'N',"AAC"},
{"Histidine",'H',"CAT"},
{"Histidine",'H',"CAC"},
{"GlutamicAcid",'E',"GAA"},
{"GlutamicAcid",'E',"GAG"},
{"AsparticAcid",'D',"GAT"},
{"AsparticAcid",'D',"GAC"},
{"Lysine",'K',"AAA"},
{"Lysine",'K',"AAG"},
{"Arginine",'R',"CGT"},
{"Arginine",'R',"CGC"},
{"Arginine",'R',"CGA"},
{"Arginine",'R',"CGG"},
{"Arginine",'R',"AGA"},
{"Arginine",'R',"AGG"},
{"Stop",'*',"TAA"},
{"Stop",'*',"TAG"},
{"Stop",'*',"TGA"}
};

char codonReader::table[4][4][4];
int codonReader::baseLookup[256];
char *codonReader::nameTable[256];

char compBases[5][2]=
{
	{'c','g'},
	{'a','t'},
	{'C','G'},
	{'A','T'},
	{'*','*'}
};

char getCompBase(char b)
{
	int i,j;
	for (i=0;i<5;++i)
		for (j=0;j<2;++j)
			if (b==compBases[i][j])
				return compBases[i][(j+1)%2];
	return '?';
}

char codonReader::translate(char codon[3],char strand)
{
	int i;
	char code[3];
	if (strand=='+')
		for (i=0;i<3;++i)
			code[i]=codon[i];
	else
	{
		for (i=0;i<3;++i)
		{
			switch(codon[2-i])
			{
			case 'A':
				code[i]='T';
				break;
			case 'T':
				code[i]='A';
				break;
			case 'C':
				code[i]='G';
				break;
			case 'G':
				code[i]='C';
				break;
			case 'a':
				code[i]='t';
				break;
			case 't':
				code[i]='a';
				break;
			case 'c':
				code[i]='g';
				break;
			case 'g':
				code[i]='c';
				break;
			}
		}
		for (i=0;i<3;++i)
			codon[i]=code[i];
	}
	return table[baseLookup[code[0]]][baseLookup[code[1]]][baseLookup[code[2]]];
}

codonReader::codonReader()
{
	int i,j,k;
	aminoAcid *aa;
	baseLookup['a']=baseLookup['A']=0;
	baseLookup['c']=baseLookup['C']=1;
	baseLookup['g']=baseLookup['G']=2;
	baseLookup['t']=baseLookup['T']=3;
	for (i=0;i<4;++i)
		for (j=0;j<4;++j)
			for (k=0;k<4;++k)
			{
				aa=&allAminoAcids[i*16+j*4+k];
				table[baseLookup[aa->code[0]]][baseLookup[aa->code[1]]][baseLookup[aa->code[2]]]=aa->letter;
				nameTable[aa->letter]=aa->name;
			}

}

int faSequenceFile::init(char *fn)
{
	int c;
	// idea is to use ftell to get position so it does not matter whether fgetc counts \n as 1 or 2 characters
	long lineEnd,lineStart;
	if (fp!=0)
		fclose (fp);
	fp=fopen(fn,"rb");
	if (fp==0)
	{
		dcerror(1,"Could not open file %s\n",fn);
		return 0;
	}
	while (!isspace(c=fgetc(fp))) ; // skip first line
	do {
		startPos=ftell(fp);
	} while (isspace(c=fgetc(fp)));
	lineStart=startPos;
	do {
		lineEnd=ftell(fp);
	} while (!isspace(c=fgetc(fp)));
	lineLength=lineEnd-lineStart;
	do {
		lineStart=ftell(fp);
	} while (isspace(c=fgetc(fp)));
	gapLength=lineStart-lineEnd;
	// first base has position 1
	return 1;
}

int faSequenceFile::getSequence(char *s,long ePos,long len)
{
	long inLine,lineNo;
	inLine=(ePos-1L)%lineLength;
	lineNo=(ePos-1L)/lineLength;
	fseek(fp,startPos+lineNo*(lineLength+gapLength)+inLine,SEEK_SET);
	for (;len>0;--len)
	{
		*s++=fgetc(fp);
		if (++inLine==lineLength)
		{
			inLine=0;
			++lineNo;
			fseek(fp,startPos+lineNo*(lineLength+gapLength)+inLine,SEEK_SET);
		}
	}
	*s='\0';
	return 1;
}

faSequenceFile::faSequenceFile()
{
	fp=0;
}

void getComplementarySequence(char *t,char *s)
{
	int i,j;
	char ch;
	for (i=0,j=strlen(s)-1;i<strlen(s);++i,--j)
	{
		switch (s[i])
			{
			case 'A':
				ch='T';
				break;
			case 'T':
				ch='A';
				break;
			case 'C':
				ch='G';
				break;
			case 'G':
				ch='C';
				break;
			case 'a':
				ch='t';
				break;
			case 't':
				ch='a';
				break;
			case 'c':
				ch='g';
				break;
			case 'g':
				ch='c';
				break;
			}
		t[j]=ch;
	}
t[i]='\0';
}

char *RNA2DNA(char *seq)
{
	char *ptr;
	for (ptr=seq;*ptr;++ptr)
	{
		switch (*ptr)
		{
			case 'u':
				*ptr='t';
				break;
			case 'U':
				*ptr='T';
				break;
			default:
				break;
		}
	}
return seq;
}
