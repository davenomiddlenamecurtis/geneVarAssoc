#include <stdio.h>
#include <string.h>
#include "dcerror.hpp"

float two_x_two_chisq(float tab[2][2])
{
float ex[2][2],col_tot[2],row_tot[2],N;
int r,c;
float dchi;
N=dchi=0;
for (r=0;r<2;++r)
  row_tot[r]=0;
for (c=0;c<2;++c)
  col_tot[c]=0;
for (r=0;r<2;++r)
  for (c=0;c<2;++c)
    {
    row_tot[r]+=tab[r][c];
    col_tot[c]+=tab[r][c];
    N+=tab[r][c];
    }
if (row_tot==0 || col_tot==0)
  return 0;
for (r=0;r<2;++r)
  for (c=0;c<2;++c)
    {
    ex[r][c]=row_tot[r]*col_tot[c]/N;
    dchi+=(tab[r][c]-ex[r][c])*(tab[r][c]-ex[r][c])/ex[r][c];
    }
return dchi;
}

int main(int argc,char *argv[])
{
	FILE *fi,*fo,*fCc[2][10];
	int gc[2][3],chr,whichChr,nCc[2],oneGc[3],g,cc,f;
	long pos,oldPos;
	float chiThresh,chisq,tab[2][2];
	char fn[30],line[200],geneName[30],rest[200],consequence[50];
	fi=fopen(argv[1],"r");
	fo=fopen(argv[2],"w");
	for (cc=0;cc<2;++cc)
	{
		fgets(line,199,fi);
		nCc[cc]=0;
		*rest='\0';
		while (sscanf(line,"%s %[^\n]",fn,rest)>=1)
		{
			fCc[cc][nCc[cc]]=fopen(fn,"r");
			++nCc[cc];
			strcpy(line,rest);
			*rest='\0';
		}
	}
	fgets(line,199,fi);
	sscanf(line,"%f",&chiThresh);
	fclose(fi);
	while (1)
	{
		for (cc=0;cc<2;++cc)
			for (g=0;g<3;++g)
				gc[cc][g]=0;
		oldPos=-1L;
		for (cc=0;cc<2;++cc)
		{
			for (f=0;f<nCc[cc];++f)
			{
				if (!fgets(line,199,fCc[cc][f]))
				{
					if (oldPos==-1L)
						goto done;
					else
						dcerror(1,"Could not read line from fCc[%d][%d]\n",cc,f);
				}
			sscanf(line,"%*ld %s %d %ld %d %d %d %[^\n]",geneName,&chr,&pos,&oneGc[0],&oneGc[1],&oneGc[2],consequence);
			if (oldPos==-1L)
				oldPos=pos;
			else
				if (pos!=oldPos)
					dcerror(2,"Positions do not match: %ld and %ld in line from fCc[%d][%d]\n",pos,oldPos,cc,f);
			for (g=0;g<3;++g)
				gc[cc][g]+=oneGc[g];
			}
			if (chr==24) //Y
			{
				tab[cc][0]=gc[cc][0];
				tab[cc][1]=gc[cc][2];
			}
			else
			{
				tab[cc][0]=2*gc[cc][0]+gc[cc][1];
				tab[cc][1]=2*gc[cc][2]+gc[cc][1];
			}
		}
		chisq=two_x_two_chisq(tab);
		if (chisq>chiThresh)
		{
		fprintf(fo," %40s %2d %9ld ",geneName,chr,pos);
		for (cc=0;cc<2;++cc)
			for (g=0;g<3;++g)
				fprintf(fo,"%5d ",gc[cc][g]);
		fprintf(fo,"%5.1f %s\n",chisq,consequence);
		}

	}
done:
	fclose(fo);
	return 0;
}
