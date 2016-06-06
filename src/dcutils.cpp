/* DCUTILS.CXX */ 
/* Copyright Dave Curtis 1994 */ 
/* dcurtis@hgmp.mrc.ac.uk */ 
/* no warranty or liability of any kind is accepted, expressed or implied */ 
 
#if !_WINDOWS
#include <string.h>
#include <ctype.h>
#include "dcstr.hpp"

// !!! WATCH OUT !!! These are quick and dirty and can break!!!

extern "C" { 

#if 0
#if !defined SUNCC && !defined __DECCXX && ! defined HGMPSOLARIS
void *memmove(void *d,void *s,unsigned c)
{
size_t i;
char *dc=(char*)d,*sc=(char*)s;
if (s>d)
  for (i=0;i<c;++i) *dc++=*sc++;
else
  for (i=c-1;i>=0;++i) dc[i]=sc[i];
return d;
}
#endif
#endif

#ifdef HGMPSOLARIS

char *strupr(char *s)
{
char *ptr;
ptr=s;
while (*ptr)
  {
  if (isalpha(*ptr) && islower(*ptr)) *ptr=toupper(*ptr);
  ++ptr;
  }
return s;
}

char *strlwr(char *s)
{
char *ptr;
ptr=s;
while (*ptr)
  {
  if (isalpha(*ptr) && isupper(*ptr)) *ptr=tolower(*ptr);
  ++ptr;
  }
return s;
}

int stricmp(char *s,char *t)
{
dcstring sbuff(s),tbuff(t);
strupr(sbuff);
strupr(tbuff);
return strcmp(sbuff,tbuff);
}

int strnicmp(char *s,char *t,int l)
{
dcstring sbuff(s),tbuff(t);
strupr(sbuff);
strupr(tbuff);
return strncmp(sbuff,tbuff,l);
}

#endif

};

#endif

