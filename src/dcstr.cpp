#include "dcstr.hpp"

int dcstring::assign(char *s)
{
if (p) delete [] p, p=0; 
if (s) return (p=dupstr(s))!=0; 
else return 1; 
}

int dcstring::assign(const dcstring &old)
{ 
if (p)
  delete [] p, p=0;
if (old.p)
  return ((p=dupstr(old.p))!=0);
else return 1;
}

char *dcstring::dupstr(const char *s)
{
char *dupstrptr;
if (s==0)
  dupstrptr=0;
else
  {
  dupstrptr=new char[strlen(s)+1];
  if (dupstrptr!=NULL)
    strcpy(dupstrptr,s);
  }
return(dupstrptr);
}

dcstring::dcstring(const char *s)
{
p=0; 
if (s)
  p=dupstr(s); 
}

dcstring::dcstring(const dcstring &old) 
{ 
p=0; 
if (old.p) 
  p=dupstr(old.p); 
}

dcstring::dcstring() 
{ 
p=0; 
}

