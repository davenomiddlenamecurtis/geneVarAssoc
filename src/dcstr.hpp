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

#ifndef DCSTRHPP
#define DCSTRHPP 1

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef save
#undef save
#endif

class dcstring {
protected:
char *dupstr(const char *s);
char *p;
public:
dcstring(const char *s);
dcstring(const dcstring &old);
dcstring();
~dcstring() { if (p) delete (p); }
char& operator[](int i)
{ 
return p[i]; 
}
dcstring& operator=(const char *s)
{
if (p) 
  delete (p),p=0; 
if (s)
  p=dupstr(s);
return *this; 
}
dcstring& operator=(const dcstring &old)
{
if (p)
  delete (p),p=0;
if (old.p)
  p=dupstr(old.p); 
return *this;
}
int assign(char *s);
int assign(const dcstring &old);
int operator==(const dcstring &old)
  { return (p==0)!=(old.p==0)?0:p==0?1:!strcmp(p,old.p); }
int operator!=(const dcstring &old)
  { return (p==0)!=(old.p==0)?1:p==0?0:strcmp(p,old.p); }
void strfree() { if (p) { delete (p); p=0; } }
operator char* () const { return p; }
void move(dcstring &old) { if (p) delete (p); p=old.p; old.p=0; }
int is_null() { return p==0; }
void set_null() { p=0; }
int save(FILE *fp) { short len=is_null()?0:strlen(p)+1; 
return fwrite((char*)&len,sizeof(short),1,fp) && 
(is_null() || fwrite(p,len,1,fp)); }
int load(FILE *fp) { short len; strfree();
return fread((char*)&len,sizeof(short),1,fp) && 
(len==0 || ((p=new char[len])!=0 && fread(p,len,1,fp))); }
};

#endif 


