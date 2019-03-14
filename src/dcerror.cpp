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

#include "dcerror.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef USEDCASSERT
// catch assertion failures so Windows does not do something clever with them
#ifdef  __cplusplus
extern "C" {
#endif

_CRTIMP void __cdecl _assert(void *str, void *fn, unsigned line)
{
fprintf(stderr,"Assertion failed: %s, file %s, line %d\n",(char *)str,(char *)fn,line);
exit(1);
}

#ifdef  __cplusplus
}
#endif
#endif

int default_error_func(const char *format,va_list arg_ptr)
{
#if 0
char s[2000];
vsprintf(s,format,arg_ptr);
fprintf(stderr,"%s",s);
#else
vfprintf(stderr,format,arg_ptr);
#endif
// system("pause");
return 0;
}
 
int error_object::operator()(int e, const char *format,...)
  {
  va_list arg_ptr;
  if (!format)
#ifndef SUN_CC
    return operator()(e,"%s\n",strerror(e));
#else
    return operator()(e,"unknown type of error\n");
#endif
  status=e;
  va_start(arg_ptr,format);
  if (showing)
    {
    if (display_func) 
      display_func(format,arg_ptr);
    else
      {
      default_error_func(format,arg_ptr);
      }
    }
  va_end(arg_ptr);
  if (fatal)
	  exit(1);
  return 0;
  }

error_object::error_object() { showing=fatal=1; status=0; display_func=0;}

error_object::~error_object() {;}
void error_object::set_display(void (*df)(const char *,va_list)) { display_func=df; }
int error_object::stat() { return status; }
void error_object::clear() { status=0; }
void error_object::hide() { showing=0; }
void error_object::show() { showing=1; }
void error_object::warn() { fatal=0; }
void error_object::kill() { fatal=1; }

error_object dcerror;

void default_error(int n,char *s)
  { dcerror(n,s); }

