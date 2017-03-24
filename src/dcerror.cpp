/* DCERROR.CXX */ 
/* Copyright Dave Curtis 1994-2012 */ 
/* david.curtis@qmul.ac.uk */ 
/* no warranty or liability of any kind is accepted, expressed or implied */ 
 
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

int default_error_func(char *format,va_list arg_ptr)
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

int error_object::operator()(int e,char *format,...)
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
void error_object::set_display(void (*df)(char *,va_list)) { display_func=df; }
int error_object::stat() { return status; }
void error_object::clear() { status=0; }
void error_object::hide() { showing=0; }
void error_object::show() { showing=1; }
void error_object::warn() { fatal=0; }
void error_object::kill() { fatal=1; }

error_object dcerror;

void default_error(int n,char *s)
  { dcerror(n,s); }

