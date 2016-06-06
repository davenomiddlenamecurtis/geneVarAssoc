/* DCUTILS.HPP */ 
/* Copyright Dave Curtis 1994 */ 
/* dcurtis@hgmp.mrc.ac.uk */ 
/* no warranty or liability of any kind is accepted, expressed or implied */ 
 
#ifndef DCUTILSHPP
#define DCUTILSHPP

/* You may not need some, or indeed any, of these, depending on what
your compiler's standard library supplies. If you do #ifdef any out, then
you _must_ make sure to #ifdef out the corresponding function
definitions in dcutils.cxx, otherwise the linking will be all messed up.
*/


#include <ctype.h>

extern "C" {
#ifdef HGMPSOLARIS
extern int stricmp(char *s,char *t);   
extern int strnicmp(char *s,char *t,int l);   
extern char *strupr(char *s);
extern char *strlwr(char *s);
#endif
#if 0
#if !defined SUNCC && !defined __DECCXX && ! defined HGMPSOLARIS
extern void *memmove(void *d,void *s,unsigned int c);
#endif
#endif
};

#endif
