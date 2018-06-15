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
