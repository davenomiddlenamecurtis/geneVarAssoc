/* DCINDEX.CXX */ 
/* Copyright Dave Curtis 1994-2012 */ 
/* david.curtis@qmul.ac.uk */ 
/* no warranty or liability of any kind is accepted, expressed or implied */ 

/* this index will return 0L if it does not find a matching record */
 
#include <fcntl.h> // for  O_CREAT etc.
#ifndef MSDOS
#include <unistd.h>
#else
#include <io.h> // for unlink()
#include <dos.h> // for  O_CREAT etc.
#endif

#ifndef _O_BINARY
#define _O_BINARY 0
// if not defined yet then hopefully we do not need it - makes sure index files are binary
#endif

#include <string.h>
#include "dcerror.hpp"
#include "dcindex.hpp"

extern "C" {
struct btree *btopen(char *path, int flags, int mode);
int btclose(struct btree *bt);
int btinsert(char *argkey, long recnbr, struct btree *bt);
int btdelete(long node_nbr, struct btree *bt);
int btnext(long *node_nbr,struct btnode *cno,struct	btree *bt);
int btprevious(long *node_nbr,struct btnode *cno,struct btree *bt);
int bttail(long *node_nbr,struct btnode *cno,struct btree *bt);
int bthead(long *node_nbr,struct btnode *cno,struct btree *bt);
int btfind(char *key1,long *node_nbr,struct btnode *cno,struct btree *bt);
};

int dc_index::add(char *key,long rec)
{
if (!h1) { dcerror(1,"index file is not open"); return 0; }
char add_buff[BT_KSIZ];
strncpy(add_buff,key, BT_KSIZ-1);
add_buff[BT_KSIZ - 1] = '\0';
int rv=btinsert(add_buff,rec,h1);
if (rv<0) { dcerror(1,"btree error in add()"); return 0; }
else return 1;
}

long dc_index::get_last()
{
if (!h1) { dcerror(1,"index file is not open"); return 0L; }
int rv=bttail(&node_nbr,&cnode,h1);
if (rv==-1) { dcerror(1,"btree() error in get_last()"); return 0L; }
else return rv==0?cnode.recno:0L;
}

int dc_index::remove()
{
if (!h1) { dcerror(1,"index file is not open"); return 0L; }
int rv=btdelete(node_nbr,h1);
if (rv<0) { dcerror(1,"btree() error in remove()"); return 0L; }
else return 1;
}

long dc_index::get_next()
{
if (!h1) { dcerror(1,"index file is not open"); return 0L; }
int rv=btnext(&node_nbr,&cnode,h1);
if (rv==-1) 
  { 
  dcerror(1,"btree() error in get_next()"); 
  return 0L; 
  }
else if (rv==BT_EOF) return 0l; // at end
// there seems to be an error in the btree routines so that when > 50,000 returns recPos=-1
// need to fix this sometime
else if (cnode.recno<0L) return 0l; // pretend we got to the end even though really an error
else return cnode.recno;
}

long dc_index::get_prev()
{
int rv=btprevious(&node_nbr,&cnode,h1);
if (rv==-1) { dcerror(1,"btree() error in get_prev()"); return 0L; }
else if (rv==BT_EOF) return 0l; // at start
else return cnode.recno;
}

long dc_index::get_first()
{
int rv=bthead(&node_nbr,&cnode,h1);
if (rv==-1) { dcerror(1,"btree() error in get_first()"); return 0L; }
else return rv==0?cnode.recno:0L;
}

int dc_index::find_matching_node(char *key,long rec)
{
// find first exact match, then try going backwards and forwards from it
// till hit right record number
int rv;
rv=btfind(key,&node_nbr,&cnode,h1);
if (rv==BT_NREC) 
  { return 0; }
else if (rv==-1) { dcerror(1,"btree error in find_matching_node()"); return 0; }
do
   { 
   if (cnode.recno==rec) return 1;
   int rrv=btnext(&node_nbr,&cnode,h1);
   if (rrv==-1) { dcerror(1,"btree() error in find_matching_node()"); return 0; }
   else if (rrv==BT_EOF) break; // at end
   } while (!strcmp(key,cnode.key));
rv=btfind(key,&node_nbr,&cnode,h1);
if (rv==BT_NREC) 
  { return 0; }
else if (rv==-1) { dcerror(1,"btree error in find_matching_node()"); return 0; }
do
   { 
   if (cnode.recno==rec) return 1;
   int rrv=btprevious(&node_nbr,&cnode,h1);
   if (rrv==-1) { dcerror(1,"btree() error in find_matching_node()"); return 0; }
   else if (rrv==BT_EOF) break; // at start
   } while (!strcmp(key,cnode.key));
return 0;
}


long dc_index::exact_find(char *key)
{
if (!h1) { dcerror(1,"index file is not open"); return 0L; }
int rv=btfind(key,&node_nbr,&cnode,h1);
if (rv==BT_NREC) 
  { 
#if 0
  dcerror(1,"key \"%s\" not found in exact_find()",key); 
  // not necessarily a real error
#endif
  return 0L; 
  }
else if (rv==-1) { dcerror(1,"btree error in exact_find()"); return 0L; }
else return cnode.recno;
}

long dc_index::near_find(char *key)
{
if (!h1) { dcerror(1,"index file is not open"); return 0L; }
int rv=btfind(key,&node_nbr,&cnode,h1);
if (rv==-1) { dcerror(1,"btree error in near_find()"); return 0L; }
else return cnode.recno;
}

int dc_index::make_new(char *name)
{
if (h1) close();
::unlink(name);
node_nbr=0L;
if ((h1 = btopen(name,O_CREAT|O_RDWR|_O_BINARY,0600)) ==NULL)
  { dcerror(1,"Could not create index file %s",name); return 0; }
else return 1;
}

int dc_index::open_old(char *name)
{
if (h1) close();
node_nbr=0L;
if ((h1 = btopen(name,O_RDWR|_O_BINARY,0600)) ==NULL)
  { 
  dcerror(1,"Could not open index file %s",name); 
  return 0; 
  }
else
  {
  get_first();
  return 1;
  }
}

dc_index::~dc_index()
{
close();
}

void dc_index::close()
{
if (h1) btclose(h1);
h1=NULL;
}

dc_index::dc_index()
{ 
h1=NULL; 
cnode.key[0]='\0';
}

