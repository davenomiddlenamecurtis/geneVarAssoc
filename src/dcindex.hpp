/* DCINDEX.HPP */ 
/* Copyright Dave Curtis 1994 */ 
/* dcurtis@hgmp.mrc.ac.uk */ 
/* no warranty or liability of any kind is accepted, expressed or implied */ 
 
#ifndef DCINDEXHPP
#define DCINDEXHPP 1

#include "btree.h"

class dc_index {
BTREE *h1;
long node_nbr;
BTNODE cnode;
public:
int dump(char *fn);
int is_open() { return h1!=0; }
int open_old(char *name);
int make_new(char *name);
void close();
int add(char *key,long rec);
int remove(); // can only remove current item from btree
long get_first();
long get_last();
long get_prev();
long get_next();
long exact_find(char *key);
long near_find(char *key);
int find_matching_node(char *key,long rec);
const char *current_key() { return cnode.key; }
dc_index();
~dc_index();
};

#endif
