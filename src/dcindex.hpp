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
