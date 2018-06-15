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

#ifndef DCEXPRESSHPP
#define DCEXPRESSHPP 1

#include "dcerror.hpp"
#include "dcstr.hpp"
#include <math.h>
#if 0
#include "statval.h"
#endif
#include <stdio.h>

#define TOKENMAXLEN 100

class dcexpr_val {
public:
virtual int is_string_really()=0;
virtual operator char*()=0;
virtual operator double()=0;
virtual dcexpr_val & operator= (dcexpr_val &v)=0;
// virtual ~dcexpr_val()=0;
// above doesn't seem to work
virtual ~dcexpr_val();
};

class dcexpr_double:public dcexpr_val {
double val;
char buff[50];
public:
virtual int is_string_really();
virtual operator char*();
virtual operator double();
virtual dcexpr_double & operator= (dcexpr_val &v);
dcexpr_double & operator=(double v) { val=v; return *this; }
// can make this virtual one day in dcexpr_val if ever needed
~dcexpr_double();
dcexpr_double(double v=0.0);
};

class dcexpr_string:public dcexpr_val {
char *buff;
public:
virtual int is_string_really();
virtual operator char*();
virtual operator double();
virtual dcexpr_string & operator= (dcexpr_val &val);
~dcexpr_string();
dcexpr_string(const char *t,int len=0);
};

class dcvnode {
char *str;
public:
int nbranches;
dcvnode *branch[4];
dcvnode(int nb);
virtual ~dcvnode(); // destroy correct object
dcvnode& operator=(dcvnode &old);
int add(dcvnode *vn,int b);
void wipe();
virtual dcexpr_val *eval()=0;
virtual dcvnode *copy()=0;
int matches(const char *s);
};

extern void add_un_op(char *lab,dcexpr_val *(*f)(dcvnode *));
extern void add_bin_op_same(char *lab,dcexpr_val *(*f)(dcvnode *,dcvnode *));
extern void add_bin_op_next(char *lab,dcexpr_val *(*f)(dcvnode *,dcvnode *));

class express {
protected:
dcvnode *head;
char token[TOKENMAXLEN ];
const char *vbin_op(const char *s,dcvnode **br,int level);
const char *vun_op(const char*s,dcvnode **br);
const char *vbracket(const char*s,dcvnode **br);
virtual const char *vprimitive(const char *s,dcvnode **br);
virtual const char *get_next(const char *s);
public:
void wipe();
express();
virtual ~express(); // only virtual so no warning
dcexpr_val *eval();
int parse(const char *s);
void debug(FILE *dbf) { debugFile=dbf; }
static FILE *debugFile;
};

class variable : public dcvnode {
public:
variable();
virtual ~variable();
dcexpr_val *eval()=0;
};

class vconstant : public dcvnode {
double value;
public:
vconstant(double v):dcvnode(0){ value=v;}
dcvnode *copy();
~vconstant();
dcexpr_val *eval();
};

class vstrconstant : public dcvnode {
char *value;
public:
vstrconstant(char *v);
dcvnode *copy();
~vstrconstant();
dcexpr_val *eval();
};

struct new_op_t { char *str; dcvnode *inst; };

extern char dc_expr_buff[];

#define MAX_N_OPS 30
#define MAX_OP_LEVEL 15
#define EVAL_R1 \
if ((r1=b1->eval())==NULL) return NULL; 

#define EVAL_R2 \
if ((r2=b2->eval())==NULL) return NULL; 

#define EVAL_BOTH \
if ((r1=b1->eval())==NULL) return NULL; \
if ((r2=b2->eval())==NULL) { delete r1; return NULL; }


#endif

