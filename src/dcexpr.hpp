/* DCEXPR.HPP */
/* Copyright Dave Curtis 1994 */
/* dcurtis@hgmp.mrc.ac.uk */
/* no warranty or liability of any kind is accepted, expressed or implied */

#ifndef DCEXPRESSHPP
#define DCEXPRESSHPP 1

#include "dcerror.hpp"
#include "dcstr.hpp"
#include <math.h>
#if 0
#include "statval.h"
#endif
#include <stdio.h>

class dcexpr_val {
public:
virtual int is_string_really()=0;
virtual operator char*()=0;
virtual operator double()=0;
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
~dcexpr_double();
dcexpr_double(double v=0.0);
};

class dcexpr_string:public dcexpr_val {
char *buff;
public:
virtual int is_string_really();
virtual operator char*();
virtual operator double();
~dcexpr_string();
dcexpr_string(char *t,int len=0);
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

extern char *vprimitive(char*,dcvnode **);
extern char *vun_op(char*,dcvnode **);
extern char *vbin_op(char*,dcvnode **,int);
extern char *vbracket(char*,dcvnode **);
extern int vcompute(char,dcvnode **,dcvnode **);
extern int un_vcompute(char,dcvnode **);
extern void add_un_op(char *lab,dcexpr_val *(*f)(dcvnode *));
extern void add_bin_op_same(char *lab,dcexpr_val *(*f)(dcvnode *,dcvnode *));
extern void add_bin_op_next(char *lab,dcexpr_val *(*f)(dcvnode *,dcvnode *));

class express {
protected:
dcvnode *head;
dcstring token;
char *vbin_op(char *s,dcvnode **br,int level);
char *vun_op(char*s,dcvnode **br);
char *vbracket(char*s,dcvnode **br);
virtual char *vprimitive(char *s,dcvnode **br);
int vcompute(char o,dcvnode **i,dcvnode **j);
int un_vcompute(char op,dcvnode **i);
virtual char *get_next(char *s);
public:
void wipe();
express();
virtual ~express(); // only virtual so no warning
dcexpr_val *eval();
int parse(char *s);
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

#endif

