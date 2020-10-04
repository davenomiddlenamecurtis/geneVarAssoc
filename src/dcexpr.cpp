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

#include "dcexpr.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#ifndef MSDOS
#define stricmp strcasecmp
#define strnicmp strncasecmp
#define strupr mystrupr
#include <ctype.h>
static char *mystrupr(char *s)
{
	char *t;
	t=s;
	while (*t)
	{
		if (isalpha(*t))
			*t=toupper(*t);
		++t;
	}
	return s;
}
#endif

struct new_op_t un_ch[MAX_N_OPS],bin_ch[MAX_N_OPS][MAX_OP_LEVEL];
char *all_op[MAX_N_OPS*2];
int n_un_ops,n_bin_ops[MAX_OP_LEVEL],n_levels,n_all_ops,last_op_level;
// for initialisers to work, all these must be zero
// hopefully they should be

FILE *express::debugFile;

express::express()
{ head=NULL; token[0]='\0';  }

dcexpr_val *express::eval()
{ 
dcexpr_val *rv;
rv=head->eval(); 
if (rv==NULL) dcerror(ENOMEM,"Error evaluating expression - probably out of memory");
return rv;
}

dcexpr_val::~dcexpr_val() {;}

int dcexpr_double::is_string_really() { return 0; }
dcexpr_double::operator char*() 
{ 
	char *ptr;
	sprintf(buff,"%f",val); // default is to 6 decimal places
	ptr=buff;
	while (*ptr)
		++ptr;
	while (*(--ptr)=='0') // stripping trailing zeros and decimal point
		*ptr='\0';
	if (*ptr=='.')
		*ptr='\0';
	return buff;
}
dcexpr_double::operator double() { return val; }
dcexpr_double::~dcexpr_double() {;}
dcexpr_double::dcexpr_double(double v) { val=v; }
dcexpr_double & dcexpr_double::operator= (dcexpr_val &v)
{ val=(double)v; return *this; }


int dcexpr_string::is_string_really() { return 1; }
dcexpr_string::operator char*() { return buff; }
dcexpr_string::operator double() { return atof(buff); }
dcexpr_string::~dcexpr_string() { if (buff) delete buff; }
dcexpr_string::dcexpr_string(const char *t,int len) 
{ 
if (len)
  {
  if ((buff=new char[len+1])!=NULL)
    {
    strncpy(buff,t,len);
    buff[len]='\0';
    }
  else dcerror(ENOMEM);
  }
else
  {
  if ((buff=new char[strlen(t)+1])!=NULL)
    strcpy(buff,t);
  else dcerror(ENOMEM); 
  }
}

dcexpr_string & dcexpr_string::operator= (dcexpr_val &val)
{
	char *s=(char*)val;
	if (buff)
		delete buff;
	buff=new char[strlen(s)+1];
	strcpy(buff,s);
	return *this;
}

char dc_expr_buff[2000];

class vnvbin_op:public dcvnode {
dcexpr_val *(*func)(dcvnode *,dcvnode *);
public:
vnvbin_op(dcexpr_val *(*f)(dcvnode *,dcvnode *));
~vnvbin_op(){;}
dcexpr_val *eval();
dcvnode *copy();
};

dcvnode *vnvbin_op::copy()
{
return new vnvbin_op(func);
}

vnvbin_op::vnvbin_op(dcexpr_val *(*f)(dcvnode *,dcvnode *)) : dcvnode(2)
{
func=f;
}

dcexpr_val *vnvbin_op::eval()
{
if (branch[0]==NULL || branch[1]==NULL)
  return 0;
dcexpr_val *rv=func(branch[0],branch[1]);
if (express::debugFile)
	fprintf(express::debugFile,"%s\t",(char*)(*rv));
return rv;
}

class vnvun_op:public dcvnode {
dcexpr_val *(*func)(dcvnode *);
public:
vnvun_op(dcexpr_val *(*f)(dcvnode *));
~vnvun_op(){;}
dcexpr_val *eval();
dcvnode *copy();
};

dcvnode *vnvun_op::copy()
{
return new vnvun_op(func);
}

vnvun_op::vnvun_op(dcexpr_val *(*f)(dcvnode *)) : dcvnode(1)
{
func=f;
}

dcexpr_val *vnvun_op::eval()
{
if (branch[0]==NULL) return NULL;
dcexpr_val *rv=func(branch[0]);
if (express::debugFile)
	fprintf(express::debugFile,"%s\t",(char*)(*rv));
return rv;
}

variable::~variable() {;}

variable::variable() : dcvnode(0) {;}

dcvnode *vconstant::copy() { return new vconstant(value); }
vconstant::~vconstant(){;}
dcexpr_val *vconstant::eval() { return new dcexpr_double(value); }

dcvnode *vstrconstant::copy() { return new vstrconstant(value); }
vstrconstant::vstrconstant(char *s) : dcvnode(0)
  { if ((value=new char[strlen(s)+1])!=NULL) strcpy(value,s); }
vstrconstant::~vstrconstant() { if (value) delete value;}

dcexpr_val *vstrconstant::eval() { return new dcexpr_string(value); }

dcexpr_val *iftrue_op(dcvnode *b1,dcvnode *b2)
{
dcexpr_val *r1;
EVAL_R1;
double test=double(*r1);
delete r1;
if (test) return b2->eval();
else return new dcexpr_string("");
}

dcexpr_val *mult_op(dcvnode *b1,dcvnode *b2)
{
dcexpr_val *r1,*r2;
EVAL_BOTH;
double test=double(*r1)* double(*r2);
delete r1;
delete r2;
return new dcexpr_double(test);
}

dcexpr_val *or_op(dcvnode *b1,dcvnode *b2)
{
dcexpr_val *r1,*r2;
EVAL_R1;
double test=double(*r1)!=0;
delete r1;
if (!test) 
  {
  EVAL_R2;
  test=double(*r2)!=0;
  delete r2;
  }
return new dcexpr_double(test);
}

dcexpr_val *and_op(dcvnode *b1,dcvnode *b2)
{
dcexpr_val *r1,*r2;
EVAL_R1;
double test=double(*r1)!=0;
delete r1;
if (test)
  {
  EVAL_R2;
  test=double(*r2)!=0;
  }
return new dcexpr_double(test);
}

dcexpr_val *div_op(dcvnode *b1,dcvnode *b2)
{
dcexpr_val *r1,*r2;
EVAL_BOTH;
double rv=double(*r1);
if (rv)
	rv=rv/double(*r2); // so 0/0 is 0 rather than NaN
delete r1; delete r2;
return new dcexpr_double(rv);
}

dcexpr_val *comma_op(dcvnode *b1,dcvnode *b2)
{
dcexpr_val *r1,*r2;
EVAL_BOTH;
delete r1;
return r2;
}

dcexpr_val *mod_op(dcvnode *b1,dcvnode *b2)
{
dcexpr_val *r1,*r2;
EVAL_BOTH;
int t1=(int)double(*r1),t2=(int)double(*r2),m;
delete r1; delete r2;
if (t2<0) t2= -t2;
return new dcexpr_double(t2?((m=t1%t2)>=0)?m:m+t2 : 0.0);
}
// x mod 0 = 0

dcexpr_val *add_op(dcvnode *b1,dcvnode *b2)
{
dcexpr_val *r1,*r2;
EVAL_BOTH;
double rv=double(*r1) + double(*r2);
delete r1; delete r2;
return new dcexpr_double(rv);
}


dcexpr_val *sub_op(dcvnode *b1,dcvnode *b2)
{
dcexpr_val *r1,*r2;
EVAL_BOTH;
double rv=double(*r1) - double(*r2);
delete r1; delete r2;
return new dcexpr_double(rv);
}

dcexpr_val *lt_op(dcvnode *b1,dcvnode *b2)
{
dcexpr_val *r1,*r2;
EVAL_BOTH;
double rv=double(*r1) < double(*r2);
delete r1; delete r2;
return new dcexpr_double(rv);
}

dcexpr_val *lte_op(dcvnode *b1,dcvnode *b2)
{
dcexpr_val *r1,*r2;
EVAL_BOTH;
double rv=double(*r1) <= double(*r2);
delete r1; delete r2;
return new dcexpr_double(rv);
}

dcexpr_val *gt_op(dcvnode *b1,dcvnode *b2)
{
dcexpr_val *r1,*r2;
EVAL_BOTH;
double rv=double(*r1) > double(*r2);
delete r1; delete r2;
return new dcexpr_double(rv);
}

dcexpr_val *gte_op(dcvnode *b1,dcvnode *b2)
{
dcexpr_val *r1,*r2;
EVAL_BOTH;
double rv=double(*r1) >= double(*r2);
delete r1; delete r2;
return new dcexpr_double(rv);
}

dcexpr_val *equ_op(dcvnode *b1,dcvnode *b2)
{ 
dcexpr_val *r1,*r2;
EVAL_BOTH;
double rv,rv1,rv2;
if (r1->is_string_really() && r2->is_string_really())
  rv=!strcmp((char*)(*r1),(char*)(*r2));
else
  rv=(rv1=double(*r1)) == (rv2=double(*r2));
delete r1; delete r2;
return new dcexpr_double(rv);
}

dcexpr_val *neq_op(dcvnode *b1,dcvnode *b2)
{ 
double rv;
dcexpr_val *r1;
if ((r1=equ_op(b1,b2))==NULL) return NULL;
rv=double(*r1)==0.0;
delete r1;
return new dcexpr_double(rv);
}

dcexpr_val *pow_op(dcvnode *b1,dcvnode *b2)
{
dcexpr_val *r1,*r2;
EVAL_BOTH;
double rv=pow(double(*r1),double(*r2));
delete r1; delete r2;
return new dcexpr_double(rv);
}

dcexpr_val *abs_op(dcvnode *b1)
{
dcexpr_val *r1;
EVAL_R1;
double rv=fabs(double(*r1));
delete r1;
return new dcexpr_double(rv);
}

dcexpr_val *negate_op(dcvnode *b1)
{
dcexpr_val *r1;
EVAL_R1;
double rv=0.0-double(*r1);
delete r1;
return new dcexpr_double(rv);
}

dcexpr_val *exp_op(dcvnode *b1)
{
dcexpr_val *r1;
EVAL_R1;
double rv=exp(double(*r1));
delete r1;
return new dcexpr_double(rv);
}

dcexpr_val *not_op(dcvnode *b1)
{
dcexpr_val *r1;
EVAL_R1;
double rv=(double(*r1)==0.0);
delete r1;
return new dcexpr_double(rv);
}

dcexpr_val *ln_op(dcvnode *b1)
{
dcexpr_val *r1;
EVAL_R1;
double rv=log(double(*r1));
delete r1;
return new dcexpr_double(rv);
}

dcexpr_val *log_op(dcvnode *b1)
{
dcexpr_val *r1;
EVAL_R1;
double rv=log10(double(*r1));
delete r1;
return new dcexpr_double(rv);
}

dcexpr_val *sin_op(dcvnode *b1)
{
dcexpr_val *r1;
EVAL_R1;
double rv=sin(double(*r1));
delete r1;
return new dcexpr_double(rv);
}

dcexpr_val *arcsin_op(dcvnode *b1)
{
dcexpr_val *r1;
EVAL_R1;
double rv=asin(double(*r1));
delete r1;
return new dcexpr_double(rv);
}

dcexpr_val *cos_op(dcvnode *b1)
{
dcexpr_val *r1;
EVAL_R1;
double rv=cos(double(*r1));
delete r1;
return new dcexpr_double(rv);
}

dcexpr_val *arccos_op(dcvnode *b1)
{
dcexpr_val *r1;
EVAL_R1;
double rv=acos(double(*r1));
delete r1;
return new dcexpr_double(rv);
}

dcexpr_val *tan_op(dcvnode *b1)
{
dcexpr_val *r1;
EVAL_R1;
double rv=tan(double(*r1));
delete r1;
return new dcexpr_double(rv);
}

dcexpr_val *arctan_op(dcvnode *b1)
{
dcexpr_val *r1;
EVAL_R1;
double rv=atan(double(*r1));
delete r1;
return new dcexpr_double(rv);
}

dcvnode& dcvnode::operator=(dcvnode &old)
{
for (int i=0;i<nbranches;++i)
  if (old.branch[i])
    {
    if ((branch[i]=old.branch[i]->copy())!=NULL)
      *branch[i]=*old.branch[i];
    }
return *this;
}

dcvnode::dcvnode(int nb)
 {
 nbranches=nb;
 for (int i=0;i<nbranches;++i) branch[i]=NULL;
 }

dcvnode::~dcvnode()
 {
 wipe();
 }

void dcvnode::wipe()
 {
 for (int i=0;i<nbranches;++i) if (branch[i]!=NULL) delete branch[i];
 }

int dcvnode::add(dcvnode *vn,int b)
 {
 if (b>=nbranches)
   return dcerror(1,"Syntax error"); // Was error 3, blank line ??!?
 branch[b]=vn;
 return 1;
 }


express::~express()
 {
 wipe();
 }
 
void express::wipe()
{
 if (head!=NULL)
   {
   delete head;
   head=NULL;
   }
 }

int express::parse(const char *s)
 {
 const char *ptr;
 if ((s=get_next(s))==NULL) return 0;
 if (token[0]=='\0') 
  {
  if (*s) return dcerror(2,"No recognisable tokens in string: %s",s);
  else return dcerror(2,"Blank line");
  }
 if (head!=NULL)
   {
   delete head;
   head=NULL;
   }
 if ((s=vbin_op(s,&head,0))==NULL) return 0;
 if (token[0]!='\0') return dcerror(3,"Syntax error: %s",(char*)token);
 else for (ptr=s;*ptr;++ptr)
  if (!isspace(*ptr))
    return dcerror(3,"Syntax error: %s",s);
 return 1;
 }

const char *express::vbin_op(const char *s,dcvnode **br,int level)
 {
 dcvnode *tempbr;
 if (level<=n_levels)
      {
      if ((s=vbin_op(s,br,level+1))==NULL) return NULL;
      for (int i=0;i<n_bin_ops[level];++i)
         if (!stricmp(bin_ch[level][i].str,token))
           {
           if ((s=get_next(s))==NULL) return NULL;
           if ((s=vbin_op(s,&tempbr,level+1))==NULL) return NULL;
           dcvnode *temp;
           if ((temp=bin_ch[level][i].inst->copy())==NULL) return NULL;
           temp->add(*br,0);
           temp->add(tempbr,1);
           *br=temp;
           i=-1;
           }
       }
  else
     s=vun_op(s,br);
  return s;
  }

const char *express::vun_op(const char*s,dcvnode **br)
    {
    int i;
    int op=0;
    for (i=0;i<n_un_ops;++i)
     if (!stricmp(un_ch[i].str,token))
       {
       op=1;
       if ((s=get_next(s))==NULL) return NULL;
       if ((s=vun_op(s,br))==NULL) return NULL;
       break;
       }
    if (!op) 
	  if ((s=vbracket(s,br))==NULL) 
	    return NULL;  // deepest level
    if (op)
       {
       dcvnode *temp;
       if ((temp=un_ch[i].inst->copy())==NULL) return NULL;
       temp->add(*br,0);
       *br=temp;
       }
    return s;
    }

const char *express::vbracket(const char*s,dcvnode **br)
    {
    if (!strcmp("(",token))
     {
     if ((s=get_next(s))==NULL) return NULL;
     if ((s=vbin_op(s,br,0))==NULL) return NULL;
     if (strcmp(")",token))
	 { dcerror(1,"expected \')\': %s",(char *)token); return NULL; }
     if ((s=get_next(s))==NULL) return NULL;
     }
    else s=vprimitive(s,br);
    return s;
    }

const char *express::vprimitive(const char *s,dcvnode **br)
 {
 int c;
 char *ptr;
 if (s==NULL) return NULL;
 if (token[0]=='\"')
 {
	 ptr=strchr(token+1,'\"');
	 if (!ptr)
	 {
		 dcerror(1,"Syntax error: %s",(char*)token);
		 return NULL;
	 }
	 else
	 {
		 *ptr='\0';
		 *br=new vstrconstant(token+1);
	 }
 }
 else
 {
	 *br=new vconstant(atof(token));
	 c=token[0];
	 if((c!='.')&&!isdigit(c))
	 {
		 dcerror(1,"Syntax error in \"%s\": %s",s,(char*)token);
		 return NULL;
	 }
	 if(*br==NULL) { dcerror(ENOMEM); return NULL; }
 }
 return get_next(s);
 }

const char *express::get_next(const char *s)
 {
 int i,len;
 double dummy;
 const char *ptr;
 if (s==NULL) return NULL;
 while (isspace(*s)) ++s;
 strcpy(token,"");
 for (i=0;i<n_all_ops;++i)
   if (!strnicmp(s,all_op[i],strlen(all_op[i])))
    {
    strcpy(token,all_op[i]);
    return s+strlen(all_op[i]);
    }
 if (*s=='(' || *s==')')
   {
   token[0]=*s;
   token[1]='\0';
   return s+1;
   }
 else if (*s=='\"')
 {
	 ptr=strchr(s+1,'\"');
	 if (!ptr)
		 return s; // should produce syntax error
	 else
	 {
		 strncpy(token,s,ptr-s+1);
		 token[ptr-s+1]='\0';
		 return ptr+1;
	 }
 }
 else if (len=0,(i=sscanf(s,"%lf%n",&dummy,&len))>=1)
   {
   sprintf(dc_expr_buff,"%g",dummy);
   strcpy(token,dc_expr_buff);
#if 0   
   return s+(i==2?len-1:strlen(s));  // i is 1 at end of string
#else
   return s+(len!=0?len:strlen(s));
#endif
   }
 else return s;
 }

void add_un_op(char *lab,dcexpr_val *(*f)(dcvnode *)) 
{ 
un_ch[n_un_ops].str=lab;
un_ch[n_un_ops++].inst=new vnvun_op(f); 
all_op[n_all_ops++]=lab;
}
 
void add_bin_op_same(char *lab,dcexpr_val *(*f)(dcvnode *,dcvnode *)) 
{ 
bin_ch[last_op_level][n_bin_ops[last_op_level]].str=lab;
bin_ch[last_op_level][n_bin_ops[last_op_level]++].inst=new vnvbin_op(f); 
all_op[n_all_ops++]=lab;
}

void add_bin_op_next(char *lab,dcexpr_val *(*f)(dcvnode *,dcvnode *)) 
{ 
++last_op_level;
if (n_levels<last_op_level)
  n_levels=last_op_level;
bin_ch[last_op_level][n_bin_ops[last_op_level]].str=lab;
bin_ch[last_op_level][n_bin_ops[last_op_level]++].inst=new vnvbin_op(f); 
all_op[n_all_ops++]=lab;
}
	
void add_bin_op(char *lab,dcexpr_val *(*f)(dcvnode *,dcvnode *),int level) 
{ 
last_op_level=level;
if (n_levels<last_op_level)
  n_levels=last_op_level;
bin_ch[last_op_level][n_bin_ops[last_op_level]].str=lab;
bin_ch[last_op_level][n_bin_ops[last_op_level]++].inst=new vnvbin_op(f); 
all_op[n_all_ops++]=lab;
}

class op_initialiser
{
public:
op_initialiser();
~op_initialiser();
};

op_initialiser::~op_initialiser()
{
int i,j;
for (i=0;i<n_un_ops;++i) delete un_ch[i].inst;
for (i=0;i<n_levels;++i) for (j=0;bin_ch[i][j].str;++j) delete bin_ch[i][j].inst;
}

op_initialiser op_initialiser_instance;

op_initialiser::op_initialiser()
{		
add_bin_op(",",comma_op,0);
add_bin_op_next("|",or_op);
add_bin_op_same("OR",or_op);
add_bin_op_next("&",and_op);
add_bin_op_same("AND",and_op);
add_bin_op_next("IFTRUE",iftrue_op);
add_bin_op_next("!=",neq_op);
add_bin_op_same("=",equ_op);
add_bin_op_same("IS",equ_op);
add_bin_op_next("<=",lte_op); //  these pairs must be in this order
add_bin_op_same("<",lt_op);   // otherwise e.g. > will be found before >=
add_bin_op_same(">=",gte_op);
add_bin_op_same(">",gt_op);
add_bin_op_next("-",sub_op);
add_bin_op_same("+",add_op);
add_bin_op_next("MOD",mod_op);
add_bin_op_next("/",div_op);
add_bin_op_same("*",mult_op);
add_bin_op_next("POW",pow_op);

add_un_op("-",negate_op);
add_un_op("!",not_op);
add_un_op("NOT",not_op);
add_un_op("ABS",abs_op);
add_un_op("LOG",log_op);
add_un_op("LN",ln_op);
add_un_op("EXP",exp_op);
add_un_op("SIN",sin_op);
add_un_op("ARCSIN",arcsin_op);
add_un_op("COS",cos_op);
add_un_op("ARCCOS",arccos_op);
add_un_op("TAN",tan_op);
add_un_op("ARCTAN",arctan_op);
}

