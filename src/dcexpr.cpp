/* DCEXPR.CXX */ 
/* Copyright Dave Curtis 1994 */ 
/* dcurtis@hgmp.mrc.ac.uk */ 
/* no warranty or liability of any kind is accepted, expressed or implied */ 
 
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
char *mystrupr(char *s)
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

#define EVAL_R1 \
if ((r1=b1->eval())==NULL) return NULL; 

#define EVAL_R2 \
if ((r2=b2->eval())==NULL) return NULL; 

#define EVAL_BOTH \
if ((r1=b1->eval())==NULL) return NULL; \
if ((r2=b2->eval())==NULL) { delete r1; return NULL; }

express::express()
{ head=NULL; token[0]='\0'; }

dcexpr_val *express::eval()
{ 
dcexpr_val *rv;
rv=head->eval(); 
if (rv==NULL) dcerror(ENOMEM,"Error evaluating expression - probably out of memory");
return rv;
}

dcexpr_val::~dcexpr_val() {;}

int dcexpr_double::is_string_really() { return 0; }
dcexpr_double::operator char*() { sprintf(buff,"%f",val); return buff; }
dcexpr_double::operator double() { return val; }
dcexpr_double::~dcexpr_double() {;}
dcexpr_double::dcexpr_double(double v) { val=v; }

int dcexpr_string::is_string_really() { return 1; }
dcexpr_string::operator char*() { return buff; }
dcexpr_string::operator double() { return atof(buff); }
dcexpr_string::~dcexpr_string() { if (buff) delete buff; }
dcexpr_string::dcexpr_string(char *t,int len) 
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

char dc_expr_buff[2000];

class vnvbin_op:public vnode {
dcexpr_val *(*func)(vnode *,vnode *);
public:
vnvbin_op(dcexpr_val *(*f)(vnode *,vnode *));
~vnvbin_op(){;}
dcexpr_val *eval();
vnode *copy();
};

vnode *vnvbin_op::copy()
{
return new vnvbin_op(func);
}

vnvbin_op::vnvbin_op(dcexpr_val *(*f)(vnode *,vnode *)) : vnode(2)
{
func=f;
}

dcexpr_val *vnvbin_op::eval()
{
if (branch[0]==NULL || branch[1]==NULL)
  return 0;
return func(branch[0],branch[1]);
}

class vnvun_op:public vnode {
dcexpr_val *(*func)(vnode *);
public:
vnvun_op(dcexpr_val *(*f)(vnode *));
~vnvun_op(){;}
dcexpr_val *eval();
vnode *copy();
};

vnode *vnvun_op::copy()
{
return new vnvun_op(func);
}

vnvun_op::vnvun_op(dcexpr_val *(*f)(vnode *)) : vnode(1)
{
func=f;
}

dcexpr_val *vnvun_op::eval()
{
if (branch[0]==NULL) return NULL;
return func(branch[0]);
}

variable::~variable() {;}

variable::variable() : vnode(0) {;}

vnode *vconstant::copy() { return new vconstant(value); }
vconstant::~vconstant(){;}
dcexpr_val *vconstant::eval() { return new dcexpr_double(value); }

vnode *vstrconstant::copy() { return new vstrconstant(value); }
vstrconstant::vstrconstant(char *s) : vnode(0)
  { if ((value=new char[strlen(s)+1])!=NULL) strcpy(value,s); }
vstrconstant::~vstrconstant() { if (value) delete value;}

dcexpr_val *vstrconstant::eval() { return new dcexpr_string(value); }

dcexpr_val *iftrue_op(vnode *b1,vnode *b2)
{
dcexpr_val *r1;
EVAL_R1;
double test=double(*r1);
delete r1;
if (test) return b2->eval();
else return new dcexpr_string("");
}

dcexpr_val *mult_op(vnode *b1,vnode *b2)
{
dcexpr_val *r1,*r2;
EVAL_BOTH;
double test=double(*r1)* double(*r2);
delete r1;
delete r2;
return new dcexpr_double(test);
}

dcexpr_val *or_op(vnode *b1,vnode *b2)
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

dcexpr_val *and_op(vnode *b1,vnode *b2)
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

dcexpr_val *div_op(vnode *b1,vnode *b2)
{
dcexpr_val *r1,*r2;
EVAL_BOTH;
double rv=double(*r1) / double(*r2);
delete r1; delete r2;
return new dcexpr_double(rv);
}

dcexpr_val *comma_op(vnode *b1,vnode *b2)
{
dcexpr_val *r1,*r2;
EVAL_BOTH;
delete r1;
return r2;
}

dcexpr_val *mod_op(vnode *b1,vnode *b2)
{
dcexpr_val *r1,*r2;
EVAL_BOTH;
int t1=(int)double(*r1),t2=(int)double(*r2),m;
delete r1; delete r2;
if (t2<0) t2= -t2;
return new dcexpr_double(t2?((m=t1%t2)>=0)?m:m+t2 : 0.0);
}
// x mod 0 = 0

dcexpr_val *add_op(vnode *b1,vnode *b2)
{
dcexpr_val *r1,*r2;
EVAL_BOTH;
double rv=double(*r1) + double(*r2);
delete r1; delete r2;
return new dcexpr_double(rv);
}


dcexpr_val *sub_op(vnode *b1,vnode *b2)
{
dcexpr_val *r1,*r2;
EVAL_BOTH;
double rv=double(*r1) - double(*r2);
delete r1; delete r2;
return new dcexpr_double(rv);
}

dcexpr_val *lt_op(vnode *b1,vnode *b2)
{
dcexpr_val *r1,*r2;
EVAL_BOTH;
double rv=double(*r1) < double(*r2);
delete r1; delete r2;
return new dcexpr_double(rv);
}

dcexpr_val *lte_op(vnode *b1,vnode *b2)
{
dcexpr_val *r1,*r2;
EVAL_BOTH;
double rv=double(*r1) <= double(*r2);
delete r1; delete r2;
return new dcexpr_double(rv);
}

dcexpr_val *gt_op(vnode *b1,vnode *b2)
{
dcexpr_val *r1,*r2;
EVAL_BOTH;
double rv=double(*r1) > double(*r2);
delete r1; delete r2;
return new dcexpr_double(rv);
}

dcexpr_val *gte_op(vnode *b1,vnode *b2)
{
dcexpr_val *r1,*r2;
EVAL_BOTH;
double rv=double(*r1) >= double(*r2);
delete r1; delete r2;
return new dcexpr_double(rv);
}

dcexpr_val *equ_op(vnode *b1,vnode *b2)
{ 
dcexpr_val *r1,*r2;
EVAL_BOTH;
double rv;
if (r1->is_string_really() && r2->is_string_really())
  rv=!strncmp(strupr((char*)(*r1)),strupr((char*)(*r2)),strlen((char*)(*r2)));
else
  rv=double(*r1) == double(*r2);
delete r1; delete r2;
return new dcexpr_double(rv);
}

dcexpr_val *neq_op(vnode *b1,vnode *b2)
{ 
double rv;
dcexpr_val *r1;
if ((r1=equ_op(b1,b2))==NULL) return NULL;
rv=double(*r1)==0.0;
delete r1;
return new dcexpr_double(rv);
}

dcexpr_val *pow_op(vnode *b1,vnode *b2)
{
dcexpr_val *r1,*r2;
EVAL_BOTH;
double rv=pow(double(*r1),double(*r2));
delete r1; delete r2;
return new dcexpr_double(rv);
}

dcexpr_val *abs_op(vnode *b1)
{
dcexpr_val *r1;
EVAL_R1;
double rv=fabs(double(*r1));
delete r1;
return new dcexpr_double(rv);
}

dcexpr_val *negate_op(vnode *b1)
{
dcexpr_val *r1;
EVAL_R1;
double rv=0.0-double(*r1);
delete r1;
return new dcexpr_double(rv);
}

dcexpr_val *exp_op(vnode *b1)
{
dcexpr_val *r1;
EVAL_R1;
double rv=exp(double(*r1));
delete r1;
return new dcexpr_double(rv);
}

dcexpr_val *not_op(vnode *b1)
{
dcexpr_val *r1;
EVAL_R1;
double rv=(double(*r1)==0.0);
delete r1;
return new dcexpr_double(rv);
}

dcexpr_val *ln_op(vnode *b1)
{
dcexpr_val *r1;
EVAL_R1;
double rv=log(double(*r1));
delete r1;
return new dcexpr_double(rv);
}

dcexpr_val *log_op(vnode *b1)
{
dcexpr_val *r1;
EVAL_R1;
double rv=log10(double(*r1));
delete r1;
return new dcexpr_double(rv);
}

dcexpr_val *sin_op(vnode *b1)
{
dcexpr_val *r1;
EVAL_R1;
double rv=sin(double(*r1));
delete r1;
return new dcexpr_double(rv);
}

dcexpr_val *arcsin_op(vnode *b1)
{
dcexpr_val *r1;
EVAL_R1;
double rv=asin(double(*r1));
delete r1;
return new dcexpr_double(rv);
}

dcexpr_val *cos_op(vnode *b1)
{
dcexpr_val *r1;
EVAL_R1;
double rv=cos(double(*r1));
delete r1;
return new dcexpr_double(rv);
}

dcexpr_val *arccos_op(vnode *b1)
{
dcexpr_val *r1;
EVAL_R1;
double rv=acos(double(*r1));
delete r1;
return new dcexpr_double(rv);
}

dcexpr_val *tan_op(vnode *b1)
{
dcexpr_val *r1;
EVAL_R1;
double rv=tan(double(*r1));
delete r1;
return new dcexpr_double(rv);
}

dcexpr_val *arctan_op(vnode *b1)
{
dcexpr_val *r1;
EVAL_R1;
double rv=atan(double(*r1));
delete r1;
return new dcexpr_double(rv);
}

vnode& vnode::operator=(vnode &old)
{
for (int i=0;i<nbranches;++i)
  if (old.branch[i])
    {
    if ((branch[i]=old.branch[i]->copy())!=NULL)
      *branch[i]=*old.branch[i];
    }
return *this;
}

vnode::vnode(int nb)
 {
 nbranches=nb;
 for (int i=0;i<nbranches;++i) branch[i]=NULL;
 }

vnode::~vnode()
 {
 wipe();
 }

void vnode::wipe()
 {
 for (int i=0;i<nbranches;++i) if (branch[i]!=NULL) delete branch[i];
 }

int vnode::add(vnode *vn,int b)
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

int express::parse(char *s)
 {
 char *ptr;
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

char *express::vbin_op(char *s,vnode **br,int level)
 {
 vnode *tempbr;
 if (level<=n_levels)
      {
      if ((s=vbin_op(s,br,level+1))==NULL) return NULL;
      for (int i=0;i<n_bin_ops[level];++i)
         if (!stricmp(bin_ch[level][i].str,token))
           {
           if ((s=get_next(s))==NULL) return NULL;
           if ((s=vbin_op(s,&tempbr,level+1))==NULL) return NULL;
           vnode *temp;
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

char *express::vun_op(char*s,vnode **br)
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
       vnode *temp;
       if ((temp=un_ch[i].inst->copy())==NULL) return NULL;
       temp->add(*br,0);
       *br=temp;
       }
    return s;
    }

char *express::vbracket(char*s,vnode **br)
    {
    if (!strcmp("(",token))
     {
     if ((s=get_next(s))==NULL) return NULL;
     if ((s=vbin_op(s,br,0))==NULL) return NULL;
     if (strcmp(")",token)) { dcerror(1,"expected \')\': %s",(char *)token); return NULL; }
     if ((s=get_next(s))==NULL) return NULL;
     }
    else s=vprimitive(s,br);
    return s;
    }

char *express::vprimitive(char *s,vnode **br)
 {
 int c;
 if (s==NULL) return NULL;
 *br=new vconstant(atof(token));
 c=token[0];
 if ((c!='.')&&!isdigit(c)) 
   {
   dcerror(1,"Syntax error: %s",(char*)token); 
   return NULL; 
   }
 if (*br==NULL) { dcerror(ENOMEM); return NULL; }
 return get_next(s);
 }

char *express::get_next(char *s)
 {
 int i,len;
 double dummy;
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

void add_un_op(char *lab,dcexpr_val *(*f)(vnode *)) 
{ 
un_ch[n_un_ops].str=lab;
un_ch[n_un_ops++].inst=new vnvun_op(f); 
all_op[n_all_ops++]=lab;
}
 
void add_bin_op_same(char *lab,dcexpr_val *(*f)(vnode *,vnode *)) 
{ 
bin_ch[last_op_level][n_bin_ops[last_op_level]].str=lab;
bin_ch[last_op_level][n_bin_ops[last_op_level]++].inst=new vnvbin_op(f); 
all_op[n_all_ops++]=lab;
}

void add_bin_op_next(char *lab,dcexpr_val *(*f)(vnode *,vnode *)) 
{ 
++last_op_level;
if (n_levels<last_op_level)
  n_levels=last_op_level;
bin_ch[last_op_level][n_bin_ops[last_op_level]].str=lab;
bin_ch[last_op_level][n_bin_ops[last_op_level]++].inst=new vnvbin_op(f); 
all_op[n_all_ops++]=lab;
}
	
void add_bin_op(char *lab,dcexpr_val *(*f)(vnode *,vnode *),int level) 
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
add_bin_op_same("<",lt_op);   // otherwise ie.g. > will be found before >=
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

