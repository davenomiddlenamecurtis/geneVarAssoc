/*
 * Copyright (C) 1988, Marcus J. Ranum, William Welch Medical Library
 * $Author: mjr $
 */

/* Modified by Dave Curtis 1994 */

#ifndef NODCURTISCHANGES
#ifndef MSDOS
#include <unistd.h>
#else
#include <io.h>
#endif

#include <stdlib.h>

#ifndef DEC
#include <string.h>
#define bcopy(s,d,n) memcpy(d,s,n)
#endif

int errno; /* DC 26/6/02 I am putting this here because MSVC does not find it */
#endif

/*
 * $Log:	btree.c,v $
 * Revision 1.1  88/06/01  21:35:07  mjr
 * Initial revision
 * 
 */

#include	<stdio.h>
#include	"btree.h"

#ifndef NODCURTISCHANGES
static long popl(struct btree *bt);
static long popr(struct btree *bt);
int bthead(long *node_nbr,struct btnode *cno,struct btree *bt);
int btprevious(long	*node_nbr,struct btnode *cno,struct btree *bt);
int btnext(long *node_nbr,struct btnode *cno,struct btree *bt);
#endif


/* if we wish to store our disk data in network byte order */
#ifdef BYTEORDER
#include	<sys/types.h>
#include	<netinet/in.h>
#endif

#define	BT_NSIZ	(sizeof(struct btnode))
#define	BT_SSIZ	(sizeof(struct btsuper))


/* the errno for btree problems. we use negative # - */
/* so btperror can use the real UNIX errno */
int	bterrno = 0;
char	*bterrs[] = {
	"No btree error",
	"Bad index magic number",
	"History stack overflow",
	"Cannot delete node zero",
	0
};

#ifdef vax
extern	void	bcopy();
extern	void	free();
extern	void	exit();
extern	void	perror();
#endif


/* write the btree superblock to disk */
static int
wsuper(bt)
struct	btree	*bt;
{
#ifdef BYTEORDER
	struct	btsuper	boge;
#endif
#ifndef NODCURTISCHANGES
#ifndef _SYS_UNISTD_H 
#ifndef _UNISTD_H 
	extern	long	lseek();
#endif
#endif
#endif

	if (lseek(bt->fd, 0L, 0) < 0)
		return (-1);

#ifdef BYTEORDER
	boge.magic = htonl(bt->sblk.magic);
	boge.free = htonl(bt->sblk.free);
	boge.root = htonl(bt->sblk.root);
	boge.list = htonl(bt->sblk.list);

	if (write(bt->fd, (char *) &boge, BT_SSIZ) != BT_SSIZ)
		return (-1);
#else
	if (write(bt->fd, (char *) &bt->sblk, BT_SSIZ) != BT_SSIZ)
		return (-1);
#endif

	return (0);
}



/* dynamically allocate a control structure for an open btree */
struct btree	*
btopen(path, flags, mode)
char	*path;
int	flags;
int	mode;
{
	struct	btree	*bt;
	int	r;
#ifdef NODCURTISCHANGES
	extern	char	*malloc();
#endif

	/* lets be dynamic, shall we ? */
#ifndef lint
	/* this to avoid the possible pointer alignment lint message */
	if ((bt = (struct btree *) malloc(sizeof(struct btree))) == NULL)
		return (NULL);
#else
	bt = (struct btree *)0;
#endif

	if ((bt->fd = open(path, flags, mode)) > -1) {

		r = read(bt->fd, (char *) &bt->sblk, BT_SSIZ);

		/* if read nothing, must be a new guy, right ? */
		if (r == 0) {
			bt->sblk.magic = BT_MAGIC;
			bt->sblk.free = 1L;
			bt->sblk.root = 0L;
			bt->sblk.list = 0L;

			if (wsuper(bt) == 0)
				r = BT_SSIZ;
		}
#ifdef BYTEORDER
		else {
			/* read something, decode the numbers */
			bt->sblk.magic = ntohl(bt->sblk.magic);
			bt->sblk.free = ntohl(bt->sblk.free);
			bt->sblk.root = ntohl(bt->sblk.root);
			bt->sblk.list = ntohl(bt->sblk.list);
		}
#endif

		/* cleverly check ret value from either read or write */
		if (r != BT_SSIZ) {
			(void) close(bt->fd);
			(void) free((char *) bt);
			return (NULL);
		}

		/* check that ole magic number */
		if (bt->sblk.magic != BT_MAGIC) {
			bterrno = BT_BAD_MAGIC;
			(void) close(bt->fd);
			(void) free((char *) bt);
			return (NULL);
		}
	} else {
		/* couldnt even open the bloody file */
		(void) free((char *) bt);
		return (NULL);
	}

#ifndef NODCURTISCHANGES
bt->slev=0;
#endif
	/* zero the cache - record numbers will never be -1, */
	/* so the cache will load as activity takes place */
	for(r = 0; r < BT_CSIZ; r++)
#ifndef NODCURTISCHANGES
		bt->cache[r].nodeno = -1L;
#else
		bt->cache[r].recno = -1L;
#endif

	return (bt);
}



/* close and deallocate the control structure */
btclose(bt)
struct	btree	*bt;
{
	int	t;

	t = wsuper(bt);
	(void) close(bt->fd);
	(void) free((char *) bt);
	return (t);
}



/* write a node to disk */
static	int
wnode(nbr, npt, bt)
long	nbr;
struct	btnode	*npt;
struct	btree	*bt;
{
#ifdef NODCURTISCHANGES
	extern	long	lseek();
	extern	char	*strncpy();
#endif
#ifdef	BYTEORDER
	struct	btnode	boge;
#endif

	if (lseek(bt->fd, (long) ((nbr * BT_NSIZ) + BT_SSIZ), 0) == -1)
		return (-1);

#ifndef NODCURTISCHANGES
	npt->nodeno = nbr;
#endif

#ifdef	BYTEORDER
#ifndef NODCURTISCHANGES
	boge.nodeno = htonl(npt->nodeno);
#endif
	boge.recno = htonl(npt->recno);
	boge.lptr = htonl(npt->lptr);
	boge.rptr = htonl(npt->rptr);
	boge.deleted = htons(npt->deleted);
	boge.balance = htons(npt->balance);
	(void) strncpy(boge.key, npt->key, BT_KSIZ - 1);
	if (write(bt->fd, (char *) &boge, BT_NSIZ) != BT_NSIZ)
		return (-1);
#else
	if (write(bt->fd, (char *) npt, BT_NSIZ) != BT_NSIZ)
		return (-1);
#endif

	/* update cache -  if the write succeeded */
	(void)bcopy((char *)npt,(char *)&bt->cache[(int)nbr % BT_CSIZ],BT_NSIZ);
	return (0);
}



/* read a node from disk */
static	int
rnode(nbr, npt, bt)
long	nbr;
struct	btnode	*npt;
struct	btree	*bt;
{
#ifdef NODCURTISCHANGES
	extern	long	lseek();
	extern	char	*strncpy();
#endif
	int	hash;
#ifdef	BYTEORDER
	struct	btnode	boge;
#endif
	hash = (int)nbr % BT_CSIZ;

#ifndef NODCURTISCHANGES
/* !! this seemed wrong - assumes recno==node number !! */
	if(bt->cache[hash].nodeno != nbr) {
#else

	/* check for cache hit - simple hash - braindead, really */
	if(bt->cache[hash].recno != nbr) {
#endif

		/* if no cache hit, load from disk */
		if (lseek(bt->fd, (long) ((nbr * BT_NSIZ) + BT_SSIZ), 0) == -1)
			return (-1);
#ifndef BYTEORDER
		if (read(bt->fd, (char *) &bt->cache[hash], BT_NSIZ) != BT_NSIZ)
			return (-1);
#else
		if (read(bt->fd, (char *) &boge, BT_NSIZ) != BT_NSIZ)
			return (-1);

#ifndef NODCURTISCHANGES
		bt->cache[hash].nodeno = ntohl(boge.nodeno);
#endif
		bt->cache[hash].recno = ntohl(boge.recno);
		bt->cache[hash].lptr = ntohl(boge.lptr);
		bt->cache[hash].rptr = ntohl(boge.rptr);
		bt->cache[hash].deleted = ntohs(boge.deleted);
		bt->cache[hash].balance = ntohs(boge.balance);
		(void) strncpy(bt->cache[hash].key, boge.key, BT_KSIZ - 1);
#endif
	}

	/* loaded OK, now copy the updated cached data to the target */
	(void)bcopy((char *) &bt->cache[hash],(char *)npt,BT_NSIZ);

	return (0);
}



static	int
pushl(bt, nbr)
struct	btree	*bt;
long	nbr;
{
	if (++(bt->lstak.sptr) >= STACK_LENGTH) {
		bterrno = BT_BAD_STACK;
		exit(0);
	}
	bt->lstak.ele[bt->lstak.sptr] = nbr;
	bt->lstak.lev[bt->lstak.sptr] = ++bt->slev;
#ifndef NODCURTISCHANGES
	return (0);
#else
	return;
#endif
}


static	int
pushr(bt, nbr)
struct	btree	*bt;
long	nbr;
{
	if (++(bt->rstak.sptr) >= STACK_LENGTH) {
		bterrno = BT_BAD_STACK;
		exit(0);
	}
	bt->rstak.ele[bt->rstak.sptr] = nbr;
	bt->rstak.lev[bt->rstak.sptr] = ++bt->slev;
#ifndef NODCURTISCHANGES
	return (0);
#else
	return;
#endif
}



static	long
popr(bt)
struct	btree	*bt;
{

	bt->slev = bt->rstak.lev[bt->rstak.sptr];

	while (bt->lstak.lev[bt->lstak.sptr] > bt->slev)
		(bt->lstak.sptr)--;

	if (bt->rstak.sptr == 0)
		return (0);

	bt->slev--;
	return (bt->rstak.ele[(bt->rstak.sptr)--]);
}



static	long
popl(bt)
struct	btree	*bt;
{

	bt->slev = bt->lstak.lev[bt->lstak.sptr];

	while (bt->rstak.lev[bt->rstak.sptr] > bt->slev)
		(bt->rstak.sptr)--;

	if (bt->lstak.sptr == 0)
		return (0);

	bt->slev--;
	return (bt->lstak.ele[(bt->lstak.sptr)--]);
}



static	void
bt_link(alpha1, node1, alpha2, node2)
int	alpha1;
int	alpha2;
struct	btnode	*node1;
struct	btnode	*node2;
{
	if (alpha1 == -1 && alpha2 == -1)
		node1->lptr = node2->lptr;
	else if (alpha1 == -1 && alpha2 == 1)
		node1->lptr = node2->rptr;
	else if (alpha1 == 1 && alpha2 == -1)
		node1->rptr = node2->lptr;
	else
		node1->rptr = node2->rptr;
}



/* number a link according to alpha */
static	void
nbr_link(nbr, alpha, node1)
long		*nbr;
int		alpha;
struct	btnode	*node1;
{
	*nbr = (alpha == 1) ? node1->rptr : node1->lptr;
}



/* set a link according to alpha */
static	void
link_nbr(alpha, node1, nbr)
int		alpha;
struct	btnode	*node1;
long		nbr;
{

	if (alpha == 1)
		node1->rptr = nbr;
	else
		node1->lptr = nbr;
}



static	void
node_bal(alpha, node1, node2, node3)
int	alpha;
struct	btnode	*node1;
struct	btnode	*node2;
struct	btnode	*node3;
{

	if (node1->balance == alpha) {
		node2->balance = -alpha;
		node3->balance = 0;
	} else if (node1->balance == 0)
		node2->balance = node3->balance = 0;
	else {
		node2->balance = 0;
		node3->balance = alpha;
	}
}



/* change the record number in a node */
btsetrec(nbr, newrec, bt)
long	nbr;
long	newrec;
struct	btree	*bt;
{
	struct	btnode	tmpnode;

	if(rnode(nbr, &tmpnode, bt) <0)
		return(-1);

	tmpnode.recno = newrec;

	if(wnode(nbr, &tmpnode, bt) <0)
		return(-1);


	return(0);
}



/* insert the node into the tree */
btinsert(argkey, recnbr, bt)
char	*argkey;
long	recnbr;
struct	btree	*bt;
{
	long	top;
	long	p_rec;
	long	s_rec;
	long	q_rec;
	long	r_rec;
	struct	btnode	q_node;
	struct	btnode	s_node;
	struct	btnode	r_node;
	struct	btnode	p_node;
	int	compare;
	extern	char	*strncpy();


	/* the very first node gets inserted specially */
	if (bt->sblk.root == 0) {
		bt->sblk.root = 1;
		p_node.balance = p_node.rptr = p_node.lptr = 0;
		p_node.deleted = BT_ACTIVE;
		p_node.recno = recnbr;

		(void) strncpy(p_node.key, argkey, BT_KSIZ - 1);
		p_node.key[BT_KSIZ - 1] = '\0';
		if(wnode(1L, &p_node, bt) < 0)
			return(-1);

		bt->sblk.free++;
		bt->sblk.list++;
		if(wsuper(bt) <0)
			return(-1);
		return (0);
	}
	top = -1;
	p_rec = bt->sblk.root;
	s_rec = bt->sblk.root;
	while (1) {
		if(rnode(p_rec, &p_node, bt) <0)
			return(-1);
		if ((compare = strcmp(argkey, p_node.key)) < 0) {

			/* (move left) */
			q_rec = p_node.lptr;

			if (q_rec == 0) {
				/* insert here */
				q_rec = bt->sblk.free++;
				p_node.lptr = q_rec;
				break;	/* loop exit */
			} else {
				/* look again from this node */
				if(rnode(q_rec, &q_node, bt) <0)
					return(-1);
				if (q_node.balance != 0) {
					top = p_rec;
					s_rec = q_rec;
				}
			}
			p_rec = q_rec;

		} else {
			/* (move right) */
			q_rec = p_node.rptr;

			if (q_rec == 0) {
				/* insert here */
				q_rec = bt->sblk.free++;
				p_node.rptr = q_rec;
				break;	/* loop exit */
			} else {
				/* look again from this node */
				if(rnode(q_rec, &q_node, bt) <0)
					return(-1);
				if (q_node.balance != 0) {
					top = p_rec;
					s_rec = q_rec;
				}
				p_rec = q_rec;
			}
		}
	}

	/* Step 5 (insert key at q_node) */
	q_node.lptr = q_node.rptr = 0;
	q_node.balance = 0;
	q_node.deleted = BT_ACTIVE;
	q_node.recno = recnbr;
	(void) strncpy(q_node.key, argkey, BT_KSIZ);
	q_node.key[BT_KSIZ - 1] = '\0';

	if (wnode(q_rec, &q_node, bt) == -1)
		return (-1);

	if(wnode(p_rec, &p_node, bt) <0)
		return(-1);
	if(rnode(s_rec, &s_node, bt) <0)
		return(-1);
	if(wsuper(bt) <0)
		return(-1);

	/* (adjust balance factors) */
	if (strcmp(argkey, s_node.key) < 0) {
		r_rec = p_rec = s_node.lptr;
	} else {
		r_rec = p_rec = s_node.rptr;
	}

	while (p_rec != q_rec) {
		if(rnode(p_rec, &p_node, bt) <0)
			return(-1);
		if (strcmp(argkey, p_node.key) < 0) {
			p_node.balance = -1;
			if(wnode(p_rec, &p_node, bt) < 0)
				return(-1);
			p_rec = p_node.lptr;
		} else {
			p_node.balance = 1;
			if(wnode(p_rec, &p_node, bt) < 0)
				return(-1);
			p_rec = p_node.rptr;
		}
	}

	compare = (strcmp(argkey, s_node.key) < 0) ? -1 : 1;
	if (s_node.balance == 0) {
		bt->sblk.list++;
		s_node.balance = compare;
		if(wnode(s_rec, &s_node, bt) < 0)
			return(-1);
		if(wsuper(bt) <0)
			return(-1);
		return (0);
	} else if (s_node.balance == -compare) {
		bt->sblk.list++;
		s_node.balance = 0;
		if(wnode(s_rec, &s_node, bt) < 0)
			return(-1);
		if(wsuper(bt) <0)
			return(-1);
		return (0);
	} else {
		/* (tree out of balance) */
		bt->sblk.list++;
		if(rnode(s_rec, &s_node, bt) <0)
			return(-1);
		if(rnode(r_rec, &r_node, bt) <0)
			return(-1);
		if (r_node.balance == compare) {
			/* (single rotate) */
			p_rec = r_rec;
			bt_link(compare, &s_node, -compare, &r_node);
			link_nbr(-compare, &r_node, s_rec);
			s_node.balance = r_node.balance = 0;
		} else {
			/* (double rotate) */
			nbr_link(&p_rec, -compare, &r_node);
			if(rnode(p_rec, &p_node, bt) <0)
				return(-1);
			bt_link(-compare, &r_node, compare, &p_node);
			link_nbr(compare, &p_node, r_rec);
			bt_link(compare, &s_node, -compare, &p_node);
			link_nbr(-compare, &p_node, s_rec);
			node_bal(compare, &p_node, &s_node, &r_node);
			p_node.balance = 0;
			if(wnode(p_rec, &p_node, bt) <0)
				return(-1);
		}

		if (top == -1) {
			bt->sblk.root = p_rec;
		} else {
			/* balanced at top of a sub-tree */
			if(rnode(top, &q_node, bt) < 0)
				return(-1);
			if (s_rec == q_node.rptr)
				q_node.rptr = p_rec;
			else
				q_node.lptr = p_rec;
			if(wnode(top, &q_node, bt) <0)
				return(-1);
		}
		if(wnode(s_rec, &s_node, bt) <0)
			return(-1);
		if(wnode(r_rec, &r_node, bt) <0)
			return(-1);
		if(wsuper(bt) <0)
			return(-1);
		return (0);
	}
}



/* drop a node */
btdelete(node_nbr, bt)
long	node_nbr;
struct	btree	*bt;
{
	struct	btnode	cno;

	if (node_nbr == 0) {
		bterrno = BT_BAD_ROOT;
		return (-1);
	} else {
		if (rnode(node_nbr, &cno, bt) == -1) 
			return(-1);
		cno.deleted = BT_DELETED;
		if (wnode(node_nbr, &cno, bt) == -1) {
			return (-1);
		} else {
			bt->sblk.list--;
			if(wsuper(bt) <0)
				return(-1);
		}
	}
	return (0);
}



/* find the next node */
btnext(node_nbr, cno, bt)
long	*node_nbr;
struct	btnode	*cno;
struct	btree	*bt;
{
#if defined NODCURTISCHANGES
	long	popl();
	long	popr();
#endif
	long	s_nod;

	s_nod = *node_nbr;

	if (*node_nbr == 0) {
		/* undefined current node - wind to beginning of file */
		return (bthead(node_nbr,cno,bt));
	}
	do {
		if (cno->rptr == 0) {
			/* can't move right */
			if (bt->lstak.sptr == 0) {
				/* none in stack */
				if(rnode(*node_nbr, cno, bt) <0)
					return(-1);
				return (BT_EOF);
			} else {
				/* can't go right & stack full (pop stack) */
				s_nod = popl(bt);
				if(rnode(s_nod, cno, bt) < 0)
					return(-1);
			}
		} else {
			/* move right */
			pushr(bt, s_nod);
			s_nod = cno->rptr;
			if(rnode(s_nod, cno, bt) <0 )
				return(-1);

			while (cno->lptr != 0) {
				/* bottom left */
				pushl(bt, s_nod);
				/* of this sub-tree */
				s_nod = cno->lptr;
				if(rnode(s_nod, cno, bt) <0)
					return(-1);
			}
		}
	} while (cno->deleted == BT_DELETED);

	*node_nbr = s_nod;
	return (0);
}



/* go to the tail of the file */
bttail(node_nbr, cno, bt)
struct	btree	*bt;
long	*node_nbr;
struct	btnode	*cno;
{
	long	s_nod;

	bt->rstak.sptr = 0;
	bt->lstak.sptr = 0;
	bt->rstak.ele[0] = 0;
	bt->lstak.ele[0] = 0;
#ifndef NODCURTISCHANGES
	bt->rstak.lev[0] = 0;
	bt->lstak.lev[0] = 0;
#endif
	/* begin at list head */
	s_nod = bt->sblk.root;

	if(rnode(s_nod, cno, bt) <0)
		return(-1);
	while (1) {
		/* search to right */
		if (cno->rptr != 0) {
			pushr(bt, s_nod);
			s_nod = cno->rptr;
			if(rnode(s_nod, cno, bt) <0)
				return(-1);
		} else {
			if (cno->deleted == BT_DELETED) {
				/* skip all deleted nodes */
				while (cno->deleted == BT_DELETED) {
#ifndef NODCURTISCHANGES
/* see comments for bthead() */
					if (btprevious(&s_nod, cno, bt) == BT_EOF) {
						if(btnext(&s_nod, cno, bt)<0)
							return(-1);
						*node_nbr = s_nod;
						return (BT_EOF);
#else
					if (btnext(&s_nod, cno, bt) == BT_EOF) {
						if(btprevious(&s_nod, cno, bt)<0)
							return(-1);
						*node_nbr = s_nod;
						return (0);
#endif
					}
				}
				*node_nbr = s_nod;
				return (0);
			} else {
				/* at end of a branch */
				*node_nbr = s_nod;
				return (0);
			}
		}
	}
}


/* go to the head of the file */
bthead(node_nbr, cno, bt)
struct	btree	*bt;
long	*node_nbr;
struct	btnode	*cno;
{
	long	s_nod;

	bt->rstak.sptr = 0;
	bt->lstak.sptr = 0;
	bt->rstak.ele[0] = 0;
	bt->lstak.ele[0] = 0;
#ifndef NODCURTISCHANGES
	bt->rstak.lev[0] = 0;
	bt->lstak.lev[0] = 0;
#endif

	/* begin at list head */
	s_nod = bt->sblk.root;

	if(rnode(s_nod, cno, bt) <0)
		return(-1);
	while (1) {
		/* search to left */
		if (cno->lptr != 0) {
			pushl(bt, s_nod);
			s_nod = cno->lptr;
			if(rnode(s_nod, cno, bt) <0)
				return(-1);
		} else {
			if (cno->deleted == BT_DELETED) {
				/* skip all deleted nodes */
				while (cno->deleted == BT_DELETED) {
#ifndef NODCURTISCHANGES
/* go forwards to undeleted one				
   and if all deleted position at start and
   return BT_EOF
   
   btnext() and btprevious() skip deleted records
*/
					if (btnext(&s_nod, cno, bt) == BT_EOF) {
						if(btprevious(&s_nod, cno, bt)<0)
							return(-1);
						*node_nbr = s_nod;
						return (BT_EOF);
#else
					if (btprevious(&s_nod, cno, bt) == BT_EOF) {
						if(btnext(&s_nod, cno, bt)<0)
							return(-1);
						*node_nbr = s_nod;
						return (0);
#endif
					}
				}
				*node_nbr = s_nod;
				return (0);
			} else {
				/* at end of a branch */
				*node_nbr = s_nod;
				return (0);
			}
		}
	}
}



/* find a key */
btfind(char *key1,long *node_nbr, struct btnode *cno,struct btree *bt)
{
	int	direction;
	long	s_nod;

	bt->rstak.sptr = 0;
	bt->lstak.sptr = 0;
	bt->rstak.ele[0] = 0;
	bt->lstak.ele[0] = 0;

	bt->slev = 0;		/* tree level at start of search */

	/* begin at list head */
	s_nod = bt->sblk.root;

	if(rnode(s_nod, cno, bt) <0)
		return(-1);
	while ((direction = strcmp(key1, cno->key)) != 0 ||
		cno->deleted == BT_DELETED) {

		if (direction > 0) {
			/* search to right */
			if (cno->rptr != 0) {
				pushr(bt, s_nod);
				s_nod = cno->rptr;
				if(rnode(s_nod, cno, bt) <0)
					return(-1);
			} else if (cno->deleted == BT_DELETED) {
				/* skip all deleted nodes */
				while (cno->deleted == BT_DELETED) {
					if (btnext(&s_nod, cno, bt) == BT_EOF) {
						if(btprevious(&s_nod, cno, bt)<0)
							return(-1);
						*node_nbr = s_nod;
						return (BT_NREC);
					}
				}
				*node_nbr = s_nod;
				return (BT_NREC);
			} else {
				/* at end of a branch */
				*node_nbr = s_nod;
				return (BT_NREC);
			}
#ifndef NODCURTISCHANGES
		} else if (direction<0) {
#else
		} else {
#endif
			/* search to left */
			if (cno->lptr != 0) {
				pushl(bt, s_nod);
				s_nod = cno->lptr;
				if(rnode(s_nod, cno, bt) <0)
					return(-1);
			} else if (cno->deleted == BT_DELETED) {
				while (cno->deleted == BT_DELETED) {
					if (btnext(&s_nod, cno, bt) == BT_EOF) {
						if(btprevious(&s_nod, cno, bt)< 0)
							return(-1);
						*node_nbr = s_nod;
						return (BT_NREC);
					}
				}
				*node_nbr = s_nod;
				return (BT_NREC);
			} else {
				*node_nbr = s_nod;
				return (BT_NREC);
			}
#ifndef NODCURTISCHANGES
		} else {
/* found a match, but it's deleted, so look forwards
   and backwards for an undeleted match */		
			long	kept_s_nod;
			struct	btnode	kept_cno;
			int found_one;
			kept_s_nod=s_nod;
			kept_cno=*cno;
			found_one=0;
			while (1) {
				if (btprevious(&s_nod, cno, bt) == BT_EOF)
					break; 
					/* all previous records deleted */
				if (strcmp(key1,cno->key)) break;
				if (cno->deleted != BT_DELETED)
					{ found_one=1; break; }
				}
			if (found_one==0) {
			s_nod=kept_s_nod;
			*cno=kept_cno;
			while (1) {
				if (btnext(&s_nod, cno, bt) == BT_EOF)
					{ 
					if (btprevious(&s_nod, cno, bt)<0)
						return -1;
					else break; 
					}
				if (strcmp(key1,cno->key)) break;
				if (cno->deleted != BT_DELETED)
					{ found_one=1; break; }
				}
			}
			if (found_one==0) {
				*node_nbr = s_nod;
				return (BT_NREC);
			}
#endif
		}
	}
	*node_nbr = s_nod;
	return (0);
}



/* get the previous node */
btprevious(long	*node_nbr,struct btnode	*cno,struct btree *bt)
{
#if defined NODCURTISCHANGES
	long	popr();
	long	popl();
#endif
	long	s_nod;

	s_nod = *node_nbr;

	/* if we are called without a node, wind to the end of file */
	if (*node_nbr == 0) {
		return (bttail(node_nbr,cno,bt));
	}


	do {
		if (cno->lptr == 0) {
			/* can't move left */
			if (bt->rstak.sptr == 0) {
				/* none in stack */
				if(rnode(*node_nbr, cno, bt) <0)
					return(-1);
				return (BT_EOF);
				/* don't reset node_nbr */
			} else {
				/* can't go left & stack full (pop stack) */
				s_nod = popr(bt);
				if(rnode(s_nod, cno, bt) <0)
					return(-1);
			}
		} else {
			/* left then to bottom right - is previous */
			pushl(bt, s_nod);
			s_nod = cno->lptr;
			if(rnode(s_nod, cno, bt) <0)
				return(-1);
			while (cno->rptr != 0) {
				/* bottom right */
				pushr(bt, s_nod);
				/* of this sub-tree */
				s_nod = cno->rptr;
				if(rnode(s_nod, cno, bt) <0)
					return(-1);
			}
		}
	} while (cno->deleted == BT_DELETED);

	*node_nbr = s_nod;
	return (0);
}



/* print the btree error message */
void
btperror(str)
char	*str;
{
	extern	int	errno;

	/* is it ours ?? */
	if (errno == 0) {
		if (str[strlen(str) - 1] != ':')
			(void) fprintf(stderr, "%s: %s\n", str, bterrs[bterrno]);
		else
			(void) fprintf(stderr, "%s %s\n", str, bterrs[bterrno]);
		bterrno = 0;
	} else {
		perror(str);
	}
}
