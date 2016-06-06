/*
 * Copyright (C) 1988, Marcus J. Ranum, William Welch Medical Library
 * $Author: mjr $
 *
 * $Header: btree.h,v 1.1 88/06/01 21:35:55 mjr Rel $: btree.h
 *
 * $Log:	btree.h,v $
 * Revision 1.1  88/06/01  21:35:55  mjr
 * Initial revision
 * 
 */

/* Modified by Dave Curtis 1994 */

#ifndef _INCL_BTREE_H

#define BT_NREC	1	/* no such record */
#define BT_EOF	2	/* end of the tree (either end) */
#define BT_ERR	-1	/* something went wrong */

#define BT_KSIZ	40	/* size of keys to store (or trunc) */
#define BT_CSIZ	30	/* # of nodes to cache readonly */
			/* current cache alg gives ~59% hits */

/* this indicates a node deleted */
#define	BT_DELETED	1
#define	BT_ACTIVE	0

#define	BT_MAGIC	0x72251

/* btree stack */
#define STACK_LENGTH 30		/* length of history stacks */

struct btstack {
	long	ele[STACK_LENGTH];	/* stack elements */
	int	lev[STACK_LENGTH];	/* stack levels */
	int	sptr;			/* stack pointer */
};

/* a disk resident btree super block */
struct	btsuper {
	long	magic;		/* generic magic number */
	long	free;		/* pointer to next free (basically, EOF) */
	long	root;		/* pointer to root node */
	long	list;		/* number of active records */
};

struct btnode {
#ifndef NODCURTISCHANGES
	long	nodeno;		/* different to recno !!! */
#endif
	long	recno;		/* index to external file, or whatever */
	long	lptr;		/* left-pointer */
	long	rptr;		/* right-pointer */
	char	key[BT_KSIZ];	/* datum */
	short	deleted;	/* deleted flag */
	short	balance;	/* balance flag */
};
#define	BTNODE	struct	btnode

/* a memory resident btree super block */
/* including room to hold a disk super block */
struct btree {
	int	fd;		/* file descriptor */
	int	slev;		/* history stack level */
	struct	btsuper	sblk;	/* copy of superblock */
	struct	btstack	rstak;	/* history stack */
	struct	btstack	lstak;	/* history stack */
	struct	btnode	cache[BT_CSIZ];	/* read-only cache */
};
#define	BTREE	struct	btree

#if defined NODCURTISCHANGES || ! defined DCINDEXHPP
BTREE	*btopen();
#endif

void	btperror();

extern	int	bterrno;	/* btree error number/flag */
extern	char	*bterrs[];	/* error message list */
/* error codes - match bterrs */
#define	BT_BAD_MAGIC	1	/* bad index file magic number */
#define	BT_BAD_STACK	2	/* history stack overflow */
#define	BT_BAD_ROOT	3	/* failed attempt to delete node 0 */

#define _INCL_BTREE_H
#endif
