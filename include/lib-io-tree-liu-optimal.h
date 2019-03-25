/*
 *  lib-io-tree-liu-optimal.h
 *  lib-io-tree
 *
 *  Created by defbond on 8/10/10.
 *  Copyright 2010 LIP/ENS-Lyon. All rights reserved.
 *
 */

#ifndef LIB_IO_TREE_LIU_OPTIMAL_H
#define LIB_IO_TREE_LIU_OPTIMAL_H

#ifdef __cplusplus


#include <list>
#include "lib-io-tree-utils.h"


struct iter_node_t{
	unsigned int index;
	double cw;
	iter_node_t * pNext;
};

struct iter_seg_t {
	iter_node_t * pStart;/*this is the node after the next valley*/
	iter_node_t * pEnd; /*this is the next valley*/
	unsigned int root;
	double value;
	double Hi;
	double standalone_Hi;
	int seq_index;
	iter_node_t * hi;
};

void PebbleOrderingRecur(Ctree * tree,unsigned int sub_root, OrdoLiu_t & SubSchedule, int quiet, int & count);
double PebbleOrderingIter(const int N,const int *prnts,const double *nwghts,const double *ewghts,const int * chstart, const int * chend, int * children, const int root, int *schedule);

#endif

#endif
