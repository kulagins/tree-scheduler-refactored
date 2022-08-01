/*
 *  lib-io-tree-minmem.h
 *  lib-io-tree
 *
 *  Created by defbond on 8/10/10.
 *  Copyright 2010 LIP/ENS-Lyon. All rights reserved.
 *
 */
#ifndef LIB_IO_TREE_MINMEM_H
#define LIB_IO_TREE_MINMEM_H

#include <list>
#include "tree.h"

#define uint64_t int

struct node_peak_t{
  int index;
  double mpeak;
  double mavail;
  node_peak_t(){
    index = -1;
    mpeak = 0;
    mavail = 0;
  }
  node_peak_t(int index_){
    index = index_;
    mpeak = 0.0;
    mavail = 0;
  }
};

struct iter_node_t{
  int index;
  double cost;
  double mpeak;
  iter_node_t(){
    index = -1;
    mpeak = 0;
    cost = -1;
  }

  iter_node_t * pSCbegin;
  iter_node_t * pSCend;
  iter_node_t * pNext;
  iter_node_t * pPrev;
};


void MinMem(Tree *tree, double MaxOutDeg, double &Required_memory, schedule_traversal &Schedule,
            int quiet);//added by Changjiang

void MinMemArray(int N, int root, double *nwghts, double *ewghts, int *chstart, int *children, double MaxOutDeg,
                 double &Required_memory, int *Schedule, int quiet, int &count);

void GreedyMinMem(Tree *tree, double &Required_memory);

#endif
