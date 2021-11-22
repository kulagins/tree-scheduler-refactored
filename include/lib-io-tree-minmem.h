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

class MinMemDLL {
  protected:
  iter_node_t * pSentinel;
    bool custom_sentinel;
    uint64_t ui_size;
  public:
    MinMemDLL(){
      pSentinel = new iter_node_t;
      pSentinel->pNext = pSentinel;
      pSentinel->pPrev = pSentinel;
      custom_sentinel = false;
      ui_size = 0;
    };
    MinMemDLL(iter_node_t * pSentinel){
      this->pSentinel = pSentinel;
      pSentinel->pNext = pSentinel;
      pSentinel->pPrev = pSentinel;
      custom_sentinel = true;
      ui_size = 0;
    };

    ~MinMemDLL(){
      if(!custom_sentinel){
        delete pSentinel;
      }
    };


    void clear(){
      pSentinel->pNext = pSentinel;
      pSentinel->pPrev = pSentinel;
      ui_size = 0;
    }

void push_back(iter_node_t * pNewNode){
      pSentinel->pPrev->pNext = pNewNode;
      pNewNode->pPrev = pSentinel->pPrev;
      pNewNode->pNext = pSentinel;
      pSentinel->pPrev = pNewNode;
      ui_size++;
    };


void insert_after(iter_node_t * pInsertAfterNode, iter_node_t * pNewNode){
      pInsertAfterNode->pNext->pPrev = pNewNode;
      pNewNode->pNext = pInsertAfterNode->pNext;
      pNewNode->pPrev = pInsertAfterNode;
      pInsertAfterNode->pNext = pNewNode;
      ui_size++;
    };



void erase(iter_node_t * pOldNode){
      pOldNode->pPrev->pNext = pOldNode->pNext;
      pOldNode->pNext->pPrev = pOldNode->pPrev;
      ui_size--;
    };

    void splice(iter_node_t * pNewLast,uint64_t removed_elements){
      if(removed_elements>0){	
        pSentinel->pPrev = pNewLast;
        pNewLast->pNext = pSentinel;

        ui_size-= removed_elements;
      }
    };

    uint64_t splice(iter_node_t * pNewLast,iter_node_t * pLastRemoved){
      iter_node_t * pCurNode = pNewLast->pNext;
      uint64_t ui_removed_count = 1;
      while(pCurNode != this->end() && pCurNode !=pLastRemoved){pCurNode = pCurNode->pNext; ui_removed_count++;}

      pNewLast->pNext = pLastRemoved->pNext;
      pLastRemoved->pNext->pPrev = pNewLast;


      ui_size-= ui_removed_count;
      return ui_removed_count;
    };


void splice(iter_node_t * pNewLast,iter_node_t * pLastRemoved, uint64_t removed_elements){
      if(removed_elements>0){
        pNewLast->pNext = pLastRemoved->pNext;
        pLastRemoved->pNext->pPrev = pNewLast;

        ui_size-= removed_elements;
      }
    };

    iter_node_t * begin(){
      return pSentinel->pNext;
    };
    iter_node_t * last(){
      return pSentinel->pPrev;
    };
    iter_node_t * end(){
      return pSentinel;
    };

uint64_t size(){
      return ui_size;
    };
};

        
void MinMem(Tree * tree, double MaxOutDeg, double & Required_memory, schedule_traversal & Schedule, int quiet);//added by Changjiang

void MinMemArray(int N,int root, double * nwghts, double * ewghts, int * chstart, int * children,double MaxOutDeg, double & Required_memory, int * Schedule, int quiet, int & count);

#endif
