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

#ifdef __cplusplus


#include <list>
#include "lib-io-tree-utils.h"

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
#ifdef DEBUG_USING_MINMEM
  double mavail;
  double subcut_value;
#endif
#ifdef DEBUG_LIU
  int liu_rank;
#endif
  iter_node_t(){
    index = -1;
    mpeak = 0;
    cost = -1;
#ifdef DEBUG_USING_MINMEM
    subcut_value = numeric_limits<double>::infinity();
    mavail = 0;
#endif
#ifdef DEBUG_LIU
    liu_rank = -1;
#endif
  }

  iter_node_t * pSCbegin;
  iter_node_t * pSCend;
#ifndef USE_DICHOTOMY
  iter_node_t * pNext;
  iter_node_t * pPrev;
#endif
};

#ifdef USE_DICHOTOMY
typedef struct s_list_item_t {
	iter_node_t * node;
	s_list_item_t * pNext;
	s_list_item_t * pPrev;
} s_list_item_t;
#endif


class MinMemDLL {
  protected:

#ifdef USE_DICHOTOMY
    s_list_item_t * items;
    s_list_item_t * pSentinel;
    uint64_t ui_max_count;
#else
    iter_node_t * pSentinel;
#endif
    bool custom_sentinel;
    uint64_t ui_size;
  public:
#ifdef USE_DICHOTOMY
    MinMemDLL(uint64_t max_count){
      pSentinel = new s_list_item_t;
      pSentinel->pNext = pSentinel;
      pSentinel->pPrev = pSentinel;
      custom_sentinel = false;
      ui_size = 0;
      items = new s_list_item_t[max_count];
      ui_max_count = max_count;
    };

#else
    MinMemDLL(){
      pSentinel = new iter_node_t;
      pSentinel->pNext = pSentinel;
      pSentinel->pPrev = pSentinel;
      custom_sentinel = false;
      ui_size = 0;
    };
#endif

#ifdef USE_DICHOTOMY
    MinMemDLL(iter_node_t * pSentinel,uint64_t max_count){
	items = new s_list_item_t[max_count];
      this->pSentinel = &items[pSentinel->index];
      this->pSentinel->node = pSentinel;
      this->pSentinel->pNext = this->pSentinel;
      this->pSentinel->pPrev = this->pSentinel;
      custom_sentinel = true;
      ui_size = 0;
      ui_max_count = max_count;
    };
#else
    MinMemDLL(iter_node_t * pSentinel){
      this->pSentinel = pSentinel;
      pSentinel->pNext = pSentinel;
      pSentinel->pPrev = pSentinel;
      custom_sentinel = true;
      ui_size = 0;
    };
#endif

    ~MinMemDLL(){
      if(!custom_sentinel){
        delete pSentinel;
      }

#ifdef USE_DICHOTOMY
	delete[] items;
#endif
    };


    void clear(){
      pSentinel->pNext = pSentinel;
      pSentinel->pPrev = pSentinel;
      ui_size = 0;
    }

void push_back(iter_node_t * pNewNode){
#ifdef USE_DICHOTOMY
      items[pNewNode->index].node = pNewNode;
      pSentinel->pPrev->pNext = &items[pNewNode->index];
      items[pNewNode->index].pPrev = pSentinel->pPrev;
      items[pNewNode->index].pNext = pSentinel;
      pSentinel->pPrev = &items[pNewNode->index];
      ui_size++;
#else
      pSentinel->pPrev->pNext = pNewNode;
      pNewNode->pPrev = pSentinel->pPrev;
      pNewNode->pNext = pSentinel;
      pSentinel->pPrev = pNewNode;
      ui_size++;
#endif
    };


void insert_after(iter_node_t * pInsertAfterNode, iter_node_t * pNewNode){
#ifdef USE_DICHOTOMY
      items[pNewNode->index].node = pNewNode;
      items[pInsertAfterNode->index].pNext->pPrev = &items[pNewNode->index];
      items[pNewNode->index].pNext = items[pInsertAfterNode->index].pNext;
      items[pNewNode->index].pPrev = &items[pInsertAfterNode->index];
      items[pInsertAfterNode->index].pNext = &items[pNewNode->index];
#else
      pInsertAfterNode->pNext->pPrev = pNewNode;
      pNewNode->pNext = pInsertAfterNode->pNext;
      pNewNode->pPrev = pInsertAfterNode;
      pInsertAfterNode->pNext = pNewNode;
#endif
      ui_size++;
    };



void erase(iter_node_t * pOldNode){
#ifdef USE_DICHOTOMY
      items[pOldNode->index].pPrev->pNext = items[pOldNode->index].pNext;
      items[pOldNode->index].pNext->pPrev = items[pOldNode->index].pPrev;
#else
      pOldNode->pPrev->pNext = pOldNode->pNext;
      pOldNode->pNext->pPrev = pOldNode->pPrev;
#endif
      ui_size--;
    };

#ifndef USE_DICHOTOMY
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
#endif


void splice(iter_node_t * pNewLast,iter_node_t * pLastRemoved, uint64_t removed_elements){
      if(removed_elements>0){
#ifdef USE_DICHOTOMY
        items[pNewLast->index].pNext = items[pLastRemoved->index].pNext;
        items[pLastRemoved->index].pNext->pPrev = &items[pNewLast->index];
#else
#if STRONG_ASSERT        
        iter_node_t * pCurNode = pNewLast->pNext;
        uint64_t ui_removed_count = 1;
        while(pCurNode != this->end() && pCurNode !=pLastRemoved){
          ui_removed_count++;


          pCurNode = pCurNode->pNext;


        }

        assert(ui_removed_count == removed_elements);
#endif
        pNewLast->pNext = pLastRemoved->pNext;
        pLastRemoved->pNext->pPrev = pNewLast;
#endif

        ui_size-= removed_elements;
      }
    };

#ifdef USE_DICHOTOMY
    s_list_item_t * begin(){
#else
    iter_node_t * begin(){
#endif
      return pSentinel->pNext;
    };

#ifdef USE_DICHOTOMY
    s_list_item_t * last(){
#else
    iter_node_t * last(){
#endif
      return pSentinel->pPrev;
    };

#ifdef USE_DICHOTOMY
    s_list_item_t * end(){
#else
    iter_node_t * end(){
#endif
      return pSentinel;
    };

uint64_t size(){
      return ui_size;
    };

#ifdef USE_DICHOTOMY
void copy(MinMemDLL * otherDLL){
//  this->clear();
	this->ui_size = otherDLL->ui_size;
//  if(otherDLL->ui_max_count>this->ui_max_count){
//    this->ui_max_count = otherDLL->ui_max_count;
//    this->items = (s_list_item_t*)realloc(this->items, this->ui_max_count*sizeof(s_list_item_t));
//  }

  this->items[otherDLL->pSentinel->node->index].node = otherDLL->pSentinel->node;
  this->items[otherDLL->pSentinel->node->index].pNext = &this->items[otherDLL->pSentinel->pNext->node->index];
  this->items[otherDLL->pSentinel->node->index].pPrev = &this->items[otherDLL->pSentinel->pPrev->node->index];
  this->pSentinel = &this->items[otherDLL->pSentinel->node->index];

  s_list_item_t * cur_item = otherDLL->begin(); 
  while(cur_item!=otherDLL->end()){
    this->items[cur_item->node->index].node = cur_item->node;
    this->items[cur_item->node->index].pNext = &this->items[cur_item->pNext->node->index];
    this->items[cur_item->node->index].pPrev = &this->items[cur_item->pPrev->node->index];
    cur_item = cur_item->pNext;
  }


//  cur_item = otherDLL->begin(); 
//  cerr<<endl;
//  while(cur_item!=otherDLL->end()){
//    cerr<<" "<<this->items[cur_item->node->index].node->index<<"("<<this->items[cur_item->node->index].pPrev->node->index<<"|"<<this->items[cur_item->node->index].pNext->node->index<<")";
//    cur_item = cur_item->pNext;
//  }
//  cerr<<endl;


}

s_list_item_t * getItem(iter_node_t * node){
	return &items[node->index];
}
 
#endif
};



//void MinMem(Ctree * tree, double MaxOutDeg , double & Required_memory, list<unsigned int> & Schedule, int quiet, int & count);
        
void MinMem(Ctree * tree,double MaxOutDeg, double & Required_memory, schedule_t & Schedule, int quiet, uint64_t & count);//added by Changjiang

void MinMemArray(int N,int root, double * nwghts, double * ewghts, int * chstart, int * children,double MaxOutDeg, double & Required_memory, int * Schedule, int quiet, int & count);

#endif

#endif
