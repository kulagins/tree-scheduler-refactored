/*
 *  lib-io-tree-liu-optimal.cpp
 *  lib-io-tree
 *
 *  Created by defbond on 8/10/10.
 *  Copyright 2010 LIP/ENS-Lyon. All rights reserved.
 *
 */
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <sys/time.h>

#include <iostream>
#include <fstream>
#include <limits>


#include "lib-io-tree.h"
#include "lib-io-tree-liu-optimal.h"

/*** LIU algorithm *****/
void H(Ctree * tree, schedule_t::iterator vi, double Vi,  OrdoLiu_t & ordo,schedule_t::iterator & hi,int & offset, double & Hi,int quiet,int & count){
  /*compute initial Hi value*/
  Hi = Vi;
  hi = vi;
  double prev_value = Vi;

  /*start exploration after the previous valley*/
  vi++;
  for (schedule_t::iterator current_item = vi; current_item!=ordo.schedule.end(); current_item++){
    offset++;
    //		count++;
    double current_value =prev_value + tree->GetNode(*current_item)->GetNW();
    for(unsigned int i=0;i<tree->GetNode(*current_item)->GetChildren()->size();i++){
      current_value -= tree->GetNode(*current_item)->GetChild(i)->GetNW();

      /*to adapt to our model consider peaknode (fi + sum(fj) +ni == getcost()) and then edge node (fi = getEW())*/ 
    }

    if(current_value<0){current_value = 0;}

    if(current_value>=Hi){
      Hi = current_value;
      hi = current_item;
      //			if(!quiet){cerr<<"hi = "<<*hi<<"("<<*current_item<<") | Hi = "<<Hi<<"| Prev = "<<prev_value<<endl;}
    }

    prev_value = current_value;
  }
}

void V(Ctree * tree, schedule_t::iterator hi, double Hi, OrdoLiu_t & ordo, schedule_t::iterator &vi, int & offset,double & Vi, int quiet,int & count){
  Vi = Hi;
  vi = hi;
  double prev_value = Hi;

  /*start exploration after the previous valley*/
  hi++;

  for (schedule_t::iterator current_item = hi; current_item!=ordo.schedule.end(); current_item++){
    offset++;
    double current_value =prev_value + tree->GetNode(*current_item)->GetNW();
    for(unsigned int i=0;i<tree->GetNode(*current_item)->GetChildren()->size();i++){
      current_value -= tree->GetNode(*current_item)->GetChild(i)->GetNW();
    }
    if(current_value<0){current_value = 0;}

    if(current_value<=Vi){
      Vi = current_value;
      vi = current_item;
      //			if(!quiet){cerr<<"vi = "<<*vi<<"("<<*current_item<<") | Vi = "<<Vi<<"| Prev_value = "<<prev_value<<endl;}
    }
    prev_value = current_value;
  }
}


void PCost(Ctree * tree,int root,OrdoLiu_t & ordo, int quiet, int & count){
  ordo.val_seg.clear();

  schedule_t::iterator hi,vi,prev_valley;
  double Hi,Vi;


  if(ordo.schedule.size()==1){
    ordo.max_pebble_cost = tree->GetNode(*(ordo.schedule.begin()))->GetNW();


    val_seg_t segment;
    segment.begin = ordo.schedule.begin();
    segment.begin_index = 0;
    segment.end_index = 0;
    segment.end = ordo.schedule.begin();
    segment.orig_ordo = & ordo;
    segment.value = 0;
    ordo.val_seg.push_back(segment);

  }
  else{


    /* compute Pcost*/
    /********** working version **********/

    if(!quiet){cerr<<endl;}
    // compute Pcost
    vi = ordo.schedule.begin();
    Vi = tree->GetNode(*vi)->GetNW();
    ordo.max_pebble_cost = 0;
    int i = 0;
    int offset = 0;
    while(*vi!=root){
      val_seg_t segment;
      prev_valley = vi;
      segment.begin_index = offset;
      H(tree, vi, Vi, ordo, hi, offset, Hi,quiet,count);
      V(tree, hi, Hi, ordo, vi, offset, Vi,quiet,count);

      segment.begin = prev_valley;
      if(i>0){
        segment.begin_index++;
        segment.begin++;
      }


      segment.end = vi;
      segment.end_index = (int)offset;
      segment.orig_ordo = & ordo;
      segment.value = Hi - Vi;
      ordo.val_seg.push_back(segment);

      if(Hi>ordo.max_pebble_cost){
        ordo.max_pebble_cost = Hi;
      }

      if(!quiet){cerr<<"\t i ="<<i<<" | vi = "<<*vi<<"  | Vi = "<<Vi<<" | hi = "<<*hi<<" Hi = |  "<<Hi<<endl;}
      i++;
    }

    /********** End of working version of PCost **************/
  }

  if(!quiet){	

    cerr<<"Segments of T"<<root<<endl;
    for (list<val_seg_t>::iterator cur_seg = ordo.val_seg.begin();cur_seg!=ordo.val_seg.end();cur_seg++){
      cerr<<"| "<<*(cur_seg->begin)<<"["<<cur_seg->begin_index<<"]("<<cur_seg->value<<")"<<*(cur_seg->end)<<"["<<cur_seg->end_index<<"]";
    }
    cerr<<endl;

  }

}

bool desc (val_seg_t i,val_seg_t j) {
  if(i.orig_ordo == j.orig_ordo){
    return (i.begin_index < j.begin_index);
    //return 0;
  }

  return (i.value>j.value); 

}


void Combine(Ctree * tree, unsigned int sub_root, vector<OrdoLiu_t> * ordo_subtrees, OrdoLiu_t & ordo, int quiet, int & count){
  ordo.schedule.clear();
  ordo.val_seg.clear();	
  /*compute pcosts of subtrees*/
  for(unsigned int i=0;i<tree->GetNode(sub_root)->GetChildren()->size();i++){
    if(!quiet){cerr<<"Computing PCost of T"<<tree->GetNode(sub_root)->GetChild(i)->GetId()<<endl;}
    //		PCost(tree, tree->GetNode(sub_root)->GetChild(i)->GetId(),ordo_subtrees->at(i),quiet, count);
    //ordo.val_seg.insert(ordo.val_seg.end(),ordo_subtrees->at(i).val_seg.begin(),ordo_subtrees->at(i).val_seg.end());
    ordo.val_seg.splice(ordo.val_seg.end(),ordo_subtrees->at(i).val_seg);
  }

  if(!quiet){
    cerr<<"Unsorted segments"<<endl;
    for (list<val_seg_t>::iterator cur_seg = ordo.val_seg.begin();cur_seg!=ordo.val_seg.end();cur_seg++){
      if(!quiet){cerr<<"| "<<*(cur_seg->begin)<<"("<<cur_seg->value<<")"<<*(cur_seg->end);}
    }
    cerr<<endl;
  }

  /*sort the segment according to their value in descending order*/
  ordo.val_seg.sort(desc);

  //	count += ordo.val_seg.size() * log(ordo.val_seg.size());

  /*merge*/
  if(!quiet){cerr<<endl<<"Merging Segments of T"<<sub_root<<endl;}
  for (list<val_seg_t>::iterator cur_seg = ordo.val_seg.begin();cur_seg!=ordo.val_seg.end();cur_seg++){
    //OrdoLiu_t * cur_ordo = ordo.val_seg[j].orig_ordo;
    if(!quiet){cerr<<"| "<<*(cur_seg->begin)<<"("<<cur_seg->value<<")"<<*(cur_seg->end);}

    /*explore until we reach the last valley*/
    schedule_t::iterator end = cur_seg->end ;
    end++;

    ordo.schedule.insert(ordo.schedule.end(),cur_seg->begin,end);
  }

  //val_seg_t subroot_segment;
  //subroot_segment->pStart	= 
  ordo.schedule.push_back(sub_root);


  if(!quiet){cerr<<endl;}

}


//void PebbleOrderingRecur(Ctree * tree,unsigned int sub_root, OrdoLiu_t & SubSchedule, int quiet, int & count){
//
//  count++;
//
//  SubSchedule.schedule.clear();
//
//  if (tree->GetNode(sub_root)->IsLeaf()){
//    if(!quiet){cerr<<"Treating T"<<sub_root<<endl;}
//    if(!quiet){cerr<<"Leaf found : "<<sub_root<<endl;}
//    SubSchedule.schedule.push_back(sub_root);
//  }
//  else{
//    vector<OrdoLiu_t> * ordo_subtrees = new vector<OrdoLiu_t>(tree->GetNode(sub_root)->GetChildren()->size());
//
//    for(unsigned int i=0;i<tree->GetNode(sub_root)->GetChildren()->size();i++){
//      PebbleOrderingRecur(tree,tree->GetNode(sub_root)->GetChild(i)->GetId(), ordo_subtrees->at(i), quiet, count);
//      if(!quiet){
//        cerr<<"sub schedule of T"<<tree->GetNode(sub_root)->GetChild(i)->GetId()<<endl;
//        for(schedule_t::iterator cur_item=ordo_subtrees->at(i).schedule.begin();cur_item!=ordo_subtrees->at(i).schedule.end();cur_item++){
//          cerr<<*cur_item<<" ";
//        }
//        cerr<<endl;
//      }
//    }
//
//    if(!quiet){cerr<<"Treating T"<<sub_root<<endl;}
//    Combine(tree, sub_root, ordo_subtrees, SubSchedule, quiet, count);
//    delete ordo_subtrees;
//  }
//  //Compute pebble cost
//  PCost(tree, sub_root,SubSchedule,quiet, count);
//  //    PCostSeg(tree, sub_root,ordo.val_seg,quiet, count);
//
//  if(!quiet){cerr<<"End of treatment of T"<<sub_root<<endl;}
//}







//double PebbleOrderingRecurAlgorithm(int N, int *prnts, double *nwghts, double *ewghts, int *schedule){
//  Ctree * tree = new Ctree(N,prnts,nwghts,ewghts);
//  Ctree * tree_liu = new Ctree();
//
//  ConvertToLiu(tree,tree_liu);
//  delete tree;
//
//  double Mr;
//  int count = 0;
//  OrdoLiu_t ordo;
//  ordo.max_pebble_cost = 0;
//  ordo.fi = 0;
//
//  PebbleOrderingRecur(tree_liu,tree_liu->GetRootId(), ordo, true, count);
//  Mr = ordo.max_pebble_cost;
//
//  unsigned int i = 0;
//  for (schedule_t::iterator last=ordo.schedule.begin(); last!=ordo.schedule.end(); ++last){
//    if((int)(*last)<=N){
//      schedule[i++] = (int)(*last);
//    }
//  }
//
//  delete tree_liu;
//
//  return Mr;
//}


//double PebbleOrderingRecurAlgorithm_timed(int N, int *prnts, double *nwghts, double *ewghts, int *schedule,double * usec,int quiet){
//  Ctree * tree = new Ctree(N,prnts,nwghts,ewghts);
//  Ctree * tree_liu = new Ctree();
//
//  ConvertToLiu(tree,tree_liu);
//
//
//  delete tree;
//
//  double Mr;
//  int count = 0;
//  OrdoLiu_t ordo;
//  ordo.max_pebble_cost = 0;
//  ordo.fi = 0;
//
//  *usec = -u_wseconds();
//  PebbleOrderingRecur(tree_liu,tree_liu->GetRootId(), ordo, quiet, count);
//  *usec += u_wseconds();
//  Mr = ordo.max_pebble_cost;
//
//
//  unsigned int i = 0;
//  for (schedule_t::iterator last=ordo.schedule.begin(); last!=ordo.schedule.end(); ++last){
//    if((int)(*last)<=N){
//      schedule[i++] = (int)(*last);
//    }
//  }
//
//  delete tree_liu;
//
//  return Mr;
//}



int ascend_liu_comp(const void * a, const void * b){
  const iter_seg_t * pa = (iter_seg_t *)a;
  const iter_seg_t * pb = (iter_seg_t *)b;
  double diff = pb->value - pa->value;


  if(pa->root == pb->root){
    return (pa->seq_index-pb->seq_index);
  }
  else{

    if(diff==0.0){
      return 0;
    }
    else if (diff>0.0){
      return 1;
    }
    else{
      return -1;
    }

  }
}

void HIter(iter_node_t * vi,double Vi, iter_node_t ** phi, double & Hi,const double * nwghts, int root){
  /*compute initial Hi value*/
  Hi = Vi;
  *phi = vi;
  double peb = Vi;

  /*start exploration after the previous valley*/
  iter_node_t * cur_node = vi;

  do{
    cur_node = cur_node->pNext;
    peb += nwghts[cur_node->index] - cur_node->cw;

    if(peb<0){peb = 0;}

    if(peb>=Hi){
      Hi = peb;
      *phi = cur_node;
    }
  }
  while(cur_node->index != root);


#if VERBOSE	
  cerr<<"New Hill found "<<(*phi)->index<<endl;
#endif

}

void VIter(iter_node_t * hi,double Hi, iter_node_t ** pvi, double & Vi, const double * nwghts, int root,int quiet){
  /*compute initial Vi value*/
  /*start exploration at the previous hill*/
  Vi = Hi;
  *pvi = hi;
  double peb = Hi;

  iter_node_t * cur_node = hi;

  do{
    cur_node = cur_node->pNext;
    peb += nwghts[cur_node->index] - cur_node->cw;

    if(peb<0){peb = 0;}

    if(peb<=Vi){
      Vi = peb;
      *pvi = cur_node;
    }
  }
  while(cur_node->index != root);


#if VERBOSE	
  cerr<<"New Valley found "<<(*pvi)->index<<endl;
#endif
}

void H1Iterseg(int segVi_index, int & segHi_index,iter_seg_t * segments, int segment_count,int root){
  /*compute initial Hi value*/
  segHi_index = segVi_index;


  double peb = 0;

  int prevroot = -1;


  for(int i=segVi_index;i<segment_count;i++){
    iter_seg_t * cur_segment = &segments[i];

    //if this is the first exploration, we need to update hills
    peb += cur_segment->standalone_Hi;		
    if(i==(segment_count-1)){			
      peb -= cur_segment->pStart->cw;
    }

    if(peb<0){peb = 0;}
    cur_segment->Hi = peb;
#if VERBOSE	
    cerr<<"Candidate for next hill "<<cur_segment->hi->index<<" peb = "<<peb<<endl;
#endif

    if(peb>=segments[segHi_index].Hi){
      segHi_index= i;
    }

    /*go to valley*/
    peb -= cur_segment->value;
    if(peb<0){peb = 0;}

    //go to next segment (they are organized as a stack so its easy)
    prevroot = cur_segment->root;
  }

#if VERBOSE	
  cerr<<"SEG New Hill found : "<<segments[segHi_index].hi->index<<endl;
#endif
}

void HIterseg(int segVi_index, int & segHi_index,iter_seg_t * segments, int segment_count,int root){
  /*compute initial Hi value*/
  segHi_index = segVi_index;


  double peb = 0;//(*psegHi)->Hi - (*psegHi)->value;

  for(int i=segVi_index;i<segment_count;i++){
    iter_seg_t * cur_segment = &segments[i];

    peb = cur_segment->Hi;


#if VERBOSE	
    cerr<<"Candidate for next hill "<<cur_segment->hi->index<<" peb = "<<peb<<endl;
#endif

    if(peb>=segments[segHi_index].Hi){
      segHi_index= i;
    }

  }

#if VERBOSE	
  cerr<<"SEG New Hill found : "<<segments[segHi_index].hi->index<<endl;
#endif
}

void VIterseg(int segHi_index, int & segVi_index,iter_seg_t * segments, int segment_count,int root){
  /*compute initial Vi value*/
  segVi_index = segHi_index;


  double peb = 0;//(*psegHi)->Hi - (*psegHi)->value;

  for(int i=segHi_index;i<segment_count;i++){
    iter_seg_t * cur_segment = &segments[i];

    peb = cur_segment->Hi-cur_segment->value;


#if VERBOSE	
    cerr<<"Candidate for next valley "<<cur_segment->pEnd->index<<" peb = "<<peb<<endl;
#endif

    if(peb<=segments[segVi_index].Hi-segments[segVi_index].value){
      segVi_index= i;
    }

  }

#if VERBOSE	
  cerr<<"SEG New Valley found : "<<segments[segVi_index].pEnd->index<<endl;
#endif
}

int PCostIterative(const int * prnts, const double * nwghts, iter_node_t * schedlist, int N,int subroot,iter_seg_t *segments,int segment_count){

#if VERBOSE	
  cerr<<"\tComputing PCost of node "<<subroot<<" and merging "<<segment_count<<" segments"<<endl;
#endif
  //update the schedule so as to correspond to the segment order

  for(int is = 0; is<segment_count-1;is++){
    iter_seg_t * cur_segment = &segments[is];
    iter_seg_t * next_segment = &segments[is+1];
#if VERBOSE	
    cerr<<"\t\t "<<is<<"th segment = ("<<cur_segment->pStart->index<<" - "<<cur_segment->pEnd->index<<"), val = "<<cur_segment->value<<endl;
#endif
    cur_segment->pEnd->pNext = next_segment->pStart;
  }

  assert(segments[segment_count-1].pEnd->index == subroot);


#if VERBOSE	
  cerr<<"\t\t "<<segment_count-1<<"th segment = ("<<segments[segment_count-1].pStart->index<<" - "<<segments[segment_count-1].pEnd->index<<"), val = "<<segments[segment_count-1].value<<endl;
#endif

#if VERBOSE	
  cerr<<endl<<"\tNew schedule : {";
  for(int is = 0; is<segment_count;is++){
    iter_seg_t * cur_segment = &segments[is];
    iter_node_t * cur_item = cur_segment->pStart;			
    while(cur_item != cur_segment->pEnd){
      cerr<<cur_item->index<<"("<<nwghts[cur_item->index]<<") ["<<cur_segment->root<<"] ";
      cur_item = cur_item->pNext;	
    }
    cerr<<cur_item->index<<"("<<nwghts[cur_item->index]<<") ["<<cur_segment->root<<"] ";
  }
  cerr<<"}"<<endl;
#endif


  //compute the segments of the new schedule
  iter_node_t *prev_valley;

  /* compute Pcost*/

  /********** working version **********/

#if VERBOSE	
  cerr<<endl;
#endif
  // compute Pcost
  iter_seg_t * segVi = &segments[0];
  iter_seg_t * segHi = &segments[0];
  prev_valley = segVi->pStart;

  iter_seg_t new_segment;

  int i = 0;
  int new_local_count = 0;


  //search for first segment
  int segVi_index = 0;
  int segHi_index =0;

  H1Iterseg(segVi_index, segHi_index,segments, segment_count,subroot);
  segHi = &segments[segHi_index];
  VIterseg(segHi_index, segVi_index,segments, segment_count,subroot);
  segVi = &segments[segVi_index];

  new_segment.pStart = prev_valley;	
  new_segment.pEnd = segVi->pEnd;
  new_segment.root = subroot;
  new_segment.value = segHi->Hi - segVi->Hi + segVi->value;
  new_segment.Hi = segHi->Hi;
  new_segment.standalone_Hi = segHi->Hi;
  new_segment.hi = segHi->hi;
  new_segment.seq_index = new_local_count;

  if(new_segment.pEnd->index!=subroot){
    prev_valley = (segVi+1)->pStart;
  }	

  //we start from the valley
  segVi->Hi -= segVi->value;
  segVi->hi = segVi->pEnd;
  segVi->value = 0;
  i++;

  while(new_segment.pEnd->index!=subroot){	
    HIterseg(segVi_index, segHi_index,segments, segment_count,subroot);
    segHi = &segments[segHi_index];
    VIterseg(segHi_index, segVi_index,segments, segment_count,subroot);
    segVi = &segments[segVi_index];

    segments[new_local_count++] = new_segment;

    new_segment.pStart = prev_valley;
    new_segment.pEnd = segVi->pEnd;
    new_segment.root = subroot;
    new_segment.value = segHi->Hi - segVi->Hi + segVi->value;
    new_segment.Hi = segHi->Hi;
    new_segment.standalone_Hi = segHi->Hi - segments[new_local_count-1].Hi + segments[new_local_count-1].value;
    new_segment.hi = segHi->hi;		
    new_segment.seq_index = new_local_count;

#if VERBOSE	
    cerr<<"\t i ="<<i<<" | vi = "<<segVi->pEnd->index<<"  | Vi = "<<segVi->Hi-segVi->value<<" | hi = "<<segHi->hi->index<<" Hi = |  "<<segHi->Hi<<endl;
#endif

    //node next to the previous valley		
    if(new_segment.pEnd->index!=subroot){
      prev_valley = (segVi+1)->pStart;
    }

    //we start from the valley
    segVi->Hi -= segVi->value;
    segVi->hi = segVi->pEnd;
    segVi->value = 0;

    i++;
  }


  segments[new_local_count++] = new_segment;

  /********** End of working version of PCost **************/
#if VERBOSE	
  cerr<<"New segment count = "<<new_local_count<<endl;
#endif

  return new_local_count-segment_count;
}



//double PebbleOrderingIter(const int N,const int *prnts,const double *nwghts,const double *ewghts,const int * chstart, const int * chend, int * children, const int root, int *schedule){
//
//  /*chstart and chend are indexed from 1 to N inclusive*/
//
//  int * po = new int[N+1];
//  memset ( (void *) po, 0, (N+1)*sizeof(*po) );
//
//  int label =1;
//  poaux(chstart, children, N, root, po, &label);
//
//  //now that we have a PO, we can process the tree from the leaves up to the root
//
//
//  iter_seg_t * segments = new iter_seg_t[N+1];
//
//  //    int first_segment = 0;
//  int total_segment_count = 0;
//
//
//  iter_node_t * schedlist = new iter_node_t[N+1];
//  int * first_segment = new int[N+1];
//  memset ( (void *) first_segment, 0, (N+1)*sizeof(*first_segment) );
//
//
//  /*Post order is done, compute weight then sort children*/
//  for(int i = 1;i<N+1;i++){
//    int nd = po[i];
//#if VERBOSE
//    cerr<<endl<<endl<<"processing node "<<nd<<endl;
//#endif
//
//    /*if nd is a leaf*/
//    if(chstart[nd]==chstart[nd+1]){
//      //we know that we add one segment
//      schedlist[nd].index = nd;
//      schedlist[nd].pNext = &schedlist[nd];
//      iter_seg_t * my_segment = &segments[total_segment_count];
//      first_segment[nd]=total_segment_count;
//#if VERBOSE
//      cerr<<"\tFirst Segment is in "<<first_segment[nd]<<"th position"<<endl;
//#endif
//      total_segment_count++;
//      schedlist[nd].cw = 0;
//      my_segment->pStart = &schedlist[nd];
//      my_segment->pEnd = &schedlist[nd];
//      my_segment->Hi = nwghts[nd];
//      my_segment->standalone_Hi = my_segment->Hi;
//      my_segment->hi = &schedlist[nd];
//      my_segment->root = nd;
//      my_segment->value = 0;
//      my_segment->seq_index = 0;
//
//#if VERBOSE
//      cerr<<"node "<<nd<<" is a leaf of w="<<nwghts[nd]<<endl;
//#endif
//    }
//    else{
//
//      //the first segment is the first segment of its first child
//      first_segment[nd]=first_segment[children[chstart[nd]]];
//#if VERBOSE
//      cerr<<"\tFirst Segment is in "<<first_segment[nd]<<"th position"<<endl;
//#endif
//
//      //Combine the segments of my children (there
//      //sort the segment according to their value in descending order
//      qsort ( &segments[first_segment[nd]], (total_segment_count - first_segment[nd]), sizeof(iter_seg_t), ascend_liu_comp );
//
//      //add a last segment containing the root
//      schedlist[nd].index = nd;
//      schedlist[nd].pNext = &schedlist[nd];
//      schedlist[nd].cw = 0;
//      for(int j = chstart[nd]; j<chstart[nd+1];j++){
//        int ch = children[j];
//#if VERBOSE
//        cerr<<"\t"<<ch<<" is a child of "<<nd<<endl;
//#endif
//        schedlist[nd].cw += nwghts[ch];
//      }
//
//
//      iter_seg_t * my_segment = &segments[total_segment_count];
//      total_segment_count++;
//      my_segment->pStart = &schedlist[nd];
//      my_segment->pEnd = &schedlist[nd];
//      my_segment->Hi = nwghts[nd];
//      my_segment->standalone_Hi = my_segment->Hi;
//      my_segment->hi = &schedlist[nd];
//      my_segment->root = nd;
//      my_segment->value = 0;
//      my_segment->seq_index = 0;
//
//
//      //NB : PCost will construct the schedule from the segment list and return the new segments count
//      total_segment_count+= PCostIterative(prnts, nwghts,schedlist,N, nd,&segments[first_segment[nd]],(total_segment_count-first_segment[nd]));
//
//#if VERBOSE
//      cerr<<"node "<<nd<<" is not a leaf of w="<<nwghts[nd]<<endl;
//#endif
//
//    }
//
//  }
//
//
//  double Mr = segments[0].Hi;
//
//
//  int offset = 0;
//  iter_node_t * cur_item = segments[0].pStart;
//
//  if(cur_item->index<=0.5*N){
//    schedule[offset++] = cur_item->index;
//  }
//  while(cur_item != cur_item->pNext){
//    cur_item = cur_item->pNext;
//    if(cur_item->index<=0.5*N){
//      schedule[offset++] = cur_item->index;
//    }
//  }
//
//
//  delete[] first_segment;
//  delete[] segments;
//  delete[] schedlist;
//  delete[] po;
//  return Mr;
//
//}

//double PebbleOrderingIterAlgorithm(int N, int *prnts, double *nwghts, double *ewghts, int *schedule){
//  double Mr;
//  //int count = 0;
//  int * chstart,*chend,*children;
//  int root;
//
//
//  po_construct(N, prnts, &chstart,&chend,&children, &root);
//
//
//  int * prnts_liu;
//  double * nwghts_liu;
//  double * ewghts_liu;
//
//  ConvertToLiu(prnts,nwghts,ewghts,N,chstart,children,&prnts_liu,&nwghts_liu,&ewghts_liu);
//
//  delete[] chstart;
//  delete[] chend;
//  delete[] children;
//
//  //update chstart et children
//  po_construct(2*N, prnts_liu, &chstart,&chend,&children, &root);
//
//
//  Mr = PebbleOrderingIter(2*N,prnts_liu,nwghts_liu,ewghts_liu,chstart,chend,children,root,schedule);
//
//
//  delete[] chstart;
//  delete[] chend;
//  delete[] children;
//
//  delete[] prnts_liu;
//  delete[] nwghts_liu;
//  delete[] ewghts_liu;
//
//  return Mr;
//}

//double PebbleOrderingIterAlgorithm_timed(int N, int *prnts, double *nwghts, double *ewghts, int *schedule, double * usec,int quiet){
//  double Mr;
//  //int count = 0;
//  int * chstart,*chend,*children;
//  int root;
//
//
//  po_construct(N, prnts, &chstart,&chend,&children, &root);
//
//  int * prnts_liu;
//  double * nwghts_liu;
//  double * ewghts_liu;
//
//  ConvertToLiu(prnts,nwghts,ewghts,N,chstart,children,&prnts_liu,&nwghts_liu,&ewghts_liu);
//
//  delete[] chstart;
//  delete[] chend;
//  delete[] children;
//
//  //update chstart et children
//  po_construct(2*N, prnts_liu, &chstart,&chend,&children, &root);
//
//  *usec = -u_wseconds();
//  Mr = PebbleOrderingIter(2*N,prnts_liu,nwghts_liu,ewghts_liu,chstart,chend,children,root,schedule);
//  *usec += u_wseconds();
//
//
//  delete[] chstart;
//  delete[] chend;
//  delete[] children;
//  delete[] prnts_liu;
//  delete[] nwghts_liu;
//  delete[] ewghts_liu;
//
//
//  return Mr;
//}





