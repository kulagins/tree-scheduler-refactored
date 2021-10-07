/*
 *  lib-io-tree-minmem.cpp
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
#include <sys/resource.h> 


#include <iostream>
#include <fstream>
#include <limits>


#include "lib-io-tree.h"
#include "lib-io-tree-minmem.h"



#ifdef DEBUG_USING_MINMEM
void explore(Node * node, double available_memory,list<Node*> * L_init, schedule_t * S_init,  double & cut_value, list<Node*> & min_sub_cut, schedule_t & sub_schedule, double & Mpeak, int quiet, int depth,uint64_t & count,iter_node_t * minmem_trace,uint64_t N)
#else
void explore(Node * node, double available_memory,list<Node*> * L_init, schedule_t * S_init,  double & cut_value, list<Node*> & min_sub_cut, schedule_t & sub_schedule, double & Mpeak, int quiet, int depth,uint64_t & count)
#endif
{

#ifdef DEBUG_USING_MINMEM
  if(count<2*N){
    minmem_trace[count].index = node->GetId();
    minmem_trace[count].mavail = available_memory;
  }
  uint64_t mycount = count;
#endif
  count++;


#if VERBOSE
  char * spacing;
  spacing = (char *)malloc((depth*2+1)*sizeof(char));
  for(int i=0;i<depth*2;++i){spacing[i]=' ';}
  spacing[depth*2]='\0';
#endif
#if VERBOSE
  cerr<<spacing<<"[ node "<<node->GetId()<<" ] is being explored with memory " <<available_memory-node->GetEW()<<" while its cost is "<<node->GetCost()<<endl;
#endif


  /* if node is unreachable, return +infty */
  if (node->GetCost() > available_memory) {
#if VERBOSE
    cerr<<spacing<<"[ node "<<node->GetId()<<" ] not enough memory"<<endl;
#endif
    Mpeak = node->GetCost();
    cut_value = numeric_limits<double>::infinity( );
#ifdef DEBUG_USING_MINMEM
    if(mycount<2*N){
      minmem_trace[mycount].subcut_value = cut_value;
      minmem_trace[mycount].mpeak = Mpeak;
    }
#endif
    return;
  }


  /* if this is a leaf, return 0 */
  if (node->IsLeaf()){
#if VERBOSE
    cerr<<mycount<<" "<<"[ node "<<node->GetId()<<" ] is a leaf"<<endl;
#endif
    sub_schedule.push_back(node->GetId());
    cut_value = 0;
    Mpeak = numeric_limits<double>::infinity( );
#ifdef DEBUG_USING_MINMEM
    if(mycount<2*N){
      minmem_trace[mycount].subcut_value = cut_value;
      minmem_trace[mycount].mpeak = Mpeak;
    }
#endif
    return;
  }

  if(L_init != NULL){
    if(L_init->size()>0){
      min_sub_cut.assign(L_init->begin(),L_init->end());
      sub_schedule.assign(S_init->begin(),S_init->end());
    }
    else{
      sub_schedule.push_back(node->GetId());
      /* place every child in the candidate nodes*/
      min_sub_cut.assign(node->GetChildren()->begin(),node->GetChildren()->end());
    }
  }
  else{
    sub_schedule.push_back(node->GetId());
    /* place every child in the candidate nodes*/
    min_sub_cut.assign(node->GetChildren()->begin(),node->GetChildren()->end());
  }

  list<Node*> * candidates = new list<Node*>(min_sub_cut);


#if  VERBOSE
  if (node->GetChildren()->size()>1) {
    cerr<<spacing<<"[ node "<<node->GetId()<<" ] initial subcut is [";
    for (list<Node*>::iterator current_node=candidates->begin(); current_node!=candidates->end(); ++current_node){
      cerr<<" "<<(*current_node)->GetId()<<"("<<(*current_node)->Mpeak<<")";
    }
    cerr<<" ]"<<endl;
  }
#endif


  cut_value = 0;
  for (list<Node*>::iterator current_node=candidates->begin(); current_node!=candidates->end(); ++current_node){
    (*current_node)->Mavail = available_memory;
    for (list<Node*>::iterator other_nodes=min_sub_cut.begin(); other_nodes!=min_sub_cut.end(); ++other_nodes){
      if((*other_nodes)->GetId() != (*current_node)->GetId()){
        (*current_node)->Mavail -= (*other_nodes)->GetEW();
      }
    }
    cut_value +=(*current_node)->GetEW();
  }

  while (!candidates->empty()) {

#if VERBOSE
    cerr<<spacing<<"*****************************************************************"<<endl;
    cerr<<spacing<<"[ node "<<node->GetId()<<" ] candidates are [";
    for (list<Node*>::iterator current_node=candidates->begin(); current_node!=candidates->end(); ++current_node){
      cerr<<" "<<(*current_node)->GetId()<<"("<<(*current_node)->Mpeak<<"|"<<(*current_node)->Mavail <<")";
    }
    cerr<<" ]"<<endl;
#endif


    for (list<Node*>::iterator current_node=candidates->begin(); current_node!=candidates->end(); ++current_node){
      double m_j;
      list<Node *> Lj;
      schedule_t Sj;


#ifdef DEBUG_USING_MINMEM
      explore(*current_node, (*current_node)->Mavail,NULL,NULL, m_j, Lj,Sj,(*current_node)->Mpeak,quiet,depth+1,count,minmem_trace,N);
#else
      explore(*current_node, (*current_node)->Mavail,NULL,NULL, m_j, Lj,Sj,(*current_node)->Mpeak,quiet,depth+1,count);
#endif

      if (m_j <= (*current_node)->GetEW()){
#if VERBOSE
        cerr<<spacing<<"  [ node "<<(*current_node)->GetId()<<" ] has been removed ("<<(*current_node)->GetEW()<<" vs "<<m_j<<"), new cut size is "<<min_sub_cut.size()-1<<endl;
#endif
        min_sub_cut.remove(*current_node);
        min_sub_cut.splice(min_sub_cut.end(),Lj);
        sub_schedule.splice(sub_schedule.end(),Sj);
      }
      else{
#if VERBOSE
        cerr<<spacing<<"  [ node "<<(*current_node)->GetId()<<" ] is better than its subtree ("<<(*current_node)->GetEW()<<" vs "<<m_j<<")"<<endl;
#endif
      }
    }

    candidates->clear();
    for (list<Node*>::iterator current_node=min_sub_cut.begin(); current_node!=min_sub_cut.end(); ++current_node){
      (*current_node)->Mavail = available_memory;
      for (list<Node*>::iterator other_nodes=min_sub_cut.begin(); other_nodes!=min_sub_cut.end(); ++other_nodes){
        if((*other_nodes)->GetId() != (*current_node)->GetId() ){
          (*current_node)->Mavail -= (*other_nodes)->GetEW();
        }
      }

      //add this node to candidates
      //			if(!quiet){cerr<<spacing<<"Mpeak = "<<(*current_node)->Mpeak<<" Mavail = "<<(*current_node)->Mavail<<endl;}

      if((*current_node)->Mavail >= (*current_node)->Mpeak)
      {
        candidates->push_back(*current_node);
        //				if(!quiet){cerr<<spacing<<"node "<<(*current_node)->GetId()<<" kept"<<endl;}
      }
    }

#if VERBOSE
    cerr<<spacing<<"[ node "<<node->GetId()<<" ] new subcut value is "<<cut_value<<endl;
    cerr<<spacing<<"[ node "<<node->GetId()<<" ] new subcut is [";

    for (list<Node*>::iterator current_node=min_sub_cut.begin(); current_node!=min_sub_cut.end(); ++current_node){
      double available_memory_after_subroot = available_memory - cut_value + node->GetEW();
      cerr<<" "<<(*current_node)->GetId()<<"("<<(*current_node)->Mpeak<<"|"<<available_memory_after_subroot  + (*current_node)->GetEW()<<")";
    }
    cerr<<"]"<<endl;
    cerr<<spacing<<"*****************************************************************"<<endl;
    cerr<<spacing<<candidates->size()<<" candidates left"<<endl;
#endif
  }



  cut_value = 0;
  Mpeak = numeric_limits<double>::infinity( );
  for (list<Node*>::iterator current_node=min_sub_cut.begin(); current_node!=min_sub_cut.end(); ++current_node){
    cut_value += (*current_node)->GetEW();


    if(Mpeak>(*current_node)->Mpeak + available_memory - (*current_node)->Mavail){
      Mpeak = (*current_node)->Mpeak + available_memory - (*current_node)->Mavail;
    }
  }

#ifdef DEBUG_USING_MINMEM
  if(mycount<2*N){
    minmem_trace[mycount].subcut_value = cut_value;
    minmem_trace[mycount].mpeak = Mpeak;
  }
#endif


#if VERBOSE
  double available_memory_after_subroot = available_memory - cut_value + node->GetEW();
  cerr<<spacing<<"[ node "<<node->GetId()<<" ] final subcut is [";
  for (list<Node*>::iterator last=min_sub_cut.begin(); last!=min_sub_cut.end(); ++last){
    cerr<<" "<<(*last)->GetId()<<"("<<(*last)->Mpeak<<"|"<<available_memory_after_subroot + (*last)->GetEW()<<")";
  }
  cerr<<"]"<<endl;

  cerr<<spacing<<"[ node "<<node->GetId()<<" ] final subcut value is "<<cut_value<<" subroot->mpeak is "<<Mpeak<<endl;
#endif

  delete candidates;
  return;
}

#ifdef DEBUG_USING_MINMEM
void MinMem(Tree * tree,double MaxOutDeg, double & Required_memory, schedule_t & Schedule, int quiet, uint64_t & count,iter_node_t * minmem_trace)
#else
void MinMem(Tree * tree,double MaxOutDeg, double & Required_memory, schedule_t & Schedule, int quiet, uint64_t & count)
#endif
{
  double M = numeric_limits<double>::infinity();
  double Mpeak = MaxOutDeg;

  struct rlimit lim;
  lim.rlim_cur = RLIM_INFINITY;
  lim.rlim_max = RLIM_INFINITY;

  setrlimit(RLIMIT_STACK,&lim);

  Node * root=tree->GetRoot();

  if(!quiet){
    cerr<<"Max out deg = "<<Mpeak<<endl;
  }
  count = 0;
  list<Node*> L;
  //	list<Node*> * L = new list<Node*>();
  //	list<Node*> * prevL = new list<Node*>();
  Schedule.clear();
  while(Mpeak<numeric_limits<double>::infinity()){
    Required_memory = Mpeak;
#ifdef DEBUG_USING_MINMEM
    explore(root, Required_memory, &L, &Schedule,  M, L, Schedule, Mpeak, quiet,0,count,minmem_trace,tree->GetNodes()->size());
#else
    explore(root, Required_memory, &L, &Schedule,  M, L, Schedule, Mpeak, quiet,0,count);
#endif

    //		list<Node*> * tmp = prevL;
    //		prevL = L;
    //		L = tmp;
#if VERBOSE
    cerr<<"[MinMemO] Available memory = "<<Required_memory<<"  Mpeak returned by explore = "<<Mpeak<<"  cutvalue = "<<M<<" cutsize = "<<L.size()<<endl;
#endif
  }

  //	tree->Print(cerr);	
  //	delete L;
  //	delete prevL;
}

#ifdef DEBUG_USING_MINMEM
double MinMemRecurAlgorithm(int N, int *prnts, double *nwghts, double *ewghts,double *mswghts, int *schedule,iter_node_t * minmem_trace)
#else
  double MinMemRecurAlgorithm(int N, int *prnts, double *nwghts, double *ewghts,double *mswghts, int *schedule)
#endif
    {
      Tree * tree = new Tree(N,prnts,nwghts,ewghts,mswghts);

      double Mr;
      double Mp = MaxOutDegree(tree,true);
      uint64_t count = 0;
      schedule_t * sub_sched = new schedule_t();

#ifdef DEBUG_USING_MINMEM
      MinMem(tree,Mp, Mr, *sub_sched, true, count,minmem_trace);
#else
      MinMem(tree,Mp, Mr, *sub_sched, true, count);
#endif

      unsigned int i = 0;
      for (schedule_t::reverse_iterator last=sub_sched->rbegin(); last!=sub_sched->rend(); ++last){
        schedule[i++] = (int)(*last);
      }

      delete sub_sched;
      delete tree;

      return Mr;
    }

//#ifdef DEBUG_USING_MINMEM
//    double MinMemRecurAlgorithm_timed(int N, int *prnts, double *nwghts, double *ewghts, int *schedule,double * usec,int quiet,iter_node_t * minmem_trace)
//#else
//      double MinMemRecurAlgorithm_timed(int N, int *prnts, double *nwghts, double *ewghts, int *schedule,double * usec,int quiet)
//#endif
//      {
//        Tree * tree = new Tree(N,prnts,nwghts,ewghts);
//
//        double Mr;
//        double Mp = MaxOutDegree(tree,quiet);
//        uint64_t count = 0;
//        schedule_t * sub_sched = new schedule_t();
//
//        *usec = -u_wseconds();
//#ifdef DEBUG_USING_MINMEM
//        MinMem(tree,Mp, Mr, *sub_sched, quiet, count,minmem_trace);
//#else
//        MinMem(tree,Mp, Mr, *sub_sched, quiet, count);
//#endif
//        *usec += u_wseconds();
//
//
//
//        unsigned int i = 0;
//        for (schedule_t::reverse_iterator last=sub_sched->rbegin(); last!=sub_sched->rend(); ++last){
//          schedule[i++] = (int)(*last);
//        }
//
//        delete sub_sched;
//        delete tree;
//
//        return Mr;
//      }

    void exploreArray(int N,double * nwghts, double * ewghts, int * chstart, int *children,int nd, double available_memory,list<node_peak_t> * L_init,int * S_init,  double & cut_value, list<node_peak_t> & min_sub_cut, int * sub_schedule,int & sched_head, double & Mpeak, int quiet, int depth,uint64_t & count){
      count++;

      char * spacing; if(!quiet){spacing = (char *)malloc((depth*2+1)*sizeof(char)); for(int i=0;i<depth*2;++i){spacing[i]='\t';} spacing[depth*2]='\0'; }

      double cost = ewghts[nd]+nwghts[nd];
      for(int j = chstart[nd]; j<chstart[nd+1];j++){ cost += ewghts[children[j]]; }


      /* if node is unreachable, return +infty */
      if (cost > available_memory + MIN_ACCURACY) {
        if(!quiet){cerr<<spacing<<"not enough mem for "<<nd<<endl;}
        Mpeak = cost;
        cut_value = numeric_limits<double>::infinity( );
        return;
      }

      /* if this is a leaf, return 0 */
      if (chstart[nd]==chstart[nd+1]){
        if(!quiet){fprintf(stderr,"%s%d is leaf\n",spacing,nd);}
        sub_schedule[sched_head--] = nd;
        cut_value = 0;
        Mpeak = numeric_limits<double>::infinity( );
        return;
      }

      cut_value = 0;

      list<node_peak_t>::iterator begin;

      int candidates_count = 0;


      if(L_init != NULL){
        if(L_init->size()>0){
          min_sub_cut.assign(L_init->begin(),L_init->end());
          candidates_count = min_sub_cut.size();
          begin = min_sub_cut.begin();
          for (list<node_peak_t>::iterator node=min_sub_cut.begin(); node!=min_sub_cut.end(); ++node){cut_value += ewghts[node->index];}
          memcpy(&sub_schedule[sched_head],&S_init[sched_head],(N +1 - sched_head)*sizeof(int));
        }
        else{
          sub_schedule[sched_head--] = nd;
          /* place every child in the candidate nodes*/
          candidates_count = chstart[nd+1]-chstart[nd];
          for(int j = chstart[nd]; j<chstart[nd+1];j++){ min_sub_cut.push_back(node_peak_t(children[j])); cut_value += ewghts[children[j]];}
          begin = min_sub_cut.begin();
        }
      }
      else{
        sub_schedule[sched_head--] = nd;
        /* place every child in the candidate nodes*/
        candidates_count = chstart[nd+1]-chstart[nd];
        begin = min_sub_cut.end();
        begin--;
        for(int j = chstart[nd]; j<chstart[nd+1];j++){ 
          min_sub_cut.push_back(node_peak_t(children[j])); 
          cut_value += ewghts[children[j]];
        }
        begin++;
        //		begin = min_sub_cut.begin();		
      }

      int sub_cut_size = distance(begin,min_sub_cut.end());

      while (candidates_count>0) {
        Mpeak = numeric_limits<double>::infinity( );

        list<node_peak_t>::iterator current_node=begin;
        int local_count = 0;
        while(local_count++<sub_cut_size){
          double m_j = numeric_limits<double>::infinity();
          //			list<node_peak_t> Lj;
          //schedule_t Sj;
          int tmp_sched_head = sched_head;
          double mp = current_node->mpeak;
          if(available_memory - cut_value + ewghts[current_node->index] >= current_node->mpeak ){

            double l_mavail = available_memory - cut_value + ewghts[current_node->index];

            if(!quiet){fprintf(stderr,"%ssubroot %d subtree avail mem %lf\n",spacing,current_node->index, l_mavail);}
            //				exploreArray(nwghts,ewghts,chstart,children,current_node->index, l_mavail,NULL,NULL, m_j, Lj,Sj,current_node->mpeak,quiet,depth+1,count);
            //				exploreArray(N,nwghts,ewghts,chstart,children,current_node->index, l_mavail,NULL,NULL, m_j, Lj,sub_schedule,tmp_sched_head,current_node->mpeak,quiet,depth+1,count);				

            list<node_peak_t>::iterator current_subcut_end = min_sub_cut.end();
            current_subcut_end--;
            exploreArray(N,nwghts,ewghts,chstart,children,current_node->index, l_mavail,NULL,NULL, m_j, min_sub_cut,sub_schedule,tmp_sched_head,current_node->mpeak,quiet,depth+1,count);				
            current_subcut_end++;
            mp = current_node->mpeak;

            if(!quiet){
              fprintf(stderr,"%sm_j %lf f_j %lf \n",spacing,m_j,ewghts[current_node->index]);
              cerr<<spacing<<"Sj : {"; for (int j = tmp_sched_head+1;j<=sched_head;j++){ cerr<<sub_schedule[j]<<" "; } cerr<<"}"<<endl;
            }


            if(!quiet){cerr<<spacing<<"m_j = "<<m_j<<" ew "<< ewghts[current_node->index]<<endl;}


            if (m_j > ewghts[current_node->index] ){
              //pop the end of the list
              min_sub_cut.erase(current_subcut_end,min_sub_cut.end());
            }
          }

          if (m_j <= ewghts[current_node->index]){
            //cut_value = cut_value - ewghts[current_node->index] + m_j;
            list<node_peak_t>::iterator to_erase = current_node;
            current_node--;
            local_count--;

            begin--;
            min_sub_cut.erase(to_erase);
            begin++;


            sub_cut_size= distance(begin,min_sub_cut.end());

            sched_head = tmp_sched_head;

          }

          current_node++;
        }


        cut_value = 0;

        candidates_count=0;
        //		if(min_sub_cut.size()>0){
        list<node_peak_t>::iterator node=begin;
        for (local_count=0; local_count<sub_cut_size;++local_count){
          cut_value += ewghts[node->index];
          ++node;
        }

        node=begin;
        for (local_count=0; local_count<sub_cut_size;++local_count){
          if(available_memory - cut_value + ewghts[node->index] >= node->mpeak ){
            candidates_count++;
          }

          //				cerr<<"Mpeak "<<Mpeak<<" vs "<<new_node->mpeak + cut_value - ewghts[new_node->index]<<endl;
          Mpeak = min(Mpeak,node->mpeak + cut_value - ewghts[node->index]);
          ++node;
        }
        //		}
        //		Mpeak += cut_value;
        //		cerr<<endl<<Mpeak<<endl;
        if(!quiet){ cerr<<spacing<<candidates_count<<" candidates left"<<endl; }		
      }

      //	if(!quiet){
      //		//		if(min_sub_cut.size()>0){
      //		int local_count;
      //		cerr<<spacing<<"MSC : {";
      //		list<node_peak_t>::iterator last=begin;
      //		for (local_count=0; local_count<sub_cut_size;++local_count){
      //			cerr<<last->index<<" ";
      //			++last;
      //		}
      //		cerr<<"}"<<endl;
      //		//		}
      //		
      //		cerr<<spacing<<"SUB SCHED : {";
      //		for (int j = sched_head+1;j<=prev_sched_head;j++){ cerr<<sub_schedule[j]<<" "; }
      //		//		for (list<int>::iterator last=sub_schedule.begin(); last!=sub_schedule.end(); ++last){
      //		//			cerr<<*last<<" ";
      //		//		}
      //		cerr<<"}"<<endl;
      //		
      //		
      //		cerr<<spacing<<"TOTAL SCHED : {";
      //		for (int j = sched_head+1;j<N;j++){ cerr<<sub_schedule[j]<<" "; }
      //		//		for (list<int>::iterator last=sub_schedule.begin(); last!=sub_schedule.end(); ++last){
      //		//			cerr<<*last<<" ";
      //		//		}
      //		cerr<<"}"<<endl;
      //		
      //		cerr<<spacing<<"CUTVALUE :"<<cut_value<<endl;
      //		cerr<<spacing<<"MPEAK :"<<Mpeak<<endl;
      //	}

      return;
    }

#if DEBUG_MEMUSAGE
    void exploreArray2(int N,double * nwghts, double * ewghts, int * chstart, int *children,iter_node_t * subroot,double available_memory,bool Linit,iter_node_t * NodeArray, double & cut_value, MinMemDLL * L,int * sub_schedule,int & sched_head, int depth,uint64_t & count, double * memsched_debug)
#else
#ifdef DEBUG_USING_MINMEM
      void exploreArray2(int N,double * nwghts, double * ewghts, int * chstart, int *children,iter_node_t * subroot,double available_memory,bool Linit,iter_node_t * NodeArray, double & cut_value, MinMemDLL * L,int * sub_schedule,int & sched_head, int depth,uint64_t & count, iter_node_t * minmem_trace, iter_node_t * minmema_trace)
#else
      void exploreArray2(int N,double * nwghts, double * ewghts, int * chstart, int *children,iter_node_t * subroot,double available_memory,bool Linit,iter_node_t * NodeArray, double & cut_value, MinMemDLL * L,int * sub_schedule,int & sched_head, int depth,uint64_t & count)
#endif
#endif
      {
#if DEBUG_MEMUSAGE
        if(!Linit || depth>0){
          assert(available_memory==memsched_debug[sched_head+1]);
          //  cerr<<"[ node "<<subroot->index<<" ] assert passed"<<endl;
        }
#endif
#ifdef DEBUG_USING_MINMEM
        if(count<2*N){
          minmema_trace[count].index = subroot->index;
          minmema_trace[count].mavail = available_memory + ewghts[subroot->index];
          subroot->mavail = available_memory + ewghts[subroot->index];
          assert(subroot->index ==  minmem_trace[count].index);
          assert(available_memory + ewghts[subroot->index] ==  minmem_trace[count].mavail);
        }
        uint64_t mycount = count;
#endif
        count++;
#if STRONG_ASSERT
        assert(count<(uint64_t)(2*N*N));	
#endif

#if VERBOSE
        char * spacing; 
        spacing = (char *)malloc((depth*2+1)*sizeof(char)); 
        for(int i=0;i<depth*2;++i){spacing[i]=' ';} spacing[depth*2]='\0'; 
#endif
#if VERBOSE
        cerr<<spacing<<"[ node "<<subroot->index<<" ] is being explored with memory " <<available_memory<<" while its cost is "<<subroot->cost<<endl;
#endif


        /* if node is unreachable, return +infty */
        if (subroot->cost - ewghts[subroot->index]> available_memory) {
#if VERBOSE
          cerr<<spacing<<"[ node "<<subroot->index<<" ] not enough memory"<<endl;
#endif
          subroot->mpeak = subroot->cost;
          cut_value = numeric_limits<double>::infinity( );
#ifdef DEBUG_USING_MINMEM
          subroot->subcut_value = cut_value;
          if(mycount<2*N){
            minmema_trace[mycount].subcut_value = cut_value;
            minmema_trace[mycount].mpeak = subroot->mpeak;

            assert(subroot->mpeak ==  minmem_trace[mycount].mpeak);
            assert(subroot->subcut_value ==  minmem_trace[mycount].subcut_value);
          }
#endif
          return;
        }

        //#define DEBUG_USING_LIU

        /* if this is a leaf, return 0 */
        if (chstart[subroot->index]==chstart[subroot->index+1]){
#if VERBOSE
          cerr<<spacing<<"[ node "<<subroot->index<<" ] is a leaf"<<endl;
#endif
#if DEBUG_MEMUSAGE
          memsched_debug[sched_head] = memsched_debug[sched_head+1 ] + ewghts[subroot->index];
#endif
          sub_schedule[sched_head--] = subroot->index;
          cut_value = 0;
          subroot->mpeak = numeric_limits<double>::infinity( );
#ifdef DEBUG_USING_MINMEM
          subroot->subcut_value = cut_value;
          if(mycount<2*N){
            minmema_trace[mycount].subcut_value = cut_value;
            minmema_trace[mycount].mpeak = subroot->mpeak;
            assert(subroot->mpeak ==  minmem_trace[mycount].mpeak);
            assert(subroot->subcut_value ==  minmem_trace[mycount].subcut_value);
          }
#endif
          return;
        }

        cut_value = 0;
        int candidates_count = 0;
        if(Linit && L->size() >0){
          candidates_count = L->size();
#ifndef USE_DICHOTOMY
          iter_node_t * cur_node = L->begin();
          while (cur_node != L->end()){
            cut_value += ewghts[cur_node->index];
            cur_node = cur_node->pNext;
          }
          subroot->pSCend = L->last();
          subroot->pSCbegin = L->begin();
#else
          s_list_item_t * cur_node = L->begin();
          while (cur_node != L->end()){
            cut_value += ewghts[cur_node->node->index];
            cur_node = cur_node->pNext;
          }
          subroot->pSCend = L->last()->node;
          subroot->pSCbegin = L->begin()->node;
#endif
        }
        else{
#if DEBUG_MEMUSAGE
          memsched_debug[sched_head] = memsched_debug[sched_head+1 ] + 2*ewghts[subroot->index] + nwghts[subroot->index] - subroot->cost;
#endif
          sub_schedule[sched_head--] = subroot->index;
          /* place every child in the candidate nodes*/
          candidates_count = chstart[subroot->index+1]-chstart[subroot->index];
          iter_node_t * cur_node;

#ifndef DEBUG_LIU    
          iter_node_t * prev_node = subroot;
          for(int j = chstart[subroot->index]; j<chstart[subroot->index+1];j++){
            cur_node = &NodeArray[children[j]];
            L->insert_after(prev_node,cur_node);
            cut_value += ewghts[cur_node->index];
            prev_node = cur_node;
          }
          subroot->pSCbegin = &NodeArray[children[chstart[subroot->index]]];
          subroot->pSCend = &NodeArray[children[chstart[subroot->index+1]-1]];
#else
          int * seen = new int [chstart[subroot->index+1] - chstart[subroot->index]];

          for(int k=0; k<chstart[subroot->index+1] - chstart[subroot->index]; k++) {seen[k] = 0;}
#ifdef DEBUG_LIU_RANK    
          if (chstart[subroot->index] != chstart[subroot->index+1]-1) {
            cerr<<"children of node "<<subroot->index<< "(liu rank)";
            for(int j = chstart[subroot->index]; j<chstart[subroot->index+1];j++) {
              cerr<<" "<<NodeArray[children[j]].index<<"("<<NodeArray[children[j]].liu_rank<<")";
            }
            cerr<<endl;
          }
#endif

          for(int j = chstart[subroot->index]; j<chstart[subroot->index+1];j++) {

            int max_liu_rank_index = -1;
            for(int k = chstart[subroot->index]; k<chstart[subroot->index+1];k++){
              if (seen[k - chstart[subroot->index]])
                continue;
              if ((max_liu_rank_index == -1) || (NodeArray[children[k]].liu_rank > NodeArray[children[max_liu_rank_index]].liu_rank)){
                max_liu_rank_index = k;
              }
            }
            assert(max_liu_rank_index != -1);

            seen[max_liu_rank_index - chstart[subroot->index]] = 1;
            cur_node = &NodeArray[children[max_liu_rank_index]];

#ifdef DEBUG_LIU_RANK      
            if (chstart[subroot->index] != chstart[subroot->index+1]-1) { cerr<<" "<<cur_node->index<<"("<<cur_node->liu_rank<<")"; }
#endif      
            L->insert_after(subroot,cur_node);
            cut_value += ewghts[cur_node->index];

            //first node explored
            if(j==chstart[subroot->index])
              subroot->pSCend = cur_node;
            //last node explored
            if(j==chstart[subroot->index+1]-1)
              subroot->pSCbegin = cur_node;
          }
          delete[] seen;
#ifdef DEBUG_LIU_RANK
          if (chstart[subroot->index] != chstart[subroot->index+1]-1) { cerr<<endl; }
#endif
#endif

        }

//TODO CONTINUE HERE
        int cut_size = candidates_count;
#if defined(DEBUG_LIU_RANK) or VERBOSE
        if (chstart[subroot->index] != chstart[subroot->index+1]-1) {
          spacing[0]='\0';
          cerr<<spacing<<"[ node "<<subroot->index<<" ] initial subcut is [";
#ifndef USE_DICHOTOMY
          iter_node_t * cur_node = subroot->pSCbegin;
          while(cut_size>0 && cur_node != subroot->pSCend->pNext){
            cerr<<" "<<cur_node->index<<"("<<cur_node->mpeak<<")";
            cur_node = cur_node->pNext;
          }
          cerr<<" ]"<<endl;
#else
          s_list_item_t * cur_node = L->getItem(subroot->pSCbegin);
	  s_list_item_t * end_node = L->getItem(subroot->pSCend);
          while(cut_size>0 && cur_node != end_node->pNext){
            cerr<<" "<<cur_node->node->index<<"("<<cur_node->node->mpeak<<")";
            cur_node = cur_node->pNext;
          }
          cerr<<" ]"<<endl;
#endif
        }
#endif


        bool first_loop = true;
        while (candidates_count>0) {

#ifndef USE_DICHOTOMY
          iter_node_t * current_node = subroot->pSCbegin;
#if VERBOSE
          cerr<<spacing<<"*****************************************************************"<<endl;
          cerr<<spacing<<"[ node "<<subroot->index<<" ] candidates are [";
          while(cut_size>0 && current_node != subroot->pSCend->pNext){
            double available_memory_after_subroot = available_memory - cut_value + ewghts[subroot->index];
            if (available_memory_after_subroot + ewghts[current_node->index]>= current_node->mpeak)
            {
              cerr<<" "<<current_node->index<<"("<<current_node->mpeak<<"|"<<available_memory_after_subroot + ewghts[current_node->index] <<")";
            }
            else{
              cerr<<" BAD{"<<current_node->index<<"("<<current_node->mpeak<<"|"<<available_memory_after_subroot + ewghts[current_node->index]<<"])";

            }
            current_node = current_node->pNext;
          }
          cerr<<" ]"<<endl;
          current_node = subroot->pSCbegin;
#endif
#else
          s_list_item_t * current_item = L->getItem(subroot->pSCbegin);
	  s_list_item_t * end_item = L->getItem(subroot->pSCend);
#endif

#ifndef USE_DICHOTOMY
          while(cut_size>0 && current_node != subroot->pSCend->pNext){
#else
          while(cut_size>0 && current_item != end_item->pNext){
		  iter_node_t * current_node = current_item->node;
#endif
            double m_j = numeric_limits<double>::infinity();
            int tmp_sched_head = sched_head;
            uint64_t prevCutSize = L->size();
            uint64_t subcut_size = 0;
            double available_memory_after_subroot = available_memory - cut_value + ewghts[subroot->index];
#if STRONG_ASSERT
            iter_node_t * next_node = current_node->pNext;
#endif

            if (available_memory_after_subroot + ewghts[current_node->index]>= current_node->mpeak || first_loop)
              {
#if DEBUG_MEMUSAGE
                exploreArray2(N,nwghts,ewghts,chstart,children,current_node,available_memory_after_subroot,false,NodeArray, m_j,L,sub_schedule,tmp_sched_head,depth+1,count,memsched_debug);
#else

#ifdef DEBUG_USING_MINMEM
                exploreArray2(N,nwghts,ewghts,chstart,children,current_node,available_memory_after_subroot,false,NodeArray, m_j,L,sub_schedule,tmp_sched_head,depth+1,count,minmem_trace,minmema_trace);
#else
                exploreArray2(N,nwghts,ewghts,chstart,children,current_node,available_memory_after_subroot,false,NodeArray, m_j,L,sub_schedule,tmp_sched_head,depth+1,count);
#endif

#endif
                subcut_size = L->size() - prevCutSize;

                if (m_j > ewghts[current_node->index] ){
                  //pop the end of the list
#if VERBOSE
                  cerr<<spacing<<"  [ node "<<current_node->index<<" ] is better than its subtree ("<<ewghts[current_node->index]<<" vs "<<m_j<<")"<<endl;
#endif
                L->splice(current_node,current_node->pSCend, subcut_size);
#ifndef NOASSERT
                if(subcut_size>0){
                  assert(current_node->pNext == current_node->pSCend->pNext);
                }
#endif

                //          for(int i = tmp_sched_head+1;i<=sched_head;i++){sub_schedule[i]=-100000;}
              }
              else if (m_j <= ewghts[current_node->index] ){
                //erase node
                if( subcut_size>0){
                  // if there are nodes in the subcut
                  if(current_node==subroot->pSCbegin){
                    subroot->pSCbegin = current_node->pSCbegin;
                    //              cerr<<"[ node "<<subroot->index<<" ] subcut new begin "<<subroot->pSCbegin->index<<endl;
                  }
                  if(current_node==subroot->pSCend){
                    subroot->pSCend = current_node->pSCend;
#ifdef USE_DICHOTOMY
	  	    end_item = L->getItem(subroot->pSCend);
#endif
                    //              cerr<<"[ node "<<subroot->index<<" ] subcut new end "<<subroot->pSCend->index<<endl;
                  }
                }
                else{
                  // if the subcut of the children is empty but subroot's cut has other nodes than the one we are going to remove
                  if(cut_size>1){
                    if(current_node==subroot->pSCbegin && !(current_node==subroot->pSCend)){
#ifndef USE_DICHOTOMY
                      subroot->pSCbegin = current_node->pNext;
#else
                      subroot->pSCbegin = current_item->pNext->node;
#endif
                      //              cerr<<"[ node "<<subroot->index<<" ] subcut new begin "<<subroot->pSCbegin->index<<endl;
                    }
                    else if(current_node==subroot->pSCend && !(current_node==subroot->pSCbegin)){
#ifndef USE_DICHOTOMY
                      subroot->pSCend = current_node->pPrev;
#else
                      subroot->pSCend = current_item->pPrev->node;
	  	      end_item = L->getItem(subroot->pSCend);
#endif
                      //              cerr<<"[ node "<<subroot->index<<" ] subcut new end "<<subroot->pSCend->index<<endl;
                    }
                  }
                }

                L->erase(current_node);
#if STRONG_ASSERT
                assert(current_node->pPrev==current_node->pNext->pPrev);
#endif

#ifndef DEBUG_USING_MINMEM
                cut_value = cut_value - ewghts[current_node->index] + m_j;
#endif

                cut_size += subcut_size -1;
                sched_head = tmp_sched_head;

#if VERBOSE
                cerr<<spacing<<cut_value<<"  [ node "<<current_node->index<<" ] has been removed ("<<ewghts[current_node->index]<<" vs "<<m_j<<"), new cut size is "<<cut_size<<endl;
#endif

                //advance the pointer to skip the new nodes in the cut
                if(subcut_size>0){
                  current_node = current_node->pSCend;
#ifdef USE_DICHOTOMY
		  current_item = L->getItem(current_node); 
#endif
                }
#ifndef USE_DICHOTOMY
#if STRONG_ASSERT
                assert(current_node->pNext == next_node);
#endif
#endif
              }

          }

#ifndef USE_DICHOTOMY
          current_node = current_node->pNext;
#else
          current_item = current_item->pNext;
#endif
        }

        first_loop=false;

#ifndef USE_DICHOTOMY
        cut_value = 0;
        current_node = subroot->pSCbegin;
        while(cut_size>0 && subroot->pSCend->pNext != current_node){
          cut_value += ewghts[current_node->index];
          current_node = current_node->pNext;
        }

        double available_memory_after_subroot = available_memory - cut_value + ewghts[subroot->index];
        candidates_count=0;
        current_node = subroot->pSCbegin;
        while(cut_size>0 && subroot->pSCend->pNext != current_node){
          if (available_memory_after_subroot  + ewghts[current_node->index]>= current_node->mpeak){ 
            candidates_count++;
          }
          current_node = current_node->pNext;
        }

#if VERBOSE
        cerr<<spacing<<"[ node "<<subroot->index<<" ] new subcut value is "<<cut_value<<endl;
        cerr<<spacing<<"[ node "<<subroot->index<<" ] new subcut is [";
        current_node = subroot->pSCbegin;
        while(cut_size>0 && subroot->pSCend->pNext != current_node){
          cerr<<" "<<current_node->index<<"("<<current_node->mpeak<<"|"<<available_memory_after_subroot  + ewghts[current_node->index]<<")";
          current_node = current_node->pNext;
        }
        cerr<<"]"<<endl;
        cerr<<spacing<<"*****************************************************************"<<endl;
        cerr<<spacing<<candidates_count<<" candidates left"<<endl;
#endif
#else
        current_item = L->getItem(subroot->pSCbegin);
	end_item = L->getItem(subroot->pSCend);
        cut_value = 0;
        while(cut_size>0 && current_item != end_item->pNext){
          cut_value += ewghts[current_item->node->index];
          current_item = current_item->pNext;
        }

        double available_memory_after_subroot = available_memory - cut_value + ewghts[subroot->index];
        current_item = L->getItem(subroot->pSCbegin);
        candidates_count=0;
        while(cut_size>0 && current_item != end_item->pNext){
          if (available_memory_after_subroot  + ewghts[current_item->node->index]>= current_item->node->mpeak){ 
            candidates_count++;
          }
          current_item = current_item->pNext;
        }
#endif
      }


#ifndef USE_DICHOTOMY
    cut_value = 0;
    iter_node_t * current_node = subroot->pSCbegin;
    while(cut_size>0 && subroot->pSCend->pNext != current_node){
      cut_value += ewghts[current_node->index];
      current_node = current_node->pNext;
    }
#ifdef DEBUG_USING_MINMEM
    subroot->subcut_value = cut_value;
#endif

#if VERBOSE
    double available_memory_after_subroot = available_memory - cut_value + ewghts[subroot->index];
    current_node = subroot->pSCbegin;
    cerr<<spacing<<"[ node "<<subroot->index<<" ] final subcut is [";
    while(cut_size>0 && subroot->pSCend->pNext != current_node){
      cerr<<" "<<current_node->index<<"("<<current_node->mpeak<<"|"<<available_memory_after_subroot + ewghts[current_node->index]<<")";
      current_node = current_node->pNext;
    }
    cerr<<"]"<<endl;
#endif


    current_node = subroot->pSCbegin;
    subroot->mpeak = numeric_limits<double>::infinity();
    while(cut_size>0 && subroot->pSCend->pNext != current_node){
#if VERBOSE
      cerr<<spacing<<"[ node "<<subroot->index<<" ] current Mpeak is "<<subroot->mpeak<<" while [ node "<<current_node->index<<" ] Mpeak is "<<current_node->mpeak<<" which gives "<<current_node->mpeak + cut_value - ewghts[current_node->index]<<" in the current cut"<<endl;
#endif
      subroot->mpeak = min(subroot->mpeak,current_node->mpeak + cut_value - ewghts[current_node->index]);




      current_node = current_node->pNext;
    }
#ifdef DEBUG_USING_MINMEM
    if(mycount<2*N){
      minmema_trace[mycount].subcut_value = cut_value;
      minmema_trace[mycount].mpeak = subroot->mpeak;
      assert(subroot->mpeak ==  minmem_trace[mycount].mpeak);
      assert(subroot->subcut_value ==  minmem_trace[mycount].subcut_value);
    }
#endif

#if VERBOSE
    cerr<<spacing<<"[ node "<<subroot->index<<" ] final subcut value is "<<cut_value<<" subroot->mpeak is "<<subroot->mpeak<<endl;
    //  cerr<<spacing<<"[ node "<<subroot->index<<" ] whole schedule : {"; for (int j = sched_head + 1; j<N;j++){ cerr<<sub_schedule[j]<<" "; } cerr<<"}"<<endl;
#endif
#else
    s_list_item_t * current_item = L->getItem(subroot->pSCbegin);
    s_list_item_t * end_item = L->getItem(subroot->pSCend);
    cut_value = 0;
    while(cut_size>0 && current_item != end_item->pNext){
      cut_value += ewghts[current_item->node->index];
      current_item = current_item->pNext;
    }

    double available_memory_after_subroot = available_memory - cut_value + ewghts[subroot->index];
    current_item = L->getItem(subroot->pSCbegin);
    subroot->mpeak = numeric_limits<double>::infinity();
    while(cut_size>0 && current_item != end_item->pNext){
      subroot->mpeak = min(subroot->mpeak,current_item->node->mpeak + cut_value - ewghts[current_item->node->index]);
      current_item = current_item->pNext;
    }
#endif
    return;
  }




void MinMemArray(int N, int root, double * nwghts, double * ewghts, int * chstart, int * children,double MaxOutDeg, double & Required_memory, int * Schedule, int quiet, uint64_t & count,int * prnts){
  double Cut_value = 0;

  struct rlimit lim;
  lim.rlim_cur = RLIM_INFINITY;
  lim.rlim_max = RLIM_INFINITY;

  setrlimit(RLIMIT_STACK,&lim);


#ifdef DEBUG_LIU
  double time_liu = 0;
  //double mem_liu = PebbleOrderingIterAlgorithm_timed( N, prnts, nwghts, ewghts, Schedule, &time_liu,1);
  double mem_liu = PostOrderIterAlgorithm_timed( N, prnts, nwghts, ewghts, Schedule, &time_liu,1);
  //schedule now contains liu's schedule
#endif



  iter_node_t * NodeArray = new iter_node_t[N+1];
  NodeArray[root].mpeak = MaxOutDeg;

#if DEBUG_MEMUSAGE
  double * memsched_debug = new double [N+1];
#endif
#ifdef USE_DICHOTOMY
  double max_memory = 0;
#endif
  for(int i=1;i<N+1;i++){
    NodeArray[i].index = i;
    NodeArray[i].cost = ewghts[i]+nwghts[i];
    for(int j = chstart[i]; j<chstart[i+1];j++){ NodeArray[i].cost += ewghts[children[j]]; }
#ifdef USE_DICHOTOMY
    max_memory += ewghts[i]+nwghts[i];
#endif
  }

#ifdef DEBUG_LIU
  for(int i=0;i<N;i++){
    int scheduled_node_index = Schedule[i];
    NodeArray[scheduled_node_index].liu_rank = N-1-i;
  }
#endif


#ifdef DEBUG_USING_MINMEM
  iter_node_t * minmema_trace = new iter_node_t[(uint64_t)((uint64_t)2*(uint64_t)N)];
  iter_node_t * minmem_trace = new iter_node_t[(uint64_t)((uint64_t)2*(uint64_t)N)];
  double stub = MinMemRecurAlgorithm(N, prnts, nwghts, ewghts,Schedule, minmem_trace);
#endif

#ifndef USE_DICHOTOMY
  MinMemDLL * MinCut = new MinMemDLL(&NodeArray[root]);
#else
  MinMemDLL * MinCut = new MinMemDLL(&NodeArray[root],N+1);
  MinMemDLL * prev_MinCut = new MinMemDLL(&NodeArray[root],N+1);
#endif

  int sched_head = N-1;
  bool initialize_cut = false;
  count = 0;
#ifndef USE_DICHOTOMY
  while(NodeArray[root].mpeak<numeric_limits<double>::infinity()){
    Required_memory = NodeArray[root].mpeak;
#if DEBUG_MEMUSAGE
    memsched_debug[sched_head+1] = Required_memory - Cut_value;
#endif
#if DEBUG_MEMUSAGE
    exploreArray2(N,nwghts,ewghts,chstart,children,&NodeArray[root],Required_memory, initialize_cut,NodeArray, Cut_value, MinCut, Schedule, sched_head,0,count,memsched_debug);
#else

#ifdef DEBUG_USING_MINMEM
    exploreArray2(N,nwghts,ewghts,chstart,children,&NodeArray[root],Required_memory, initialize_cut,NodeArray, Cut_value, MinCut, Schedule, sched_head,0,count,minmem_trace,minmema_trace);
#else
    exploreArray2(N,nwghts,ewghts,chstart,children,&NodeArray[root],Required_memory, initialize_cut,NodeArray, Cut_value, MinCut, Schedule, sched_head,0,count);
#endif

#endif
    initialize_cut = true;
#if STRONG_ASSERT
    assert(Required_memory<=NodeArray[root].mpeak);
#endif
#if VERBOSE
    cerr<<"[MinMem] Available memory = "<<Required_memory<<"  Mpeak returned by explore = "<<NodeArray[root].mpeak<<"  cutvalue = "<<Cut_value<<" cutsize = "<<MinCut->size()<<endl;
#endif
  }
#else
  uint64_t prev_sched_head = N-1;
  double min_memory = MaxOutDeg;
  max_memory = 2*MaxOutDeg; /* tmp need to check if it is too low*/
  /*start with max out deg : good heuristic */
  double half_mem = min_memory;//min_memory + floor((max_memory-min_memory)/2);
  double prev_peak = -1;
  double prev_valid_peak = -1;
  MinMemDLL * current_cut = MinCut;
  MinMemDLL * prev_cut = prev_MinCut;
  while(max_memory != min_memory){

#if VERBOSE
    cerr<<"Input cut size = "<<current_cut->size()<<endl;
#endif 
    exploreArray2(N,nwghts,ewghts,chstart,children,&NodeArray[root],half_mem, initialize_cut,NodeArray, Cut_value, current_cut, Schedule, sched_head,0,count);


//  s_list_item_t * cur_item = current_cut->begin(); 
//  cerr<<endl;
//  while(cur_item!=current_cut->end()){
//    cerr<<" "<<cur_item->node->index;
//    cur_item = cur_item->pNext;
//  }
//  cerr<<endl;





#if VERBOSE
    cerr<<"[MinMem] Available memory = "<<half_mem<<"  Mpeak returned by explore = "<<NodeArray[root].mpeak<<"  cutvalue = "<<Cut_value<<" cutsize = "<<MinCut->size()<<" prev_peak = "<<prev_peak<<" prev_valid_peak = "<<prev_valid_peak<<" Min memory = "<<min_memory<<" Max memory = "<<max_memory<<endl;
    cerr<<"Output cut size = "<<current_cut->size()<<endl;
#endif

    if(NodeArray[root].mpeak==numeric_limits<double>::infinity() && prev_peak == -1) {
      break;
    }

    if(prev_valid_peak == half_mem /*&& NodeArray[root].mpeak==numeric_limits<double>::infinity()*/){
      break;
    }

    if(NodeArray[root].mpeak==numeric_limits<double>::infinity()){
      //      sched_head = N-1;
      //      initialize_cut = false;
      initialize_cut = true;
      sched_head = prev_sched_head;

      current_cut->copy(prev_cut);
      //prev_cut->copy(current_cut);
      //          MinMemDLL * tmp = current_cut;
      //        current_cut = prev_cut;
      //      prev_cut = tmp;
      max_memory = half_mem;
      half_mem = min_memory + ceil((max_memory-min_memory)/2);
      //      if(half_mem< prev_valid_peak){
      //        half_mem = prev_valid_peak;
      //      }

    }
    else{
      prev_cut->copy(current_cut);
      prev_sched_head = sched_head;
      //      min_memory = NodeArray[root].mpeak;
      min_memory = half_mem;
      half_mem = min_memory + ceil((max_memory-min_memory)/2);
//      if(half_mem< NodeArray[root].mpeak){
//        half_mem = NodeArray[root].mpeak;
//      }

      initialize_cut = true;

      prev_valid_peak = NodeArray[root].mpeak;
    }

    prev_peak = NodeArray[root].mpeak;

  }

//  delete tmpCut;
  Required_memory = half_mem;

#endif


#if DEBUG_MEMUSAGE
  delete[] memsched_debug;
#endif
  delete MinCut;
#ifdef USE_DICHOTOMY
  delete prev_MinCut;
#endif
  delete[] NodeArray;

#ifdef DEBUG_USING_MINMEM
  delete[] minmem_trace;
  delete[] minmema_trace;
#endif
}


double MinMemArrayRecurAlgorithm(int N, int *prnts, double *nwghts, double *ewghts, int *schedule){
  double Mr;
  double Mp;
  uint64_t count = 0;
  int * chstart,*chend,*children;
  int root;

  po_construct(N, prnts, &chstart,&chend,&children, &root);
  Mp = MaxOutDegree(N, nwghts, ewghts, chstart,children);


  MinMemArray(N,root,nwghts,ewghts,chstart,children,Mp, Mr, schedule, true, count,prnts);

  delete[] chstart;
  delete[] chend;
  delete[] children;

  return Mr;
}

double MinMemArrayRecurAlgorithm_timed(int N, int *prnts, double *nwghts, double *ewghts, int *schedule,double * usec,int quiet){
  double Mr;
  double Mp;
  uint64_t count = 0;
  int * chstart,*chend,*children;
  int root;

  po_construct(N, prnts, &chstart,&chend,&children, &root);
  Mp = MaxOutDegree(N, nwghts, ewghts, chstart,children);

  *usec = -u_wseconds();
  MinMemArray(N,root,nwghts,ewghts,chstart,children,Mp, Mr, schedule, quiet, count,prnts);
  *usec += u_wseconds();

  delete[] chstart;
  delete[] chend;
  delete[] children;
  return Mr;
}
