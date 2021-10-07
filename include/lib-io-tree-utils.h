/*
 *  lib-io-tree-utils.h
 *  lib-io-tree
 *
 *  Created by defbond on 8/10/10.
 *  Copyright 2010 LIP/ENS-Lyon. All rights reserved.
 *
 */
#ifndef LIB_IO_TREE_UTILS_H
#define LIB_IO_TREE_UTILS_H

#ifdef __cplusplus

#include <ostream>
#include <iostream>
#include <stdio.h>
#include <list>

#include <vector>
#include <forward_list>
#include <assert.h>
#include <map>


#ifndef DEBUG_MEMUSAGE
#define DEBUG_MEMUSAGE 0
#endif
#ifndef VERBOSE
#define VERBOSE 0
#endif
#ifndef STRONG_ASSERT
#define STRONG_ASSERT 0
#endif

using namespace std;

#ifndef MAX_COMBI_SIZE
#define MAX_COMBI_SIZE 5
#endif

extern double BANDWIDTH;

typedef enum
{
  FURTHEST_NODE = 1,
  BEST_K_COMBI,
  BEST_FIT_ABS,
  FIRST_FIT_ABS,
  BEST_FIT,
  FIRST_FIT,
  BEST_INC_COMBI,
  BEST_COMBI,
  LARGEST_FIT,
  IMMEDIATELY
} io_method_t;

double u_wseconds(void);

      }

      nodes->resize(new_node_count);

      unsigned int i = 0;
      for(vector<Task*>::iterator iter = nodes->begin();iter!=nodes->end();iter++){
        *iter = new Task();
        (*iter)->SetId(i++);
      }
    }

    void AddNode(Task * newNode){
      nodes->push_back(newNode);
    }

    void AddRoot(Task * newNode){
      assert(root_count == 1);
      nodes->push_back(newNode);
      root_index = nodes->size()-1;
    }

    Task * GetRoot() const{
      return nodes->at(root_index-1);
    }

  void SetothersideID(unsigned int qtreeID)
  {
    Qtree_id = qtreeID;
  }

  unsigned int GetothersideID()
  {
    return Qtree_id;
  }

    void SetTreeId(unsigned int _id){
        tree_id = _id;
    }
    
    Task * GetNode(unsigned int node_id) const{
      return nodes->at(node_id-1);
    }
  }

  void MergetoParent()
  {
    this->GetParent()->SetMSW(this->GetMSW() + this->GetParent()->GetMSW());
    this->GetParent()->RemoveChild(this->id);
    this->GetParent()->GetChildren()->insert(this->GetParent()->GetChildren()->end(), this->children->begin(), this->children->end());
    //cout<<", children: ";
    for (vector<Cnode *>::iterator iter = this->children->begin(); iter != this->children->end(); ++iter)
    {
      //cout<<(*iter)->GetId()<<" ";
      (*iter)->SetParent(this->GetParent());
      (*iter)->SetParentId(this->GetParent()->GetId());
    }
    //cout<<endl;
    this->children->clear();
    this->~Cnode();
  }

  // bool operator < (const Cnode& other) const
  // {
  //     double thiscost = GetMSCost();
  //      return ( thiscost < other.GetMSCost());
  //  }

    Task * GetNodeByPos(unsigned int node_idx) const{
      return nodes->at(node_idx);
    }

    const vector<Task*> * GetNodes() const{
      return nodes;
    }
    
    void addNode(Task* newnode){
        this->nodes->push_back(newnode);
    }
    
    void setOriginalTree(Tree* origTree){
      if (!this->originalTreeInitialized){
        this->originalTree = origTree;
        this->originalTreeInitialized = true;
      }
      else{
        cout << "Original Tree can only be set once!"<<endl;
        exit (EXIT_FAILURE);
      }
    }

    Tree * getOriginalTree(){
      return this->originalTree;
    }
};

struct s_node_t
{
  int parent;
  vector<int> children;
  double edge_weight;
  double node_weight;
  double Mpeak;
  double Mavail;
};

struct s_io_t
{
  unsigned int node;
  double unloaded_data;
};

typedef list<int> schedule_t;
typedef list<Task*> cut_t;

typedef pair<unsigned int, double> io_t;
typedef map<unsigned int, double> io_map;
typedef pair<unsigned int, unsigned int> node_sche;
typedef pair<unsigned int, double> node_ew;

struct OrdoLiu_t;
struct val_seg_t
{
  schedule_t::iterator begin;
  unsigned int begin_index;
  schedule_t::iterator end;
  unsigned int end_index;
  OrdoLiu_t *orig_ordo;
  double value;
};

struct OrdoLiu_t
{
  double max_pebble_cost;
  double fi;
  list<val_seg_t> val_seg;
  schedule_t schedule;
}; 


void ConvertToLiu(const Tree * tree_us, Tree * tree_liu) ;
void ConvertToLiu(const int * oldprnts,const double * oldnwghts,const double * oldewghts, int N,const int* chstart,const int * children, int ** pprnts, double ** pnwghts, double ** pewghts);
void parse_tree(const char *filename,Tree * tree);
void parse_tree(const char *filename,int * N ,int **prnts,double **nwghts,double **ewghts, double **mswghts);

extern "C"
{
#endif
  void po_construct(const int N, const int *prnts, int **chstart, int **chend, int **children, int *root);
  void poaux(const int *chstart, const int *children, int N, int r, int *por, int *label);
#ifdef __cplusplus
} /* closing brace for extern "C" */

bool check_schedule(int * prnts,int * sched,int N);

double MaxOutDegree(Tree * tree,int quiet);
double MaxOutDegree(int N, double * nwghts, double * ewghts, int * chstart,int * children);

void NextValley(Task * node, double available_memory,  double & cut_value, list<Task*> & min_sub_cut, list<unsigned int> & sub_schedule, double & Inc, int quiet, int depth,int & count);
//double IOCounter(Tree & tree, schedule_t & sub_schedule, double available_memory,bool divisible,int quiet);
//double IOCounter(int N, double * nwghts, double * ewghts, int * chstart,int * children, int * schedule, double available_memory,bool divisible,int quiet, io_method_t method=FURTHEST_NODE);
double IOCounter(Tree* tree, int N, double * nwghts, double * ewghts, int * chstart,int * children, int * schedule, double available_memory,bool divisible,int quiet,unsigned int & com_freq, vector<unsigned int>* brokenEdges, io_method_t method);
double IOCounterWithVariableMem(Tree* tree, int N, double * nwghts, double * ewghts, int * chstart,int * children, int * schedule, vector<double> availableMemorySizesA2, int &currentProcessor,
                                         std::map<int, int> &taskToPrc, std::map<int, bool> &isProcBusy, bool divisible,int quiet,unsigned int & com_freq, vector<unsigned int>* brokenEdges, io_method_t method);
Tree* BuildSubtree(Tree* tree, Task* SubtreeRoot, unsigned int new_tree_size, int** prnts, double** ewghts, double** timewghts, double** spacewghts, int * chstart, int * children);

#endif
#endif
