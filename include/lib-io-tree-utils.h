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
#include "cluster.h"

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

class Cnode
{
protected:
  bool cost_computed;
  double cost;
  double edge_weight;
  double node_weight;
  double MS_weight = 0; //assume execution time for any node is larger than 0
  double makespan_nocommu;
  bool makespan_computed = false;
  vector<Cnode *> *children;
  Cnode *parent;
  unsigned int parent_id;
  unsigned int id;
  bool broken = false;
  int label;
  double MS_sequentialPart, MS_parallelPart;
  double makespan_difference;
  unsigned int Qtree_id;

public:
  Cnode()
  {
    id = 0;
    Mpeak = 0;
    parent_id = 0;
    Mavail = 0;
    parent = 0;
    cost_computed = false;
    children = new vector<Cnode *>();
  }
  Cnode(double nw, double ew, double mw)
  {
    id = 0;
    Mpeak = 0;
    parent_id = 0;
    Mavail = 0;
    parent = 0;
    cost_computed = false;
    children = new vector<Cnode *>();

    edge_weight = ew;
    node_weight = nw;
    MS_weight = mw;
    makespan_nocommu = mw;
  }
  Cnode(unsigned int pparent_id, double nw, double ew, double mw)
  {
    id = 0;
    Mpeak = 0;
    Mavail = 0;
    parent = 0;
    cost_computed = false;
    makespan_computed = false;
    children = new vector<Cnode *>();

    edge_weight = ew;
    node_weight = nw;
    MS_weight = mw;
    makespan_nocommu = mw;
    parent_id = pparent_id;
  }

  ~Cnode()
  {
    for (vector<Cnode *>::iterator iter = children->begin(); iter != children->end(); iter++)
    {
      delete *iter;
    }
    delete children;
  }

  void SetMSDiff(double slack)
  {
    makespan_difference = slack;
  }

  double GetMSDiff()
  {
    return makespan_difference;
  }

  void SetParent(Cnode *pparent)
  {
    this->parent = pparent;
  }

  void AddChild(Cnode *pchild)
  {
    this->children->push_back(pchild);
    cost_computed = false;
  }

  vector<Cnode *> *GetChildren()
  {
    return children;
  }

  Cnode *GetChild(unsigned int node_id)
  {
    return children->at(node_id);
  }

  Cnode *GetParent()
  {
    return parent;
  }

  bool IsLeaf() const
  {
    return children->size() == 0;
  }

  bool IsRoot() const
  {
    return parent_id == 0;
  }

  double GetCost()
  {
    if (!cost_computed)
    {
      cost = edge_weight + node_weight;
      for (vector<Cnode *>::iterator iter = children->begin(); iter != children->end(); iter++)
      {
        cost += (*iter)->GetEW();
      }
      cost_computed = true;
    }
    return cost;
  }

  void SetParentId(unsigned int pparent_id)
  {
    parent_id = pparent_id;
  }

  void SetId(unsigned int pid)
  {
    id = pid;
  }

  void SetLabel(int pid)
  {
    label = pid;
  }

  void SetEW(double ew)
  {
    edge_weight = ew;
  }

  void SetNW(double nw)
  {
    node_weight = nw;
  }

  void SetMSW(double mw)
  {
    MS_weight = mw;
  }

  unsigned int GetParentId() const
  {
    return parent_id;
  }

  double GetEW() const
  {
    return edge_weight;
  }

  double GetNW() const
  {
    return node_weight;
  }

  double GetMSW() const
  {
    return MS_weight;
  }

  unsigned int GetId() const
  {
    return id;
  }

  int GetLabel() const
  {
    return label;
  }

  void Print(ostream &out) const
  {
    out << max((unsigned int)0, GetParentId()) << " " << GetNW() << " " << GetEW() << endl;
    for (vector<Cnode *>::iterator iter = children->begin(); iter != children->end(); iter++)
    {
      (*iter)->Print(out);
    }
  }

  void BreakEdge()
  {
    broken = true; //break this edge
  }

  void RestoreEdge()
  {
    broken = false; //resotre this edge
  }

  bool IsBorken()
  {
    if (broken == true)
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  void updateMSCost()
  {
    makespan_nocommu = MS_sequentialPart + MS_parallelPart;
  }

  double GetSequentialPart()
  {
    return MS_sequentialPart;
  }

  double GetParallelPart()
  {
    return MS_parallelPart;
  }

  double GetMSsequential(bool updateEnforce, double &MS_parallel)
  {
    if ((makespan_computed == true) & (updateEnforce == false))
    {
      return MS_sequentialPart;
    }

    MS_sequentialPart = MS_weight;
    MS_parallelPart = 0;
    double temp;
    for (vector<Cnode *>::iterator iter = this->GetChildren()->begin(); iter != this->GetChildren()->end(); ++iter)
    {
      if ((*iter)->IsBorken())
      {
        //cout<<"edge "<<(*iter)->GetId()<<" broken"<<endl;
        temp = (*iter)->GetMSCost(true, updateEnforce);
        if (temp > MS_parallelPart)
        {
          MS_parallelPart = temp;
        }
      }
      else
      {
        MS_sequentialPart += (*iter)->GetMSsequential(updateEnforce, MS_parallelPart);
        if (updateEnforce == true)
        {
          (*iter)->updateMSCost();
        }
      }
    }

    if (MS_parallelPart > MS_parallel)
    {
      MS_parallel = MS_parallelPart;
    }

    return MS_sequentialPart;
  }

  double GetMSminusComu()
  {
    //return(this->GetMSCost(false,false)-this->GetEW()/BANDWIDTH);
    return (makespan_nocommu - edge_weight / BANDWIDTH);
  }

  double GetMSminusW()
  {
    return (makespan_nocommu + edge_weight / BANDWIDTH - MS_weight);
  }

  void SetMSUncomputed()
  {
    makespan_computed = false;
  }

  double GetMSCost(bool commulication = false, bool updateEnforce = false)
  {
    if ((makespan_computed == true) & (updateEnforce == false))
    {
      if (commulication == true)
      {
        return makespan_nocommu + edge_weight / BANDWIDTH;
      }
      else
      {
        return makespan_nocommu;
      }
    }

    MS_parallelPart = 0;
    MS_sequentialPart = this->GetMSsequential(updateEnforce, MS_parallelPart); //MS_parallelPart will be update here.
    makespan_nocommu = MS_sequentialPart + MS_parallelPart;

    makespan_computed = true;
    if (commulication == true)
    {
      //cout<<id<<"-"<<makespan_nocommu<<endl;//test
      return makespan_nocommu + edge_weight / BANDWIDTH;
    }

    //cout<<id<<"-"<<makespan_nocommu<<endl; //test
    return makespan_nocommu;
  }

  void SetothersideID(unsigned int qtreeID)
  {
    Qtree_id = qtreeID;
  }

  unsigned int GetothersideID()
  {
    return Qtree_id;
  }

  void RemoveChild(unsigned int childId)
  {
    for (vector<Cnode *>::iterator iter = this->children->begin(); iter != this->children->end(); ++iter)
    {
      if ((*iter)->GetId() == childId)
      {
        this->children->erase(iter);
        break;
      }
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

  unsigned int Ci;
  double Mpeak;
  double Mavail;
};

class Ctree
{
protected:
  vector<Cnode *> *nodes;
  unsigned int root_index;
  unsigned int root_count;
  unsigned int offset_id;
  unsigned int tree_id;

public:
  Ctree()
  {
    root_index = 0;
    root_count = 0;
    offset_id = 0;
    tree_id = 1;
    nodes = new vector<Cnode *>();
  }

  Ctree(int N, int *prnts, double *nwghts, double *ewghts, double *mswghts)
  {
    root_index = 1;
    root_count = 0;
    offset_id = 0;
    tree_id = 1;
    nodes = new vector<Cnode *>();

    this->AllocateNodes(N);

    for (int i = 1; i < N + 1; i++)
    {
      Cnode *cur_node = this->GetNode(i);
      cur_node->GetChildren()->clear();
      cur_node->SetEW(ewghts[i]);
      cur_node->SetNW(nwghts[i]);
      cur_node->SetMSW(mswghts[i]);
      cur_node->SetId(i);
      cur_node->SetLabel(i);
    }

    for (int i = 1; i < N + 1; i++)
    {
      Cnode *cur_node = this->GetNode(i);

      if (prnts[i] > 0)
      {
        cur_node->SetParentId(prnts[i]);
        cur_node->SetParent(this->GetNode(prnts[i]));
        this->GetNode(prnts[i])->AddChild(cur_node);
      }
      else
      {
        cur_node->SetParentId(0);
        this->SetRootId(i);
        this->SetTreeId(i);
      }
    }
    // test
    //    Cnode* currentNode;
    //    for (unsigned int i=2; i<=tree_size; ++i) {
    //        currentNode=treeobj->GetNode(i);
    //        std::cout<<currentNode->GetId()<<" "<<currentNode->GetParent()->GetId()<<"\n";
    //        for (vector<Cnode*>::iterator iter=currentNode->GetChildren()->begin(); iter!=currentNode->GetChildren()->end();
    //             ++iter) {
    //            std::cout<<"   "<<(*iter)->GetId()<<"\n";
    //        }
    //    }
  }

  ~Ctree()
  {
    if (root_index != 0 && nodes->size() > 0)
    {
      delete GetRoot();
    }

    delete nodes;
  }

  //TODO recoder ca en recursif
  void Print(ostream &out) const
  {
    //		if(root_index!=0 && nodes->size()>0){
    //			out<<nodes->size()<<endl;
    //
    //			GetRoot()->Print(out);
    //		}
    //

    out << nodes->size() << endl;

    for (vector<Cnode *>::iterator iter = nodes->begin(); iter != nodes->end(); iter++)
    {
      out << max((unsigned int)0, (*iter)->GetParentId() /*+1-offset_id*/) << " " << (*iter)->GetNW() << " " << (*iter)->GetEW() << endl;
    }
  }

  void AllocateNodes(int new_node_count)
  {
    if (root_count > 0 && nodes->size() > 0)
    {
      delete GetRoot();
    }

    nodes->resize(new_node_count);

    unsigned int i = 0;
    for (vector<Cnode *>::iterator iter = nodes->begin(); iter != nodes->end(); iter++)
    {
      *iter = new Cnode();
      (*iter)->SetId(i++);
    }

    offset_id = nodes->front()->GetId();
  }

  void AddNode(Cnode *newNode)
  {
    nodes->push_back(newNode);
  }

  void AddRoot(Cnode *newNode)
  {
    root_count++;
    assert(root_count == 1);
    nodes->push_back(newNode);
    root_index = nodes->size() - 1;
  }

  Cnode *GetRoot() const
  {
    return nodes->at(root_index - 1);
  }

  unsigned int GetRootId() const
  {
    return root_index;
  }

  void SetRootId(unsigned int root_id)
  {
    root_index = root_id;
  }

  void SetTreeId(unsigned int _id)
  {
    tree_id = _id;
  }

  Cnode *GetNode(unsigned int node_id) const
  {
    return nodes->at(node_id - 1);
  }

  Cnode *GetNodeByPos(unsigned int node_idx) const
  {
    return nodes->at(node_idx);
  }

  const vector<Cnode *> *GetNodes() const
  {
    return nodes;
  }

  void addNode(Cnode *newnode)
  {
    this->nodes->push_back(newnode);
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
typedef list<Cnode *> cut_t;

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

void ConvertToLiu(const Ctree *tree_us, Ctree *tree_liu);
void ConvertToLiu(const int *oldprnts, const double *oldnwghts, const double *oldewghts, int N, const int *chstart, const int *children, int **pprnts, double **pnwghts, double **pewghts);
void parse_tree(const char *filename, Ctree *tree);
void parse_tree(const char *filename, int *N, int **prnts, double **nwghts, double **ewghts, double **mswghts);

extern "C"
{
#endif
  void po_construct(const int N, const int *prnts, int **chstart, int **chend, int **children, int *root);
  void poaux(const int *chstart, const int *children, int N, int r, int *por, int *label);
#ifdef __cplusplus
} /* closing brace for extern "C" */

bool check_schedule(int *prnts, int *sched, int N);

double MaxOutDegree(Ctree *tree, int quiet);
double MaxOutDegree(int N, double *nwghts, double *ewghts, int *chstart, int *children);

void NextValley(Cnode *node, double available_memory, double &cut_value, list<Cnode *> &min_sub_cut, list<unsigned int> &sub_schedule, double &Inc, int quiet, int depth, int &count);
//double IOCounter(Ctree & tree, schedule_t & sub_schedule, double available_memory,bool divisible,int quiet);
//double IOCounter(int N, double * nwghts, double * ewghts, int * chstart,int * children, int * schedule, double available_memory,bool divisible,int quiet, io_method_t method=FURTHEST_NODE);
double IOCounter(Ctree *tree, int N, double *nwghts, double *ewghts, int *chstart, int *children, int *schedule, double available_memory, bool divisible, int quiet, unsigned int &com_freq, vector<unsigned int> *brokenEdges, io_method_t method);
double IOCounterWithVariableMem(Ctree *tree, int N, double *nwghts, double *ewghts, int *chstart, int *children, int *schedule,
                                Cluster *cluster, bool divisible, int quiet, unsigned int &com_freq, vector<unsigned int> *brokenEdges, io_method_t method);
Ctree *BuildSubtree(Ctree *tree, Cnode *SubtreeRoot, unsigned int new_tree_size, int **prnts, double **ewghts, double **timewghts, double **spacewghts, int *chstart, int *children);

#endif
#endif
