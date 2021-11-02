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
#include <algorithm>
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

typedef enum {
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

class Task {
protected:
    bool cost_computed;
    double cost;
    double edge_weight;
    double node_weight;
    double MS_weight = 0;//assume execution time for any node is larger than 0
    double makespan_nocommu;
    bool makespan_computed = false;
    vector<Task *> *children;
    Task *parent;
    unsigned int parent_id;
    unsigned int id;
    bool broken = false;
    int label;
    double MS_sequentialPart, MS_parallelPart;
    double makespan_difference;
    unsigned int Qtree_id;

public :
    Task() {
        id = 0;
        Mpeak = 0;
        parent_id = 0;
        Mavail = 0;
        parent = 0;
        cost_computed = false;
        children = new vector<Task *>();
    }

    Task(double nw, double ew, double mw) {
        id = 0;
        Mpeak = 0;
        parent_id = 0;
        Mavail = 0;
        parent = 0;
        cost_computed = false;
        children = new vector<Task *>();

        edge_weight = ew;
        node_weight = nw;
        MS_weight = mw;
        makespan_nocommu = mw;
    }

    Task(unsigned int pparent_id, double nw, double ew, double mw) {
        id = 0;
        Mpeak = 0;
        Mavail = 0;
        parent = 0;
        cost_computed = false;
        makespan_computed = false;
        children = new vector<Task *>();

        edge_weight = ew;
        node_weight = nw;
        MS_weight = mw;
        makespan_nocommu = mw;
        parent_id = pparent_id;
    }

    ~Task() {
        for (vector<Task *>::iterator iter = children->begin(); iter != children->end(); iter++) {
            delete *iter;
        }
        delete children;
    }

    void SetMSDiff(double slack) {
        makespan_difference = slack;
    }

    double GetMSDiff() {
        return makespan_difference;
    }

    void SetParent(Task *pparent) {
        this->parent = pparent;
    }

    void AddChild(Task *pchild) {
        this->children->push_back(pchild);
        cost_computed = false;
    }

    vector<Task *> *GetChildren() {
        return children;
    }

    Task *GetChild(unsigned int node_id) {
        return children->at(node_id);
    }

    Task *GetParent() {
        return parent;
    }

    bool IsLeaf() const {
        return children->size() == 0;
    }

    bool IsRoot() const {
        return parent_id == 0;
    }

    double GetCost() {
        if (!cost_computed) {
            cost = edge_weight + node_weight;
            for (vector<Task *>::iterator iter = children->begin(); iter != children->end(); iter++) {
                cost += (*iter)->GetEW();
            }
            cost_computed = true;
        }
        return cost;
    }

    void SetParentId(unsigned int pparent_id) {
        parent_id = pparent_id;
    }

    void SetId(unsigned int pid) {
        id = pid;
    }

    void SetLabel(int pid) {
        label = pid;
    }

    void SetEW(double ew) {
        edge_weight = ew;
    }

    void SetNW(double nw) {
        node_weight = nw;
    }

    void SetMSW(double mw) {
        MS_weight = mw;
    }

    unsigned int GetParentId() const {
        return parent_id;
    }

    double GetEW() const {
        return edge_weight;
    }

    double GetNW() const {
        return node_weight;
    }

    double GetMSW() const {
        return MS_weight;
    }

    unsigned int GetId() const {
        return id;
    }

    int GetLabel() const {
        return label;
    }

    void Print(ostream &out) const {
        out << max((unsigned int) 0, GetParentId()) << " " << GetNW() << " " << GetEW() << endl;
        for (vector<Task *>::iterator iter = children->begin(); iter != children->end(); iter++) {
            (*iter)->Print(out);
        }

    }

    void BreakEdge() {
        broken = true;//break this edge
    }

    void RestoreEdge() {
        broken = false;//resotre this edge
    }

    bool IsBroken() {
        if (broken == true) {
            return true;
        } else {
            return false;
        }
    }

    void updateMSCost() {
        makespan_nocommu = MS_sequentialPart + MS_parallelPart;
    }

    double GetSequentialPart() {
        return MS_sequentialPart;
    }

    double GetParallelPart() {
        return MS_parallelPart;
    }

    double GetMSsequential(bool updateEnforce, double &MS_parallel) {
        if ((makespan_computed == true) & (updateEnforce == false)) {
            return MS_sequentialPart;
        }

        MS_sequentialPart = MS_weight;
        MS_parallelPart = 0;
        double temp;
        for (vector<Task *>::iterator iter = this->GetChildren()->begin(); iter != this->GetChildren()->end(); ++iter) {
            if ((*iter)->IsBroken()) {
                //cout<<"edge "<<(*iter)->GetId()<<" broken"<<endl;
                temp = (*iter)->GetMSCost(true, updateEnforce);
                if (temp > MS_parallelPart) {
                    MS_parallelPart = temp;
                }
            } else {
                MS_sequentialPart += (*iter)->GetMSsequential(updateEnforce, MS_parallelPart);
                if (updateEnforce == true) {
                    (*iter)->updateMSCost();
                }
            }
        }

        if (MS_parallelPart > MS_parallel) {
            MS_parallel = MS_parallelPart;
        }

        return MS_sequentialPart;
    }

    double GetMSminusComu() {
        if (Cluster::getFixedCluster()->isHomogeneous()) {
            return (makespan_nocommu - edge_weight / Cluster::getFixedCluster()->getHomogeneousBandwidth());
        } else throw "Cluster not homogeneous";
    }

    double GetMSminusW() {
        if (Cluster::getFixedCluster()->isHomogeneous()) {
            return (makespan_nocommu + edge_weight / Cluster::getFixedCluster()->getHomogeneousBandwidth() - MS_weight);
        } else throw "Cluster not homogeneous";
    }

    void SetMSUncomputed() {
        makespan_computed = false;
    }

    double GetMSCost(bool commulication = false, bool updateEnforce = false) {
        if (!Cluster::getFixedCluster()->isHomogeneous()) throw "Cluster not homogeneous";

        if ((makespan_computed == true) & (updateEnforce == false)) {
            if (commulication == true) {
                return makespan_nocommu + edge_weight / Cluster::getFixedCluster()->getHomogeneousBandwidth();
            } else {
                return makespan_nocommu;
            }
        }

        MS_parallelPart = 0;
        MS_sequentialPart = this->GetMSsequential(updateEnforce, MS_parallelPart);//MS_parallelPart will be update here.
        makespan_nocommu = MS_sequentialPart + MS_parallelPart;

        makespan_computed = true;
        if (commulication == true) {
            //cout<<id<<"-"<<makespan_nocommu<<endl;//test
            return makespan_nocommu + edge_weight / Cluster::getFixedCluster()->getHomogeneousBandwidth();
        }

        //cout<<id<<"-"<<makespan_nocommu<<endl; //test
        return makespan_nocommu;
    }

    void SetothersideID(unsigned int qtreeID) {
        Qtree_id = qtreeID;
    }

    unsigned int GetothersideID() {
        return Qtree_id;
    }

    void RemoveChild(unsigned int childId) {
        for (vector<Task *>::iterator iter = this->children->begin(); iter != this->children->end(); ++iter) {
            if ((*iter)->GetId() == childId) {
                this->children->erase(iter);
                break;
            }
        }
    }

    void MergetoParent() {
        this->GetParent()->SetMSW(this->GetMSW() + this->GetParent()->GetMSW());
        this->GetParent()->RemoveChild(this->id);
        this->GetParent()->GetChildren()->insert(this->GetParent()->GetChildren()->end(), this->children->begin(),
                                                 this->children->end());
        //cout<<", children: ";
        for (vector<Task *>::iterator iter = this->children->begin(); iter != this->children->end(); ++iter) {
            //cout<<(*iter)->GetId()<<" ";
            (*iter)->SetParent(this->GetParent());
            (*iter)->SetParentId(this->GetParent()->GetId());
        }
        //cout<<endl;
        this->children->clear();
        this->~Task();
    }

    unsigned int Ci;
    double Mpeak;
    double Mavail;

    double Sequence();

    double SplitSubtrees(bool twolevel, list<Task *> &parallelRoots, unsigned long &sequentialLength);

    list<Task *> fillParallelRootsUntilBestMakespan(vector<double> &makespansOfSplittings,
                                                    unsigned long stepsUntilMinimalMakespan) const;


};

class Tree {
protected:
    vector<Task *> *nodes;
    unsigned int root_index;
    unsigned int root_count;
    unsigned int offset_id;
    unsigned int tree_id;

    static Tree *originalTree;

public:

    Tree() {
        root_index = 0;
        root_count = 0;
        offset_id = 0;
        tree_id = 1;
        nodes = new vector<Task *>();
    }

    Tree(int N, int *prnts, double *nwghts, double *ewghts, double *mswghts) {
        root_index = 1;
        root_count = 0;
        offset_id = 0;
        tree_id = 1;
        nodes = new vector<Task *>();

        this->AllocateNodes(N);

        for (int i = 1; i < N + 1; i++) {
            //cout << "node id: " << i<<endl;
            Task *cur_node = this->GetNode(i);
            cur_node->GetChildren()->clear();
            cur_node->SetEW(ewghts[i]);
            cur_node->SetNW(nwghts[i]);
            cur_node->SetMSW(mswghts[i]);
            cur_node->SetId(i);
            cur_node->SetLabel(i);
        }

        for (int i = 1; i < N + 1; i++) {
            Task *cur_node = this->GetNode(i);

            if (prnts[i] > 0) {
                cur_node->SetParentId(prnts[i]);
                cur_node->SetParent(this->GetNode(prnts[i]));
                this->GetNode(prnts[i])->AddChild(cur_node);
            } else {
                cur_node->SetParentId(0);
                this->SetRootId(i);
                this->SetTreeId(i);
            }
        }
    }

    Tree(vector<Task *> nodes, Tree *originalTree) {
        root_index = 1;
        root_count = 0;
        offset_id = 0;
        tree_id = 1;

        *(this->nodes) = nodes;
        this->originalTree = originalTree;
    }


    ~Tree() {
        if (root_index != 0 && nodes->size() > 0) {
            delete GetRoot();
        }

        delete nodes;
    }


    void Print(ostream &out) const {
        out << nodes->size() << endl;

        for (vector<Task *>::iterator iter = nodes->begin(); iter != nodes->end(); iter++) {
            out << max((unsigned int) 0, (*iter)->GetParentId()/*+1-offset_id*/) << " " << (*iter)->GetNW() << " "
                << (*iter)->GetEW() << endl;
        }

    }


    void AllocateNodes(int new_node_count) {
        if (root_count > 0 && nodes->size() > 0) {
            delete GetRoot();
        }

        nodes->resize(new_node_count);

        unsigned int i = 0;
        for (vector<Task *>::iterator iter = nodes->begin(); iter != nodes->end(); iter++) {
            *iter = new Task();
            (*iter)->SetId(i++);
        }

        offset_id = nodes->front()->GetId();
    }
    void reverse_vector(){
        reverse(nodes->begin(),nodes->end());
    }

    void AddNode(Task *newNode) {
        nodes->push_back(newNode);
    }

    void AddRoot(Task *newNode) {
        root_count++;
        assert(root_count == 1);
        nodes->push_back(newNode);
        root_index = nodes->size() - 1;
    }

    Task *GetRoot() const {
        return nodes->at(root_index - 1);
    }

    unsigned int GetRootId() const {
        return root_index;
    }

    void SetRootId(unsigned int root_id) {
        root_index = root_id;
    }

    void SetTreeId(unsigned int _id) {
        tree_id = _id;
    }

    Task *GetNode(unsigned int node_id) const {
// in most cases, the task with a specific ID is found on the position that is exactly before the ID
// however, if this is not the case we do a linear search to find where the task is

        unsigned int i = 0;
        Task * task =nodes->at(node_id - 1);
        if (task->GetId() == node_id){
                return task;
        }else{
            
            unsigned long treeSize = this->GetNodes()->size();
            for (i=0 ;i<treeSize;i++) {
                task = this->GetNodeByPos(i);
                if (task->GetId() == node_id){
                    return task;
                }
            }
        }
        cout << this->GetNodeByPos(i-1)->GetId();
        throw "Task not Found!";
    }

    Task *GetNodeByPos(unsigned int node_idx) const {
        assert(node_idx<this->GetNodes()->size());
        return nodes->at(node_idx);
    }

    const vector<Task *> *GetNodes() const {
        return nodes;
    }

    void addNode(Task *newnode) {
        this->nodes->push_back(newnode);
    }

    static void setOriginalTree(Tree *origTree) {
        Tree::originalTree = origTree;
    }

    static Tree *getOriginalTree() {
        return Tree::originalTree;
    }

    void printBrokenEdges() {
        cout << "Print broken edges" << endl;
        unsigned long treeSize = this->GetNodes()->size();
        for (unsigned int i = treeSize; i >= 1; --i) {
            Task *currentnode = this->GetNode(i);
            if (currentnode->IsBroken()) {
                cout << i << " ";
            }
        }
        cout << "End" << endl;


    }

    Tree *BuildQtree();

    unsigned int HowmanySubtrees(bool quiet);

    bool MemoryEnough(Task *Qrootone, Task *Qroottwo, bool leaf, double memory_size, int *chstart, int *children);

    double ImprovedSplit(unsigned int number_processor, int *chstart, int *childrenID);

    double Merge(unsigned int num_subtrees, int *chstart, int *childrenID, bool CheckMemory);

    double MergeV2(unsigned int num_subtrees, unsigned int processor_number, double const memory_size, int *chstart,
                   int *childrenID, bool CheckMemory);

    double SplitAgain();

    double SplitAgainV2(unsigned int processor_number, unsigned int num_subtrees, std::map<int, int> &taskToPrc,
                        std::map<int, bool> &isProcBusy);

    vector<Task *> buildCriticalPath();
};


typedef list<int> schedule_t;

typedef map<unsigned int, double> io_map;
typedef pair<unsigned int, unsigned int> node_sche;
typedef pair<unsigned int, double> node_ew;


void parse_tree(const char *filename, int *N, int **prnts, double **nwghts, double **ewghts, double **mswghts);
Tree * read_tree(const char *filename);

extern "C"
{
void po_construct(const int N, const int *prnts, int **chstart, int **chend, int **children, int *root);
void poaux(const int *chstart, const int *children, int N, int r, int *por, int *label);
} /* closing brace for extern "C" */

double MaxOutDegree(Tree *tree, int quiet);

double MaxOutDegree(int N, double *nwghts, double *ewghts, int *chstart, int *children);

double IOCounter(Tree *tree, int N, double *nwghts, double *ewghts, int *chstart, int *children, int *schedule,
                 double available_memory, bool divisible, int quiet, unsigned int &com_freq,
                 vector<unsigned int> *brokenEdges, io_method_t method);

double
IOCounterWithVariableMem(Tree *tree, int N, double *nwghts, double *ewghts, int *chstart, int *children, int *schedule,
                         vector<double> availableMemorySizesA2, int &currentProcessor,
                         std::map<int, int> &taskToPrc, std::map<int, bool> &isProcBusy, bool divisible, int quiet,
                         unsigned int &com_freq, vector<unsigned int> *brokenEdges, io_method_t method);

Tree *BuildSubtree(Tree *tree, Task *SubtreeRoot, unsigned int new_tree_size, int **prnts, double **ewghts,
                   double **timewghts, double **spacewghts, int *chstart, int *children);

void popSmallestRootsToFitToCluster(list<Task *> &parallelRoots, unsigned long amountSubtrees);

void breakPreparedEdges(Task *root, list<Task *> &parallelRoots);

double getWeightPQ(list<Task *> &parallelRoots, Task *currentNode);

double getWeightSurplusFromSmallestNodes(list<Task *> &parallelRoots);

#endif
#endif
