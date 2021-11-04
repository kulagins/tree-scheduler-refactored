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
    bool __root{};

public :
    unsigned int Ci;
    double Mpeak;
    double Mavail;

    Task() {
        id = 0;
        Mpeak = 0;
        parent_id = 0;
        Mavail = 0;
        parent = 0;
        cost_computed = false;
        children = new vector<Task *>();
    }

    Task(double nw, double ew, double mw, bool root = false) {
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
        if (root) {
            this->__root = true;
        } else {
            this->__root = false;
        }
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

    Task(const Task &otherTask, const unsigned int newId, Task *newParent)
            : cost_computed{otherTask.cost_computed},
              cost{otherTask.cost},
              edge_weight{otherTask.edge_weight},
              node_weight{otherTask.node_weight},
              MS_weight{otherTask.MS_weight},
              makespan_nocommu{otherTask.makespan_nocommu},
              makespan_computed{otherTask.makespan_computed},
              parent{newParent ? newParent : 0},
              parent_id{newParent ? newParent->getId() : 0},
              id{newId},
              broken{otherTask.broken},
              label{otherTask.label},
              MS_sequentialPart{otherTask.MS_sequentialPart},
              MS_parallelPart{otherTask.MS_parallelPart},
              makespan_difference{otherTask.makespan_difference},
              Qtree_id{otherTask.Qtree_id},
              Ci{otherTask.Ci},
              Mpeak{otherTask.Mpeak},
              Mavail{otherTask.Mavail} {
        children = new vector<Task *>();
        __root = otherTask.isRoot();
    }

    ~Task() {
        for (vector<Task *>::iterator iter = children->begin(); iter != children->end(); iter++) {
            delete *iter;
        }
        delete children;
    }

    void setMakespanDiff(double slack) {
        makespan_difference = slack;
    }

    double getMakespanDiff() {
        return makespan_difference;
    }

    void setParent(Task *pparent) {
        this->parent = pparent;
    }

    void addChild(Task *pchild) {
        this->children->push_back(pchild);
        cost_computed = false;
    }

    vector<Task *> *getChildren() {
        return children;
    }

    Task *getChild(unsigned int node_id) {
        return children->at(node_id);
    }

    Task *getParent() {
        return parent;
    }

    bool isLeaf() const {
        return children->size() == 0;
    }

    bool isRoot() const {
        return this->__root;
    }

    double getCost() {
        if (!cost_computed) {
            cost = edge_weight + node_weight;
            for (vector<Task *>::iterator iter = children->begin(); iter != children->end(); iter++) {
                cost += (*iter)->getEdgeWeight();
            }
            cost_computed = true;
        }
        return cost;
    }

    void setParentId(unsigned int pparent_id) {
        parent_id = pparent_id;
    }

    void setId(unsigned int pid) {
        id = pid;
    }

    void setLabel(int pid) {
        label = pid;
    }

    void setEdgeWeight(double ew) {
        edge_weight = ew;
    }

    void setNodeWeight(double nw) {
        node_weight = nw;
    }

    void setMakespanWeight(double mw) {
        MS_weight = mw;
    }

    unsigned int getParentId() const {
        return parent_id;
    }

    double getEdgeWeight() const {
        return edge_weight;
    }

    double getNodeWeight() const {
        return node_weight;
    }

    double getMakespanWeight() const {
        return MS_weight;
    }

    unsigned int getId() const {
        return id;
    }

    int getLabel() const {
        return label;
    }

    void Print(ostream &out) const {
        out << max((unsigned int) 0, getParentId()) << " " << getNodeWeight() << " " << getEdgeWeight() << endl;
        for (vector<Task *>::iterator iter = children->begin(); iter != children->end(); iter++) {
            (*iter)->Print(out);
        }

    }

    void breakEdge() {
        broken = true;//break this edge
    }

    void restoreEdge() {
        broken = false;//resotre this edge
    }

    bool isBroken() {
        if (broken == true) {
            return true;
        } else {
            return false;
        }
    }

    void updateMakespanCost() {
        makespan_nocommu = MS_sequentialPart + MS_parallelPart;
    }

    double getSequentialPart() {
        return MS_sequentialPart;
    }

    double getParallelPart() {
        return MS_parallelPart;
    }

    double getMakespanSequential(bool updateEnforce, double &MS_parallel) {
        if ((makespan_computed == true) & (updateEnforce == false)) {
            return MS_sequentialPart;
        }

        MS_sequentialPart = MS_weight;
        MS_parallelPart = 0;
        double temp;
        for (vector<Task *>::iterator iter = this->getChildren()->begin(); iter != this->getChildren()->end(); ++iter) {
            if ((*iter)->isBroken()) {
                //cout<<"edge "<<(*iter)->getId()<<" broken"<<endl;
                temp = (*iter)->getMakespanCost(true, updateEnforce);
                if (temp > MS_parallelPart) {
                    MS_parallelPart = temp;
                }
            } else {
                MS_sequentialPart += (*iter)->getMakespanSequential(updateEnforce, MS_parallelPart);
                if (updateEnforce == true) {
                    (*iter)->updateMakespanCost();
                }
            }
        }

        if (MS_parallelPart > MS_parallel) {
            MS_parallel = MS_parallelPart;
        }

        return MS_sequentialPart;
    }

    double getMakespanMinusComu() {
        if (Cluster::getFixedCluster()->isHomogeneous()) {
            return (makespan_nocommu - edge_weight / Cluster::getFixedCluster()->getHomogeneousBandwidth());
        } else throw "Cluster not homogeneous";
    }

    double getMakespanMinusW() {
        if (Cluster::getFixedCluster()->isHomogeneous()) {
            return (makespan_nocommu + edge_weight / Cluster::getFixedCluster()->getHomogeneousBandwidth() - MS_weight);
        } else throw "Cluster not homogeneous";
    }

    void setMakespanUncomputed() {
        makespan_computed = false;
    }

    double getMakespanCost(bool commulication = false, bool updateEnforce = false) {
        if (!Cluster::getFixedCluster()->isHomogeneous()) throw "Cluster not homogeneous";

        if ((makespan_computed == true) & (updateEnforce == false)) {
            if (commulication == true) {
                return makespan_nocommu + edge_weight / Cluster::getFixedCluster()->getHomogeneousBandwidth();
            } else {
                return makespan_nocommu;
            }
        }

        MS_parallelPart = 0;
        MS_sequentialPart = this->getMakespanSequential(updateEnforce,
                                                        MS_parallelPart);//MS_parallelPart will be update here.
        makespan_nocommu = MS_sequentialPart + MS_parallelPart;

        makespan_computed = true;
        if (commulication == true) {
            //cout<<id<<"-"<<makespan_nocommu<<endl;//test
            return makespan_nocommu + edge_weight / Cluster::getFixedCluster()->getHomogeneousBandwidth();
        }

        //cout<<id<<"-"<<makespan_nocommu<<endl; //test
        return makespan_nocommu;
    }

    void setOtherSideId(unsigned int qtreeID) {
        Qtree_id = qtreeID;
    }

    unsigned int getOtherSideId() {
        return Qtree_id;
    }

    void removeChild(unsigned int childId) {
        for (vector<Task *>::iterator iter = this->children->begin(); iter != this->children->end(); ++iter) {
            if ((*iter)->getId() == childId) {
                this->children->erase(iter);
                break;
            }
        }
    }

    void mergeToParent() {
        this->getParent()->setMakespanWeight(this->getMakespanWeight() + this->getParent()->getMakespanWeight());
        this->getParent()->removeChild(this->id);
        this->getParent()->getChildren()->insert(this->getParent()->getChildren()->end(), this->children->begin(),
                                                 this->children->end());
        //cout<<", children: ";
        for (vector<Task *>::iterator iter = this->children->begin(); iter != this->children->end(); ++iter) {
            //cout<<(*iter)->getId()<<" ";
            (*iter)->setParent(this->getParent());
            (*iter)->setParentId(this->getParent()->getId());
        }
        //cout<<endl;
        this->children->clear();
        this->~Task();
    }

    double Sequence();

    double SplitSubtrees(bool twolevel, list<Task *> &parallelRoots, unsigned long &sequentialLength);

    list<Task *> fillParallelRootsUntilBestMakespan(vector<double> &makespansOfSplittings,
                                                    unsigned long stepsUntilMinimalMakespan) const;


};

class Tree {
protected:
    vector<Task *> *tasks;
    unsigned int root_count;
    unsigned int offset_id;
    unsigned int tree_id;

    Task *root;
    static Tree *originalTree;

public:

    Tree() {
        root_count = 0;
        offset_id = 0;
        tree_id = 1;
        tasks = new vector<Task *>();
    }

    Tree(int N, int *prnts, double *nwghts, double *ewghts, double *mswghts) {
        root_count = 0;
        offset_id = 0;
        tree_id = 1;
        tasks = new vector<Task *>();

        this->allocateTasks(N);

        for (int i = 1; i < N + 1; i++) {
            //cout << "node id: " << i<<endl;
            Task *cur_node = this->getTask(i);
            cur_node->getChildren()->clear();
            cur_node->setEdgeWeight(ewghts[i]);
            cur_node->setNodeWeight(nwghts[i]);
            cur_node->setMakespanWeight(mswghts[i]);
            cur_node->setId(i);
            cur_node->setLabel(i);
        }

        for (int i = 1; i < N + 1; i++) {
            Task *cur_node = this->getTask(i);

            if (prnts[i] > 0) {
                cur_node->setParentId(prnts[i]);
                cur_node->setParent(this->getTask(prnts[i]));
                this->getTask(prnts[i])->addChild(cur_node);
            } else {
                cur_node->setParentId(0);
                this->setRootId(i);
                this->setTreeId(i);
            }
        }
    }

    Tree(vector<Task *> *nodes, Task * root, Tree *originalTree) {
        root_count = 1;
        offset_id = 0;
        tree_id = 1;
        root = root;
        setRootId(root->getId());
        this->tasks = nodes;
        this->originalTree = originalTree;
    }


    ~Tree() {
        if (this->getRootId() != 0 && tasks->size() > 0) {
            delete getRoot();
        }

        delete tasks;
    }


    void Print(ostream &out) const {
        out << tasks->size() << endl;

        for (vector<Task *>::iterator iter = tasks->begin(); iter != tasks->end(); iter++) {
            out << max((unsigned int) 0, (*iter)->getParentId()/*+1-offset_id*/) << " " << (*iter)->getNodeWeight()
                << " "
                << (*iter)->getEdgeWeight() << endl;
        }

    }


    void allocateTasks(int new_node_count) {
        if (root_count > 0 && tasks->size() > 0) {
            delete getRoot();
        }

        tasks->resize(new_node_count);

        unsigned int i = 0;
        for (vector<Task *>::iterator iter = tasks->begin(); iter != tasks->end(); iter++) {
            *iter = new Task();
            (*iter)->setId(i++);
        }

        offset_id = tasks->front()->getId();
    }

    void reverseVector() {
        reverse(tasks->begin(), tasks->end());
    }

    void addTask(Task *newNode) {
        tasks->push_back(newNode);
    }

    void addRoot(Task *newNode) {
        root_count++;
        assert(root_count == 1);
        assert(newNode->isRoot());
        tasks->push_back(newNode);
        this->root = newNode;
    }

    Task *getRoot() const {
        assert(root_count == 1);
        return this->root;
    }

    unsigned int getRootId() const {
        assert(root_count == 1);
        return this->getRoot()->getId();
    }

    void setRootId(unsigned int root_id) {
        assert(root_count == 1);
        this->root->setId(root_id);

    }

    void setTreeId(unsigned int _id) {
        tree_id = _id;
    }

    Task *getTask(unsigned int node_id) const {
        Task *task = tasks->at(node_id - 1);
        if (task->getId() == node_id) {
            return task;
        } else {
            for(Task * taskSequential: * this->getTasks()){
                if(taskSequential->getId() == node_id){}
                return task;
            }
            throw runtime_error("Task not found for id " + to_string(node_id));
        }
    }

    Task *getTaskByPos(unsigned int node_idx) const {
        assert(node_idx < this->getTasks()->size());
        return tasks->at(node_idx);
    }

    const vector<Task *> *getTasks() const {
        return tasks;
    }

    static void setOriginalTree(Tree *origTree) {
        Tree::originalTree = origTree;
    }

    static Tree *getOriginalTree() {
        return Tree::originalTree;
    }

    void printBrokenEdges() {
        cout << "Print broken edges" << endl;
        unsigned long treeSize = this->getTasks()->size();
        for (unsigned int i = treeSize; i >= 1; --i) {
            Task *currentnode = this->getTask(i);
            if (currentnode->isBroken()) {
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

Tree *read_tree(const char *filename);

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

Tree *BuildSubtree(Tree *tree, Task *subtreeRoot);

Tree *BuildSubtreeOld(Tree *tree, Task *SubtreeRoot, unsigned int new_tree_size, int **prnts, double **ewghts,
                      double **timewghts, double **spacewghts, int *chstart, int *children);

void popSmallestRootsToFitToCluster(list<Task *> &parallelRoots, unsigned long amountSubtrees);

void breakPreparedEdges(Task *root, list<Task *> &parallelRoots);

double getWeightPQ(list<Task *> &parallelRoots, Task *currentNode);

double getWeightSurplusFromSmallestNodes(list<Task *> &parallelRoots);

#endif
#endif
