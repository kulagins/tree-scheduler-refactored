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
#include <float.h>
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

void freeProcessorIfAvailable(Task *task);

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
    Processor *assignedProcessor;
    vector<Processor *> *feasibleProcessors;
    double minMemUnderlying;

    double tMax;

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
        assignedProcessor = nullptr;
        Qtree_id = 0;
        feasibleProcessors = new vector<Processor *>();
        minMemUnderlying = 0;
        tMax = 0;
    }

    Task(double nw, double ew, double mw, bool root = false) {
        id = 0;
        Mpeak = 0;
        parent_id = 0;
        Mavail = 0;
        parent = 0;
        cost_computed = false;
        children = new vector<Task *>();
        assignedProcessor = nullptr;
        edge_weight = ew;
        node_weight = nw;
        MS_weight = mw;
        makespan_nocommu = mw;
        Qtree_id = 0;
        if (root) {
            this->__root = true;
        } else {
            this->__root = false;
        }
        feasibleProcessors = new vector<Processor *>();
        minMemUnderlying = 0;
        tMax = 0;
    }

    Task(unsigned int pparent_id, double nw, double ew, double mw) {
        id = 0;
        Mpeak = 0;
        Mavail = 0;
        parent = 0;
        cost_computed = false;
        makespan_computed = false;
        children = new vector<Task *>();
        assignedProcessor = nullptr;
        edge_weight = ew;
        node_weight = nw;
        MS_weight = mw;
        makespan_nocommu = mw;
        parent_id = pparent_id;
        Qtree_id = 0;
        feasibleProcessors = new vector<Processor *>();
        tMax = 0;
        minMemUnderlying = 0;
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
        __root = false;
        assignedProcessor = nullptr;
        feasibleProcessors = new vector<Processor *>();
        minMemUnderlying = otherTask.minMemUnderlying;
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
        this->parent_id = pparent->getId();
    }

    void addChild(Task *pchild) {
        //   cout << "add child " << pchild->getId() << " with parent " << pchild->getParentId() <<" to "<< this->getId()<< " and children ";
        //    for(Task* child: *pchild->getChildren()){
        //       cout<<child->getId()<<" ";
        //   }
        //   cout<<endl;
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

    Processor *getAssignedProcessor() const {
        return assignedProcessor;
    }

    void setAssignedProcessor(Processor *assignedProcessor) {
        this->assignedProcessor = assignedProcessor;
    }

    double getAssignedProcessorSpeed() {
        if (this->assignedProcessor == NULL) throw "NO PROCESSOR";
        else return this->assignedProcessor->getProcessorSpeed();
    }

    void toggleRootStatus(bool newStatus) {
        if (this->isRoot() != newStatus) {
            __root = !__root;
        } else {
            throw "Trying to set root to a position that it is alrady in";
        }
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

    vector<Processor *> *getFeasibleProcessors() {
        return this->feasibleProcessors;
    }

    int getLabel() const {
        return label;
    }

    /*  vector<Processor *> *getFeasibleProcessors(){
          return &this->feasibleProcessors;
      }
  */
    //Todo: sort?
    void addFeasibleProcessor(Processor *proc) {
        auto position_it = find(feasibleProcessors->begin(), feasibleProcessors->end(), proc);
        if (position_it == feasibleProcessors->end()) {
            this->feasibleProcessors->push_back(proc);
        }
    }

    void setMinMemUnderlying(double minMem) {
        this->minMemUnderlying = minMem;
    }

    double getMinMemUnderlying() {
        return this->minMemUnderlying;
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
        this->broken = false;//resotre this edge
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
        for (Task *child: *this->getChildren()) {
            if (child->isBroken()) {
                //cout<<"edge "<<(*iter)->getId()<<" broken"<<endl;
                temp = child->getMakespanCost(true, updateEnforce);
                if (temp > MS_parallelPart) {
                    MS_parallelPart = temp;
                }
            } else {
                MS_sequentialPart += child->getMakespanSequential(updateEnforce, MS_parallelPart);
                if (updateEnforce == true) {
                    child->updateMakespanCost();
                }
            }
        }

        if (MS_parallelPart > MS_parallel) {
            MS_parallel = MS_parallelPart;
        }

        return MS_sequentialPart;
    }

    double getMakespanSequentialWithSpeeds(bool updateEnforce, double &MS_parallel) {

        double assignedProcSpeed = this->getAssignedProcessorSpeed();
        if ((makespan_computed == true) & (updateEnforce == false)) {
            return MS_sequentialPart;
        }

        MS_sequentialPart = MS_weight / assignedProcSpeed;
        MS_parallelPart = 0;
        double temp;
        for (Task *child: *this->getChildren()) {
            if (child->isBroken()) {
                //cout<<"edge "<<(*iter)->getId()<<" broken"<<endl;
                temp = child->getMakespanCostWithSpeeds(true, updateEnforce);
                if (temp > MS_parallelPart) {
                    MS_parallelPart = temp;
                }
            } else {
                MS_sequentialPart += child->getMakespanSequentialWithSpeeds(updateEnforce, MS_parallelPart);
                if (updateEnforce == true) {
                    child->updateMakespanCost();
                }
            }
        }

        if (MS_parallelPart > MS_parallel) {
            MS_parallel = MS_parallelPart;
        }

        return MS_sequentialPart;
    }

    void setMakespanUncomputed() {
        makespan_computed = false;
    }

    double getMakespanCost(bool commulication = false, bool updateEnforce = false) {
        if (!(Cluster::getFixedCluster())->isBandwidthHomogeneous()) throw "Cluster not homogeneous";

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
            // cout<<id<<"-"<<makespan_nocommu<<endl;//test
            // cout<< (edge_weight / Cluster::getFixedCluster()->getHomogeneousBandwidth())<<endl;
            return makespan_nocommu + edge_weight / Cluster::getFixedCluster()->getHomogeneousBandwidth();
        }

        //  cout<<id<<"-"<<makespan_nocommu<<endl; //test
        return makespan_nocommu;
    }

    double getMakespanCostWithSpeeds(bool commulication = false, bool updateEnforce = false) {
        if (!(Cluster::getFixedCluster())->isBandwidthHomogeneous()) throw "Cluster not homogeneous";
        // double assignedRootProcessorSpeed = this->getAssignedProcessorSpeed();

        if ((makespan_computed == true) & (updateEnforce == false)) {
            if (commulication == true) {
                return makespan_nocommu + edge_weight / Cluster::getFixedCluster()->getHomogeneousBandwidth();
            } else {
                return makespan_nocommu;
            }
        }

        MS_parallelPart = 0;
        MS_sequentialPart = this->getMakespanSequentialWithSpeeds(updateEnforce,
                                                                  MS_parallelPart);//MS_parallelPart will be update here.
        makespan_nocommu = MS_sequentialPart + MS_parallelPart;

        makespan_computed = true;
        if (commulication == true) {
            return makespan_nocommu + edge_weight / Cluster::getFixedCluster()->getHomogeneousBandwidth();
        }
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

    double SplitSubtrees(bool twolevel, list<Task *> &parallelRoots, unsigned long &sequentialLength, int limit);

    list<Task *> fillParallelRootsUntilBestMakespan(vector<double> &makespansOfSplittings,
                                                    unsigned long stepsUntilMinimalMakespan) const;

    void precomputeMinMems(Tree *tree, bool greedy = false);

    void updateTMax() {
        if (this->feasibleProcessors->size() == 0) {
            //cout << "Task" << to_string(this->getId()) << " has 0 feasible processors during iterations. "
            //       << Cluster::getFixedCluster()->getNumberFreeProcessors() << " are still free.";
            throw "No Schedule Possible";
        }
        if (this->feasibleProcessors->size() == 1) this->tMax = DBL_MAX;
        else {
            double s_max = this->getFastestFeasibleProcessor()->getProcessorSpeed();
            double beta = Cluster::getFixedCluster()->getHomogeneousBandwidth();
            this->tMax = this->edge_weight / beta + this->MS_weight / s_max;
        }
    }

    double getTMax() { return tMax; }

    void setTMax(double tMax) { this->tMax = tMax; }


    // Asserts that feasibleProcessors is ordered!
    Processor *getFastestFeasibleProcessor() {
        sort(feasibleProcessors->begin(), feasibleProcessors->end(),
             [](Processor *a, Processor *b) { return a->getProcessorSpeed() > b->getProcessorSpeed(); });
        return feasibleProcessors->front();
    }

    void deleteFeasible(Processor *proc) {
        auto position_it = find(feasibleProcessors->begin(), feasibleProcessors->end(), proc);
        if (position_it != feasibleProcessors->end()) {
            feasibleProcessors->erase(position_it);
        }
    }

    void assignFeasibleProcessorsToSubtree(double minMem);

    double computeMinMemUnderlying(Tree *tree, bool greedy);

    vector<Task *> tasksInSubtreeRootedHere() {
        vector<Task *> result;
        result.push_back(this);
        vector<Task *> candidates;
        for (Task *child: *this->getChildren()) {
            if (!child->isBroken())
                candidates.push_back(child);
        }

        while (!candidates.empty()) {
            Task *candidate = candidates.back();
            candidates.pop_back();

            result.push_back(candidate);
            for (Task *child: *candidate->getChildren()) {
                if (!child->isBroken())
                    candidates.push_back(child);
            }

        }
        return result;

    }


    vector<Task *> breakNBiggestChildren(int n) {
        vector<Task *> newlyBroken;

        struct {
            bool operator()(Task *a, Task *b) const {
                return a->getMakespanWeight() < b->getMakespanWeight();
            }
        } weightLess;
        std::sort(this->getChildren()->begin(), this->getChildren()->end(), weightLess);
        int border = n < this->getChildren()->size() ? n : this->getChildren()->size();
        //if (n < this->getChildren()->size()) {
        //    cout << "n smaller" << n << " " << this->getChildren()->size() << endl;
        //}
        for (int i = 0; i < border; i++) {
            (this->getChildren()->at(i))->breakEdge();
            newlyBroken.push_back(this->getChildren()->at(i));
        }
        return newlyBroken;
    }

    void restoreBrokenChildren() {
        for (Task *child: *this->getChildren()) {
            child->restoreEdge();
        }
    }

    bool isAnyChildBroken() {
        for (Task *child: *this->getChildren()) {
            if (child->isBroken()) return true;
        }
        return false;
    }
};


class Tree {
protected:
    vector<Task *> *tasks;
    unsigned int root_count;
    unsigned int offset_id;
    unsigned int tree_id;
    unsigned int size;
    Task *root;
    Task *taskMaxMakespan;
    Task *taskMaxMemRequirement;
    static Tree *originalTree;


public:

    Tree() {
        root_count = 0;
        offset_id = 0;
        tree_id = 1;
        tasks = new vector<Task *>();
        size = 0;
        taskMaxMakespan = NULL;
        taskMaxMemRequirement = NULL;
    }

    Tree(int N, int *prnts, double *nwghts, double *ewghts, double *mswghts) {
        root_count = 0;
        offset_id = 0;
        tree_id = 1;
        tasks = new vector<Task *>();
        taskMaxMakespan = NULL;
        taskMaxMemRequirement = NULL;
        Task *cur_node;

        for (int i = 1; i < N + 1; i++) {
            if (prnts[i] > 0) {
                cur_node = new Task(prnts[i], nwghts[i], ewghts[i], mswghts[i]);
                cur_node->setLabel(i);
                this->addTask(cur_node);
            } else {
                cur_node = new Task(nwghts[i], ewghts[i], mswghts[i], true);
                this->addRoot(cur_node);
                this->setTreeId(i);
            }
            cur_node->setId(i);
        }

        size = getTasks()->size();
    }

    Tree(vector<Task *> *nodes, Task *root, Tree *originalTree) {
        root_count = 1;
        offset_id = 0;
        tree_id = 1;
        this->root = root;
        root->toggleRootStatus(true);
        setRootId(root->getId());
        this->tasks = nodes;
        this->size = nodes->size();
        this->originalTree = originalTree;
        taskMaxMakespan = NULL;
        taskMaxMemRequirement = NULL;
    }


    ~Tree() {
        // cout<<"get root count "<<root_count<<endl;
        if (this->getRootId() != 0 && tasks->size() > 0) {
            delete getRoot();
        }

        delete tasks;
        delete taskMaxMakespan;
        delete taskMaxMemRequirement;
    }


    void Print(ostream &out) const {
        out << tasks->size() << endl;

        for (vector<Task *>::iterator iter = tasks->begin(); iter != tasks->end(); iter++) {
            out << (*iter)->getId() << ", parent: " << max((unsigned int) 0, (*iter)->getParentId()/*+1-offset_id*/)
                << " " << (*iter)->getNodeWeight()
                << " "
                << (*iter)->getEdgeWeight() << "; " << endl;
            out << "children: " << (*iter)->getChildren()->size() << endl;
            for (Task *child: *(*iter)->getChildren()) {
                out << "\t" << child->getId() << endl;
            }

        }
        cout << endl;
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
        size = new_node_count;
    }

    void reverseVector() {
        reverse(tasks->begin(), tasks->end());
    }

    void addTask(Task *newNode) {
        // cout << "add task " << newNode->getId() << " with parent " << newNode->getParentId() << "and children ";
        //  for(Task* child: *newNode->getChildren()){
        //     cout<<child->getId()<<" ";
        //  }
        //  cout<<endl;
        tasks->push_back(newNode);
        size++;
    }

    void addRoot(Task *newNode) {
        root_count++;
        assert(root_count == 1);
        assert(newNode->isRoot());
        tasks->push_back(newNode);
        this->root = newNode;
        size++;
    }

    Task *getRoot() const {
        assert(root_count == 1);
        return this->root;
    }

    //TODO: hier exception
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
        Task *task = tasks->at(node_id == 0 ? 0 : node_id - 1);
        if (task->getId() == node_id) {
            return task;
        } else {
            for (Task *taskSequential: *this->getTasks()) {
                if (taskSequential->getId() == node_id) {
                    return taskSequential;
                }
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

    void setTasks(vector<Task *> *c) {
        this->tasks = c;
    }

    static void setOriginalTree(Tree *origTree) {
        Tree::originalTree = origTree;
    }

    static Tree *getOriginalTree() {
        return Tree::originalTree;
    }

    int getSize() const {
        return size;
    }

    vector<Task *> getBrokenTasks() {
        vector<Task *> broken;

        for (Task *task: *getTasks()) {
            if (task->isBroken()) broken.push_back(task);
        }
        return broken;
    }

    Task *getTaskMaxMakespan() const {
        return taskMaxMakespan;
    }

    void setTaskMaxMakespan(Task *taskMaxMakespan) {
        Tree::taskMaxMakespan = taskMaxMakespan;
    }

    Task *getTaskMaxMemRequirement() const {
        return taskMaxMemRequirement;
    }

    void setTaskMaxMemRequirement(Task *taskMaxMemRequirement) {
        Tree::taskMaxMemRequirement = taskMaxMemRequirement;
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

    unsigned long countBrokenEdges() {
        unsigned long treeSize = this->getTasks()->size();
        unsigned long counter = 0;
        for (unsigned int i = treeSize; i >= 1; --i) {
            Task *currentnode = this->getTask(i);
            if (currentnode->isBroken()) {
                counter++;
            }
        }
        return counter;

    }

    Tree *BuildQtree();

    Tree *BuildQtreeOld();

    unsigned int HowmanySubtrees(bool quiet);

    bool
    MemoryEnough(Task *Qrootone, Task *Qroottwo, bool leaf, double available_memory_size, double &requiredMemorySize);

    double
    CheckRequiredMemoryIfMerged(Task *Qrootone, Task *Qroottwo, bool leaf);

    bool MemoryEnoughOld(Tree *tree, Task *Qrootone, Task *Qroottwo, bool leaf, double memory_size, int *chstart,
                         int *children);

    double ImprovedSplit();

    double ASAP();

    double Merge(bool CheckMemory);

    double
    MergeOld(unsigned int num_subtrees, unsigned int processor_number, double const memory_size, bool CheckMemory);

    double
    MergeV2(unsigned int num_subtrees, unsigned int processor_number, double const memory_size, bool CheckMemory);

    double SplitAgain();

    double SplitAgainV2(unsigned int processor_number, unsigned int num_subtrees, std::map<int, int> &taskToPrc,
                        std::map<int, bool> &isProcBusy);

    vector<Task *> buildCriticalPath(Tree *Qtree);

    int numberOfLeaves();

    double avgNodeWeight();

    double avgEdgeWeight();

    double avgMSWeight();

    void clearComputedValues();

    void cleanAssignedAndReassignFeasible() {
        for (Task *task: *tasks) {
            task->assignFeasibleProcessorsToSubtree(task->getMinMemUnderlying());
            freeProcessorIfAvailable(task);
        }
        for (Processor *item: Cluster::getFixedCluster()->getProcessors()) {
            if (item->isBusy) {
                item->isBusy = false;
                item->setAssignedTaskId(-1);
                item->setAssignedTask(NULL);
                item->setOccupiedMemorySize(0);
            }
        }
        assert(Cluster::getFixedCluster()->getNumberProcessors() ==
               Cluster::getFixedCluster()->getNumberFreeProcessors());
    }

};

double SplitAgainOld(Tree *tree, unsigned int processor_number, unsigned int num_subtrees);

typedef map<unsigned int, double> io_map;
typedef pair<unsigned int, unsigned int> node_sche;
typedef pair<unsigned int, double> node_ew;

Tree *read_tree(const char *filename);

extern "C"
{
void po_construct(const int N, const int *prnts, int **chstart, int **chend, int **children, int *root);
} /* closing brace for extern "C" */

double MaxOutDegree(Tree *tree, int quiet);

vector<double> maxAndAvgFanout(Tree *tree);

int maxDepth(Task *root);

double IOCounter(Tree *subtree, int *schedule,
                 bool divisible, int quiet, unsigned int &com_freq,
                 vector<unsigned int> *brokenEdges, io_method_t method);

double
IOCounterWithVariableMem(Tree *tree, int *schedule,
                         Cluster *cluster, bool divisible, int quiet, unsigned int &com_freq,
                         vector<unsigned int> *brokenEdges, io_method_t method);

Tree *BuildSubtree(Tree *tree, Task *subtreeRoot);

Tree *BuildSubtreeOld(Tree *tree, Task *SubtreeRoot, unsigned int new_tree_size, int **prnts, double **ewghts,
                      double **timewghts, double **spacewghts, int *chstart, int *children);

void popSmallestRootsToFitToCluster(list<Task *> &parallelRoots, unsigned long amountSubtrees, int limit);

void breakPreparedEdges(Task *root, list<Task *> &parallelRoots);

list<Task *> buildParallelRootsFromSequentialSet(Task *root, list<Task *> &sequentialSet);

double getWeightPQ(list<Task *> &parallelRoots, Task *currentNode);

double getWeightSurplusFromSmallestNodes(list<Task *> &parallelRoots, int limit);

int *
copyScheduleBackwards(schedule_traversal *schedule_f);

class SeqSet {
protected:
    vector<Task *> seqSet;
    vector<Task *> parallelRoots;
    double makespan;
    int numberSteps;

public:
    SeqSet(Tree *tree, double makespan) {
        Task *root = tree->getRoot();
        this->seqSet = root->tasksInSubtreeRootedHere();
        this->parallelRoots = tree->getBrokenTasks();
        this->makespan = makespan;
        numberSteps = 0;
    }

    SeqSet(Tree *tree, double makespan, int steps) {
        Task *root = tree->getRoot();
        this->seqSet = root->tasksInSubtreeRootedHere();
        this->parallelRoots = tree->getBrokenTasks();
        this->makespan = makespan;
        numberSteps = steps;
    }

    string print() {
        string result = "";
        result += "Size of SeqSet: " + to_string(seqSet.size()) + ",\t Number of Subtrees: " +
                  to_string(parallelRoots.size())
                  + ",\t Makespan: " + to_string(makespan) + ",\t Steps: " + to_string(numberSteps);
        return result;
    }

    string printDetailed() {
        string result = "";
        result += "Size of SeqSet: " + to_string(seqSet.size()) + ",\t Number of Subtrees: " +
                  to_string(parallelRoots.size())
                  + ",\t Makespan: " + to_string(makespan);
        throw "not implemented!";
        return result;
    }

    void implementSeqSet(Tree *tree) {
        for (Task *task: *tree->getTasks()) {
            task->restoreEdge();
        }
        for (Task *task: parallelRoots) {
            task->breakEdge();
        }

    }
};

#endif
#endif
