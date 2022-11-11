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
#include <set>

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

extern Tree *veryOriginalTree;

struct FastestProcessor {
    bool operator()(Processor *a, const Processor *b) const {
        return a->getProcessorSpeed() >= b->getProcessorSpeed();
    }
};

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
    set<Processor *, FastestProcessor> *feasibleProcessors;
    double minMemUnderlying;
    int level;
    double tMax;

public :
    unsigned int Ci;
    double Mpeak;
    double Mavail;
    bool needsRecomputeMemReq;

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
        feasibleProcessors = new set<Processor *, FastestProcessor>();
        minMemUnderlying = 0;
        tMax = 0;
        level = -1;
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
        feasibleProcessors = new set<Processor *, FastestProcessor>();
        minMemUnderlying = 0;
        tMax = 0;
        level = -1;
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
        feasibleProcessors = new set<Processor *, FastestProcessor>();
        tMax = 0;
        minMemUnderlying = 0;
        level = -1;
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
        feasibleProcessors = new set<Processor *, FastestProcessor>();
        for (Processor *p: *otherTask.feasibleProcessors) {
            feasibleProcessors->insert(p);
        }
        minMemUnderlying = otherTask.minMemUnderlying;
        level = -1;
    }

    ~Task() {
        this->feasibleProcessors->clear();
        delete feasibleProcessors;
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
        if (this->assignedProcessor == NULL) return -1;
        else return this->assignedProcessor->getProcessorSpeed();
    }

    void toggleRootStatus(bool newStatus) {
        if (this->isRoot() != newStatus) {
            __root = !__root;
        } else {
            throw std::runtime_error("Trying to set root to a position that it is alrady in");
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

    void setCost(double costToSet) {
        cost = costToSet;
        cost_computed = true;
    }

    void setCostComputed(bool compute) {
        this->cost_computed = compute;
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

    set<Processor *, FastestProcessor> *getFeasibleProcessors() {
        return this->feasibleProcessors;
    }

    void setFeasibleProcessors(set<Processor *, FastestProcessor> *v) {
        this->feasibleProcessors->clear();
        for (const auto &item: *v) {
            this->feasibleProcessors->insert(item);
        }

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
        this->feasibleProcessors->insert(proc);
        //throw "TODO";
    }

    void setMinMemUnderlying(double minMem) {
        this->minMemUnderlying = minMem;
    }

    double getMinMemUnderlying() {
        return this->minMemUnderlying;
    }

    int getLevel() {
        return this->level;
    }

    void setLevel(int newLevel) {
        this->level = newLevel;
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
        if (assignedProcSpeed == -1)
            throw runtime_error("No processor assigned to " + to_string(this->getId()));
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

        // cout << "on " << this->getId() << "MS P " << MS_parallelPart << " seq " << MS_sequentialPart <<

        //      " ew " << edge_weight << " BW " << Cluster::getFixedCluster()->getHomogeneousBandwidth() << endl;
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
        //cout << "on " << this->getId() << "MS P " << MS_parallelPart << " seq " << MS_sequentialPart <<

        //     " ew " << edge_weight << " BW " << Cluster::getFixedCluster()->getHomogeneousBandwidth() << endl;
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
            throw std::runtime_error("No Schedule Possible");
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
        //sort(feasibleProcessors->begin(), feasibleProcessors->end(),
        //     [](Processor *a, Processor *b) { return a->getProcessorSpeed() > b->getProcessorSpeed(); });
        assert((*feasibleProcessors->begin())->getProcessorSpeed() >=
               (*feasibleProcessors->end().operator--())->getProcessorSpeed());
        return *feasibleProcessors->begin();
    }

    void deleteFeasible(Processor *proc) {
        auto position_it = find(feasibleProcessors->begin(), feasibleProcessors->end(), proc);
        if (position_it != feasibleProcessors->end()) {
            feasibleProcessors->erase(position_it);
        }
    }

    void assignFeasibleProcessorsToSubtree(double minMem);

    double computeMinMemUnderlyingAndAssignFeasible(Tree *tree, bool greedy);

    vector<Task *> getTasksInSubtreeRootedHere() {
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
        std::sort(this->getChildren()->begin(), this->getChildren()->end(), [](Task *a, Task *b) {
            return a->getMakespanCost(true, false) > b->getMakespanCost(true, false);
        });

        vector<Task *> newlyBroken;
        int border = n < this->getChildren()->size() ? n : this->getChildren()->size();

        for (int i = 0; i < border; i++) {
            if (!this->getChildren()->at(i)->isBroken()) {
                (this->getChildren()->at(i))->breakEdge();
                newlyBroken.push_back(this->getChildren()->at(i));
            }
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
    int numberTasksWMinMem;
    int deepestLevel;

    Tree() {
        root_count = 0;
        offset_id = 0;
        tree_id = 1;
        tasks = new vector<Task *>();
        size = 0;
        taskMaxMakespan = NULL;
        taskMaxMemRequirement = NULL;
        numberTasksWMinMem = 0;
        deepestLevel = -1;
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
        numberTasksWMinMem = 0;
        deepestLevel = -1;
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
        Tree::originalTree = originalTree;
        taskMaxMakespan = NULL;
        taskMaxMemRequirement = NULL;
        numberTasksWMinMem = 0;
        deepestLevel = -1;
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


    void Print(ostream &out, int border = 1) const {
        out << tasks->size() << endl;
        int maxOut = 0;
        for (vector<Task *>::iterator iter = tasks->begin(); iter != tasks->end(); iter++) {
            maxOut++;
            if (maxOut > tasks->size() / border) break;
            out << (*iter)->getId() << ", parent: " << max((unsigned int) 0, (*iter)->getParentId()/*+1-offset_id*/)
                << " " << (*iter)->getNodeWeight()
                << " "
                << (*iter)->getEdgeWeight() << "; " << endl;
            out<<"Assigned to: "<<(*iter)->getAssignedProcessorSpeed()<<endl;
          //  out << "children: " << (*iter)->getChildren()->size() << endl;
           // for (Task *child: *(*iter)->getChildren()) {
          //      out << "\t" << child->getId() << endl;
          //  }

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

    void removeTask(int id) {
        auto iterator = std::find_if(tasks->begin(), tasks->end(),
                                     [id](Task *task) { return task->getId() == id; });
        tasks->erase(iterator);

        iterator = std::find_if(tasks->begin(), tasks->end(), [id](Task *task) { return task->getId() == id; });
        assert(iterator ==
               tasks->end());
        size--;
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
        if (node_id < tasks->size()) {
            Task *task = tasks->at(node_id == 0 ? 0 : node_id - 1);
            if (task->getId() == node_id) {
                return task;
            }
        }

        for (Task *taskSequential: *this->getTasks()) {
            if (taskSequential->getId() == node_id) {
                return taskSequential;
            }
        }
        throw runtime_error("Task not found for id " + to_string(node_id));
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

    Tree *BuildQtree(bool sumMakespans =false);

    Tree *BuildQtreeOld();

    unsigned int HowmanySubtrees(bool quiet);

    unsigned int HowmanySubtreesAndWeights(bool quiet);

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
        for (Task *task: getBrokenTasks()) {
            task->getFeasibleProcessors()->clear();
            task->assignFeasibleProcessorsToSubtree(task->getMinMemUnderlying());
            freeProcessorIfAvailable(task);
        }
        if (Cluster::getFixedCluster()->getNumberProcessors() !=
            Cluster::getFixedCluster()->getNumberFreeProcessors()) {
            Cluster::getFixedCluster()->freeAllBusyProcessors();
        }

    }

    void reassignRootProcessorToSubtree(Task *subtreeRoot) {
        Processor *processorOfRoot = subtreeRoot->getAssignedProcessor();
        vector<Task *> tasksUnderRoot = subtreeRoot->getTasksInSubtreeRootedHere();
        for (Task *task: tasksUnderRoot) {
            task->setAssignedProcessor(processorOfRoot);
        }
       // processorOfRoot->assignTask(subtreeRoot);
    }

    Task* findNextBrokenParent(Task* child) {
        Task *parent = child->getParent();
        while (parent != nullptr && !parent->isBroken()) {
            parent = parent->getParent();
        }
        if (parent != nullptr) {
            if (!parent->isBroken()) {
                parent = parent->getParent();
            }
            assert(parent->isBroken() == true);
            return parent;
        }
        return parent;
    }

    void mergeTaskToOnlyChild(Task *mergeRoot);

    void mergeTaskToAllChildren(Task *mergeRoot);

    void mergeLinearChains();

    void renumberAllTasks();

    void levelsToTasks();

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
        this->seqSet = root->getTasksInSubtreeRootedHere();
        this->parallelRoots = tree->getBrokenTasks();
        this->makespan = makespan;
        numberSteps = 0;
    }

    SeqSet(Tree *tree, double makespan, int steps) {
        Task *root = tree->getRoot();
        this->seqSet = root->getTasksInSubtreeRootedHere();
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
