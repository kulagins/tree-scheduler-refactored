//
//  heuristics.cpp
//  memCom2
//
//  Created by changjiang GOU on 11/05/2018.
//  Copyright © 2018 ROMA. All rights reserved.
//

#include <stdio.h>
#include <time.h>
#include <list>
#include <algorithm>
#include <tuple>
#include "../include/heuristics.h"
#include "../include/lib-io-tree-free-methods.h"
#include "../include/cluster.h"
//INFO: paths are necessary for Clion. In case of problems please ping Svetlana.
#include "../include/lib-io-tree-minmem.h"
#include <bits/stdc++.h>
/*#include "../include/inputParser.h" */

//#include <omp.h>
bool cmp_noincreasing(Task *a, Task *b) {
    return (a->getMakespanCost(true, false) >= b->getMakespanCost(true, false));
};

bool cmp_nodecreasing(Task *a, Task *b) { return (a->getMakespanCost(true, false) < b->getMakespanCost(true, false)); };

bool cmp_noIn_noCommu(Task *a, Task *b) {
    return (a->getMakespanCost(false, false) >= b->getMakespanCost(false, false));
};

bool cmp_Mem_nodecreasing(Task *a, Task *b) {
    return (a->getNodeWeight() >= b->getNodeWeight());
};


struct CompareMapEntries {
    int val;

    CompareMapEntries(const int &val) : val(val) {}
};

bool operator==(const std::pair<int, int> &p, const CompareMapEntries &c) {
    return c.val == p.second;
}

bool operator==(const CompareMapEntries &c, const std::pair<int, int> &p) {
    return c.val == p.second;
}

template<class T, class U>
void getTwoLargestElementTypetwo(T container, U &Largest, U &secondLargest) {
    if (container->front()->getMakespanCost(true, false) > container->at(1)->getMakespanCost(true, false)) {
        Largest = container->front();
        secondLargest = container->at(1);
    } else {
        Largest = container->at(1);
        secondLargest = container->front();
    }

    if (container->size() > 2) {
        vector<Task *>::iterator iter = container->begin();
        iter = iter + 2;
        for (; iter != container->end(); ++iter) {
            if ((*iter)->getMakespanCost(true, false) > Largest->getMakespanCost(true, false)) {
                secondLargest = Largest;
                Largest = *iter;
            } else if ((*iter)->getMakespanCost(true, false) > secondLargest->getMakespanCost(true, false)) {
                secondLargest = *iter;
            }
        }
    }
}

void getTwoLargestElementTypethree(vector<Task *> *container, vector<Task *>::iterator &Largest,
                                   vector<Task *>::iterator &secondLargest) {
    if (container->front()->getMakespanCost(false, false) >= container->back()->getMakespanCost(false, false)) {
        Largest = container->begin();
        secondLargest = Largest;
        advance(secondLargest, 1);
    } else {
        secondLargest = container->begin();
        Largest = secondLargest;
        advance(Largest, 1);
    }

    if (container->size() > 2) {
        vector<Task *>::iterator iter = container->begin();
        advance(iter, 2);
        for (; iter != container->end(); ++iter) {
            if ((*iter)->getMakespanCost(false, false) > (*Largest)->getMakespanCost(false, false)) {
                secondLargest = Largest;
                Largest = iter;
            } else if ((*iter)->getMakespanCost(false, false) > (*secondLargest)->getMakespanCost(false, false)) {
                secondLargest = iter;
            }
        }
    }
}

void getTwoSmallestElement(list<Task *> *container, list<Task *>::iterator &Smallest,
                           list<Task *>::iterator &secondSmallest) {
    if (container->front()->getMakespanWeight() <= container->back()->getMakespanWeight()) {
        Smallest = container->begin();
        secondSmallest = Smallest;
        advance(secondSmallest, 1);
    } else {
        secondSmallest = container->begin();
        Smallest = secondSmallest;
        advance(Smallest, 1);
    }

    if (container->size() > 2) {
        list<Task *>::iterator iter = container->begin();
        advance(iter, 2);
        for (; iter != container->end(); ++iter) {
            if ((*iter)->getMakespanWeight() < (*Smallest)->getMakespanWeight()) {
                secondSmallest = Smallest;
                Smallest = iter;
            } else if ((*iter)->getMakespanWeight() < (*secondSmallest)->getMakespanWeight()) {
                secondSmallest = iter;
            }
        }
    }
}


bool MiWi_nodecreasing(Task *a, Task *b) {
    return (a->getNodeWeight() * a->getMakespanWeight() < b->getNodeWeight() * b->getMakespanWeight());
};


double Task::SplitSubtrees(bool twolevel, list<Task *> &parallelRoots, unsigned long &sequentialLength, int limit) {
    parallelRoots.clear();
    parallelRoots.emplace_front(this);

    //cost from unsplitted root
    vector<double> makespansOfSplittings{this->getMakespanCost(true, true)}; // take communication cost into account

    double MS_sequential = this->getEdgeWeight(), weightSurplusFromSmallestTasks, weightsTasksPriorityQueue;
    if (Cluster::getFixedCluster()->isBandwidthHomogeneous()) {
        MS_sequential /= Cluster::getFixedCluster()->getHomogeneousBandwidth();
    }

    Task *currentNode = this;

    while (!currentNode->isLeaf()) {
        MS_sequential = MS_sequential + currentNode->getMakespanWeight();
        //hier
        weightsTasksPriorityQueue = getWeightPQ(parallelRoots, currentNode);
        if (weightsTasksPriorityQueue == -1) break;
        else {
            currentNode = *max_element(parallelRoots.begin(), parallelRoots.end(), cmp_nodecreasing); //non-decreasing
        }

        weightSurplusFromSmallestTasks = getWeightSurplusFromSmallestNodes(parallelRoots, limit);
        //cout<<"makespan "<<MS_sequential+weightSurplusFromSmallestTasks+weightsTasksPriorityQueue<<endl;
        //  cout<<MS_sequential<<" "<<weightSurplusFromSmallestTasks<<" "<<weightsTasksPriorityQueue<<"; ";
        makespansOfSplittings.push_back(MS_sequential + weightSurplusFromSmallestTasks + weightsTasksPriorityQueue);
    }
    //  cout<<"MS"<<endl;
    //  for(double d: makespansOfSplittings){
    //     cout<<d<<" ";
    //  }
    //  cout<<endl;

    if (twolevel) {
        return *std::min_element(makespansOfSplittings.begin(), makespansOfSplittings.end());
    }

    //return broken edges, i.e., root of subtrees
    auto smallestMS_iter = min_element(makespansOfSplittings.begin(), makespansOfSplittings.end());
    sequentialLength = smallestMS_iter - makespansOfSplittings.begin();
    parallelRoots = fillParallelRootsUntilBestMakespan(makespansOfSplittings, sequentialLength);

    unsigned long amountSubtrees;
    if (parallelRoots.size() > 1) {
        amountSubtrees = parallelRoots.size() + 1;
    } else {
        amountSubtrees = 1;
    }

    popSmallestRootsToFitToCluster(parallelRoots, amountSubtrees, limit);
    breakPreparedEdges(this, parallelRoots);
    //    cout<<"makespan from the tree root "<<root->getMakespanCost(true,true)<<endl;

    return *smallestMS_iter;
}

list<Task *>
Task::fillParallelRootsUntilBestMakespan(vector<double> &makespansOfSplittings,
                                         unsigned long stepsUntilMinimalMakespan) const {

    vector<Task *> *children;
    unsigned int i = 0;
    list<Task *> parallelRoots = {const_cast<Task *>(this)};
    Task *currentNode = const_cast<Task *>(this);

    while (i < stepsUntilMinimalMakespan) {
        parallelRoots.remove(currentNode);

        children = currentNode->getChildren();
        for (auto &iter: *children) {
            if (!iter->isBroken()) {
                parallelRoots.push_back(iter);
            }
        }

        currentNode = *max_element(parallelRoots.begin(), parallelRoots.end(), cmp_nodecreasing); //non-decreasing
        i++;
    }
    return parallelRoots;
}

double getWeightSurplusFromSmallestNodes(list<Task *> &parallelRoots, int limit) {
    int limitToPartitioning = limit == -1 ? Cluster::getFixedCluster()->getNumberProcessors() : limit;
    double weightSurplusFromSmallestNodes = 0;
    unsigned long surplusOfSubtreesOverProcessors = 0;
    unsigned long amountSubtrees = parallelRoots.size() + 1;

    if (amountSubtrees > limitToPartitioning) {
        parallelRoots.sort(cmp_noIn_noCommu); //non-increasing sort, computation weight, no communication
        auto iter = parallelRoots.rbegin();
        surplusOfSubtreesOverProcessors = amountSubtrees - limitToPartitioning;
        for (unsigned int i = 0; i < surplusOfSubtreesOverProcessors; ++i, ++iter) {
            weightSurplusFromSmallestNodes += (*iter)->getMakespanCost(false,
                                                                       false); // no comunication cost, ImprovedSplit never goes to here.
        }
    }
    return weightSurplusFromSmallestNodes;
}

double getWeightPQ(list<Task *> &parallelRoots, Task *currentNode) {
    double Weight_PQ = 0;
    double temp;
    parallelRoots.remove(currentNode);
    //cout<<"pop up "<<currentNode->getId()<<endl;

    vector<Task *> *children = currentNode->getChildren();
    for (vector<Task *>::iterator iter = children->begin(); iter != children->end(); iter++) {
        if ((*iter)->isBroken()) {
            temp = (*iter)->getMakespanCost(true, false);
            if (temp > Weight_PQ) {
                Weight_PQ = temp;
            }
        } else {
            parallelRoots.push_back(*iter);
            //cout<<"   insert "<<(*iter)->getId()<<endl;
        }
    }

    if (parallelRoots.empty()) {
        return -1;
    } else {
        currentNode = *max_element(parallelRoots.begin(), parallelRoots.end(), cmp_nodecreasing);//non-decreasing
    }
    temp = currentNode->getMakespanCost(true, false);
    if (temp > Weight_PQ) {
        Weight_PQ = temp;
    }
    return Weight_PQ;
}

void breakPreparedEdges(Task *root, list<Task *> &parallelRoots) {
    root->breakEdge(); //root should always be broken
    for (auto iter = parallelRoots.begin(); iter != parallelRoots.end(); ++iter) {
        (*iter)->breakEdge();
    }
}

//TODO: n^2, think of improving
list<Task *> buildParallelRootsFromSequentialSet(Task *root, list<Task *> &sequentialSet) {
    list<Task *> parallelRoots;
    for (auto iter = sequentialSet.begin(); iter != sequentialSet.end(); ++iter) {
        for (Task *child: *(*iter)->getChildren()) {
            auto it = std::find(sequentialSet.begin(), sequentialSet.end(), child);
            //If the child is not itself in the sequential set, then it is the root of a subtree
            if (it == sequentialSet.end()) {
                parallelRoots.push_back(child);
            }
        }

    }
    return parallelRoots;
}

void popSmallestRootsToFitToCluster(list<Task *> &parallelRoots, unsigned long amountSubtrees, int limit) {
    int limitToPartitioning = limit == -1 ? Cluster::getFixedCluster()->getNumberProcessors() : limit;
    if (amountSubtrees > limitToPartitioning) {
        parallelRoots.sort(cmp_noIn_noCommu); //non-increasing sort, computation weight, no communication cost
        unsigned int surplus = amountSubtrees - limitToPartitioning;
        for (unsigned int i = 0; i < surplus; ++i) {
            parallelRoots.pop_back();
        }
    }
}

void ISCore(Task *root, unsigned long num_processors,
            bool sequentialPart) { //number of processors here assumed to the same as tree'size
    list<Task *> parallelRoots;
    double MS_before;
    double MS_now;
    unsigned long SF_now; //avoid dead lock

    if (root->isLeaf()) {
        //cout<<"root is leaf, return."<<endl;
        return;
    }

    root->SplitSubtrees(false, parallelRoots,
                        SF_now,
                        -1); //SF_now will be modified in SplitSubtrees, it represents the length of sequential part, 0 means the subtree no need to partition

    if (sequentialPart == true) {
        if (SF_now == 0) {
            //cout<<"this subtree has already been fully checked, return."<<endl;
            return;
        }
    }

    parallelRoots.sort(cmp_noincreasing); //non-increasing sort, communication counted

    Task *frontNode;
    if (parallelRoots.size() > 1) { //==1 means there is no parallel part
        while (true) {
            frontNode = parallelRoots.front();
            parallelRoots.pop_front();
            MS_before = frontNode->getMakespanCost(true, false);

            //cout<<"---ISCore works on Parallel root "<<frontNode->getId()<<endl;
            ISCore(frontNode, num_processors, false);

            MS_now = frontNode->getMakespanCost(true, true); //Makespan updated enforced

            if (MS_now >= MS_before) {
                break;
            }

            if (parallelRoots.empty()) {
                break;
            }

            if (parallelRoots.front()->getMakespanCost(true, false) <= MS_now) {
                break;
            }
        }

        //cout<<"---ISCore works on Sequential root "<<root->getId()<<endl;
        ISCore(root, num_processors, true);
    }

    return;
}

double Tree::ImprovedSplit() {
    unsigned long tree_size = this->getTasks()->size();
    Task *root = this->getRoot();
    //cout<<"---ISCore works on the root"<<endl;
    ISCore(root, tree_size, false);


    double makespan = this->Merge(false);
    return makespan;
}

bool cmp_merge_smallest(const pair<double, Task *> &a, const pair<double, Task *> &b) { return a.first < b.first; };

bool estimateMS(Tree *tree, Tree *Qtree, Task *&smallestNode, Processor *&processorToMergeTo, bool CheckMemory) {
    //cout<<"   ---start compute the minimum combination"<<endl;

    Task *currentQNode;
    double increase;
    bool memoryEnough;
    Task *LargestNode;
    Task *secondLargest;
    bool leaf = false;
    const vector<Task *> *subtrees = Qtree->getTasks();
    vector<Task *> *children;
    double homogeneousBandwidth;
    if (Cluster::getFixedCluster()->isBandwidthHomogeneous()) {
        homogeneousBandwidth = Cluster::getFixedCluster()->getHomogeneousBandwidth();
    } else throw "Cluster not homogeneous bandwidth";

    if (subtrees->front()->getId() != 1) { //the root is supposed to be the first element in vector nodes
        cout << "error in function estimateMS" << endl;
        return false;
    }

    vector<Task *> tempQue;
    currentQNode = Qtree->getRoot();
    currentQNode->setMakespanDiff(0);
    children = currentQNode->getChildren();
    tempQue.insert(tempQue.end(), children->begin(), children->end());
    //cout<<"   ---compute makespan difference---"<<endl;
    while (!tempQue.empty()) {
        currentQNode = tempQue.back();
        tempQue.pop_back();
        currentQNode->setMakespanDiff(
                currentQNode->getParent()->getMakespanDiff() + currentQNode->getParent()->getParallelPart() -
                currentQNode->getMakespanCost(true, false));
        //cout<<"   subtree "<<currentQNode->getOtherSideId()<<", makespan difference: "<<currentQNode->getMakespanDiff()<<endl;
        children = currentQNode->getChildren();
        tempQue.insert(tempQue.end(), children->begin(), children->end());
    }
    //cout<<"   ----------------------------------"<<endl;

    list<pair<double, Task *>> list_increase_id;
    vector<Task *>::const_iterator iter = subtrees->begin();
    ++iter;
    unsigned long size = subtrees->size() - 1;
    //  #pragma omp parallel for
    for (unsigned int step = 0; step < size; ++step) {
        currentQNode = *(iter + step);

        if (tree->getTask(currentQNode->getOtherSideId())->isBroken()) { //this subtree has not been merged yet
            children = currentQNode->getChildren();
            if (children->empty()) { //this is a leaf node
                children = currentQNode->getParent()->getChildren();
                if (children->size() == 2) {
                    increase = children->front()->getMakespanWeight() + children->back()->getMakespanWeight() -
                               currentQNode->getParent()->getParallelPart();
                } else if (children->size() == 1) {
                    increase = -currentQNode->getEdgeWeight() / homogeneousBandwidth;
                } else {

                    getTwoLargestElementTypetwo(children, LargestNode,
                                                secondLargest); //no-increasing, communication counted

                    if (currentQNode->getId() == LargestNode->getId()) {
                        increase = currentQNode->getMakespanWeight() + secondLargest->getMakespanCost(true, false) -
                                   currentQNode->getMakespanCost(true, false);
                    } else {
                        increase = currentQNode->getMakespanWeight();
                    }
                }
            } else { //not a leaf node
                children = currentQNode->getParent()->getChildren();

                if (children->size() == 1) {
                    increase = -currentQNode->getEdgeWeight() / homogeneousBandwidth;
                } else {
                    getTwoLargestElementTypetwo(children, LargestNode,
                                                secondLargest); //no-increasing, communication counted

                    if (currentQNode->getId() == LargestNode->getId()) {
                        increase = currentQNode->getMakespanWeight() +
                                   max(secondLargest->getMakespanCost(true, false),
                                       currentQNode->getParallelPart()) -
                                   currentQNode->getMakespanCost(true, false);
                    } else {
                        increase = currentQNode->getMakespanWeight();
                    }
                }
            }

            //now consider the "slack" between the parent of current node with siblings of its parent
            increase = increase - currentQNode->getParent()->getMakespanDiff();

            //cout<<"merge, increase in MS(r) "<<increase<<endl;
            list_increase_id.push_back(pair<double, Task *>(increase, currentQNode));
        }
    }

    bool feasible = false;
    list<pair<double, Task *>>::iterator smallest_iter;
    while (feasible == false && list_increase_id.empty() == false) {
        smallest_iter = min_element(list_increase_id.begin(), list_increase_id.end(), cmp_merge_smallest);
        currentQNode = (*smallest_iter).second;
        // cout << "   increase in MS(r) estimated: " << (*smallest_iter).first << endl;

        children = currentQNode->getChildren();
        if (children->empty()) {
            leaf = true;
        }
        double requiredMemory = 0;
        requiredMemory = tree->CheckRequiredMemoryIfMerged(currentQNode->getParent(), currentQNode, leaf);
        //cout << "required memory " << requiredMemory << endl;

        if (Cluster::getFixedCluster()->areAllUnassigned(currentQNode, tree) &&
            Cluster::getFixedCluster()->hasFreeProcessor()) {
            processorToMergeTo = Cluster::getFixedCluster()->smallestFreeProcessorFitting(requiredMemory);
        } else
            processorToMergeTo = Cluster::getFixedCluster()->findSmallestFittingProcessorForMerge(currentQNode,
                                                                                                  tree,
                                                                                                  requiredMemory);

        if (CheckMemory) {
            memoryEnough = processorToMergeTo != nullptr;
            //memoryEnough =MemoryEnough(currentQNode->getParent(), currentQNode, leaf, memorySize,      requiredMemory);
        } else {
            memoryEnough = true;
        }

        if (memoryEnough) {
            feasible = true;
            smallestNode = currentQNode;
            if (processorToMergeTo != nullptr) {
                processorToMergeTo->setOccupiedMemorySize(requiredMemory);
            }
        } else {
            list_increase_id.erase(smallest_iter);
        }
    }

    //cout<<"   ---end compute the minimum combination"<<endl;
    //  cout<<"processor"<<endl;
    //  cout<<processorToMergeTo->getMemorySize()<<endl;
    return feasible;
}

bool estimateMSOld(Tree *tree, Tree *Qtree, Task *&smallestNode, double memory_size, bool CheckMemory) {
    //cout<<"   ---start compute the minimum combination"<<endl;

    Task *currentQNode;
    double increase;
    bool memoryEnough;
    Task *LargestNode;
    Task *secondLargest;
    bool leaf = false;
    const vector<Task *> *subtrees = Qtree->getTasks();
    vector<Task *> *children;

    if (subtrees->front()->getId() != 1) {//the root is supposed to be the first element in vector nodes
        cout << "error in function estimateMS" << endl;
        return false;
    }

    vector<Task *> tempQue;
    currentQNode = Qtree->getRoot();
    currentQNode->setMakespanDiff(0);
    children = currentQNode->getChildren();
    tempQue.insert(tempQue.end(), children->begin(), children->end());
    //cout<<"   ---compute makespan difference---"<<endl;
    while (!tempQue.empty()) {
        currentQNode = tempQue.back();
        tempQue.pop_back();
        currentQNode->setMakespanDiff(
                currentQNode->getParent()->getMakespanDiff() + currentQNode->getParent()->getParallelPart() -
                currentQNode->getMakespanCost(true, false));
        //cout<<"   subtree "<<currentQNode->getothersideID()<<", makespan difference: "<<currentQNode->getMakespanDiff()<<endl;
        children = currentQNode->getChildren();
        tempQue.insert(tempQue.end(), children->begin(), children->end());
    }
    //cout<<"   ----------------------------------"<<endl;


    list<pair<double, Task *>> list_increase_id;
    vector<Task *>::const_iterator iter = subtrees->begin();
    ++iter;
    unsigned long size = subtrees->size() - 1;
//  #pragma omp parallel for
    for (unsigned int step = 0; step < size; ++step) {
        currentQNode = *(iter + step);

        if (tree->getTask(currentQNode->getOtherSideId())->isBroken() == true) {//this subtree has not been merged yet
            children = currentQNode->getChildren();
            if (children->empty()) {//this is a leaf node
                children = currentQNode->getParent()->getChildren();
                if (children->size() == 2) {
                    increase = children->front()->getMakespanWeight() + children->back()->getMakespanWeight() -
                               currentQNode->getParent()->getParallelPart();
                } else if (children->size() == 1) {
                    increase = -currentQNode->getEdgeWeight() / Cluster::getFixedCluster()->getHomogeneousBandwidth();;
                } else {

                    getTwoLargestElementTypetwo(children, LargestNode,
                                                secondLargest);//no-increasing, communication counted

                    if (currentQNode->getId() == LargestNode->getId()) {
                        increase = currentQNode->getMakespanWeight() + secondLargest->getMakespanCost(true, false) -
                                   currentQNode->getMakespanCost(true, false);
                    } else {
                        increase = currentQNode->getMakespanWeight();
                    }
                }
            } else {//not a leaf node
                children = currentQNode->getParent()->getChildren();

                if (children->size() == 1) {
                    increase = -currentQNode->getEdgeWeight() / Cluster::getFixedCluster()->getHomogeneousBandwidth();;
                } else {
                    getTwoLargestElementTypetwo(children, LargestNode,
                                                secondLargest);//no-increasing, communication counted

                    if (currentQNode->getId() == LargestNode->getId()) {
                        increase = currentQNode->getMakespanWeight() +
                                   max(secondLargest->getMakespanCost(true, false), currentQNode->getParallelPart()) -
                                   currentQNode->getMakespanCost(true, false);
                    } else {
                        increase = currentQNode->getMakespanWeight();
                    }
                }
            }

            //now consider the "slack" between the parent of current node with siblings of its parent
            increase = increase - currentQNode->getParent()->getMakespanDiff();

            //cout<<"merge, increase in MS(r) "<<increase<<endl;
            list_increase_id.push_back(pair<double, Task *>(increase, currentQNode));
        }
    }

    bool feasible = false;
    list<pair<double, Task *>>::iterator smallest_iter;
    while (feasible == false && list_increase_id.empty() == false) {
        smallest_iter = min_element(list_increase_id.begin(), list_increase_id.end(), cmp_merge_smallest);
        currentQNode = (*smallest_iter).second;
        //cout << "   increase in MS(r) estimated: " << (*smallest_iter).first << endl;

        children = currentQNode->getChildren();
        if (children->empty()) {
            leaf = true;
        }

        if (CheckMemory == true) {
            double requiredMemorySize;
            memoryEnough = tree->MemoryEnough(currentQNode->getParent(), currentQNode, leaf, memory_size,
                                              requiredMemorySize);
        } else {
            memoryEnough = true;
        }

        if (memoryEnough == true) {
            feasible = true;
            smallestNode = currentQNode;
        } else {
            list_increase_id.erase(smallest_iter);
        }
    }

    //cout<<"   ---end compute the minimum combination"<<endl;
    return feasible;
}


void freeProcessorIfAvailable(Task *task) {
    if (task->getAssignedProcessor() != NULL) {
        task->getAssignedProcessor()->setAssignedTask(NULL);
        task->getAssignedProcessor()->isBusy = false;
        task->getAssignedProcessor()->setOccupiedMemorySize(0);
        task->getAssignedProcessor()->setAssignedTaskId(-1);
        task->setAssignedProcessor(NULL);
    }
}

double Tree::Merge(bool CheckMemory) {
    Task *root = this->getRoot();
    unsigned int num_subtrees = this->HowmanySubtrees(true);

    if (Cluster::getFixedCluster()->getNumberProcessors() >= num_subtrees) {
        return root->getMakespanCost(true, true);
    }

    Tree *Qtreeobj = this->BuildQtree();

    Task *node_smallest_increase;
    Task *parent;
    int shortage = num_subtrees - Cluster::getFixedCluster()->getNumberProcessors();
    double temp;
    Task *nodeone;
    Task *nodetwo;
    bool memoryEnough;

    while (shortage > 0) { //merge subtree
        //cout << "shortage " << shortage << endl;
        temp = Qtreeobj->getRoot()->getMakespanCost(true, true); //initilize ms
        temp = this->getRoot()->getMakespanCost(true, true);     //update ms

        Processor *processorToMergeTo;
        memoryEnough = estimateMS(this, Qtreeobj, node_smallest_increase, processorToMergeTo, CheckMemory);

        //when parameter checkMemory is false, memoryEnough will always be true;
        if (memoryEnough) {
            //merge currentNode (or and its sibling) to its parent
            if (node_smallest_increase->isLeaf()) {
                parent = node_smallest_increase->getParent();
                if (parent->getChildren()->size() == 2) {
                    nodeone = parent->getChildren()->front();
                    nodetwo = parent->getChildren()->back();
                    //cout<<"Merge node "<<nodeone->getOtherSideId()<<" and its sibling "<<nodetwo->getOtherSideId()<<endl;
                    nodeone->mergeToParent();
                    nodetwo->mergeToParent();
                    shortage = shortage - 2;
                    this->getTask(nodeone->getOtherSideId())->restoreEdge();
                    this->getTask(nodetwo->getOtherSideId())->restoreEdge();
                    freeProcessorIfAvailable(this->getTask(nodeone->getOtherSideId()));
                    freeProcessorIfAvailable(this->getTask(nodetwo->getOtherSideId()));
                    freeProcessorIfAvailable(this->getTask(parent->getOtherSideId()));
                    processorToMergeTo->assignTask(this->getTask(parent->getOtherSideId()));

                } else {
                    //cout<<"Merge node "<<node_smallest_increase->getOtherSideId()<<endl;
                    node_smallest_increase->mergeToParent();
                    shortage--;
                    this->getTask(node_smallest_increase->getOtherSideId())->restoreEdge();
                    freeProcessorIfAvailable(this->getTask(node_smallest_increase->getOtherSideId()));
                    freeProcessorIfAvailable(this->getTask(parent->getOtherSideId()));
                    processorToMergeTo->assignTask(this->getTask(parent->getOtherSideId()));
                }
            } else {
                //cout<<"Merge node "<<node_smallest_increase->getOtherSideId()<<endl;
                node_smallest_increase->mergeToParent();
                shortage--;
                this->getTask(node_smallest_increase->getOtherSideId())->restoreEdge();
                freeProcessorIfAvailable(this->getTask(node_smallest_increase->getOtherSideId()));
                freeProcessorIfAvailable(this->getTask(node_smallest_increase->getParent()->getOtherSideId()));
                processorToMergeTo->assignTask(
                        this->getTask(node_smallest_increase->getParent()->getOtherSideId()));
            }
            //cout<<"------------------------"<<endl;
        } else {
            break;
        }
    }

    if (shortage > 0) { //failure
        temp = -1;
    } else {
        temp = root->getMakespanCost(true, true);
    }
    delete Qtreeobj;

    return temp;
}

double
Tree::MergeOld(unsigned int num_subtrees, unsigned int processor_number, double const memory_size, bool CheckMemory) {
    Task *root = this->getRoot();

    if (processor_number >= num_subtrees) {
        return root->getMakespanCost(true, true);
    }

    Tree *Qtreeobj = this->BuildQtree();

    Task *node_smallest_increase;
    Task *parent;
    int shortage = num_subtrees - processor_number;
    double temp;
    Task *nodeone;
    Task *nodetwo;
    bool memoryEnough;

    while (shortage > 0) {//merge subtree
        //cout<<"shortage "<<shortage<<endl;
        temp = Qtreeobj->getRoot()->getMakespanCost(true, true);//initilize ms
        temp = this->getRoot()->getMakespanCost(true, true);//update ms

        //memoryEnough=increaseMS(tree, Qtreeobj, node_smallest_increase, chstart, childrenID, memory_size, CheckMemory);
        memoryEnough = estimateMSOld(this, Qtreeobj, node_smallest_increase, memory_size, CheckMemory);

        //when parameter checkMemory is false, memoryEnough will always be true;
        if (memoryEnough == true) {
            //merge currentNode (or and its sibling) to its parent
            if (node_smallest_increase->isLeaf()) {
                parent = node_smallest_increase->getParent();
                if (parent->getChildren()->size() == 2) {
                    nodeone = parent->getChildren()->front();
                    nodetwo = parent->getChildren()->back();
                    //cout<<"Merge node "<<nodeone->getothersideID()<<" and its sibling "<<nodetwo->getothersideID()<<endl;
                    nodeone->mergeToParent();
                    nodetwo->mergeToParent();
                    shortage = shortage - 2;
                    this->getTask(nodeone->getOtherSideId())->restoreEdge();
                    this->getTask(nodetwo->getOtherSideId())->restoreEdge();
                } else {
                    //cout<<"Merge node "<<node_smallest_increase->getothersideID()<<endl;
                    node_smallest_increase->mergeToParent();
                    shortage--;
                    this->getTask(node_smallest_increase->getOtherSideId())->restoreEdge();
                }
            } else {
                //cout<<"Merge node "<<node_smallest_increase->getothersideID()<<endl;
                node_smallest_increase->mergeToParent();
                shortage--;
                this->getTask(node_smallest_increase->getOtherSideId())->restoreEdge();
            }
            //cout<<"------------------------"<<endl;
        } else {
            break;
        }
    }

    if (shortage > 0) {//failure
        temp = -1;
    } else {
        temp = root->getMakespanCost(true, true);
    }
    delete Qtreeobj;

    return temp;
}


double
Tree::MergeV2(unsigned int num_subtrees, unsigned int processor_number, double const memory_size,
              bool CheckMemory) {
    if (processor_number >= num_subtrees) {
        return this->getRoot()->getMakespanCost(true, true);
    }

    this->getRoot()->getMakespanCost(true, true); //update makespan

    Tree *Qtreeobj = this->BuildQtree();
    // cout<<"qtree"<<endl;
    //  Qtreeobj->Print(cout);
    Task *currentNode;
    Task *Qroot = Qtreeobj->getRoot();
    int shortage = num_subtrees - processor_number;
    list<Task *> Llist;
    vector<unsigned int> CriticalPath;
    double temp;
    Task *largestNode;
    list<Task *>::iterator smallest;
    list<Task *>::iterator secondSmallest;
    vector<Task *> *Children;
    long pathlength;
    vector<Task *> queue;
    bool memoryCheckPass = false, leaf = false;
    Task *nodeone;
    Task *nodetwo;
    bool DeadBreak, firstTime;

    while (shortage > 0) {
        DeadBreak = true;
        firstTime = true;
        Llist.clear();
        CriticalPath.clear();
        temp = Qroot->getMakespanCost(true, true);           //update ms
        temp = this->getRoot()->getMakespanCost(true, true); //update ms

        CriticalPath.push_back(1);
        largestNode = Qroot;
        Children = largestNode->getChildren();

        while (!Children->empty()) { //initialize critical path
            temp = largestNode->getParallelPart();
            for (vector<Task *>::iterator iter = Children->begin(); iter != Children->end(); ++iter) {
                if ((*iter)->getMakespanCost(true, false) == temp) {
                    largestNode = (*iter);
                    break;
                }
            }
            CriticalPath.push_back(largestNode->getId());
            Children = largestNode->getChildren();
        }

        Children = Qroot->getChildren();
        pathlength = CriticalPath.size();

        for (unsigned int i = 1; i < pathlength; ++i) { //initialize vector L
            for (vector<Task *>::iterator iter = Children->begin(); iter != Children->end(); ++iter) {
                if ((*iter)->getId() != CriticalPath[i]) {
                    queue.push_back(*iter);
                }
            }
            Children = Qtreeobj->getTask(CriticalPath[i])->getChildren();
        }

        while (!queue.empty()) { //initialize vector L
            currentNode = queue.back();
            queue.pop_back();
            Children = currentNode->getChildren();
            for (vector<Task *>::iterator iter = Children->begin(); iter != Children->end(); ++iter) {
                Llist.push_back(*iter);
                queue.push_back(*iter);
            }
        }

        CheckOnCritical:
        if (Llist.empty()) {
            queue.push_back(Qroot);
            while (!queue.empty()) { //initialize vector L
                currentNode = queue.back();
                queue.pop_back();
                Children = currentNode->getChildren();
                for (vector<Task *>::iterator iter = Children->begin(); iter != Children->end(); ++iter) {
                    Llist.push_back(*iter);
                    queue.push_back(*iter);
                }
            }
        }

        do {
            if (Llist.size() == 1) {
                secondSmallest = Llist.begin();
                smallest = Llist.begin();
            } else {
                //   cout << "llist size "<<Llist.size() << endl;
                //for (Task *task: Llist) {
                //     cout << "task id " << task->getId() << endl;
                // }
                getTwoSmallestElement(&Llist, smallest, secondSmallest);
            }

            if ((*smallest)->isLeaf()) {
                leaf = true;
            } else {
                leaf = false;
            }
            double memoryRequired;
            if (CheckMemory) {
                memoryCheckPass = this->MemoryEnough((*smallest)->getParent(), (*smallest), leaf, memory_size,
                                                     memoryRequired);
            } else {
                memoryCheckPass = true;
            }

            if (!memoryCheckPass) {
                if ((*secondSmallest)->isLeaf()) {
                    leaf = true;
                } else {
                    leaf = false;
                }

                memoryCheckPass = this->MemoryEnough((*secondSmallest)->getParent(), *secondSmallest, leaf,
                                                     memory_size, memoryRequired);
                if (memoryCheckPass) {
                    currentNode = *secondSmallest;
                    DeadBreak = false;
                } else {
                    Llist.erase(secondSmallest);
                }
            } else {
                currentNode = *smallest;
                DeadBreak = false;
            }
        } while (!memoryCheckPass && (!Llist.empty()));


        if (DeadBreak && firstTime) {
            Llist.clear();
            firstTime = false;
            goto CheckOnCritical;
        }

        if (DeadBreak) {
            delete Qtreeobj;
            return -1; //which means failure
        }

        //merge currentNode (or and its sibling) to its parent
        if (currentNode->isLeaf()) {
            if (currentNode->getParent()->getChildren()->size() == 2) {
                nodeone = currentNode->getParent()->getChildren()->front();
                nodetwo = currentNode->getParent()->getChildren()->back();
                //cout<<"Merge node "<<nodeone->getId()<<"-"<<" and its sibling "<<nodetwo->getId()<<"-"<<endl;
                nodeone->mergeToParent();
                nodetwo->mergeToParent();
                shortage = shortage - 2;
                this->getTask(nodeone->getOtherSideId())->restoreEdge();
                this->getTask(nodetwo->getOtherSideId())->restoreEdge();
            } else {
                //cout<<"Merge node "<<currentNode->getId()<<"-"<<endl;
                currentNode->mergeToParent();
                shortage--;
                this->getTask(currentNode->getOtherSideId())->restoreEdge();
            }
        } else {
            //cout<<"Merge node "<<currentNode->getId()<<"-"<<endl;
            currentNode->mergeToParent();
            shortage--;
            this->getTask(currentNode->getOtherSideId())->restoreEdge();
        }
        //cout<<"------------------------"<<endl;
    }

    temp = this->getRoot()->getMakespanCost(true, true);
    delete Qtreeobj;

    return temp;
}

bool cmp_asap(Task *a, Task *b) { return (a->getMakespanCost(false, false) < b->getMakespanCost(false, false)); };

double Tree::ASAP() {
    unsigned int num_processors = Cluster::getFixedCluster()->getProcessors().size();
    list<Task *> PriorityQue;
    vector<Task *> BrokenEdges;
    unsigned long step_minimumMS = 0;
    double minimumMS = this->getRoot()->getMakespanCost(true, true);
    //cout<<"Excuting sequentially, makespan "<<minimumMS<<endl;
    double temp;
    Task *LargestNode;
    list<Task *>::iterator node_position;

    vector<Task *> *children = this->getRoot()->getChildren();
    while (children->size() == 1) { //avoid the linear chain
        children = children->front()->getChildren();
    }

    PriorityQue.insert(PriorityQue.end(), children->begin(), children->end());

    while (num_processors > 1) { //Breaking an edge a time
        if (PriorityQue.empty()) { //when having more processors than nodes
            break;
        }

        if (PriorityQue.size() == 1) {
            LargestNode = PriorityQue.front();
        } else {
            LargestNode = *max_element(PriorityQue.begin(), PriorityQue.end(),
                                       cmp_asap); //computation weight, no communication
        }

        node_position = find(PriorityQue.begin(), PriorityQue.end(), LargestNode);
        PriorityQue.erase(node_position);

        if (LargestNode->getParent()->getChildren()->size() > 1) {
            LargestNode->breakEdge(); //break edge
            BrokenEdges.push_back(LargestNode);
            temp = this->getRoot()->getMakespanCost(true, true);
            //cout<<"Break edge "<<LargestNode->getId()<<", makespan now: "<<temp;
            num_processors--;
            if (temp < minimumMS) {
                minimumMS = temp;
                step_minimumMS = BrokenEdges.size();
                //cout<<", makespan decreased";
            }
            //cout<<endl;
        }
        //cout<<"   pop up node "<<LargestNode->getId()<<endl;
        children = LargestNode->getChildren();
        PriorityQue.insert(PriorityQue.end(), children->begin(), children->end());
    }

    //cout<<"resotre edge ";
    unsigned long restore_index = BrokenEdges.size();
    while (restore_index > step_minimumMS) {
        BrokenEdges[restore_index - 1]->restoreEdge();
        //cout<<BrokenEdges[restore_index-1]->getId()<<" ";
        restore_index--;
    }
    //cout<<endl;

    return minimumMS;
}

bool checkRequiredMemSize(Tree *tree, Task *SubtreeRoot) {
    Tree *subtree = BuildSubtree(tree, SubtreeRoot);
    double maxout, requiredMemorySize;
    schedule_traversal *schedule_f = new schedule_traversal();
    maxout = MaxOutDegree(subtree, true);
    MinMem(subtree, maxout, requiredMemorySize, *schedule_f, true);
    //delete subtree;
    return Cluster::getFixedCluster()->smallestFreeProcessorFitting(requiredMemorySize) != nullptr;

}

bool
EstimateDecrease(int idleP, Tree *tree, vector<Task *> *criticalPath, bool *lastsubtree, Task **node_i, Task **node_j) {
    //cout<<"   --------------estimate decrease in makespan-----------------"<<endl;
    *lastsubtree = false;
    bool MSdecreased = false;
    vector<Task *> *children;
    vector<double> decreaseSequence;
    double temp, decrease = -1;
    vector<Task *> tempQue;
    Task *lastSubtreeRoot = tree->getTask(criticalPath->back()->getOtherSideId());

    double homogeneousBandwidth;
    if (Cluster::getFixedCluster()->isBandwidthHomogeneous()) {
        homogeneousBandwidth = Cluster::getFixedCluster()->getHomogeneousBandwidth();
    } else throw "Cluster not homogeneous";

    //cout<<"   Last subtree root "<<lastSubtreeRoot->getId()<<endl;
    //nodes on the last subtree of critical path
    if (idleP > 1) { //has at least 2 idle processor
        tempQue.push_back(lastSubtreeRoot);
        vector<Task *>::iterator largestNode, secondLargest;
        //cout<<"   work on the last subtree "<<criticalPath->back()->getOtherSideId()<<endl;
        while (!tempQue.empty()) {
            children = tempQue.back()->getChildren();
            tempQue.pop_back();
            if (children->size() >
                1) {                                                                        //has at least 2 children
                getTwoLargestElementTypethree(children, largestNode,
                                              secondLargest); //node_i is the largest, in terms of W
                temp = min(
                        (*largestNode)->getSequentialPart() - (*secondLargest)->getEdgeWeight() / homogeneousBandwidth,
                        (*secondLargest)->getSequentialPart() - (*largestNode)->getEdgeWeight() / homogeneousBandwidth);
                if (temp > decrease && checkRequiredMemSize(tree, *largestNode) &&
                    checkRequiredMemSize(tree, *secondLargest)) {
                    decrease = temp;
                    *node_i = *largestNode;
                    *node_j = *secondLargest;
                }

                for (vector<Task *>::iterator it = children->begin(); it != children->end(); ++it) {
                    tempQue.push_back(*it);
                    if (it != largestNode) {
                        temp = min(
                                (*largestNode)->getSequentialPart() -
                                (*secondLargest)->getEdgeWeight() / homogeneousBandwidth,
                                (*secondLargest)->getSequentialPart() -
                                (*largestNode)->getEdgeWeight() / homogeneousBandwidth);
                        if (temp > decrease && checkRequiredMemSize(tree, *largestNode) &&
                            checkRequiredMemSize(tree, *secondLargest)) {
                            decrease = temp;
                            *node_i = *largestNode;
                            *node_j = *it;
                        }
                    }
                }
            } else {
                for (vector<Task *>::iterator it = children->begin(); it != children->end(); ++it) {
                    tempQue.push_back(*it);
                }
            }
        }
    }

    if (criticalPath->size() == 1) { //only one subtree on the critical path
        if (decrease >= 0) {
            *lastsubtree = true;
            MSdecreased = true;
        }
        return MSdecreased;
    }

    //nodes on other subtrees
    double decrease_othersubtrees = -1;
    Task *output_node;
    Task *subtreeRoot;
    Task *currentNode; //current node is on the path composed of critial path nodes
    Task *nodeOnPath;
    Task *SubtreeT = criticalPath->back();
    double MS_t, W_t;

    //cout<<"other subtrees"<<endl;
    // cout<<"   working on subtree ";
    do {
        currentNode = tree->getTask(SubtreeT->getOtherSideId());
        SubtreeT = SubtreeT->getParent();
        //  cout<<"   "<<SubtreeT->getOtherSideId()<<"{ "<<endl;
        subtreeRoot = tree->getTask(SubtreeT->getOtherSideId());
        MS_t = SubtreeT->getMakespanCost(true, false);
        W_t = SubtreeT->getMakespanWeight();

        do {
            nodeOnPath = currentNode;
            currentNode = currentNode->getParent();
            tempQue.push_back(currentNode);
            while (!tempQue.empty()) {
                children = tempQue.back()->getChildren();
                tempQue.pop_back();
                for (vector<Task *>::iterator it = children->begin(); it != children->end(); ++it) {
                    if ((*it)->getId() != nodeOnPath->getId() && (!(*it)->isBroken())) {
                        //  cout<<"    "<<(*it)->getId()<<" W_i "<<(*it)->getSequentialPart()<<", MS(t) "<<MS_t<<", W_t "<<W_t<<", MS_tj "<<(*it)->getParallelPart()<<endl;
                        tempQue.push_back((*it));
                        temp = min((*it)->getSequentialPart(),
                                   MS_t - W_t - (*it)->getEdgeWeight() / homogeneousBandwidth -
                                   (*it)->getParallelPart());
                        if (temp > decrease_othersubtrees && checkRequiredMemSize(tree, *it)) {
                            decrease_othersubtrees = temp;
                            output_node = (*it);
                        }
                    }
                }
            }
        } while (!currentNode->isBroken());
        //   cout<<"   }"<<endl;
    } while (subtreeRoot->getId() != tree->getRootId());
    //  cout<<endl;
    // cout<<" decrease: "<<decrease<<" other subtrees decrease: "<<decrease_othersubtrees;
    if (decrease_othersubtrees >= 0) {
        //     cout<<"node "<<output_node->getId()<<endl;
        MSdecreased = true;
        if (decrease_othersubtrees < decrease) {
            *lastsubtree = true;
        } else {
            *node_i = output_node;
        }
    } else {
        if (decrease >= 0) {
            *lastsubtree = true;
            MSdecreased = true;
        }
    }

    return MSdecreased;
}

vector<Task *> buildCriticalPath(Task *root) {
    vector<Task *> CriticalPath;
    Task *largestNode;
    vector<Task *> *Children;
    double temp;

    CriticalPath.clear();
    CriticalPath.push_back(root);
    root->getMakespanCost(true, true);         //update critical path
    largestNode = root;
    Children = root->getChildren();
    //cout<<"critical path (subtres' roots){1 ";
    while (!Children->empty()) { //initialize critical path
        temp = largestNode->getParallelPart();
        for (vector<Task *>::iterator iter = Children->begin(); iter != Children->end(); ++iter) {
            if ((*iter)->getMakespanCost(true, false) == temp) {
                largestNode = (*iter);
                break;
            }
        }
        //cout<<largestNode->getOtherSideId()<<" ";
        CriticalPath.push_back(largestNode);
        Children = largestNode->getChildren();
    }
    return CriticalPath;
}

double Tree::SplitAgainV2(unsigned int processor_number, unsigned int num_subtrees, std::map<int, int> &taskToPrc,
                          std::map<int, bool> &isProcBusy) {
    double MS_now;
    Task *root = this->getRoot();
    Tree *Qtreeobj = this->BuildQtree();

    vector<Task *> CriticalPath; //Q nodes on Critical Path

    Task *Qroot = Qtreeobj->getRoot();
    Task *largestNode;
    Task *node_i;
    Task *node_j;
    Task *parent;
    double temp;
    vector<Task *> *Children;
    bool MSReduced, onLastSubtree;
    vector<Task *> tempVector;

    int idleProcessors = processor_number - num_subtrees;
    int currentIdleProcessor = isProcBusy[num_subtrees];
    while (idleProcessors > 0) {

        //cout<<"******** root id "<<tree->getRootId()<<" ********"<<endl;
        CriticalPath.clear();
        CriticalPath.push_back(Qroot);
        MS_now = root->getMakespanCost(true, true); //update makespan
        Qroot->getMakespanCost(true, true);         //update critical path
        largestNode = Qroot;
        Children = Qroot->getChildren();
        //cout<<"critical path (subtres' roots){1 ";
        while (!Children->empty()) { //initialize critical path
            temp = largestNode->getParallelPart();
            for (vector<Task *>::iterator iter = Children->begin(); iter != Children->end(); ++iter) {
                if ((*iter)->getMakespanCost(true, false) == temp) {
                    largestNode = (*iter);
                    break;
                }
            }
            //cout<<largestNode->getOtherSideId()<<" ";
            CriticalPath.push_back(largestNode);
            Children = largestNode->getChildren();
        }
        //cout<<"}"<<endl;

        //cout<<"Idle processor now: "<<idleProcessors<<endl;
        MSReduced = EstimateDecrease(idleProcessors, this, &CriticalPath, &onLastSubtree, &node_i, &node_j);

        if (MSReduced == true) {
            if (onLastSubtree == false) {
                // cout<<"split again cut edge "<<node_i->getId()<<endl;
                node_i->breakEdge(); //C<-C\cup C_k
                idleProcessors--;
                taskToPrc.at(node_i->getId()) = currentIdleProcessor;
                isProcBusy.at(currentIdleProcessor) = true;
                //   cout<<"is busy? "<< (isProcBusy.at(currentIdleProcessor)? "true": "false")<<endl;
                currentIdleProcessor++;

                node_i->setOtherSideId(Qtreeobj->getTasks()->size() + 1);
                parent = node_i->getParent();
                while (!parent->isBroken()) {
                    parent = parent->getParent();
                }
                Task *Qparent = Qtreeobj->getTask(parent->getOtherSideId());
                Task *Qchild;
                Task *newNode = new Task(parent->getOtherSideId(), 0, node_i->getEdgeWeight(),
                                         node_i->getSequentialPart());
                newNode->setId(Qtreeobj->getTasks()->size() + 1);
                newNode->setParent(Qparent);
                newNode->breakEdge();
                newNode->setOtherSideId(node_i->getId());
                Qparent->addChild(newNode);
                Qtreeobj->addTask(newNode);
                temp = Qparent->getMakespanWeight();
                Qparent->setMakespanWeight(temp - newNode->getMakespanWeight());
                //cout<<"create new Q node "<<newNode->getId()<<", msw "<<newNode->getMakespanWeight()<<", its parent "<<Qparent->getId()<<", msw "<<Qparent->getMakespanWeight()<<endl;

                newNode->getChildren()->clear();
                if (node_i->getParallelPart() > 0) {
                    //cout<<"went to here1."<<endl;
                    tempVector.push_back(node_i);
                    while (!tempVector.empty()) {
                        Children = tempVector.back()->getChildren();
                        tempVector.pop_back();
                        for (vector<Task *>::iterator iter = Children->begin(); iter != Children->end(); ++iter) {
                            if ((*iter)->isBroken()) {
                                //cout<<"went to here2."<<endl;
                                Qchild = Qtreeobj->getTask((*iter)->getOtherSideId());
                                newNode->addChild(Qchild);
                                Qchild->setParent(newNode);
                                Qchild->setParentId(newNode->getId());
                                Qparent->removeChild((*iter)->getOtherSideId());
                            } else {
                                tempVector.push_back((*iter));
                            }
                        }
                    }
                }
            } else {
                //  cout<<"split again cut edge "<<node_i->getId()<<" and edge "<<node_j->getId()<<endl;
                node_i->breakEdge(); //C<-C\cup C_k
                node_j->breakEdge(); //C<-C\cup C_k
                idleProcessors = idleProcessors - 2;

                taskToPrc.at(node_i->getId()) = currentIdleProcessor;
                isProcBusy.at(currentIdleProcessor) = true;
                currentIdleProcessor++;
                //    cout<<"is busy? "<<  (isProcBusy.at(currentIdleProcessor)? "true": "false")<<endl;

                taskToPrc.at(node_j->getId()) = currentIdleProcessor;
                isProcBusy.at(currentIdleProcessor) = true;
                currentIdleProcessor++;
                //      cout<<"is busy? "<<  (isProcBusy.at(currentIdleProcessor)? "true": "false")<<endl;

                node_i->setOtherSideId(Qtreeobj->getTasks()->size() + 1);
                node_j->setOtherSideId(Qtreeobj->getTasks()->size() + 2);

                Task *newNodeone = new Task(CriticalPath.back()->getId(), 0, node_i->getEdgeWeight(),
                                            node_i->getSequentialPart());
                newNodeone->setId(Qtreeobj->getTasks()->size() + 1);
                newNodeone->getChildren()->clear();
                newNodeone->setParent(CriticalPath.back());
                newNodeone->breakEdge();
                newNodeone->setOtherSideId(node_i->getId());
                CriticalPath.back()->addChild(newNodeone);
                Qtreeobj->addTask(newNodeone);
                temp = CriticalPath.back()->getMakespanWeight();
                temp = temp - newNodeone->getMakespanWeight();

                Task *newNodetwo = new Task(CriticalPath.back()->getId(), 0, node_j->getEdgeWeight(),
                                            node_j->getSequentialPart());
                newNodetwo->setId(Qtreeobj->getTasks()->size() + 1);
                newNodetwo->getChildren()->clear();
                newNodetwo->setParent(CriticalPath.back());
                newNodetwo->breakEdge();
                newNodetwo->setOtherSideId(node_j->getId());
                CriticalPath.back()->addChild(newNodetwo);
                Qtreeobj->addTask(newNodetwo);
                temp = temp - newNodetwo->getMakespanWeight();
                CriticalPath.back()->setMakespanWeight(temp);
                //cout<<"create new Q node "<<newNodetwo->getId()<<", msw "<<newNodetwo->getMakespanWeight()<<" and new node "<<newNodeone->getId()<<", msw "<<newNodeone->getMakespanWeight()<<", their parent "<<CriticalPath.back()->getId()<<", msw "<<CriticalPath.back()->getMakespanWeight()<<endl;
            }
        } else {
            break;
        }
    }

    delete Qtreeobj;

    MS_now = root->getMakespanCost(true, true);
    return MS_now;
}

double Tree::SplitAgain() {
    double MS_now;
    Tree *Qtreeobj = this->BuildQtree();
    //cout<<"qtree initially "<<endl;
    //Qtreeobj->Print(cout);
    vector<Task *> CriticalPath;//Q nodes on Critical Path
    Task *node_i;
    Task *node_j;
    Task *parent;
    double temp;

    vector<Task *> *ChildrenQTree;
    bool MSReduced, onLastSubtree;
    vector<Task *> tempVector;

    unsigned int number_subtrees = this->HowmanySubtrees(true);

    unsigned int idleProcessors = Cluster::getFixedCluster()->getNumberProcessors() - number_subtrees;
    while (idleProcessors > 0) {
        MS_now = this->getRoot()->getMakespanCost(true, true); //update makespan
        CriticalPath = buildCriticalPath(Qtreeobj);
        //cout<<"}"<<endl;

        //cout << "Idle processor now: " << idleProcessors << endl;
        MSReduced = EstimateDecrease(idleProcessors, this, &CriticalPath, &onLastSubtree, &node_i, &node_j);

        if (MSReduced == true) {
            if (onLastSubtree == false) {
                //cout << "cut edge " << node_i->getId() << endl;
                node_i->breakEdge(); //C<-C\cup C_k
                Processor *mostSuitableProcessor = Cluster::getFixedCluster()->getBiggestFreeProcessor();
                mostSuitableProcessor->assignTask(node_i);
                idleProcessors--;

                node_i->setOtherSideId(Qtreeobj->getTasks()->size() + 1);
                parent = node_i->getParent();
                while (!parent->isBroken()) {
                    parent = parent->getParent();
                }
                Task *Qparent = Qtreeobj->getTask(parent->getOtherSideId());
                Task *Qchild;
                Task *newNode = new Task(parent->getOtherSideId(), 0, node_i->getEdgeWeight(),
                                         node_i->getSequentialPart());
                newNode->setId(Qtreeobj->getTasks()->size() + 1);
                newNode->setParent(Qparent);
                newNode->breakEdge();
                newNode->setOtherSideId(node_i->getId());
                Qparent->addChild(newNode);
                Qtreeobj->addTask(newNode);
                temp = Qparent->getMakespanWeight();
                Qparent->setMakespanWeight(temp - newNode->getMakespanWeight());
                //cout<<"create new Q node "<<newNode->getId()<<", msw "<<newNode->getMakespanWeight()<<", its parent "<<Qparent->getId()<<", msw "<<Qparent->getMakespanWeight()<<endl;

                newNode->getChildren()->clear();
                if (node_i->getParallelPart() > 0) {
                    //cout<<"went to here1."<<endl;
                    tempVector.push_back(node_i);
                    while (!tempVector.empty()) {
                        ChildrenQTree = tempVector.back()->getChildren();
                        tempVector.pop_back();
                        for (vector<Task *>::iterator iter = ChildrenQTree->begin();
                             iter != ChildrenQTree->end(); ++iter) {
                            if ((*iter)->isBroken()) {
                                //cout<<"went to here2."<<endl;
                                Qchild = Qtreeobj->getTask((*iter)->getOtherSideId());
//                                cout<<"qchild "<<Qchild->getId()<<endl;
                                newNode->addChild(Qchild);
                                Qchild->setParent(newNode);
                                Qchild->setParentId(newNode->getId());
                                Qparent->removeChild((*iter)->getOtherSideId());
                            } else {
                                tempVector.push_back((*iter));
                            }
                        }
                    }
                }
            } else {
                // cout << "cut edge " << node_i->getId() << " and edge " << node_j->getId() << endl;
                node_i->breakEdge();//C<-C\cup C_k
                node_j->breakEdge();//C<-C\cup C_k
                Processor *mostSuitableProcessor = Cluster::getFixedCluster()->getBiggestFreeProcessor();
                mostSuitableProcessor->assignTask(node_i);
                mostSuitableProcessor = Cluster::getFixedCluster()->getBiggestFreeProcessor();
                mostSuitableProcessor->assignTask(node_j);
                idleProcessors = idleProcessors - 2;

                node_i->setOtherSideId(Qtreeobj->getTasks()->size() + 1);
                node_j->setOtherSideId(Qtreeobj->getTasks()->size() + 2);

                Task *newNodeone = new Task(CriticalPath.back()->getId(), 0, node_i->getEdgeWeight(),
                                            node_i->getSequentialPart());
                newNodeone->setId(Qtreeobj->getTasks()->size() + 1);
                newNodeone->getChildren()->clear();
                newNodeone->setParent(CriticalPath.back());
                newNodeone->breakEdge();
                newNodeone->setOtherSideId(node_i->getId());
                CriticalPath.back()->addChild(newNodeone);
                Qtreeobj->addTask(newNodeone);
                temp = CriticalPath.back()->getMakespanWeight();
                temp = temp - newNodeone->getMakespanWeight();

                Task *newNodetwo = new Task(CriticalPath.back()->getId(), 0, node_j->getEdgeWeight(),
                                            node_j->getSequentialPart());
                newNodetwo->setId(Qtreeobj->getTasks()->size() + 1);
                newNodetwo->getChildren()->clear();
                newNodetwo->setParent(CriticalPath.back());
                newNodetwo->breakEdge();
                newNodetwo->setOtherSideId(node_j->getId());
                CriticalPath.back()->addChild(newNodetwo);
                Qtreeobj->addTask(newNodetwo);
                temp = temp - newNodetwo->getMakespanWeight();
                CriticalPath.back()->setMakespanWeight(temp);
                //cout<<"create new Q node "<<newNodetwo->getId()<<", msw "<<newNodetwo->getMakespanWeight()<<" and new node "<<newNodeone->getId()<<", msw "<<newNodeone->getMakespanWeight()<<", their parent "<<CriticalPath.back()->getId()<<", msw "<<CriticalPath.back()->getMakespanWeight()<<endl;
            }
        } else {
            break;
        }
    }
    //  cout<<"qtree after> "<<endl;
    //   Qtreeobj->Print(cout);

    delete Qtreeobj;

    MS_now = this->getRoot()->getMakespanCost(true, true);
    return MS_now;
}

double SplitAgainOld(Tree *tree, unsigned int processor_number, unsigned int num_subtrees) {
    double MS_now;
    Task *root = tree->getRoot();
    Tree *Qtreeobj = tree->BuildQtreeOld();
    //Qtreeobj->Print(cout);
    // Qtreeobj = tree->BuildQtreeOld();
    //Qtreeobj->Print(cout);
    vector<Task *> CriticalPath;//Q nodes on Critical Path

    Task *Qroot = Qtreeobj->getRoot();
    Task *largestNode;
    Task *node_i;
    Task *node_j;
    Task *parent;
    double temp;
    vector<Task *> *Children;
    bool MSReduced, onLastSubtree;
    vector<Task *> tempVector;

    int idleProcessors = processor_number - num_subtrees;
    while (idleProcessors > 0) {
        //cout<<"******** root id "<<tree->getRootId()<<" ********"<<endl;
        CriticalPath.clear();
        CriticalPath.push_back(Qroot);
        MS_now = root->getMakespanCost(true, true);//update makespan
        Qroot->getMakespanCost(true, true);//update critical path
        largestNode = Qroot;
        Children = Qroot->getChildren();
        //cout<<"critical path (subtres' roots){1 ";
        while (!Children->empty()) {//initialize critical path
            temp = largestNode->getParallelPart();
            for (vector<Task *>::iterator iter = Children->begin(); iter != Children->end(); ++iter) {
                if ((*iter)->getMakespanCost(true, false) == temp) {
                    largestNode = (*iter);
                    break;
                }
            }
            //cout<<largestNode->getothersideID()<<" ";
            CriticalPath.push_back(largestNode);
            Children = largestNode->getChildren();
        }
        //cout<<"}"<<endl;

        //cout<<"Idle processor now: "<<idleProcessors<<endl;
        MSReduced = EstimateDecrease(idleProcessors, tree, &CriticalPath, &onLastSubtree, &node_i, &node_j);

        if (MSReduced == true) {
            if (onLastSubtree == false) {
                //cout<<"cut edge "<<node_i->getId()<<endl;
                node_i->breakEdge();//C<-C\cup C_k
                idleProcessors--;

                node_i->setOtherSideId(Qtreeobj->getTasks()->size() + 1);
                parent = node_i->getParent();
                while (!parent->isBroken()) {
                    parent = parent->getParent();
                }
                Task *Qparent = Qtreeobj->getTask(parent->getOtherSideId());
                Task *Qchild;
                Task *newNode = new Task(parent->getOtherSideId(), 0, node_i->getEdgeWeight(),
                                         node_i->getSequentialPart());
                newNode->setId(Qtreeobj->getTasks()->size() + 1);
                newNode->setParent(Qparent);
                newNode->breakEdge();
                newNode->setOtherSideId(node_i->getId());
                Qparent->addChild(newNode);
                Qtreeobj->addTask(newNode);
                temp = Qparent->getMakespanWeight();
                Qparent->setMakespanWeight(temp - newNode->getMakespanWeight());
                //cout<<"create new Q node "<<newNode->getId()<<", msw "<<newNode->getMakespanWeight()<<", its parent "<<Qparent->getId()<<", msw "<<Qparent->getMakespanWeight()<<endl;

                newNode->getChildren()->clear();
                if (node_i->getParallelPart() > 0) {
                    //cout<<"went to here1."<<endl;
                    tempVector.push_back(node_i);
                    while (!tempVector.empty()) {
                        Children = tempVector.back()->getChildren();
                        tempVector.pop_back();
                        for (vector<Task *>::iterator iter = Children->begin(); iter != Children->end(); ++iter) {
                            if ((*iter)->isBroken()) {
                                //cout<<"went to here2."<<endl;
                                Qchild = Qtreeobj->getTask((*iter)->getOtherSideId());
                                //  cout<<"qchild "<<Qchild->getId()<<endl;
                                newNode->addChild(Qchild);
                                Qchild->setParent(newNode);
                                Qchild->setParentId(newNode->getId());
                                Qparent->removeChild((*iter)->getOtherSideId());
                            } else {
                                tempVector.push_back((*iter));
                            }
                        }
                    }
                }
            } else {
                //cout<<"cut edge "<<node_i->getId()<<" and edge "<<node_j->getId()<<endl;
                node_i->breakEdge();//C<-C\cup C_k
                node_j->breakEdge();//C<-C\cup C_k
                idleProcessors = idleProcessors - 2;

                node_i->setOtherSideId(Qtreeobj->getTasks()->size() + 1);
                node_j->setOtherSideId(Qtreeobj->getTasks()->size() + 2);

                Task *newNodeone = new Task(CriticalPath.back()->getId(), 0, node_i->getEdgeWeight(),
                                            node_i->getSequentialPart());
                newNodeone->setId(Qtreeobj->getTasks()->size() + 1);
                newNodeone->getChildren()->clear();
                newNodeone->setParent(CriticalPath.back());
                newNodeone->breakEdge();
                newNodeone->setOtherSideId(node_i->getId());
                CriticalPath.back()->addChild(newNodeone);
                Qtreeobj->addTask(newNodeone);
                temp = CriticalPath.back()->getMakespanWeight();
                temp = temp - newNodeone->getMakespanWeight();

                Task *newNodetwo = new Task(CriticalPath.back()->getId(), 0, node_j->getEdgeWeight(),
                                            node_j->getSequentialPart());
                newNodetwo->setId(Qtreeobj->getTasks()->size() + 1);
                newNodetwo->getChildren()->clear();
                newNodetwo->setParent(CriticalPath.back());
                newNodetwo->breakEdge();
                newNodetwo->setOtherSideId(node_j->getId());
                CriticalPath.back()->addChild(newNodetwo);
                Qtreeobj->addTask(newNodetwo);
                temp = temp - newNodetwo->getMakespanWeight();
                CriticalPath.back()->setMakespanWeight(temp);
                //cout<<"create new Q node "<<newNodetwo->getId()<<", msw "<<newNodetwo->getMakespanWeight()<<" and new node "<<newNodeone->getId()<<", msw "<<newNodeone->getMakespanWeight()<<", their parent "<<CriticalPath.back()->getId()<<", msw "<<CriticalPath.back()->getMakespanWeight()<<endl;
            }
        } else { break; }
    }

    delete Qtreeobj;

    MS_now = root->getMakespanCost(true, true);
    return MS_now;
}


vector<Task *> Tree::buildCriticalPath(Tree *Qtreeobj) {
    Task *Qroot = Qtreeobj->getRoot();

    double temp;
    Task *largestNodeQTree;
    vector<Task *> CriticalPath;
    CriticalPath.clear();
    CriticalPath.push_back(Qroot);

    Qroot->getMakespanCost(true, true);         //update critical path
    largestNodeQTree = Qroot;
    vector<Task *> *ChildrenQTree = Qroot->getChildren();
    //cout<<"critical path (subtres' roots){1 ";
    while (!ChildrenQTree->empty()) {//initialize critical path
        temp = largestNodeQTree->getParallelPart();
        for (vector<Task *>::iterator iter = ChildrenQTree->begin(); iter != ChildrenQTree->end(); ++iter) {
            if ((*iter)->getMakespanCost(true, false) == temp) {
                largestNodeQTree = (*iter);
                break;
            }
        }
        //cout<<largestNodeQTree->getOtherSideId()<<" ";
        CriticalPath.push_back(largestNodeQTree);
        ChildrenQTree = largestNodeQTree->getChildren();
    }
    return CriticalPath;
}


list<unsigned int> markSubtreeTillBottom(Task *subtreeRoot) {
    list<unsigned int> allNodes;
    list<Task *> queue;
    allNodes.clear();
    queue.clear();
    queue.push_back(subtreeRoot);

    do {
        Task *firstTask = queue.front();
        allNodes.push_back(firstTask->getId());
        queue.pop_front();
        for (Task *childOfFirst: *firstTask->getChildren()) {
            queue.push_back(childOfFirst);
        }
    } while (!queue.empty());
    return allNodes;
}

void
Immediately(Tree *tree, int *schedule,
            double m_availble, vector<unsigned int> *brokenEdges) {
    unsigned int treeSize = tree->getTasks()->size();
    double memory_occupied = tree->getTask(schedule[treeSize - 1])->getEdgeWeight();

    list<unsigned int> allNodes;

    unsigned long subtree_size;
    list<unsigned int>::iterator iter;
    schedule_traversal *schedule_f = new schedule_traversal();
    list<int>::iterator ite_sche;
    double maxoutD, memory_required, node_cost, data_to_unload;
    int cur_task_id;
    vector<unsigned int> subtreeBrokenEdges;

    for (unsigned long rank = treeSize - 1; rank >= 1; rank--) {
        cur_task_id = schedule[rank];
        if (cur_task_id != 0) { //=0 means this node has already been moved to another processor
            //cout<<" "<<cur_task_id;
            Task *currTask = tree->getTask(cur_task_id);
            node_cost = currTask->getEdgeWeight() + currTask->getNodeWeight();
            for (auto child: *currTask->getChildren()) {
                node_cost += child->getEdgeWeight();
            }

            data_to_unload = memory_occupied + node_cost - currTask->getEdgeWeight() - m_availble;
            if (data_to_unload > 0) { // schedule the subtree that is rooted at this node onto another processor

                currTask->breakEdge(); // set it cut, used for building a quotient tree later
                brokenEdges->push_back(currTask->getOtherSideId());

                allNodes = markSubtreeTillBottom(currTask);
                subtree_size = allNodes.size();

                for (long i = rank - 1; i >= 0; i--) {
                    iter = find(allNodes.begin(), allNodes.end(), schedule[i]);
                    //If node is in subtree
                    if (iter != allNodes.end()) {
                        schedule[i] = 0; //IO counter will pass 0;
                        allNodes.erase(iter);
                    }
                    if (allNodes.size() == 1) {
                        break;
                    }
                }

                Tree *subtree = BuildSubtree(tree, currTask);

                subtree_size = subtree->getTasks()->size();

                int *schedule_copy = new int[subtree_size + 1];
                maxoutD = MaxOutDegree(subtree, true);
                schedule_f->clear();

                MinMem(subtree, maxoutD, memory_required, *schedule_f, true);
                ite_sche = schedule_f->begin();
                for (unsigned int i = subtree_size; i >= 1; --i) {
                    schedule_copy[i] = *ite_sche;
                    advance(ite_sche, 1);
                }
                schedule_copy[0] = subtree_size + 1;

                if (memory_required > m_availble) {
                    Immediately(subtree, schedule_copy,
                                m_availble, &subtreeBrokenEdges);

                    for (vector<unsigned int>::iterator iter = subtreeBrokenEdges.begin();
                         iter != subtreeBrokenEdges.end(); ++iter) {
                        brokenEdges->push_back(tree->getTask(*iter)->getOtherSideId());
                    }
                }

                delete[] schedule_copy;
                delete subtree;

                memory_occupied -= currTask->getEdgeWeight();
                memory_occupied = max(0.0, memory_occupied);
            } else { //memory is enough for executing this node
                memory_occupied += node_cost - 2 * currTask->getEdgeWeight() - currTask->getNodeWeight();
                memory_occupied = max(0.0, memory_occupied);
            }
        }
    }
    delete schedule_f;
    //cout<<endl;
}

int MemoryCheck(Tree *tree, io_method_t method, bool useMinimalAvailableProvcessor) {

    vector<std::tuple<double, Tree *, schedule_traversal *, Task *>> subtreeRoots;
    Task *currentnode;
    std::tuple<double, Tree *, schedule_traversal *, Task *> subtreeRootAndMemReq;
    tree->getRoot()->breakEdge();
    //subtreeRoots = tree->getBrokenTasks();
    for (unsigned int i = tree->getSize(); i >= 1; --i) {
        double maxoutD, memory_required;
        currentnode = tree->getTask(i);
        if (currentnode->isBroken()) {
            //cout<<i<<" ";
            Tree *subtree = BuildSubtree(tree, currentnode);
            maxoutD = MaxOutDegree(subtree, true);
            schedule_traversal *schedule_f = new schedule_traversal();
            MinMem(subtree, maxoutD, memory_required, *schedule_f, true);
            tuple<double, Tree *, schedule_traversal *, Task *> tuple{
                    memory_required,
                    subtree, schedule_f,
                    currentnode};
            subtreeRoots.push_back(tuple);
        }
    }
    sort(subtreeRoots.begin(), subtreeRoots.end(),
         [](const std::tuple<double, Tree *, schedule_traversal *, Task *> &lhs,
            const std::tuple<double, Tree *, schedule_traversal *, Task *> &rhs) {
             double lhsMem = get<0>(lhs);
             double rhsMem = get<0>(rhs);
             return lhsMem > rhsMem;
         });
    //processors are sorted when read
    unsigned int com_freq;
    vector<unsigned int> BrokenEdgesID;
    while (!subtreeRoots.empty()) {
        //   tree->HowmanySubtrees(false);
        subtreeRootAndMemReq = subtreeRoots.front();
        subtreeRoots.erase(subtreeRoots.begin());
        Processor *mostSuitableProcessor;
        double currentMemorySize;

        try {
            //this is possible, because both subtrees and processors are sorted by size
            mostSuitableProcessor = Cluster::getFixedCluster()->getBiggestFreeProcessor();
        }
        catch (exception e) {
            mostSuitableProcessor = Cluster::getFixedCluster()->getFirstFreeProcessorOrSmallest();
        }
        mostSuitableProcessor->assignTask(get<3>(subtreeRootAndMemReq));
        currentMemorySize = mostSuitableProcessor->getMemorySize();

        // cout << "using proc of memory " << mostSuitableProcessor->getMemorySize() << "for task "
        //     << get<3>(subtreeRootAndMemReq)->getId() << endl;
        //cout << "Subtree " << get<3>(subtreeRootAndMemReq)->getId() << " needs memory "
        //    << get<0>(subtreeRootAndMemReq)
        //    << endl;
        if (get<0>(subtreeRootAndMemReq) > currentMemorySize) {
            //   cout << ", larger than what is available: " << currentMemorySize << endl;

            int *schedule_copy = copyScheduleBackwards(get<2>(subtreeRootAndMemReq));

            switch (method) {
                case FIRST_FIT:
                    try {
                        Tree *&subtree = get<1>(subtreeRootAndMemReq);
                        IOCounter(subtree, schedule_copy, false, true,
                                  com_freq, &BrokenEdgesID,
                                  FIRST_FIT);
                    }
                    catch (exception e) {
                        cout << "not enough processors";
                        return -1;
                    }
                    break;
                case LARGEST_FIT:
                    IOCounter(get<1>(subtreeRootAndMemReq), schedule_copy, false, true,
                              com_freq, &BrokenEdgesID,
                              LARGEST_FIT);
                    break;
                case IMMEDIATELY:
                    Immediately(get<1>(subtreeRootAndMemReq), schedule_copy, currentMemorySize, &BrokenEdgesID);
                    break;

                default:
                    break;
            }
            delete[] schedule_copy;
        } else {
            if (Cluster::getFixedCluster()->hasFreeProcessor()) {
                mostSuitableProcessor->setOccupiedMemorySize(get<0>(subtreeRootAndMemReq));
            }
        }
        //cout<<endl;
        Cluster::getFixedCluster()->assignTasksForIds(tree);
        // cout << "deleting" << endl;
        //  delete get<1>(subtreeRootAndMemReq);
        //   delete get<2>(subtreeRootAndMemReq);
    }


    for (vector<unsigned int>::iterator iter = BrokenEdgesID.begin(); iter != BrokenEdgesID.end(); ++iter) {
        tree->getTask(*iter)->breakEdge();
    }
    // for (Processor *p: Cluster::getFixedCluster()->getProcessors()) {
    //     if(p->getAssignedTaskId()!=)
    //      p->assignTask(tree->getTask(p->getAssignedTaskId()));
    //   }
    // Cluster::getFixedCluster()->printProcessors();

    for (auto tuple: subtreeRoots) {
        delete get<1>(tuple);
        delete get<2>(tuple);
    }
    return 0;
}


int MemoryCheckHomp(Tree *tree, io_method_t method, double processor_memory_size) {
    vector<std::tuple<double, Tree *, schedule_traversal *, Task *>> subtreeRoots;
    Task *currentnode;
    std::tuple<double, Tree *, schedule_traversal *, Task *> subtreeRootAndMemReq;
    tree->getRoot()->breakEdge();
    //subtreeRoots = tree->getBrokenTasks();
    for (unsigned int i = tree->getSize(); i >= 1; --i) {
        double maxoutD, memory_required;
        currentnode = tree->getTask(i);
        if (currentnode->isBroken()) {
            //cout<<i<<" ";
            Tree *subtree = BuildSubtree(tree, currentnode);
            maxoutD = MaxOutDegree(subtree, true);
            schedule_traversal *schedule_f = new schedule_traversal();
            MinMem(subtree, maxoutD, memory_required, *schedule_f, true);
            tuple<double, Tree *, schedule_traversal *, Task *> tuple{
                    memory_required,
                    subtree, schedule_f,
                    currentnode};
            subtreeRoots.push_back(tuple);
        }
    }
    unsigned int com_freq;
    vector<unsigned int> BrokenEdgesID;
    while (!subtreeRoots.empty()) {
        //   tree->HowmanySubtrees(false);
        subtreeRootAndMemReq = subtreeRoots.front();
        subtreeRoots.erase(subtreeRoots.begin());
        cout << "Subtree " << get<3>(subtreeRootAndMemReq)->getId() << " needs memory "
             << get<0>(subtreeRootAndMemReq)
             << endl;
        if (get<0>(subtreeRootAndMemReq) > processor_memory_size) {
            cout << ", larger than what is available: " << processor_memory_size << endl;

            int *schedule_copy = copyScheduleBackwards(get<2>(subtreeRootAndMemReq));

            switch (method) {
                case FIRST_FIT:
                    try {
                        IOCounter(get<1>(subtreeRootAndMemReq), schedule_copy, false, true,
                                  com_freq, &BrokenEdgesID,
                                  FIRST_FIT);
                    }
                    catch (exception e) {
                        cout << "not enough processors";
                        return -1;
                    }
                    break;
                case LARGEST_FIT:
                    IOCounter(get<1>(subtreeRootAndMemReq), schedule_copy, false, true,
                              com_freq, &BrokenEdgesID,
                              LARGEST_FIT);
                    break;
                case IMMEDIATELY:
                    Immediately(get<1>(subtreeRootAndMemReq), schedule_copy, processor_memory_size, &BrokenEdgesID);
                    break;

                default:
                    break;
            }
            delete[] schedule_copy;
        }
    }


    for (vector<unsigned int>::iterator iter = BrokenEdgesID.begin(); iter != BrokenEdgesID.end(); ++iter) {
        tree->getTask(*iter)->breakEdge();
    }
    // for (Processor *p: Cluster::getFixedCluster()->getProcessors()) {
    //     if(p->getAssignedTaskId()!=)
    //      p->assignTask(tree->getTask(p->getAssignedTaskId()));
    //   }
    // Cluster::getFixedCluster()->printProcessors();

    for (auto tuple: subtreeRoots) {
        delete get<1>(tuple);
        delete get<2>(tuple);
    }
    return 0;
}

void MemoryCheckA2(Tree *tree, Cluster *cluster, io_method_t method, bool skipBig) {
    vector<Task *> subtreeRoots;
    vector<Task *> subtreeRootsSkipped;
    Task *currentnode;
    Task *subtreeRoot;
    tree->getRoot()->breakEdge();

    //cout<<"Subtrees' roots: ";
    unsigned long treeSize = tree->getTasks()->size();
    for (unsigned int i = treeSize; i >= 1; --i) {
        currentnode = tree->getTask(i);
        if (currentnode->isBroken()) {
            //cout<<i<<" ";
            subtreeRoots.push_back(currentnode);
            //cout << "root " << currentnode->getMakespanCost() << endl;
        }
    }
    //cout<<endl;
    sort(subtreeRoots.begin(), subtreeRoots.end(),
         [](Task *lhs, Task *rhs) { return lhs->getMakespanCost() < rhs->getMakespanCost(); });

    double maxoutD, memory_required;
    schedule_traversal *schedule_f = new schedule_traversal();
    unsigned int com_freq;
    unsigned long subtreeSize;
    list<int>::iterator ite_sche;
    vector<unsigned int> BrokenEdgesID;
    Processor *currentProcessor = cluster->getBiggestFreeProcessor();


    while (!subtreeRoots.empty()) {
        double currentMem = currentProcessor->getMemorySize();

        subtreeRoot = subtreeRoots.back();
        subtreeRoots.pop_back();
        //     cout<<"get free proc "<<currentMem<<" for task "<<subtreeRoot->getId()<<endl;

        Tree *subtree = BuildSubtree(tree, subtreeRoot);

        subtreeSize = subtree->getTasks()->size();
        int *schedule_copy = new int[subtreeSize + 1];
        maxoutD = MaxOutDegree(subtree, true);
        schedule_f->clear();

        MinMem(subtree, maxoutD, memory_required, *schedule_f, true);

        //  cout << "Subtree " << subtreeRoot->getId() << " needs memory " << memory_required;
        if (memory_required > currentMem) {
            if (!skipBig) {
                //   cout << ", larger than what is available: " << currentMem << " on proc " << currentProcessor << endl;

                ite_sche = schedule_f->begin();
                for (unsigned int i = subtreeSize; i >= 1; --i) {
                    schedule_copy[i] = *ite_sche;
                    advance(ite_sche, 1);
                }

                schedule_copy[0] = subtreeSize + 1;

                switch (method) {
                    case FIRST_FIT:
                        IOCounterWithVariableMem(subtree, schedule_copy, cluster, false, true, com_freq,
                                                 &BrokenEdgesID, FIRST_FIT);
                        break;
                    case LARGEST_FIT:
                        IOCounterWithVariableMem(subtree, schedule_copy, cluster, false, true, com_freq,
                                                 &BrokenEdgesID, LARGEST_FIT);
                        break;
                    case IMMEDIATELY:
                        Immediately(subtree, schedule_copy, currentMem, &BrokenEdgesID);
                        break;

                    default:
                        break;
                }
            } else {
                subtreeRootsSkipped.push_back(subtreeRoot);
            }
        } else {
            //  currentProcessor->assignTask(subtreeRoot);
            currentProcessor->assignTask(subtreeRoot);
            currentProcessor = cluster->getBiggestFreeProcessor();
            currentProcessor->setOccupiedMemorySize(memory_required);
            //        cout<<"got new proc in else "<<currentProcessor->getMemorySize()<<endl;
        }
        //cout<<endl;
        //
        //   cout << "Now big trees" << endl;
        while (!subtreeRootsSkipped.empty()) {
            double currentMem = currentProcessor->getMemorySize();
            subtreeRoot = subtreeRootsSkipped.back();
            subtreeRootsSkipped.pop_back();
            //     cout << ", larger than what is available: " << currentMem << " on proc " << currentProcessor << endl;

            ite_sche = schedule_f->begin();
            for (unsigned int i = subtreeSize; i >= 1; --i) {
                schedule_copy[i] = *ite_sche;
                advance(ite_sche, 1);
            }

            schedule_copy[0] = subtreeSize + 1;

            switch (method) {
                case FIRST_FIT:
                    IOCounterWithVariableMem(subtree, schedule_copy, cluster, false, true, com_freq,
                                             &BrokenEdgesID, FIRST_FIT);
                    break;
                case LARGEST_FIT:
                    IOCounterWithVariableMem(subtree, schedule_copy, cluster, false, true, com_freq,
                                             &BrokenEdgesID, LARGEST_FIT);
                    break;
                case IMMEDIATELY:
                    Immediately(subtree, schedule_copy, currentMem, &BrokenEdgesID);
                    break;

                default:
                    break;
            }
        }

        delete[] schedule_copy;
        delete subtree;
    }
    delete schedule_f;

    for (vector<unsigned int>::iterator iter = BrokenEdgesID.begin(); iter != BrokenEdgesID.end(); ++iter) {
        tree->getTask(*iter)->breakEdge();
    }
}

void Cluster::SetBandwidth(double CCR, Tree *tree) {
    double sum_edges = 0;
    double sum_weights = 0;

    for (Task *task: *tree->getTasks()) {
        sum_edges += task->getEdgeWeight();
        sum_weights += task->getMakespanWeight();
    }

    if (this->isHomogeneous()) {
        double bandwidth = sum_edges / (sum_weights * CCR);
        this->setHomogeneousBandwidth(bandwidth);
    }

}

int *
copyScheduleBackwards(schedule_traversal *schedule_f) {

    list<int>::iterator ite_sche = schedule_f->begin();
    int *schedule_copy = new int[schedule_f->size() + 1];
    for (unsigned int i = schedule_f->size(); i >= 1; --i) {
        schedule_copy[i] = *ite_sche;
        advance(ite_sche, 1);
    }
    schedule_copy[0] = schedule_f->size() + 1;

    return schedule_copy;
}

class TMaxHeap {
public:
    bool operator()(Task *t1, Task *t2) {
        return t1->getTMax() < t2->getTMax();
    }


// since we are going over the heap from top to bottom,
// everything above the idx is already compliant with the heap constraint
    static void siftUp(vector<Task *> *heap, int idx) {
        int parent_idx;
        Task *tmp;

        while (idx > 0) {
            parent_idx = (int) (idx - 1) / 2;
            tmp = heap->at(idx);
            if (heap->at(parent_idx)->getTMax() < heap->at(idx)->getTMax()) {
                heap->at(idx) = heap->at(parent_idx);
                heap->at(parent_idx) = tmp;
                idx = parent_idx;
            } else {
                break;
            }

        }

    }

};


void distributeProcessors(Tree *qTree, string chooseSubtreeAssign) {
    Task *root = qTree->getRoot();
    if (root->getFeasibleProcessors()->empty()) {
        throw "No feasible processors";
    }
    Processor *pFast = root->getFastestFeasibleProcessor();
    pFast->assignTask(root);
    removeProcessorFromAllFeasSets(pFast, qTree);

    auto tasks = qTree->getTasks();

    vector<Task *> *taskHeap = new vector<Task *>();//tasks->size());
    for (auto task: *tasks) {
        if (task->getAssignedProcessor() == NULL) {
            task->updateTMax();
            taskHeap->push_back(task);
        }
    }

    auto compareTasksForHeap = [chooseSubtreeAssign, qTree](Task *t1, Task *t2) -> bool {
        if (t1->getTMax() == DBL_MAX && t2->getTMax() == DBL_MAX) return false;
        else if (t1->getTMax() == DBL_MAX) return false;
        else if (t2->getTMax() == DBL_MAX) return true;
        else {
            //to infinity values in tmax; compare upon condition
            if (chooseSubtreeAssign == "CP") {
                vector<Task *> CriticalPath = buildCriticalPath(qTree->getRoot());
                int idt1 = t1->getId();
                const vector<Task *>::iterator &iteratort1 = std::find_if(CriticalPath.begin(), CriticalPath.end(),
                                                                          [idt1](Task *task) {
                                                                              return task->getId() == idt1;
                                                                          });
                int idt2 = t2->getId();
                const vector<Task *>::iterator &iteratort2 = std::find_if(CriticalPath.begin(), CriticalPath.end(),
                                                                          [idt2](Task *task) {
                                                                              return task->getId() == idt2;
                                                                          });
                if (iteratort1 == CriticalPath.end() && iteratort2 == CriticalPath.end())
                    return t1->getTMax() < t2->getTMax();
                else if (iteratort1 == CriticalPath.end()) return true;
                else if (iteratort2 == CriticalPath.end()) return false;
                else {
                    return t1->getTMax() < t2->getTMax();
                }
            } else if (chooseSubtreeAssign == "MW") {
                return t1->getMinMemUnderlying() < t2->getMinMemUnderlying();

            } else if (chooseSubtreeAssign == "MD") {
                return t1->getChildren()->size() < t2->getChildren()->size();
            } else return t1->getTMax() < t2->getTMax();
        }

    };
    make_heap(taskHeap->begin(), taskHeap->end(), compareTasksForHeap);

    while (taskHeap->size() > 0) {

        pop_heap(taskHeap->begin(), taskHeap->end(), compareTasksForHeap);
        Task *task = taskHeap->back();
        taskHeap->pop_back();

        Processor *pFast = task->getFastestFeasibleProcessor();
        pFast->assignTask(task);
        for (int i = 0; i < taskHeap->size(); i++) {
            taskHeap->at(i)->deleteFeasible(pFast);
            taskHeap->at(i)->updateTMax();
            //TODO: is sift up corect?
            TMaxHeap::siftUp(taskHeap, i);
            //SiftInfTmaxUpPreserveOrder(taskHeap);
        }
        //SiftInfTmaxUpPreserveOrder(taskHeap);

    }
    delete taskHeap;
}

void distributeProcessorsOld(Tree *qTree) {
    Task *root = qTree->getRoot();
    if (root->getFeasibleProcessors()->empty()) {
        throw "No feasible processors";
    }
    Processor *pFast = root->getFastestFeasibleProcessor();
    pFast->assignTask(root);
    removeProcessorFromAllFeasSets(pFast, qTree);

    auto tasks = qTree->getTasks();

    vector<Task *> *taskHeap = new vector<Task *>();//tasks->size());
    for (auto task: *tasks) {
        if (task->getAssignedProcessor() == NULL) {
            task->updateTMax();
            taskHeap->push_back(task);
        }
    }
    make_heap(taskHeap->begin(), taskHeap->end(), TMaxHeap());

    while (taskHeap->size() > 0) {
        pop_heap(taskHeap->begin(), taskHeap->end(), TMaxHeap());
        Task *task = taskHeap->back();
        taskHeap->pop_back();
        Processor *pFast = task->getFastestFeasibleProcessor();
        pFast->assignTask(task);
        for (int i = 0; i < taskHeap->size(); i++) {
            taskHeap->at(i)->deleteFeasible(pFast);
            taskHeap->at(i)->updateTMax();
            TMaxHeap::siftUp(taskHeap, i);
        }

    }
    delete taskHeap;
}

void SiftInfTmaxUpPreserveOrder(vector<Task *> *taskHeap) {
    int init_size = taskHeap->size();
    if (init_size != 0 && init_size != 1) {

        vector<Task *> withTMaxInf;
        for (const auto &item: *taskHeap) {
            if (item->getTMax() == DBL_MAX)
                withTMaxInf.push_back(item);
        }

        vector<Task *> others;
        set_difference(taskHeap->begin(), taskHeap->end(),
                       withTMaxInf.begin(), withTMaxInf.end(),
                //inserter(others, others.begin())
                       std::back_inserter(others),
                       [](auto &a, auto &b) { return a->getId() < b->getId(); });

        taskHeap->clear();
        taskHeap->insert(taskHeap->end(), withTMaxInf.begin(), withTMaxInf.end());
        taskHeap->insert(taskHeap->end(), others.begin(), others.end());
    }

    assert(taskHeap->size() == init_size);
}

string growSeqSetWhileImprovesMakespan2(list<Task *> &seqSet, Tree *tree) {

    int cntrAdditionToSS = 0;
    int cntTries = 0;
    Task *root = tree->getRoot();
    list<Task *> frontier;


    for (Task *task: *tree->getTasks()) {
        if (task->isBroken()) {
            frontier.push_back(task);
        }
    }
    frontier.sort(cmp_Mem_nodecreasing);

    double minMakespan = assignToBestProcessors(tree, {});
    tree->cleanAssignedAndReassignFeasible();
    SeqSet optimalSeqSet = SeqSet(tree, minMakespan);

    int maxNumberChildren = 0;
    while (!frontier.empty()) {
        tree->cleanAssignedAndReassignFeasible();
        // cout<<"try add new"<<endl;
        Task *potentialAddition = frontier.front();
        frontier.pop_front();

        if (maxNumberChildren < potentialAddition->getChildren()->size())
            maxNumberChildren = potentialAddition->getChildren()->size();

        if (!potentialAddition->isRoot()) potentialAddition->restoreEdge();
        for (Task *child: *potentialAddition->getChildren()) {
            child->breakEdge();
        }
        cout.precision(20);

        double potentialMakespan = assignToBestProcessors(tree, {});
        cout << "potential: " << potentialMakespan << ", current min: " << " " << minMakespan << " equals? "
             << (potentialMakespan == minMakespan ? "y" : "n");// << endl;
        cntTries++;

        if (potentialMakespan <= minMakespan) {
            cout << " add!" << endl;
            cntrAdditionToSS++;
            minMakespan = potentialMakespan;
            //no breaking back the edge, we accept this improvement
            // instead add children of newly added task
            for (Task *child: *potentialAddition->getChildren()) {
                auto it = upper_bound(frontier.begin(), frontier.end(), child,
                                      cmp_Mem_nodecreasing);
                frontier.insert(it, child);
            }
        } else {
            // cout << " no add!" << endl;
            potentialAddition->breakEdge();
            for (Task *child: *potentialAddition->getChildren()) {
                child->restoreEdge();
            }
            assignToBestProcessors(tree, {});
            // no adding children
        }

        // cout << "#trees " << tree->HowmanySubtrees(false) << " added " << cntrAdditionToSS << endl;
    }
    // cout << "tried " << cntTries << " added to SS " << cntrAdditionToSS << " max # children " << maxNumberChildren
    //<< endl;
    return "tried " + to_string(cntTries) + " added to SS " + to_string(cntrAdditionToSS);
}

string growSeqSetWhileImprovesMakespanAllowWorse(list<Task *> &seqSet, Tree *tree) {

    int cntrAdditionToSS = 0;
    int cntTries = 0;
    Task *root = tree->getRoot();
    list<Task *> frontier;


    for (Task *task: *tree->getTasks()) {
        if (task->isBroken()) {
            frontier.push_back(task);
        }
    }
    frontier.sort(cmp_Mem_nodecreasing);

    double minMakespan = assignToBestProcessors(tree, {});
    tree->cleanAssignedAndReassignFeasible();
    SeqSet optimalSeqSet = SeqSet(tree, minMakespan, 1);


    while (!frontier.empty()) {
        tree->cleanAssignedAndReassignFeasible();
        // cout<<"try add new"<<endl;
        Task *potentialAddition = frontier.front();
        frontier.pop_front();


        if (!potentialAddition->isRoot()) potentialAddition->restoreEdge();
        for (Task *child: *potentialAddition->getChildren()) {
            child->breakEdge();
        }
        cout.precision(20);

        double potentialMakespan = assignToBestProcessors(tree, {});
        // cout << "potential: " << potentialMakespan << ", current min: " << " " << minMakespan << " equals? "
        //     << (potentialMakespan == minMakespan ? "y" : "n") << endl;
        cntTries++;
        for (Task *child: *potentialAddition->getChildren()) {
            auto it = upper_bound(frontier.begin(), frontier.end(), child,
                                  cmp_Mem_nodecreasing);
            frontier.insert(it, child);
        }
        if (potentialMakespan <= minMakespan) {
            // cout << " add!" << endl;
            cntrAdditionToSS++;
            minMakespan = potentialMakespan;
            optimalSeqSet = SeqSet(tree, minMakespan, cntTries);
            //  optimalSeqSet.print();

        } else {
            //cout << " no add!" << endl;
            // no adding children
        }

        // cout << "#trees " << tree->HowmanySubtrees(false) << " added " << cntrAdditionToSS << endl;
    }
    optimalSeqSet.implementSeqSet(tree);
    assignToBestProcessors(tree, {});
    return optimalSeqSet.print();
}

string growSeqSetWhileImprovesMakespanWorseAlongChains(list<Task *> &seqSet, Tree *tree) {

    int cntrAdditionToSS = 0;
    int cntTries = 0;
    Task *root = tree->getRoot();
    list<Task *> frontier;


    for (Task *task: *tree->getTasks()) {
        if (task->isBroken()) {
            frontier.push_back(task);
        }
    }
    frontier.sort(cmp_Mem_nodecreasing);

    double minMakespan = assignToBestProcessors(tree, {});
    tree->cleanAssignedAndReassignFeasible();
    SeqSet optimalSeqSet = SeqSet(tree, minMakespan);

    int maxLengthOfChain = 0;
    while (!frontier.empty()) {
        // cout << "num free proc 4 " << Cluster::getFixedCluster()->getNumberFreeProcessors() << endl;
        tree->cleanAssignedAndReassignFeasible();
        // cout << "num free proc 5 " << Cluster::getFixedCluster()->getNumberFreeProcessors() << endl;
        // cout<<"try add new"<<endl;
        Task *potentialAddition = frontier.front();
        frontier.pop_front();
        int lengthOfChain = 0;

        while (potentialAddition->getChildren()->size() == 1) {
            //  cout<<"go along linear chain"<<endl;
            lengthOfChain++;
            potentialAddition->restoreEdge();
            potentialAddition = potentialAddition->getChildren()->at(0);
        }
        if (lengthOfChain > maxLengthOfChain) maxLengthOfChain = lengthOfChain;
        if (!potentialAddition->isRoot()) potentialAddition->restoreEdge();

        for (Task *child: *potentialAddition->getChildren()) {
            child->breakEdge();
        }

        double potentialMakespan = assignToBestProcessors(tree, {});
        cout << "potential: " << potentialMakespan << ", current min: " << " " << minMakespan << " equals? "
             << (potentialMakespan == minMakespan ? "y" : "n");// << endl;
        //  cout << "num free proc 1 " << Cluster::getFixedCluster()->getNumberFreeProcessors() << endl;
        cntTries++;

        if (potentialMakespan <= minMakespan) {
            //cout << " add!" << endl;
            cntrAdditionToSS++;
            minMakespan = potentialMakespan;
            //no breaking back the edge, we accept this improvement
            // instead add children of newly added task
            for (Task *child: *potentialAddition->getChildren()) {
                auto it = upper_bound(frontier.begin(), frontier.end(), child,
                                      cmp_Mem_nodecreasing);
                frontier.insert(it, child);
            }
        } else {
            tree->cleanAssignedAndReassignFeasible();
            // cout << " no add!" << endl;
            potentialAddition->breakEdge();
            for (Task *child: *potentialAddition->getChildren()) {
                child->restoreEdge();
            }
            //  cout << "num free proc 2 " << Cluster::getFixedCluster()->getNumberFreeProcessors() << endl;
            assignToBestProcessors(tree, {});
            //   cout << "num free proc 3 " << Cluster::getFixedCluster()->getNumberFreeProcessors() << endl;
            // no adding children
        }

        // cout << "#trees " << tree->HowmanySubtrees(false) << " added " << cntrAdditionToSS << endl;
    }
    // cout << "tried " << cntTries << " added to SS " << cntrAdditionToSS << " max # children " << maxLengthOfChain
    //<< endl;
    return "tried " + to_string(cntTries) + " added to SS " + to_string(cntrAdditionToSS);
}


string growSeqSetWhileImprovesMakespanWorseAroundRoot(list<Task *> &seqSet, Tree *tree, double part) {

    int cntrAdditionToSS = 0;
    int cntTries = 0;
    Task *root = tree->getRoot();
    list<Task *> frontier;


    for (Task *task: *tree->getTasks()) {
        if (task->isBroken()) {
            frontier.push_back(task);
        }
    }
    frontier.sort(cmp_Mem_nodecreasing);

    double minMakespan = assignToBestProcessors(tree, {});
    tree->cleanAssignedAndReassignFeasible();
    SeqSet optimalSeqSet = SeqSet(tree, minMakespan);

    int maxLengthOfChain = 0;
    while (!frontier.empty()) {
        // cout << "num free proc 4 " << Cluster::getFixedCluster()->getNumberFreeProcessors() << endl;
        tree->cleanAssignedAndReassignFeasible();
        // cout << "num free proc 5 " << Cluster::getFixedCluster()->getNumberFreeProcessors() << endl;
        // cout<<"try add new"<<endl;
        Task *potentialAddition = frontier.front();
        frontier.pop_front();

        if (!potentialAddition->isRoot()) potentialAddition->restoreEdge();

        for (Task *child: *potentialAddition->getChildren()) {
            child->breakEdge();
        }

        double potentialMakespan = assignToBestProcessors(tree, {});
        //cout << "potential: " << potentialMakespan << ", current min: " << " " << minMakespan << " equals? "
        //    << (potentialMakespan == minMakespan ? "y" : "n");// << endl;
        //  cout << "num free proc 1 " << Cluster::getFixedCluster()->getNumberFreeProcessors() << endl;
        cntTries++;
        if (potentialMakespan <= minMakespan || potentialAddition->getId() < tree->getSize() * part) {
            //cout << " add!" << endl;
            cntrAdditionToSS++;
            minMakespan = potentialMakespan;
            //no breaking back the edge, we accept this improvement
            // instead add children of newly added task
            for (Task *child: *potentialAddition->getChildren()) {
                auto it = upper_bound(frontier.begin(), frontier.end(), child,
                                      cmp_Mem_nodecreasing);
                frontier.insert(it, child);
            }
        } else {
            tree->cleanAssignedAndReassignFeasible();
            // cout << " no add!" << endl;
            potentialAddition->breakEdge();
            for (Task *child: *potentialAddition->getChildren()) {
                child->restoreEdge();
            }
            //  cout << "num free proc 2 " << Cluster::getFixedCluster()->getNumberFreeProcessors() << endl;
            assignToBestProcessors(tree, {});
            //   cout << "num free proc 3 " << Cluster::getFixedCluster()->getNumberFreeProcessors() << endl;
            // no adding children
        }

        // cout << "#trees " << tree->HowmanySubtrees(false) << " added " << cntrAdditionToSS << endl;
    }
    // cout << "tried " << cntTries << " added to SS " << cntrAdditionToSS << " max # children " << maxLengthOfChain
    //<< endl;
    return "tried " + to_string(cntTries) + " added to SS " + to_string(cntrAdditionToSS);
}

void growSeqSetWithUnfeasible(Task *task, list<Task *> &seqSet, Task *treeRoot) {
    if (task->getFeasibleProcessors()->empty()) { // && (buildParallelRootsFromSequentialSet(treeRoot, seqSet)).size() <
        //Cluster::getFixedCluster()->getNumberProcessors() - 1) {
        seqSet.push_back(task);
    }
    for (Task *child: *task->getChildren()) {
        growSeqSetWithUnfeasible(child, seqSet, treeRoot);
    }
}

string seqSetAndFeasSets(Tree *tree) {
    tree->cleanAssignedAndReassignFeasible();
    string result = "";
    list<Task *> sequentialSet;
    sequentialSet.clear();
    //sequentialSet.emplace_front(tree->getRoot());
    growSeqSetWithUnfeasible(tree->getRoot(), sequentialSet, tree->getRoot());
    auto parallelRoots = buildParallelRootsFromSequentialSet(tree->getRoot(), sequentialSet);
    breakPreparedEdges(tree->getRoot(), parallelRoots);

    if (tree->HowmanySubtrees(true) >= Cluster::getFixedCluster()->getNumberProcessors()) {
        throw "too many parallel roots after growing SeqSet: ";// + to_string(parallelRoots.size());
    }
    //result = growSeqSetWhileImprovesMakespan2(sequentialSet, tree);
    // result = growSeqSetWhileImprovesMakespanAllowWorse(sequentialSet, tree);
    //result = growSeqSetWhileImprovesMakespanWorseAlongChains(sequentialSet, tree);
    result = growSeqSetWhileImprovesMakespanWorseAroundRoot(sequentialSet, tree, 0.2);
    return result;

}

void assignCorrespondingTreeTasks(Tree *tree, Tree *qTree) {
    for (const auto &item: *tree->getTasks()) {
        assert(item->getAssignedProcessor() == NULL);
    }
    for (Task *qTask: *qTree->getTasks()) {
        Task *taskInTree = tree->getTask(qTask->getOtherSideId());
        vector<Task *> allTasksInSubtree = taskInTree->tasksInSubtreeRootedHere();
        for (Task *taskInSubtree: allTasksInSubtree) {
            taskInSubtree->setAssignedProcessor(qTask->getAssignedProcessor());
        }
    }
    for (const auto &item: *tree->getTasks()) {
        assert(item->getAssignedProcessor() != NULL);
    }
}

double assignToBestProcessors(Tree *tree, vector<Task *> newlyBroken, string chooseSubtreeAssign) {
    tree->cleanAssignedAndReassignFeasible();

    assert(Cluster::getFixedCluster()->getNumberProcessors() == Cluster::getFixedCluster()->getNumberFreeProcessors());
    clock_t time;
    time = clock();

    //TODO: is correct?
    if (newlyBroken.empty()) {
        for (auto &item: tree->getBrokenTasks()) {
            item->needsRecomputeMemReq = true;
        }
    }

    for (auto &item: newlyBroken) {
        item->needsRecomputeMemReq = true;
        Task *parent = item->getParent();
        while (parent != nullptr && !parent->isBroken()) {
            parent = parent->getParent();
        }
        if (parent != nullptr) {
            if (!parent->isBroken()) {
                parent = parent->getParent();
            }
            assert(parent->isBroken() == true);
            parent->needsRecomputeMemReq = true;
        }
    }

    for (auto &item: tree->getBrokenTasks()) {
        item->computeMinMemUnderlyingAndAssignFeasible(tree, false);
        assert(item->getMinMemUnderlying() != 0);
    }
    Tree *qTree = tree->BuildQtree();

    for (Task *task: *qTree->getTasks()) {
        if (task->getFeasibleProcessors()->empty()) {
            tree->getTask(task->getOtherSideId())->needsRecomputeMemReq = true;
            tree->getTask(task->getOtherSideId())->computeMinMemUnderlyingAndAssignFeasible(tree, false);
            task->setFeasibleProcessors(tree->getTask(task->getOtherSideId())->getFeasibleProcessors());
            if (task->getFeasibleProcessors()->empty()) {
                string exception =
                        "Task" + to_string(task->getOtherSideId()) + " has 0 feasible processors at beginning.";
                throw exception;
            }
        }
    }
    try {
        distributeProcessors(qTree, chooseSubtreeAssign);
    }
    catch (exception e) {
        Cluster::getFixedCluster()->freeAllBusyProcessors();
        return numeric_limits<double>::infinity();
    }
    assignCorrespondingTreeTasks(tree, qTree);
    double resultingMakespan = tree->getRoot()->getMakespanCostWithSpeeds(true, true);
    delete qTree;
    //  delete qTree1;

    timeForAssignment += (clock() - time);
    return resultingMakespan;
}

void removeProcessorFromAllFeasSets(Processor *processor, Tree *tree) {

    for (Task *task: *tree->getTasks()) {
        task->deleteFeasible(processor);
    }
}

Task *chooseSubtree(string subtreeChoiceCode, Tree *tree, vector<Task *> candidates) {
    // cout << "choose subtree" << endl;
    clock_t time;
    time = clock();
    int initialSize = candidates.size();
    if (subtreeChoiceCode == "LMW") {
        std::sort(candidates.begin(), candidates.end(), [](Task *a, Task *b) {
            return (a->getNodeWeight() * a->getMakespanWeight() >= b->getNodeWeight() * b->getMakespanWeight());
        });
        return candidates.front();
    } else if (subtreeChoiceCode == "CP") {
        Tree *qtree = tree->BuildQtree();
        vector<Task *> CriticalPath = buildCriticalPath(qtree->getRoot());

        auto candidatesContainOtherSideId = [candidates](Task *i) {
            //cout << "i: " << i->getId() << " " << i->getOtherSideId() << endl;
            for (Task *task: candidates) {
                //  cout << "task: " << task->getId() << " " << task->getOtherSideId() << endl;
                if (task->getId() == i->getOtherSideId()) {
                    // cout << "yes" << endl;
                    return true;
                }
            }
            return false;
        };

        Task *firstThatIsInCandidates = nullptr;

        if (CriticalPath.size() != 0) {
            vector<Task *>::iterator iterator = std::find_if(CriticalPath.begin(), CriticalPath.end(),

                                                             candidatesContainOtherSideId);
            firstThatIsInCandidates =
                    *iterator;
            if (iterator == CriticalPath.end()) {
                //cout << "not found in CP" << endl;
                firstThatIsInCandidates = *std::find_if(qtree->getTasks()->begin(),
                                                        qtree->getTasks()->end(),
                                                        candidatesContainOtherSideId);
            }

        } else {
            if (candidates.size() != 0) {
                //cout << "critical path size 0" << endl;
            }
            return candidates.front();
        }

        Task *correspondingTaskInTree = tree->getTask(firstThatIsInCandidates->getOtherSideId());
        delete qtree;
        timeChooseTree += (clock() - time);
        return correspondingTaskInTree;
    } else {
        timeChooseTree += (clock() - time);
        throw std::runtime_error("not implemented");
    }


}

vector<pair<int, vector<Task *>>> levelOrder(Task *root) {
    if (root == NULL)
        throw std::runtime_error("no root");;
    queue<Task *> parent_queue, child_queue;
    vector<pair<int, vector<Task *>>> tasksAndLevels;
    parent_queue.push(root);

    // Level 0 corresponds to the top of the tree i.e at the root level
    int level = 0;


    while (!parent_queue.empty() or !child_queue.empty()) {
        if (!parent_queue.empty()) {
            // cout << "Level " << level << ": ";
            tasksAndLevels.push_back({level, vector<Task *>()});
        }

        while (!parent_queue.empty()) {
            Task *node = parent_queue.front();
            parent_queue.pop();
            // cout << node->getId() << " ";
            tasksAndLevels.back().second.push_back(node);

            for (Task *child: *node->getChildren()) {
                if (!child->isBroken())
                    child_queue.push(child);
            }
        }
        // cout << endl;
        level++;

        if (!child_queue.empty()) {
            //  cout << "Level " << level << ": ";
            tasksAndLevels.push_back({level, vector<Task *>()});
        }


        while (!child_queue.empty()) {
            Task *node = child_queue.front();
            child_queue.pop();
            // cout << node->getId() << " ";
            tasksAndLevels.back().second.push_back(node);

            for (Task *child: *node->getChildren()) {
                if (!child->isBroken())
                    parent_queue.push(child);
            }

        }
        //  cout << endl;
        level++;
    }
    return tasksAndLevels;
}

Task *chooseTask(Task *root, Tree *tree, string nodeChoiceCode, string assignSubtreeChoiceCode) {
    clock_t time;
    time = clock();

    vector<Task *> candidates = buildCandidatesForNode(tree, nodeChoiceCode, root);

    timeChooseNode += (clock() - time);
    time = clock();
    pair<Task *, double> mins = findBestCutAmong(tree, candidates, assignSubtreeChoiceCode,
                                                 numeric_limits<double>::infinity());
    timeBestCutInNodeChoice += (clock() - time);

    return mins.first;


}

vector<Task *> buildCandidatesForNode(Tree *tree, const string &nodeChoiceCode, Task *root) {

    Tree *subtree = BuildSubtree(tree, root);
    vector<Task *> candidates;

    if (nodeChoiceCode == "M") {
        vector<Task *> eligibleTasks;
        vector<pair<int, vector<Task *>>> tasksAndLevels = levelOrder(subtree->getRoot());
        double maxLayer = tasksAndLevels.back().first;
        int middleLayer = ceil(maxLayer / 2);
        int lowerBound = maxLayer / 2;
        int upperBound = ceil(maxLayer / 2);

        for (const auto &item: tasksAndLevels) {
            if (item.first == middleLayer) {
                for (const auto &item1: item.second) {
                    if (item1->getChildren()->size() >= 2)
                        eligibleTasks.push_back(tree->getTask(item1->getOtherSideId()));
                }
            }
        }


        while (eligibleTasks.empty() && lowerBound >= 0) {
            lowerBound--;
            upperBound++;
            for (const auto &item: tasksAndLevels) {
                if (item.first <= upperBound && item.first >= lowerBound) {
                    for (const auto &item1: item.second) {
                        if (item1->getChildren()->size() >= 2)
                            eligibleTasks.push_back(tree->getTask(item1->getOtherSideId()));
                    }
                }
            }

        }
        candidates = eligibleTasks;

    } else if (nodeChoiceCode == "FFT") {
        for (Task *task: *subtree->getTasks()) {
            if (task->getChildren()->size() >= 2 && task->getId() != subtree->getRootId()) {
                candidates.push_back(tree->getTask(task->getOtherSideId()));
                break;
            }
        }
    } else if (nodeChoiceCode == "EX") {
        copy_if(subtree->getTasks()->begin(), subtree->getTasks()->end(), back_inserter(candidates),
                [](Task *i) {
                    return !i->isAnyChildBroken();
                });
    } else if (nodeChoiceCode == "FFTM") {
        int desiredSize = subtree->getTasks()->size() / 10;
        for (Task *task: *subtree->getTasks()) {
            if (task->getChildren()->size() >= 2 && candidates.size() < desiredSize) {
                candidates.push_back(tree->getTask(task->getOtherSideId()));
            }
        }
    }
    delete subtree;
    return candidates;
}

void chooseAssignSubtree(string parser, Tree *tree) {
    cout << "choose assign subtree" << endl;
    bool (*maxWi)(Task *, Task *) = [](Task *a, Task *b) {
        return a->getMakespanCost(false, true) < b->getMakespanCost(false, true);
    };

    bool (*maxDesc)(Task *, Task *) = [](Task *a, Task *b) {
        return a->getChildren()->size() < b->getChildren()->size();
    };
    vector<Task *> candidates = *tree->getTasks();
    int initialSize = tree->getTasks()->size();
    if (parser == "LMW") sort(candidates.begin(), candidates.end(), maxWi);
    else if (parser == "MD") sort(candidates.begin(), candidates.end(), maxDesc);
    else if (parser == "CP") {

        vector<Task *> CriticalPath = buildCriticalPath(tree->getRoot());
        vector<Task *> notOnCP;
        std::set_difference(candidates.begin(), candidates.end(),
                            CriticalPath.begin(), CriticalPath.end(),
                            std::inserter(notOnCP, notOnCP.begin()),
                            [](auto &a, auto &b) { return a->getId() < b->getId(); }

        );
        //std::back_inserter(notOnCP));
        candidates.clear();
        candidates.insert(candidates.end(), CriticalPath.begin(), CriticalPath.end());
        candidates.insert(candidates.end(), notOnCP.begin(), notOnCP.end());

        assert(candidates.size() == initialSize);
        assert(tree->getTasks()->size() == initialSize);
    } else
        throw std::runtime_error("not implemented");

    tree->setTasks(new vector<Task *>());
    for (Task *task: candidates) {
        tree->addTask(task);
    }

    assert(candidates.size() == initialSize);
    assert(tree->getTasks()->size() == initialSize);
    //std::for_each(candidates.begin(), candidates.end(), [](Task *n) { cout << n->getId() << ", "; });
    // cout << endl;
}


string
partitionHeuristics(Tree *tree, string subtreeChoiceCode, string nodeChoiceCode, string assignSubtreeChoiceCode,
                    int cutWhat) {
    timeForAssignment = 0;
    string result = " times: ";
    clock_t time;
    time = clock();

    Cluster::getFixedCluster()->freeAllBusyProcessors();

    tree->getRoot()->computeMinMemUnderlyingAndAssignFeasible(tree, false);

    double minMakespan = FirstCutSomeNodes(tree, assignSubtreeChoiceCode);

    if (minMakespan == std::numeric_limits<int>::max()) return "0 0 0 " + to_string(minMakespan);

    tree->cleanAssignedAndReassignFeasible();
    assert(Cluster::getFixedCluster()->getNumberFreeProcessors() ==
           Cluster::getFixedCluster()->getNumberProcessors());

    result += to_string((clock() - time) / CLOCKS_PER_SEC) + " ";
    time = clock();
    timeForAssignment = 0;
    cout << "first cut ready" << endl;

    switch (cutWhat) {
        case 0:
            cutSingleNodePerSubtreeUntilBestMakespan(tree, subtreeChoiceCode, nodeChoiceCode, assignSubtreeChoiceCode,
                                                     minMakespan, false);
            break;
        case 1:
            cutSingleNodeInAllSubtreesSimultaneously(tree, subtreeChoiceCode, nodeChoiceCode, assignSubtreeChoiceCode,
                                                     minMakespan);
            break;
        case 2:
            cutSingleNodePerSubtreeUntilBestMakespan(tree, subtreeChoiceCode, nodeChoiceCode, assignSubtreeChoiceCode,
                                                     minMakespan, true);
            break;
        default:
            throw std::runtime_error("no partitioning method");
    }

    tree->cleanAssignedAndReassignFeasible();
    for (Processor *item: Cluster::getFixedCluster()->getProcessors()) {
        if (item->isBusy) {
            item->isBusy = false;
            item->setAssignedTaskId(-1);
            item->setAssignedTask(NULL);
            item->setOccupiedMemorySize(0);
        }
    }
    double makespan = assignToBestProcessors(tree, {}, assignSubtreeChoiceCode);
    if (minMakespan != makespan) {
        cout << "MS wrong " << minMakespan << " " << makespan << endl;
        /*cout << "assignment: " << endl;
        for (const auto &item: Cluster::getFixedCluster()->getProcessors()) {
            cout << item->getAssignedTaskId() << " w mem req "
                 << (item->getAssignedTask() == NULL ? "0" : to_string(item->getAssignedTask()->getMinMemUnderlying()))
                 << " on " << item->getMemorySize() << " and " << item->getProcessorSpeed()
                 << endl;
        }
        cout << "end! res makespan " << makespan << endl; */
    }
    result += to_string((clock() - time) / CLOCKS_PER_SEC) + " ";
    return result;
}

void cutSingleNodePerSubtreeUntilBestMakespan(Tree *tree, string &subtreeChoiceCode, string &nodeChoiceCode,
                                              string &assignSubtreeChoiceCode, double &minMakespan, bool cutMultiple) {
    vector<Task *> subtreeCandidates = tree->getBrokenTasks();
    // tree->HowmanySubtrees(false);
    unsigned int numberAvailableProcessors =
            Cluster::getFixedCluster()->getNumberProcessors() - tree->HowmanySubtrees(true);
    while (!subtreeCandidates.empty() && numberAvailableProcessors != 0) {
        numberAvailableProcessors =
                Cluster::getFixedCluster()->getNumberProcessors() - tree->HowmanySubtrees(true);
        Task *subtree = chooseSubtree(subtreeChoiceCode, tree, subtreeCandidates);
        //cout << "choosen subtree " << subtree->getId() << ", cand size " << subtreeCandidates.size() << endl;
        Task *bestTask = chooseTask(subtree, tree, nodeChoiceCode, assignSubtreeChoiceCode);
        if (bestTask != NULL) {
            vector<Task *> freshlyBroken = bestTask->breakNBiggestChildren(
                    numberAvailableProcessors);
            double currentMakespan = assignToBestProcessors(tree, freshlyBroken, assignSubtreeChoiceCode);
            if (currentMakespan <= minMakespan) {
                minMakespan = currentMakespan;
                for (Task *candidate: freshlyBroken) {
                    if (candidate->getChildren()->size() >= 2) {
                        subtreeCandidates.push_back(candidate);
                    }
                }
            } else {
                bestTask->restoreBrokenChildren();
                assignToBestProcessors(tree, freshlyBroken, assignSubtreeChoiceCode);
                subtreeCandidates.erase(find(subtreeCandidates.begin(), subtreeCandidates.end(), subtree));
            }
        } else {
            const vector<Task *>::iterator &position = find(subtreeCandidates.begin(), subtreeCandidates.end(),
                                                            subtree);
            subtreeCandidates.erase(position);

        }
        numberAvailableProcessors =
                Cluster::getFixedCluster()->getNumberProcessors() - tree->HowmanySubtrees(true);
    }

}

void cutSingleNodeInAllSubtreesSimultaneously(Tree *tree, string &subtreeChoiceCode, string &nodeChoiceCode,
                                              string &assignSubtreeChoiceCode, double &minMakespan) {
    vector<Task *> subtreeCandidates = tree->getBrokenTasks();
    double currentMakespan = minMakespan;
    vector<Task *> freshlyBroken;
    vector<pair<Task *, vector<Task *>>> subtreesAndTheirCandidateNodes;
    for (auto &item: tree->getBrokenTasks()) {
        vector<Task *> candidates = buildCandidatesForNode(tree, nodeChoiceCode, item);
        subtreesAndTheirCandidateNodes.push_back(make_pair(item, candidates));
    }

    do {
        freshlyBroken.resize(0);
        for (auto &item: subtreesAndTheirCandidateNodes) {
            if (item.second.size() != 0) {
                vector<Task *> newBrokenChildren = item.second.at(0)->breakNBiggestChildren(
                        Cluster::getFixedCluster()->getNumberProcessors() - tree->HowmanySubtrees(true));
                freshlyBroken.insert(freshlyBroken.end(), newBrokenChildren.begin(), newBrokenChildren.end());
                item.second.erase(item.second.begin());
            }
        }
        if (freshlyBroken.size() == 0) return;

        currentMakespan = assignToBestProcessors(tree, freshlyBroken, assignSubtreeChoiceCode);
        if (currentMakespan < minMakespan) {
            for (auto &newlyCutSubtreeRoot: freshlyBroken) {
                auto iterator = std::find_if(subtreesAndTheirCandidateNodes.begin(),
                                             subtreesAndTheirCandidateNodes.end(),
                                             [newlyCutSubtreeRoot](const pair<Task *, vector<Task *>> pair) {
                                                 return pair.first->getId() == newlyCutSubtreeRoot->getId();
                                             });
                assert(iterator == subtreesAndTheirCandidateNodes.end());
                vector<Task *> candidates = buildCandidatesForNode(tree, nodeChoiceCode, newlyCutSubtreeRoot);
                subtreesAndTheirCandidateNodes.push_back(make_pair(newlyCutSubtreeRoot, candidates));
            }
        } else {
            for (auto &freshlyBrokenTask: freshlyBroken) {
                freshlyBrokenTask->restoreEdge();
            }
        }
    } while (true);

}

double FirstCutSomeNodes(Tree *tree, string assignSubtreeChoiceCode) {

    vector<Task *> bestFFT = buildCandidatesForNode(tree, "FFT", tree->getRoot());
    vector<Task *> bestMiddle = buildCandidatesForNode(tree, "M", tree->getRoot());

    vector<Task *> candidates(bestFFT);
    candidates.insert(candidates.end(), bestMiddle.begin(), bestMiddle.end());
    assert(candidates.size() == bestMiddle.size() + bestFFT.size());

    pair<Task *, double> mins = findBestCutAmong(tree, candidates, assignSubtreeChoiceCode);

    vector<Task *> newlyBroken;
    if (mins.first != nullptr) {
        newlyBroken = mins.first->breakNBiggestChildren(
                Cluster::getFixedCluster()->getNumberFreeProcessors());
    }
    double computedMs = assignToBestProcessors(tree, newlyBroken, assignSubtreeChoiceCode);
    assert(computedMs == mins.second);
    return computedMs;
}

Task *CutTaskWithMaxImprovement(Tree *tree, string assignSubtreeChoiceCode) {
    vector<Task *> candidates;
    std::copy_if(tree->getTasks()->begin(), tree->getTasks()->end(), std::back_inserter(candidates), [](Task *i) {
        return i->getChildren()->size() >= 2;
    });
    pair<Task *, double> mins = findBestCutAmong(tree, candidates, assignSubtreeChoiceCode);
    if (mins.first != nullptr) {
        mins.first->breakNBiggestChildren(
                Cluster::getFixedCluster()->getNumberFreeProcessors());
    }

    return mins.first;

}


pair<Task *, double>
findBestCutAmong(Tree *tree, vector<Task *> candidates, string assignSubtreeChoiceCode, double initMS) {

    tree->cleanAssignedAndReassignFeasible();
    Cluster::getFixedCluster()->freeAllBusyProcessors();

    double minMakespan = initMS == -1 ? assignToBestProcessors(tree, {}, assignSubtreeChoiceCode) : initMS;
    Task *taskMinMakespan = nullptr;
    //cout << "\t task candidates size " << candidates.size() << endl;
    for (Task *task: candidates) {
        //cout << "\t try task " << task->getId() << endl;
        task = tree->getTask(task->getId());
        if (task->getChildren()->size() >= 2 && !task->isAnyChildBroken()) {
            vector<Task *> newlyBroken = task->breakNBiggestChildren(
                    Cluster::getFixedCluster()->getNumberFreeProcessors());
            double currentMakespan = assignToBestProcessors(tree, newlyBroken, assignSubtreeChoiceCode);
            if (currentMakespan < minMakespan) {
                minMakespan = currentMakespan;
                taskMinMakespan = task;
            }
            task->restoreBrokenChildren();
            tree->cleanAssignedAndReassignFeasible();

        }
    }
    return make_pair(taskMinMakespan, minMakespan);
}

void Tree::mergeLinearChains() {
    Task *currentToBeExplored;
    list<Task *> toBeExplored;

    toBeExplored.push_back(this->getRoot());

    while (!toBeExplored.empty()) {
        currentToBeExplored = toBeExplored.front();
        toBeExplored.pop_front();
        // cout << currentToBeExplored->getId() << " w #ch " << currentToBeExplored->getChildren()->size() << endl;
        if (currentToBeExplored->getChildren()->size() == 1) {
            mergeTaskToOnlyChild(currentToBeExplored);
            toBeExplored.push_back(currentToBeExplored);
        } else {
            toBeExplored.insert(toBeExplored.end(), currentToBeExplored->getChildren()->begin(),
                                currentToBeExplored->getChildren()->end());
        }
    }

    for (const auto &item: *this->getTasks()) {
        assert(item->getChildren()->size() > 1 || item->getChildren()->size() == 0);
    }
    //TODO RENUMBER HERE OR NOT?
    // this->renumberAllTasks();
}

void Tree::renumberAllTasks() {
    int currentIdByOrder = 1;
    for (const auto &item: *this->getTasks()) {
        item->setId(currentIdByOrder);
        for (const auto &child: *item->getChildren()) {
            child->setParentId(item->getId());
        }
        currentIdByOrder++;
    }
}