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
#include "../include/heuristics.h"
#include "../include/lib-io-tree-free-methods.h"
#include "../include/cluster.h"
#include <lib-io-tree-minmem.h>

//#include <omp.h>
bool cmp_noincreasing(Task *a, Task *b) {
    return (a->getMakespanCost(true, false) >= b->getMakespanCost(true, false));
};

bool cmp_nodecreasing(Task *a, Task *b) { return (a->getMakespanCost(true, false) < b->getMakespanCost(true, false)); };

bool cmp_noIn_noCommu(Task *a, Task *b) {
    return (a->getMakespanCost(false, false) >= b->getMakespanCost(false, false));
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
void GetTwoLargestElementTypetwo(T container, U &Largest, U &secondLargest) {
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

void GetTwoLargestElementTypethree(vector<Task *> *container, vector<Task *>::iterator &Largest,
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

void GetTwoSmallestElement(list<Task *> *container, list<Task *>::iterator &Smallest,
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
    sequentialLength = smallestMS_iter - makespansOfSplittings.begin();;
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
    int limitToPartitioning = limit==-1? Cluster::getFixedCluster()->getNumberProcessors():limit;
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

void popSmallestRootsToFitToCluster(list<Task *> &parallelRoots, unsigned long amountSubtrees, int limit) {
    int limitToPartitioning = limit==-1? Cluster::getFixedCluster()->getNumberProcessors():limit;
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
                        SF_now, -1); //SF_now will be modified in SplitSubtrees, it represents the length of sequential part, 0 means the subtree no need to partition

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

        if (tree->getTask(currentQNode->getOtherSideId())->isBroken() == true) { //this subtree has not been merged yet
            children = currentQNode->getChildren();
            if (children->empty()) { //this is a leaf node
                children = currentQNode->getParent()->getChildren();
                if (children->size() == 2) {
                    increase = children->front()->getMakespanWeight() + children->back()->getMakespanWeight() -
                               currentQNode->getParent()->getParallelPart();
                } else if (children->size() == 1) {
                    increase = -currentQNode->getEdgeWeight() / homogeneousBandwidth;
                } else {

                    GetTwoLargestElementTypetwo(children, LargestNode,
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
                    GetTwoLargestElementTypetwo(children, LargestNode,
                                                secondLargest); //no-increasing, communication counted

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
        //cout<<"   increase in MS(r) estimated: "<<(*smallest_iter).first<<endl;

        children = currentQNode->getChildren();
        if (children->empty()) {
            leaf = true;
        }
        Processor *currentQNodeProcessor = tree->getTask(currentQNode->getOtherSideId())->getAssignedProcessor();
        Processor *parentProcessor = tree->getTask(currentQNode->getParent()->getOtherSideId())->getAssignedProcessor();
        double memorySize;
        if (currentQNodeProcessor != nullptr
            && parentProcessor != nullptr) {
            //both are assigned to a processor, choose a bigger one
            processorToMergeTo = currentQNodeProcessor->getMemorySize() > parentProcessor->getMemorySize() ?
                                 currentQNodeProcessor : parentProcessor;
            memorySize = processorToMergeTo->getMemorySize();
        } else if (currentQNodeProcessor == nullptr && parentProcessor != nullptr) {
            processorToMergeTo = parentProcessor;
            memorySize = parentProcessor->getMemorySize();
        } else if (currentQNodeProcessor != nullptr &&
                   parentProcessor == nullptr) {
            processorToMergeTo = currentQNodeProcessor;
            memorySize = currentQNodeProcessor->getMemorySize();
        } else if (currentQNodeProcessor == nullptr
                   || parentProcessor == nullptr) {
            //both unassigned, skip
            memorySize = 0;
        }
        double requiredMemory = 0;
        if (CheckMemory == true) {
            memoryEnough = tree->MemoryEnough(currentQNode->getParent(), currentQNode, leaf, memorySize,
                                              requiredMemory);
        } else {
            memoryEnough = true;
        }

        if (memoryEnough == true) {
            feasible = true;
            smallestNode = currentQNode;
            if(processorToMergeTo!= nullptr){
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
        //cout<<"shortage "<<shortage<<endl;
        temp = Qtreeobj->getRoot()->getMakespanCost(true, true); //initilize ms
        temp = this->getRoot()->getMakespanCost(true, true);     //update ms

        //    cout<<"merge get free proc "<<Cluster::getFixedCluster()->getFirstFreeProcessor()->getMemorySize()<<endl;
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
                    processorToMergeTo->assignTask(this->getTask(parent->getOtherSideId()));
                    this->getTask(nodeone->getOtherSideId())->setAssignedProcessor(nullptr);
                    this->getTask(nodetwo->getOtherSideId())->setAssignedProcessor(nullptr);

                } else {
                    //cout<<"Merge node "<<node_smallest_increase->getOtherSideId()<<endl;
                    node_smallest_increase->mergeToParent();
                    shortage--;
                    this->getTask(node_smallest_increase->getOtherSideId())->restoreEdge();
                    processorToMergeTo->assignTask(this->getTask(parent->getOtherSideId()));
                    this->getTask(node_smallest_increase->getOtherSideId())->setAssignedProcessor(nullptr);
                }
            } else {
                //cout<<"Merge node "<<node_smallest_increase->getOtherSideId()<<endl;
                node_smallest_increase->mergeToParent();
                shortage--;
                this->getTask(node_smallest_increase->getOtherSideId())->restoreEdge();
                processorToMergeTo->assignTask(this->getTask(node_smallest_increase->getParent()->getOtherSideId()));
                this->getTask(node_smallest_increase->getOtherSideId())->setAssignedProcessor(nullptr);
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
Tree::MergeV2(unsigned int num_subtrees, unsigned int processor_number, double const memory_size, bool CheckMemory) {
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
                GetTwoSmallestElement(&Llist, smallest, secondSmallest);
            }

            if ((*smallest)->isLeaf()) {
                leaf = true;
            } else {
                leaf = false;
            }
            double memoryRequired;
            if (CheckMemory == true) {
                memoryCheckPass = this->MemoryEnough((*smallest)->getParent(), (*smallest), leaf, memory_size,
                                                     memoryRequired);
            } else {
                memoryCheckPass = true;
            }

            if (memoryCheckPass == false) {
                if ((*secondSmallest)->isLeaf()) {
                    leaf = true;
                } else {
                    leaf = false;
                }

                memoryCheckPass = this->MemoryEnough((*secondSmallest)->getParent(), *secondSmallest, leaf,
                                                     memory_size, memoryRequired);
                if (memoryCheckPass == true) {
                    currentNode = *secondSmallest;
                    DeadBreak = false;
                } else {
                    Llist.erase(secondSmallest);
                }
            } else {
                currentNode = *smallest;
                DeadBreak = false;
            }
        } while ((memoryCheckPass == false) && (!Llist.empty()));


        if (DeadBreak == true && firstTime == true) {
            Llist.clear();
            firstTime = false;
            goto CheckOnCritical;
        }

        if (DeadBreak == true) {
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
                GetTwoLargestElementTypethree(children, largestNode,
                                              secondLargest); //node_i is the largest, in terms of W
                temp = min(
                        (*largestNode)->getSequentialPart() - (*secondLargest)->getEdgeWeight() / homogeneousBandwidth,
                        (*secondLargest)->getSequentialPart() - (*largestNode)->getEdgeWeight() / homogeneousBandwidth);
                if (temp > decrease) {
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
                        if (temp > decrease) {
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

    //cout<<"   working on subtree ";
    do {
        currentNode = tree->getTask(SubtreeT->getOtherSideId());
        SubtreeT = SubtreeT->getParent();
        //cout<<"   "<<SubtreeT->getOtherSideId()<<"{ "<<endl;
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
                        //cout<<"    "<<(*it)->getId()<<" W_i "<<(*it)->getSequentialPart()<<", MS(t) "<<MS_t<<", W_t "<<W_t<<", MS_tj "<<(*it)->getParallelPart()<<endl;
                        tempQue.push_back((*it));
                        temp = min((*it)->getSequentialPart(),
                                   MS_t - W_t - (*it)->getEdgeWeight() / homogeneousBandwidth -
                                   (*it)->getParallelPart());
                        if (temp > decrease_othersubtrees) {
                            decrease_othersubtrees = temp;
                            output_node = (*it);
                        }
                    }
                }
            }
        } while (!currentNode->isBroken());
        //cout<<"   }"<<endl;
    } while (subtreeRoot->getId() != tree->getRootId());
    //cout<<endl;

    if (decrease_othersubtrees >= 0) {
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
//    cout<<"qtree initially "<<endl;
    //  Qtreeobj->Print(cout);
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

        //cout<<"Idle processor now: "<<idleProcessors<<endl;
        MSReduced = EstimateDecrease(idleProcessors, this, &CriticalPath, &onLastSubtree, &node_i, &node_j);

        if (MSReduced == true) {
            if (onLastSubtree == false) {
                //cout<<"cut edge "<<node_i->getId()<<endl;
                node_i->breakEdge(); //C<-C\cup C_k
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
        } else {
            break;
        }
    }
    //  cout<<"qtree after> "<<endl;
    //  Qtreeobj->Print(cout);

    delete Qtreeobj;

    MS_now = this->getRoot()->getMakespanCost(true, true);
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

int MemoryCheck(Tree *tree, io_method_t method) {

    vector<Task *> subtreeRoots;
    Task *subtreeRoot, *currentnode;
    tree->getRoot()->breakEdge();


    //subtreeRoots = tree->getBrokenTasks();
    for (unsigned int i = tree->getSize(); i >= 1; --i) {
        currentnode = tree->getTask(i);
        if (currentnode->isBroken()) {
            //cout<<i<<" ";
            subtreeRoots.push_back(currentnode);
        }
    }

    double maxoutD, memory_required;
    schedule_traversal *schedule_f = new schedule_traversal();
    unsigned int com_freq;
    vector<unsigned int> BrokenEdgesID;
    while (!subtreeRoots.empty()) {
     //   tree->HowmanySubtrees(false);
        subtreeRoot = subtreeRoots.back();
        subtreeRoots.pop_back();
        Processor *biggestFreeProcessor = Cluster::getFixedCluster()->getFirstFreeProcessorOrSmallest();
        if (Cluster::getFixedCluster()->hasFreeProcessor()) {
            biggestFreeProcessor->assignTask(subtreeRoot);
        }
        double currentMemorySize = biggestFreeProcessor->getMemorySize();
        cout<<"using proc of memory "<<biggestFreeProcessor->getMemorySize()<<"for task "<<subtreeRoot->getId()<<endl;

        Tree *subtree = BuildSubtree(tree, subtreeRoot);
        // cout<<"subtree "<<subtree->getSize();
        // subtree->Print(cout);
        maxoutD = MaxOutDegree(subtree, true);
        schedule_f->clear();
        MinMem(subtree, maxoutD, memory_required, *schedule_f, true);

        cout << "Subtree " << subtreeRoot->getId() << " needs memory " << memory_required <<endl;
        if (memory_required > currentMemorySize) {
            //        cout<<", larger than what is available: "<<currentMemorySize<<endl;

            int *schedule_copy = copyScheduleBackwards(schedule_f);

            switch (method) {
                case FIRST_FIT:
                    try {
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
                    IOCounter(subtree, schedule_copy, false, true,
                              com_freq, &BrokenEdgesID,
                              LARGEST_FIT);
                    break;
                case IMMEDIATELY:
                    Immediately(subtree, schedule_copy, currentMemorySize, &BrokenEdgesID);
                    break;

                default:
                    break;
            }
            delete[] schedule_copy;
        } else {
            //cout << "we are in else in mem check" << endl;
            if (Cluster::getFixedCluster()->hasFreeProcessor()) {
                biggestFreeProcessor->setOccupiedMemorySize(memory_required);
            }
        }
        //cout<<endl;
        Cluster::getFixedCluster()->assignTasksForIds(tree);
        delete subtree;
    }
    delete schedule_f;

    for (vector<unsigned int>::iterator iter = BrokenEdgesID.begin(); iter != BrokenEdgesID.end(); ++iter) {
        tree->getTask(*iter)->breakEdge();
    }
    // for (Processor *p: Cluster::getFixedCluster()->getProcessors()) {
    //     if(p->getAssignedTaskId()!=)
    //      p->assignTask(tree->getTask(p->getAssignedTaskId()));
    //   }
    // Cluster::getFixedCluster()->printProcessors();
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
    Processor *currentProcessor = cluster->getFirstFreeProcessor();


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
            currentProcessor = cluster->getFirstFreeProcessor();
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
