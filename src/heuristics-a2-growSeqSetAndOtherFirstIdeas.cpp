//
//  heuristics.cpp
//  memCom2
//
//  Created by changjiang GOU on 11/05/2018.
//  Copyright © 2018 ROMA. All rights reserved.
//


#include <list>
#include <algorithm>
#include "../include/heuristics.h"
#include "../include/lib-io-tree-free-methods.h"

//INFO: paths are necessary for Clion. In case of problems please ping Svetlana.

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
