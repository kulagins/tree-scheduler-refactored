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

//INFO: paths are necessary for Clion. In case of problems please ping Svetlana.
#include "tMaxHeap.cpp"
#include <bits/stdc++.h>


void freeProcessorIfAvailable(Task *task) {
    if (task->getAssignedProcessor() != NULL) {
        task->getAssignedProcessor()->setAssignedTask(NULL);
        task->getAssignedProcessor()->isBusy = false;
        task->getAssignedProcessor()->setOccupiedMemorySize(0);
        task->getAssignedProcessor()->setAssignedTaskId(-1);
        task->setAssignedProcessor(NULL);
    }
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

void distributeProcessors(Tree *qTree) {
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
    // make_heap(taskHeap->begin(), taskHeap->end(), TMaxHeap());
    SiftInfTmaxUpPreserveOrder(taskHeap);
    while (taskHeap->size() > 0) {
        // pop_heap(taskHeap->begin(), taskHeap->end(), TMaxHeap());
        //TODO CHECK!
        Task *task = taskHeap->front();
        taskHeap->erase(taskHeap->begin());
        // taskHeap->pop_back();
        Processor *pFast = task->getFastestFeasibleProcessor();
        pFast->assignTask(task);
        for (int i = 0; i < taskHeap->size(); i++) {
            taskHeap->at(i)->deleteFeasible(pFast);
            taskHeap->at(i)->updateTMax();
            SiftInfTmaxUpPreserveOrder(taskHeap);
        }
        SiftInfTmaxUpPreserveOrder(taskHeap);

    }
    delete taskHeap;
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

    for (auto &item: newlyBroken) {
        Task *parent = item->getParent();
        while (parent != nullptr && !parent->isBroken()) {
            parent = parent->getParent();
        }
        if (parent != nullptr){
            if (parent->getMinMemUnderlying() != 0) {
                //cout << " needs to recompute minMem on " << parent->getId() << ", ";
            }
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
            string exception = "Task" + to_string(task->getOtherSideId()) + " has 0 feasible processors at beginning.";
            throw exception;
        }
    }
    try {
        distributeProcessors(qTree);
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
    //cout << "choose task" << endl;


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
            if (task->getChildren()->size() >= 2) {
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


double CutTaskWithMaxImprovementHeuristicChoice(Tree *tree, string assignSubtreeChoiceCode) {

    vector<Task *> bestFFT = buildCandidatesForNode(tree, "FFTM", tree->getRoot());
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
   /* for (const auto &item: newlyBroken) {
        cout << item->getId() << ", ";
    }
    cout << endl;
    cout << "try cut children of task " << mins.first->getId() << " got ms " << computedMs << " on #subtrees "
         << tree->HowmanySubtrees(true) << endl;
    //assert(computedMs == mins.second);
    */
    return computedMs;
}


pair<Task *, double>
findBestCutAmong(Tree *tree, vector<Task *> candidates, string assignSubtreeChoiceCode, double initMS) {
    tree->cleanAssignedAndReassignFeasible();
    Cluster::getFixedCluster()->freeAllBusyProcessors();
    double minMakespan = initMS == -1 ? assignToBestProcessors(tree, {}, assignSubtreeChoiceCode) : initMS;
    Task *taskMinMakespan = nullptr;

    for (Task *task: candidates) {
        task = tree->getTask(task->getId());
        if (task->getChildren()->size() >= 2 && !task->isAnyChildBroken()) {
            vector<Task *> newlyBroken = task->breakNBiggestChildren(
                    Cluster::getFixedCluster()->getNumberFreeProcessors());
            double currentMakespan = assignToBestProcessors(tree, newlyBroken, assignSubtreeChoiceCode);
              if (currentMakespan < minMakespan) {
                minMakespan = currentMakespan;
                taskMinMakespan = task;
            }
            tree->cleanAssignedAndReassignFeasible();
            task->restoreBrokenChildren();
        }
    }
    return make_pair(taskMinMakespan, minMakespan);
}



string
partitionHeuristics(Tree *tree, string subtreeChoiceCode, string nodeChoiceCode, string assignSubtreeChoiceCode) {
    timeForAssignment = 0;
    string result = "times: ";
    clock_t time;
    time = clock();
    //tree->cleanAssignedAndReassignFeasible();
    // Cluster::getFixedCluster()->sortProcessorsByMemSize();
    Cluster::getFixedCluster()->freeAllBusyProcessors();
    tree->getRoot()->computeMinMemUnderlyingAndAssignFeasible(tree, false);

    double minMakespan = CutTaskWithMaxImprovementHeuristicChoice(tree, assignSubtreeChoiceCode);

    if (minMakespan == std::numeric_limits<int>::max()) return "0 0 0 " + to_string(minMakespan);

    tree->cleanAssignedAndReassignFeasible();
    assert(Cluster::getFixedCluster()->getNumberFreeProcessors() ==
           Cluster::getFixedCluster()->getNumberProcessors());

    result += to_string((clock() - time) / CLOCKS_PER_SEC) + " ";
    time = clock();
    timeForAssignment = 0;
    cout << "first cut ready" << endl;

    vector<Task *> subtreeCandidates = tree->getBrokenTasks();
    while (!subtreeCandidates.empty() && Cluster::getFixedCluster()->getNumberFreeProcessors() != 0) {

        Task *subtree = chooseSubtree(subtreeChoiceCode, tree, subtreeCandidates);
        // cout << "subtree " << to_string(subtree->getId()) << endl;
        Task *bestTask = chooseTask(subtree, tree, nodeChoiceCode, assignSubtreeChoiceCode);
        //cout << "chose " << (bestTask == NULL ? "nothing" : to_string(bestTask->getId())) << endl;
        if (bestTask != NULL) {
            //  bestTask = tree->getTask(bestTask->getId());
            vector<Task *> freshlyBroken = bestTask->breakNBiggestChildren(
                    Cluster::getFixedCluster()->getNumberProcessors() - tree->HowmanySubtrees(true));
            double currentMakespan = assignToBestProcessors(tree, freshlyBroken, assignSubtreeChoiceCode);
            if (currentMakespan <= minMakespan) {
                //    cout << "1" << endl;
                //cout << "smaller!" << endl;
                minMakespan = currentMakespan;
                for (Task *candidate: freshlyBroken) {
                    if (candidate->getChildren()->size() >= 2) {
                        subtreeCandidates.push_back(candidate);
                    }
                }
            } else {
                //   cout << "2" << endl;
                bestTask->restoreBrokenChildren();
                subtreeCandidates.erase(find(subtreeCandidates.begin(), subtreeCandidates.end(), subtree));

            }
        } else {
            //  cout << "3" << endl;
            const vector<Task *>::iterator &position = find(subtreeCandidates.begin(), subtreeCandidates.end(),
                                                            subtree);
            subtreeCandidates.erase(position);

        }
        //   cout << "candidate ready, # trees " << tree->HowmanySubtrees(false) << "  candidates: "
        //  << subtreeCandidates.size() << endl;
        // tree->HowmanySubtrees(false);
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
    result += to_string((clock() - time) / CLOCKS_PER_SEC) + " ";
    return result + to_string(timeForAssignment / CLOCKS_PER_SEC) + " " + to_string(timeChooseTree / CLOCKS_PER_SEC) +
           " " +
           to_string(timeChooseNode / CLOCKS_PER_SEC)
           + " " + to_string(timeBestCutInNodeChoice / CLOCKS_PER_SEC) +
           // to_string(assignToBestProcessors(tree, {}, assignSubtreeChoiceCode)) +
           " ";
    // delete subtreeCandidates;
}

