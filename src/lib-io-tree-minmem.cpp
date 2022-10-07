/*
 *  lib-io-tree-minmem.cpp
 *  lib-io-tree
 *
 *  Created by defbond on 8/10/10.
 *  Copyright 2010 LIP/ENS-Lyon. All rights reserved.
 *
 */
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <set>
#include <sys/time.h>
#include <sys/resource.h>


#include <iostream>
#include <fstream>
#include <limits>


#include "../include/lib-io-tree.h"
#include "../include/lib-io-tree-minmem.h"

void fillLowestAndThreeUpperLevels(Tree *treeToBeChanged, vector<Task *> &lowestLevel,
                                   set<Task *> &firstLevel, set<Task *> &secondLevel, set<Task *> &thirdLevel);

void explore(Task *node, double available_memory, list<Task *> *L_init, schedule_traversal *S_init, double &cut_value,
             list<Task *> &min_sub_cut, schedule_traversal &sub_schedule, double &Mpeak, int quiet, int depth) {
    //cout<<"explore node id "<<node->getId()<<" avail mem "<<available_memory<<endl;
    /* if node is unreachable, return +infty */
    if (node->getCost() > available_memory) {

        Mpeak = node->getCost();
        cut_value = numeric_limits<double>::infinity();
        return;
    }


    /* if this is a leaf, return 0 */
    if (node->isLeaf()) {

        sub_schedule.push_back(node->getId());
        cut_value = 0;
        Mpeak = numeric_limits<double>::infinity();
        return;
    }

    if (L_init != NULL) {
        if (L_init->size() > 0) {
            min_sub_cut.assign(L_init->begin(), L_init->end());
            sub_schedule.assign(S_init->begin(), S_init->end());
        } else {
            sub_schedule.push_back(node->getId());
            /* place every child in the candidate nodes*/
            min_sub_cut.assign(node->getChildren()->begin(), node->getChildren()->end());
        }
    } else {
        sub_schedule.push_back(node->getId());
        /* place every child in the candidate nodes*/
        min_sub_cut.assign(node->getChildren()->begin(), node->getChildren()->end());
    }

    list<Task *> *candidates = new list<Task *>(min_sub_cut);
    cut_value = 0;
    for (list<Task *>::iterator current_node = candidates->begin(); current_node != candidates->end(); ++current_node) {
        (*current_node)->Mavail = available_memory;
        for (list<Task *>::iterator other_nodes = min_sub_cut.begin();
             other_nodes != min_sub_cut.end(); ++other_nodes) {
            if ((*other_nodes)->getId() != (*current_node)->getId()) {
                (*current_node)->Mavail -= (*other_nodes)->getEdgeWeight();
            }
        }
        cut_value += (*current_node)->getEdgeWeight();
    }

    while (!candidates->empty()) {


        for (list<Task *>::iterator current_node = candidates->begin();
             current_node != candidates->end(); ++current_node) {
            double m_j;
            list<Task *> Lj;
            schedule_traversal Sj;

            explore(*current_node, (*current_node)->Mavail, NULL, NULL, m_j, Lj, Sj, (*current_node)->Mpeak, quiet,
                    depth + 1);
            //cout<<"Task: "<<(*current_node)->getId()<<" MJ: "<<m_j<<endl;
            if (m_j <= (*current_node)->getEdgeWeight()) {

                min_sub_cut.remove(*current_node);
                min_sub_cut.splice(min_sub_cut.end(), Lj);
                sub_schedule.splice(sub_schedule.end(), Sj);
            } else {

            }
        }

        candidates->clear();
        for (list<Task *>::iterator current_node = min_sub_cut.begin();
             current_node != min_sub_cut.end(); ++current_node) {
            (*current_node)->Mavail = available_memory;
            for (list<Task *>::iterator other_nodes = min_sub_cut.begin();
                 other_nodes != min_sub_cut.end(); ++other_nodes) {
                if ((*other_nodes)->getId() != (*current_node)->getId()) {
                    (*current_node)->Mavail -= (*other_nodes)->getEdgeWeight();
                }
            }

            //add this node to candidates
            //			if(!quiet){cerr<<spacing<<"Mpeak = "<<(*current_node)->Mpeak<<" Mavail = "<<(*current_node)->Mavail<<endl;}

            if ((*current_node)->Mavail >= (*current_node)->Mpeak) {
                candidates->push_back(*current_node);
                //				if(!quiet){cerr<<spacing<<"node "<<(*current_node)->getId()<<" kept"<<endl;}
            }
        }
    }


    cut_value = 0;
    Mpeak = numeric_limits<double>::infinity();
    for (list<Task *>::iterator current_node = min_sub_cut.begin(); current_node != min_sub_cut.end(); ++current_node) {
        cut_value += (*current_node)->getEdgeWeight();


        if (Mpeak > (*current_node)->Mpeak + available_memory - (*current_node)->Mavail) {
            Mpeak = (*current_node)->Mpeak + available_memory - (*current_node)->Mavail;
        }
    }
    // if (cut_value > 0)
    //     cout << "Task: " << node->getId() << " CV: " << cut_value << endl;
    delete candidates;
    return;
}

void MinMem(Tree *tree, double MaxOutDeg, double &Required_memory, schedule_traversal &Schedule, int quiet) {
    double M = numeric_limits<double>::infinity();
    double Mpeak = MaxOutDeg;

    struct rlimit lim;
    lim.rlim_cur = RLIM_INFINITY;
    lim.rlim_max = RLIM_INFINITY;

    setrlimit(RLIMIT_STACK, &lim);

    Task *root = tree->getRoot();

    if (!quiet) {
        cerr << "Max out deg = " << Mpeak << endl;
    }

    list<Task *> L;
    Schedule.clear();
    while (Mpeak < numeric_limits<double>::infinity()) {
        Required_memory = Mpeak;
        explore(root, Required_memory, &L, &Schedule, M, L, Schedule, Mpeak, quiet, 0);


    }

}

void GreedyMinMem(Tree *tree, double &Required_memory) {
    //cout << "start computing greedy MM for root " << tree->getRoot()->getOtherSideId() << endl;
    vector<Task *> tasksWithoutMinMem;
    for (const auto &item: *tree->getTasks()) {
        if (item->getMinMemUnderlying() == 0) {
            tasksWithoutMinMem.push_back(item);
        }
    }

    for (auto it = tasksWithoutMinMem.rbegin(); it != tasksWithoutMinMem.rend(); it++) {
        vector<Task *> *children = (*it)->getChildren();
        double maxMemReqFromChildren = 0;
        for (const auto &child: *children) {
            double minMemUnderlyingFromChild = child->getMinMemUnderlying();
            if (minMemUnderlyingFromChild == 0) {
                //cout << "child has no MM" << endl;
                Tree *subtree = BuildSubtree(tree, child);
                GreedyMinMem(subtree, minMemUnderlyingFromChild);
                delete subtree;
            }
            double memReqFromChild = child->getEdgeWeight() + minMemUnderlyingFromChild;
            if (memReqFromChild > maxMemReqFromChildren) {
                maxMemReqFromChildren = memReqFromChild;
            }
        }
        Required_memory = maxMemReqFromChildren
                          + (*it)->getNodeWeight();
        // + (*it)->getEdgeWeight();
        // cout<<"result "<<Required_memory<<endl;
        // (*it)->setMinMemUnderlying(memReqForTask);
    }
}


void MinMem3Level(Tree *tree, double &Required_memory) {

    Task *root = tree->getRoot();
    if (root->getChildren()->size() == 0) {
        root->setMinMemUnderlying(root->getCost());
    }
    int originalTreeSize = tree->getSize();
    Tree *treeToBeChanged = BuildSubtree(tree, tree->getRoot());
    assert(treeToBeChanged->getSize() == tree->getSize());
    assert(treeToBeChanged->getTasks()->at(0)->getMakespanWeight() == tree->getTasks()->at(0)->getMakespanWeight());
    treeToBeChanged->levelsToTasks();

    vector<Task *> lowestLevel;
    std::set<Task *> firstLevel, secondLevel, thirdLevel;

    fillLowestAndThreeUpperLevels(treeToBeChanged, lowestLevel, firstLevel, secondLevel, thirdLevel);

    if (secondLevel.empty()) {
        double maxCost = root->getCost();
        for (const auto &item: *root->getChildren()) {
            if (item->getCost() > maxCost) maxCost = item->getCost();
        }
        root->setMinMemUnderlying(maxCost);
        treeToBeChanged->getRoot()->setMinMemUnderlying(maxCost);
    }
    while (!thirdLevel.empty()) {
        for (auto &item: firstLevel) {
            treeToBeChanged->mergeTaskToAllChildren(item);
        }
        treeToBeChanged->deepestLevel--;
        //compute MinMems for third level nodes
        for (const auto &item: thirdLevel) {
            item->setCostComputed(false);
            double maxout, requiredMemorySize;
            schedule_traversal *schedule_f = new schedule_traversal();
            Tree *subtree = BuildSubtree(treeToBeChanged, item);

            maxout = MaxOutDegree(subtree, true);
            MinMem(subtree, maxout, requiredMemorySize, *schedule_f, true);
            delete subtree;
            delete schedule_f;

            item->setMinMemUnderlying(requiredMemorySize);
            item->setCost(requiredMemorySize);
        }
        fillLowestAndThreeUpperLevels(treeToBeChanged, lowestLevel, firstLevel, secondLevel, thirdLevel);

    }

    if (!secondLevel.empty()) {
        for (auto &item: secondLevel) {
            item->setCostComputed(false);
            double maxout, requiredMemorySize;
            schedule_traversal *schedule_f = new schedule_traversal();
            Tree *subtree = BuildSubtree(treeToBeChanged, item);
            maxout = MaxOutDegree(subtree, true);
            MinMem(subtree, maxout, requiredMemorySize, *schedule_f, true);
            delete subtree;
            delete schedule_f;
            item->setMinMemUnderlying(requiredMemorySize);
        }
    }
    assert(tree->getSize() == originalTreeSize);
    assert(treeToBeChanged->getRoot()->getMinMemUnderlying() != 0);
    Required_memory = treeToBeChanged->getRoot()->getMinMemUnderlying();
}

void fillLowestAndThreeUpperLevels(Tree *treeToBeChanged, vector<Task *> &lowestLevel,
                                   set<Task *> &firstLevel, set<Task *> &secondLevel, set<Task *> &thirdLevel) {
    int deepesLevelInTree = treeToBeChanged->deepestLevel;
    lowestLevel.resize(0);
    firstLevel.clear();
    secondLevel.clear();
    thirdLevel.clear();

    //fill lowest level with leaves that are on the lowest level
    copy_if(treeToBeChanged->getTasks()->begin(), treeToBeChanged->getTasks()->end(),
            back_inserter(lowestLevel), [deepesLevelInTree](Task *a) {
                return a->getLevel() == deepesLevelInTree && a->getChildren()->size() == 0;
            });


    // first level are parents of leaves (to be merged)
    for (const auto &item: lowestLevel) {
        if (item->getParent() != NULL) {
            firstLevel.insert(item->getParent());
        }

    }

    //second, intermediate level - parents of first
    for (const auto &item: firstLevel) {
        if (item->getParent() != NULL) {
            secondLevel.insert(item->getParent());
        }
    }

    //third leve, MinMem is to be computed
    for (const auto &item: secondLevel) {
        if (item->getParent() != NULL) {
            thirdLevel.insert(item->getParent());
        }
    }
}
