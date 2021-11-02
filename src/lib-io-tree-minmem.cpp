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
#include <sys/time.h>
#include <sys/resource.h>


#include <iostream>
#include <fstream>
#include <limits>


#include "../include/lib-io-tree.h"
#include "../include/lib-io-tree-minmem.h"

void explore(Task *node, double available_memory, list<Task *> *L_init, schedule_t *S_init, double &cut_value,
             list<Task *> &min_sub_cut, schedule_t &sub_schedule, double &Mpeak, int quiet, int depth,
             uint64_t &count) {
    count++;

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
            schedule_t Sj;

            explore(*current_node, (*current_node)->Mavail, NULL, NULL, m_j, Lj, Sj, (*current_node)->Mpeak, quiet,
                    depth + 1, count);

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

    delete candidates;
    return;
}

void MinMem(Tree *tree, double MaxOutDeg, double &Required_memory, schedule_t &Schedule, int quiet, uint64_t &count) {
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
    count = 0;
    list<Task *> L;
    Schedule.clear();
    while (Mpeak < numeric_limits<double>::infinity()) {
        Required_memory = Mpeak;
        explore(root, Required_memory, &L, &Schedule, M, L, Schedule, Mpeak, quiet, 0, count);


    }
}

double MinMemRecurAlgorithm(int N, int *prnts, double *nwghts, double *ewghts, double *mswghts, int *schedule) {
    Tree *tree = new Tree(N, prnts, nwghts, ewghts, mswghts);

    double Mr;
    double Mp = MaxOutDegree(tree, true);
    uint64_t count = 0;
    schedule_t *sub_sched = new schedule_t();
    MinMem(tree, Mp, Mr, *sub_sched, true, count);

    unsigned int i = 0;
    for (schedule_t::reverse_iterator last = sub_sched->rbegin(); last != sub_sched->rend(); ++last) {
        schedule[i++] = (int) (*last);
    }

    delete sub_sched;
    delete tree;

    return Mr;
}

void exploreArray(int N, double *nwghts, double *ewghts, int *chstart, int *children, int nd, double available_memory,
                  list<node_peak_t> *L_init, int *S_init, double &cut_value, list<node_peak_t> &min_sub_cut,
                  int *sub_schedule, int &sched_head, double &Mpeak, int quiet, int depth, uint64_t &count) {
    count++;

    char *spacing;
    if (!quiet) {
        spacing = (char *) malloc((depth * 2 + 1) * sizeof(char));
        for (int i = 0; i < depth * 2; ++i) { spacing[i] = '\t'; }
        spacing[depth * 2] = '\0';
    }

    double cost = ewghts[nd] + nwghts[nd];
    for (int j = chstart[nd]; j < chstart[nd + 1]; j++) { cost += ewghts[children[j]]; }


    /* if node is unreachable, return +infty */
    if (cost > available_memory + MIN_ACCURACY) {
        if (!quiet) { cerr << spacing << "not enough mem for " << nd << endl; }
        Mpeak = cost;
        cut_value = numeric_limits<double>::infinity();
        return;
    }

    /* if this is a leaf, return 0 */
    if (chstart[nd] == chstart[nd + 1]) {
        if (!quiet) { fprintf(stderr, "%s%d is leaf\n", spacing, nd); }
        sub_schedule[sched_head--] = nd;
        cut_value = 0;
        Mpeak = numeric_limits<double>::infinity();
        return;
    }

    cut_value = 0;

    list<node_peak_t>::iterator begin;

    int candidates_count = 0;


    if (L_init != NULL) {
        if (L_init->size() > 0) {
            min_sub_cut.assign(L_init->begin(), L_init->end());
            candidates_count = min_sub_cut.size();
            begin = min_sub_cut.begin();
            for (list<node_peak_t>::iterator node = min_sub_cut.begin();
                 node != min_sub_cut.end(); ++node) { cut_value += ewghts[node->index]; }
            memcpy(&sub_schedule[sched_head], &S_init[sched_head], (N + 1 - sched_head) * sizeof(int));
        } else {
            sub_schedule[sched_head--] = nd;
            /* place every child in the candidate nodes*/
            candidates_count = chstart[nd + 1] - chstart[nd];
            for (int j = chstart[nd]; j < chstart[nd + 1]; j++) {
                min_sub_cut.push_back(node_peak_t(children[j]));
                cut_value += ewghts[children[j]];
            }
            begin = min_sub_cut.begin();
        }
    } else {
        sub_schedule[sched_head--] = nd;
        /* place every child in the candidate nodes*/
        candidates_count = chstart[nd + 1] - chstart[nd];
        begin = min_sub_cut.end();
        begin--;
        for (int j = chstart[nd]; j < chstart[nd + 1]; j++) {
            min_sub_cut.push_back(node_peak_t(children[j]));
            cut_value += ewghts[children[j]];
        }
        begin++;
    }

    int sub_cut_size = distance(begin, min_sub_cut.end());

    while (candidates_count > 0) {
        Mpeak = numeric_limits<double>::infinity();

        list<node_peak_t>::iterator current_node = begin;
        int local_count = 0;
        while (local_count++ < sub_cut_size) {
            double m_j = numeric_limits<double>::infinity();
            int tmp_sched_head = sched_head;
            double mp = current_node->mpeak;
            if (available_memory - cut_value + ewghts[current_node->index] >= current_node->mpeak) {

                double l_mavail = available_memory - cut_value + ewghts[current_node->index];

                if (!quiet) {
                    fprintf(stderr, "%ssubroot %d subtree avail mem %lf\n", spacing, current_node->index, l_mavail);
                }

                list<node_peak_t>::iterator current_subcut_end = min_sub_cut.end();
                current_subcut_end--;
                exploreArray(N, nwghts, ewghts, chstart, children, current_node->index, l_mavail, NULL, NULL, m_j,
                             min_sub_cut, sub_schedule, tmp_sched_head, current_node->mpeak, quiet, depth + 1, count);
                current_subcut_end++;
                mp = current_node->mpeak;

                if (!quiet) {
                    fprintf(stderr, "%sm_j %lf f_j %lf \n", spacing, m_j, ewghts[current_node->index]);
                    cerr << spacing << "Sj : {";
                    for (int j = tmp_sched_head + 1; j <= sched_head; j++) { cerr << sub_schedule[j] << " "; }
                    cerr << "}" << endl;
                }


                if (!quiet) { cerr << spacing << "m_j = " << m_j << " ew " << ewghts[current_node->index] << endl; }


                if (m_j > ewghts[current_node->index]) {
                    //pop the end of the list
                    min_sub_cut.erase(current_subcut_end, min_sub_cut.end());
                }
            }

            if (m_j <= ewghts[current_node->index]) {
                list<node_peak_t>::iterator to_erase = current_node;
                current_node--;
                local_count--;

                begin--;
                min_sub_cut.erase(to_erase);
                begin++;


                sub_cut_size = distance(begin, min_sub_cut.end());

                sched_head = tmp_sched_head;

            }

            current_node++;
        }


        cut_value = 0;

        candidates_count = 0;
        //		if(min_sub_cut.size()>0){
        list<node_peak_t>::iterator node = begin;
        for (local_count = 0; local_count < sub_cut_size; ++local_count) {
            cut_value += ewghts[node->index];
            ++node;
        }

        node = begin;
        for (local_count = 0; local_count < sub_cut_size; ++local_count) {
            if (available_memory - cut_value + ewghts[node->index] >= node->mpeak) {
                candidates_count++;
            }
            Mpeak = min(Mpeak, node->mpeak + cut_value - ewghts[node->index]);
            ++node;
        }
        if (!quiet) { cerr << spacing << candidates_count << " candidates left" << endl; }
    }

    return;
}

void exploreArray2(int N, double *nwghts, double *ewghts, int *chstart, int *children, iter_node_t *subroot,
                   double available_memory, bool Linit, iter_node_t *NodeArray, double &cut_value, MinMemDLL *L,
                   int *sub_schedule, int &sched_head, int depth, uint64_t &count) {
    count++;
    /* if node is unreachable, return +infty */
    if (subroot->cost - ewghts[subroot->index] > available_memory) {

        subroot->mpeak = subroot->cost;
        cut_value = numeric_limits<double>::infinity();
        return;
    }
    /* if this is a leaf, return 0 */
    if (chstart[subroot->index] == chstart[subroot->index + 1]) {
        sub_schedule[sched_head--] = subroot->index;
        cut_value = 0;
        subroot->mpeak = numeric_limits<double>::infinity();
        return;
    }

    cut_value = 0;
    int candidates_count = 0;
    if (Linit && L->size() > 0) {
        candidates_count = L->size();
        iter_node_t *cur_node = L->begin();
        while (cur_node != L->end()) {
            cut_value += ewghts[cur_node->index];
            cur_node = cur_node->pNext;
        }
        subroot->pSCend = L->last();
        subroot->pSCbegin = L->begin();
    } else {
        sub_schedule[sched_head--] = subroot->index;
        /* place every child in the candidate nodes*/
        candidates_count = chstart[subroot->index + 1] - chstart[subroot->index];
        iter_node_t *cur_node;

        iter_node_t *prev_node = subroot;
        for (int j = chstart[subroot->index]; j < chstart[subroot->index + 1]; j++) {
            cur_node = &NodeArray[children[j]];
            L->insert_after(prev_node, cur_node);
            cut_value += ewghts[cur_node->index];
            prev_node = cur_node;
        }
        subroot->pSCbegin = &NodeArray[children[chstart[subroot->index]]];
        subroot->pSCend = &NodeArray[children[chstart[subroot->index + 1] - 1]];
    }
    int cut_size = candidates_count;
    bool first_loop = true;
    while (candidates_count > 0) {

        iter_node_t *current_node = subroot->pSCbegin;

        while (cut_size > 0 && current_node != subroot->pSCend->pNext) {
            double m_j = numeric_limits<double>::infinity();
            int tmp_sched_head = sched_head;
            uint64_t prevCutSize = L->size();
            uint64_t subcut_size = 0;
            double available_memory_after_subroot = available_memory - cut_value + ewghts[subroot->index];

            if (available_memory_after_subroot + ewghts[current_node->index] >= current_node->mpeak || first_loop) {

                exploreArray2(N, nwghts, ewghts, chstart, children, current_node, available_memory_after_subroot, false,
                              NodeArray, m_j, L, sub_schedule, tmp_sched_head, depth + 1, count);
                subcut_size = L->size() - prevCutSize;

                if (m_j > ewghts[current_node->index]) {
                    //pop the end of the list

                    L->splice(current_node, current_node->pSCend, subcut_size);
                    if (subcut_size > 0) {
                        assert(current_node->pNext == current_node->pSCend->pNext);
                    }

                } else if (m_j <= ewghts[current_node->index]) {
                    //erase node
                    if (subcut_size > 0) {
                        // if there are nodes in the subcut
                        if (current_node == subroot->pSCbegin) {
                            subroot->pSCbegin = current_node->pSCbegin;
                            //              cerr<<"[ node "<<subroot->index<<" ] subcut new begin "<<subroot->pSCbegin->index<<endl;
                        }
                        if (current_node == subroot->pSCend) {
                            subroot->pSCend = current_node->pSCend;
                            //              cerr<<"[ node "<<subroot->index<<" ] subcut new end "<<subroot->pSCend->index<<endl;
                        }
                    } else {
                        // if the subcut of the children is empty but subroot's cut has other nodes than the one we are going to remove
                        if (cut_size > 1) {
                            if (current_node == subroot->pSCbegin && !(current_node == subroot->pSCend)) {
                                subroot->pSCbegin = current_node->pNext;
                                //              cerr<<"[ node "<<subroot->index<<" ] subcut new begin "<<subroot->pSCbegin->index<<endl;
                            } else if (current_node == subroot->pSCend && !(current_node == subroot->pSCbegin)) {
                                subroot->pSCend = current_node->pPrev;
                                //              cerr<<"[ node "<<subroot->index<<" ] subcut new end "<<subroot->pSCend->index<<endl;
                            }
                        }
                    }

                    L->erase(current_node);

                    cut_size += subcut_size - 1;
                    sched_head = tmp_sched_head;



                    //advance the pointer to skip the new nodes in the cut
                    if (subcut_size > 0) {
                        current_node = current_node->pSCend;
                    }
                }

            }
            current_node = current_node->pNext;
        }

        first_loop = false;
        cut_value = 0;
        current_node = subroot->pSCbegin;
        while (cut_size > 0 && subroot->pSCend->pNext != current_node) {
            cut_value += ewghts[current_node->index];
            current_node = current_node->pNext;
        }

        double available_memory_after_subroot = available_memory - cut_value + ewghts[subroot->index];
        candidates_count = 0;
        current_node = subroot->pSCbegin;
        while (cut_size > 0 && subroot->pSCend->pNext != current_node) {
            if (available_memory_after_subroot + ewghts[current_node->index] >= current_node->mpeak) {
                candidates_count++;
            }
            current_node = current_node->pNext;
        }
    }


    cut_value = 0;
    iter_node_t *current_node = subroot->pSCbegin;
    while (cut_size > 0 && subroot->pSCend->pNext != current_node) {
        cut_value += ewghts[current_node->index];
        current_node = current_node->pNext;
    }


    current_node = subroot->pSCbegin;
    subroot->mpeak = numeric_limits<double>::infinity();
    while (cut_size > 0 && subroot->pSCend->pNext != current_node) {

        subroot->mpeak = min(subroot->mpeak, current_node->mpeak + cut_value - ewghts[current_node->index]);
        current_node = current_node->pNext;
    }
    return;
}


void MinMemArray(int N, int root, double *nwghts, double *ewghts, int *chstart, int *children, double MaxOutDeg,
                 double &Required_memory, int *Schedule, int quiet, uint64_t &count, int *prnts) {
    double Cut_value = 0;

    struct rlimit lim;
    lim.rlim_cur = RLIM_INFINITY;
    lim.rlim_max = RLIM_INFINITY;

    setrlimit(RLIMIT_STACK, &lim);

    iter_node_t *NodeArray = new iter_node_t[N + 1];
    NodeArray[root].mpeak = MaxOutDeg;
    for (int i = 1; i < N + 1; i++) {
        NodeArray[i].index = i;
        NodeArray[i].cost = ewghts[i] + nwghts[i];
        for (int j = chstart[i]; j < chstart[i + 1]; j++) { NodeArray[i].cost += ewghts[children[j]]; }
    }

    MinMemDLL *MinCut = new MinMemDLL(&NodeArray[root]);

    int sched_head = N - 1;
    bool initialize_cut = false;
    count = 0;
    while (NodeArray[root].mpeak < numeric_limits<double>::infinity()) {
        Required_memory = NodeArray[root].mpeak;
        exploreArray2(N, nwghts, ewghts, chstart, children, &NodeArray[root], Required_memory, initialize_cut,
                      NodeArray, Cut_value, MinCut, Schedule, sched_head, 0, count);
        initialize_cut = true;
    }
    delete MinCut;
    delete[] NodeArray;
}