//
// Created by kulagins on 29.08.22.
//


#include "../include/lib-io-tree-minmem.h"
#include "../include/lib-io-tree-free-methods.h"
#include "../include/heuristics.h"
#include <tuple>
#include <algorithm>


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
                       back_inserter(others),
                       [](auto &a, auto &b) { return a->getId() < b->getId(); });

        taskHeap->clear();
        taskHeap->insert(taskHeap->end(), withTMaxInf.begin(), withTMaxInf.end());
        taskHeap->insert(taskHeap->end(), others.begin(), others.end());
    }

    assert(taskHeap->size() == init_size);
}