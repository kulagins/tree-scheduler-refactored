/*
 *  lib-io-tree-utils.cpp
 *  lib-io-tree
 *
 *  Created by defbond on 8/10/10.
 *  Copyright 2010 LIP/ENS-Lyon. All rights reserved.
 *
 */
#include <stdio.h>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <cmath>
#include <sstream>
#include <vector>

#include "../include/tree.h"
#include "../include/lib-io-tree.h"
#include "../include/lib-io-tree-minmem.h"
#include "../include/lib-io-tree-free-methods.h"
#include <sys/time.h>
#include <algorithm>

#ifndef CLUSTER_H
#define CLUSTER_H

using namespace std;

Tree *Tree::originalTree = NULL;

bool sort_sche(node_sche a, node_sche b) {
    return (a.second > b.second);
}

bool sort_ew(node_ew a, node_ew b) {
    return (a.second > b.second);
}

double u_wseconds(void) {
    struct timeval tp;

    gettimeofday(&tp, NULL);

    return (double) tp.tv_sec + (double) tp.tv_usec / 1000000.0;
};


// BUilds quotient tree for the whole original tree
Tree *Tree::BuildQtree() { //Qtree is for makespan side, so do not use it for space side
    root->breakEdge();
    root->getMakespanCost(true, true); //update

    Task *copy;
    Task *parent;
    Task *rootCopy;
    auto *tasksInQtree = new vector<Task *>();
    rootCopy = new Task(*root, 1, nullptr);
    rootCopy->setNodeWeight(root->getSequentialPart());
    rootCopy->setMakespanWeight(root->getSequentialPart());
    //rootCopy->toggleRootStatus(true);
    tasksInQtree->push_back(rootCopy);
    root->setOtherSideId(1);
    rootCopy->setOtherSideId(root->getId());

    Task *currentNode;
    int nodeIdCounter = 2;
    for (unsigned int i = 2; i <= this->getSize(); ++i) {
        currentNode = this->getTask(i);
        if (currentNode->isBroken()) {
            //   cout <<"cpy "<<currentNode->getId()<<" osideid "<<nodeIdCounter<<endl;
            copy = new Task(*currentNode, nodeIdCounter, nullptr);
            copy->setNodeWeight(currentNode->getSequentialPart());
            copy->setMakespanWeight(currentNode->getSequentialPart());
            currentNode->setOtherSideId(nodeIdCounter);
            copy->setOtherSideId(currentNode->getId());
            tasksInQtree->push_back(copy);
            nodeIdCounter++;
        }
    }
    nodeIdCounter = 2;

    vector<Task *> brokenTasksWithoutRoot = getBrokenTasks();
    brokenTasksWithoutRoot.erase(brokenTasksWithoutRoot.begin());

    for (Task *brokenTask: brokenTasksWithoutRoot) {
        //cout<<"broken task id " <<brokenTask->getId()<<endl;
        parent = brokenTask->getParent();
        while (parent != nullptr && !parent->isBroken()) {
            parent = parent->getParent();
        }
        if (parent != nullptr) {
            // cout<<"prnt "<<parent->getOtherSideId()<<endl;
            for (Task *childTask: *tasksInQtree) {
                if (childTask->getId() == nodeIdCounter) {
                    //   cout<<"set it on task "<<childTask->getId()<<endl;
                    childTask->setParentId(parent->getOtherSideId());
                }
            }
        }

        nodeIdCounter++;
    }

    Tree *Qtreeobj = new Tree(tasksInQtree, rootCopy,
                              this->getOriginalTree()); //Qtree only reprents makespan, not memory consumption

    for (unsigned int i = 1; i <= HowmanySubtrees(true); i++) {
        Qtreeobj->getTask(i)->breakEdge();
    }

    unsigned long treeSize = Qtreeobj->getTasks()->size();
    for (unsigned int i = 0; i < treeSize; i++) {
        Task *task = Qtreeobj->getTaskByPos(i);
        if (!task->isRoot()) {
            parent = Qtreeobj->getTask(task->getParentId());
            task->setParent(parent);
            parent->addChild(task);
        }
    }

    return Qtreeobj;
}


///Qtree corresponds to a whole original tree
Tree *Tree::BuildQtreeOld() { //Qtree is for makespan side, so do not use it for space side
    Task *root = this->getRoot();
    root->breakEdge();
    this->getRoot()->getMakespanCost(true, true); //update
    size_t tree_size = this->getTasks()->size();
    unsigned long num_subtrees = this->HowmanySubtrees(true);

    int *prnts = new int[num_subtrees + 1];
    double *ewghts = new double[num_subtrees + 1];
    double *timewghts = new double[num_subtrees + 1];
    int *brokenEdges = new int[num_subtrees + 1];

    //creat Quotient tree
    brokenEdges[1] = 1; //root node
    prnts[1] = 0;
    ewghts[1] = 0;
    timewghts[1] = root->getSequentialPart();
    unsigned int j = 2;
    root->setOtherSideId(1);

    Task *currentNode;
    for (unsigned int i = 2; i <= tree_size; ++i) {
        currentNode = this->getTask(i);
        if (currentNode->isBroken()) {
            //  cout << "broken node " << currentNode->getId() << " " << currentNode->getSequentialPart() << " i " << i
            //        << " j " << j << endl;
            currentNode->setOtherSideId(j); //corresponding node's ID on Qtree
            brokenEdges[j] = i;
            timewghts[j] = currentNode->getSequentialPart();
            ewghts[j] = currentNode->getEdgeWeight();
            ++j;
        }
    }

    for (unsigned int i = 2; i <= num_subtrees; ++i) {
        currentNode = this->getTask(brokenEdges[i])->getParent();
        while (!currentNode->isBroken()) {
            currentNode = currentNode->getParent();
        }
        prnts[i] = currentNode->getOtherSideId();
    }
    //   for (int i = 0; i < num_subtrees; i++) {
    //    cout << i << " " << prnts[i] << " " << timewghts[i] << endl;
    //   }

    Tree *Qtreeobj = new Tree(num_subtrees, prnts, timewghts, ewghts,
                              timewghts); //Qtree only reprents makespan, not memory consumption

    for (unsigned int i = 1; i <= num_subtrees; i++) {
        Qtreeobj->getTask(i)->breakEdge();                    //break edge
        Qtreeobj->getTask(i)->setOtherSideId(brokenEdges[i]); //corresponding node's ID on tree
    }

    delete[] prnts;
    delete[] ewghts;
    delete[] timewghts;
    delete[] brokenEdges;

    return Qtreeobj;
}


unsigned int Tree::HowmanySubtrees(bool quiet) {
    unsigned int number_subtrees = 0;
    this->getRoot()->breakEdge();
    const vector<Task *> *Nodes = this->getTasks();
    if (quiet == false) {
        cout << "Broken Edges { ";
    }
    for (auto it = Nodes->begin(); it != Nodes->end(); ++it) {
        if ((*it)->isBroken()) {
            number_subtrees++;
            if (quiet == false) {
                cout << (*it)->getId() << " ";
            }
        }
    }
    if (quiet == false) {
        cout << "}" << endl;
    }
    return number_subtrees;
}

bool
Tree::MemoryEnough(Task *Qrootone, Task *Qroottwo, bool leaf, double available_memory_size,
                   double &requiredMemorySize) {
    bool enough = false;

    Task *SubtreeRoot = this->getTask(Qrootone->getOtherSideId());

    vector<Task *> *childrenvector = Qrootone->getChildren();
    if ((leaf == true) & (childrenvector->size() == 2)) {
        this->getTask(childrenvector->front()->getOtherSideId())->restoreEdge();
        this->getTask(childrenvector->back()->getOtherSideId())->restoreEdge();
    } else {
        this->getTask(Qroottwo->getOtherSideId())->restoreEdge(); //restore edge temporarilly
    }


    Tree *subtree = BuildSubtree(this, SubtreeRoot);
    double maxout;
    schedule_traversal *schedule_f = new schedule_traversal();
    maxout = MaxOutDegree(subtree, true);
    MinMem(subtree, maxout, requiredMemorySize, *schedule_f, true);

    if (requiredMemorySize <= available_memory_size) {
        enough = true;
    }

    if ((leaf == true) & (childrenvector->size() == 2)) {
        this->getTask(childrenvector->front()->getOtherSideId())->breakEdge();
        this->getTask(childrenvector->back()->getOtherSideId())->breakEdge();
    } else {
        this->getTask(Qroottwo->getOtherSideId())->breakEdge();
    }

    delete subtree;
    delete schedule_f;
    return enough;
}


double Task::Sequence() {
    return this->getMakespanCost();
}

Tree *read_tree(const char *filename) {
    //  cout << filename << endl;
    ifstream OpenFile(filename);
    string line;
    stringstream line_stream;

    // count the number of tasks
    unsigned int num_tasks = 0;
    while (getline(OpenFile, line)) {
        if (!line.empty() && line[0] != '%') {
            num_tasks++;
        }
    }
    OpenFile.clear();
    OpenFile.seekg(0, ios::beg);
    Tree *tree = new Tree();

    while (getline(OpenFile, line)) {
        line_stream.clear();
        line_stream.str(line);

        if (!line_stream.str().empty()) {
            unsigned int id;
            unsigned int parent_id;
            double ew, nw, msw;
            Task *task;

            line_stream >> id >> parent_id >> nw >> msw >> ew;

            if (parent_id == 0) {
                task = new Task(nw, ew, msw, true);
                task->setId(1 + num_tasks - id);
                tree->addRoot(task);
            } else {
                task = new Task(1 + num_tasks - parent_id, nw, ew, msw);
                task->setId(1 + num_tasks - id);
                tree->addTask(task);
            }
        }

    }
    OpenFile.close();
    tree->reverseVector();
    Task *parent;
    unsigned long treeSize = tree->getTasks()->size();
    for (unsigned int i = 0; i < treeSize; i++) {
        Task *task = tree->getTaskByPos(i);
        if (!task->isRoot()) {
            parent = tree->getTask(task->getParentId());
            task->setParent(parent);
            parent->addChild(task);
        }
    }

    return tree;

}

double IOCounter(Tree &tree, schedule_traversal &sub_schedule, double available_memory, bool divisible, int quiet) {
    double memory_occupation = 0;
    double io_volume = 0;
    io_map unloaded_nodes;
    schedule_traversal loaded_nodes;

    /*iterates through the given permutation (schedule)*/
    for (schedule_traversal::iterator cur_task_id = sub_schedule.begin();
         cur_task_id != sub_schedule.end(); cur_task_id++) {
        Task *cur_node = tree.getTask(*cur_task_id);

        /*if the node was unloaded*/
        if (unloaded_nodes.find(*cur_task_id) != unloaded_nodes.end()) {
            if (!quiet) {
                cerr << "Loading " << unloaded_nodes[*cur_task_id] << "of " << *cur_task_id << " (IO)" << endl;
            }

            memory_occupation += unloaded_nodes[*cur_task_id];
            unloaded_nodes.erase(*cur_task_id);
        }

        double data_to_unload = memory_occupation + cur_node->getCost() - cur_node->getEdgeWeight() - available_memory;
        if (!quiet) {
            cerr << "min data to unload " << data_to_unload << endl;
        }
        if (data_to_unload > 0) {
            /*if we dont have enough room, unload files and update both io and occupation*/
            double unloaded_data = 0;
            /*unload furthest non unloaded node which is NOT in current_node children first*/
            schedule_traversal::reverse_iterator far_node_id = loaded_nodes.rbegin();
            while ((far_node_id != loaded_nodes.rend()) && (unloaded_data < data_to_unload)) {
                /*try to unload this node*/
                bool is_already_unloaded = false;
                double remaining_loaded_data = tree.getTask(*far_node_id)->getEdgeWeight();
                if (unloaded_nodes.find(*far_node_id) != unloaded_nodes.end()) {
                    is_already_unloaded = true;
                    remaining_loaded_data = max(remaining_loaded_data - unloaded_nodes[*far_node_id], 0.0);
                }

                double local_data_to_unload;
                if (divisible) {
                    local_data_to_unload = min(remaining_loaded_data, data_to_unload);
                } else {
                    local_data_to_unload = remaining_loaded_data;
                }

                if (!quiet) {
                    cerr << "unloading (IO) " << local_data_to_unload << " of " << *far_node_id << endl;
                }
                unloaded_data += local_data_to_unload;
                if (is_already_unloaded) {
                    unloaded_nodes[*far_node_id] += local_data_to_unload;
                } else {
                    unloaded_nodes[*far_node_id] = local_data_to_unload;
                }

                if (remaining_loaded_data == local_data_to_unload) {
                    loaded_nodes.remove(*far_node_id);
                }

                far_node_id++;
            }

            io_volume += unloaded_data;
            memory_occupation -= unloaded_data;
        }

        if (!quiet) {
            cerr << "occupation before processing " << memory_occupation << endl;
        }
        /*if we have enough memory to process the node, update occupation*/
        memory_occupation += cur_node->getCost() - 2 * cur_node->getEdgeWeight() - cur_node->getNodeWeight();
        if (!quiet) {
            cerr << "processing " << *cur_task_id << endl;
        }
        if (!quiet) {
            cerr << "unloading " << *cur_task_id << endl;
        }
        loaded_nodes.remove(*cur_task_id);

        if (!quiet) {
            cerr << "loading ";
        }
        for (vector<Task *>::iterator child = cur_node->getChildren()->begin();
             child != cur_node->getChildren()->end(); child++) {
            if (!quiet) {
                cerr << (*child)->getId() << " ";
            }
            loaded_nodes.push_back((*child)->getId());
        }
        if (!quiet) {
            cerr << endl;
        }
        if (!quiet) {
            cerr << "New occupation after processing " << memory_occupation << endl;
        }
    }

    return io_volume;
    //    cerr<<"IO Volume "<<io_volume<<endl;
}

double unload_largest_first_fit(Tree *tree, vector<unsigned int> &unloaded_nodes, list<node_ew> &loaded_nodes,
                                const double data_to_unload) {
    double unloaded_data = 0.0;

    /*loaded_nodes already sorted, unload largest nodes till there is enough space*/
    list<node_ew>::iterator largest_node = loaded_nodes.begin();
    while ((largest_node != loaded_nodes.end()) && (unloaded_data < data_to_unload)) {
        largest_node = loaded_nodes.begin();
        tree->getTask(largest_node->first)->breakEdge(); //break this edge;
        //cout<<"******Largest first: break edge "<<largest_node->first<<endl;
        double edge_weight = largest_node->second;
        unloaded_data += edge_weight;
        unloaded_nodes.push_back(largest_node->first);
        loaded_nodes.pop_front();
    }

    return unloaded_data;
}

double unload_furthest_nodes(Tree *tree, vector<unsigned int> &unloaded_nodes, list<node_sche> &loaded_nodes,
                             const double data_to_unload, bool divisible) {
    double unloaded_data = 0.0;
    // cout << "loaded nodes: size " << loaded_nodes.size();
    // list<node_sche>::iterator loaded_iterator = loaded_nodes.begin();
    // if (loaded_nodes.size() < 20)
    // {
    //     while (loaded_iterator != loaded_nodes.end())
    //     {
    //         cout << loaded_iterator->first << endl;
    //     }
    //     cout << endl;
    // }
    // else
    // {
    //     cout << "big loaded nodes";
    // }
    // cout << "unloaded nodes: ";
    // for (int i = 0; i < unloaded_nodes.size(); i++)
    // {
    //     cout << unloaded_nodes[i] << endl;
    // }
    // cout << endl;

    /*unload furthest non unloaded node which is NOT in current_node children first*/
    list<node_sche> old_loaded_nodes = loaded_nodes;
    list<node_sche>::iterator far_node = loaded_nodes.begin();
    while ((far_node != loaded_nodes.end()) && (unloaded_data < data_to_unload)) {
        // cout << far_node->first << " " << far_node->second << endl;
        /*try to unload this node*/
        far_node = loaded_nodes.begin();
        Task *taskFarNode = tree->getTask(far_node->first);
        double remaining_loaded_data = taskFarNode->getEdgeWeight();
        double local_data_to_unload;
        if (divisible) {
            local_data_to_unload = min(remaining_loaded_data, data_to_unload);
        } else {
            local_data_to_unload = remaining_loaded_data;
        }

        unloaded_data += local_data_to_unload;

        unloaded_nodes.push_back((*far_node).first);

        ////cout<<"-------LSNF remove "<<local_data_to_unload<<endl;
        ////cout<<"break edge "<<far_node->first<<endl;
        if (far_node->first == 0) {
            cout << "Problem! loaded nodes " << endl;
            list<node_sche>::iterator loaded_iterator = loaded_nodes.begin();
            if (loaded_nodes.size() < 20) {
                while (loaded_iterator != loaded_nodes.end()) {
                    cout << loaded_iterator->first << endl;
                }
            } else {
                cout << "big loaded nodes";
            }
            cout << "unloaded nodes: ";
            for (int i = 0; i < unloaded_nodes.size(); i++) {
                cout << unloaded_nodes[i] << endl;
            }

            cout << "old loaded nodes size " << old_loaded_nodes.size() << endl;
            cout << "old loaded nodes first " << old_loaded_nodes.begin()->first << endl;
            tree->getTask(far_node->first)->breakEdge(); //break this edge
            loaded_nodes.pop_front();
        }

        // cout << "after pop " << far_node->first << " " << far_node->second << endl;
    }

    return unloaded_data;
}


double unload_furthest_first_fit(Tree *tree, vector<unsigned int> &unloaded_nodes, list<node_sche> &loaded_nodes,
                                 const double data_to_unload, bool divisible) {
    double unloaded_data = 0.0;
    /*unload furthest non unloaded node which is NOT in current_node children first*/
    list<node_sche>::iterator far_node = loaded_nodes.begin();
    while ((far_node != loaded_nodes.end()) && (unloaded_data < data_to_unload)) {
        /*try to unload this node*/
        Task *farNodeTask = tree->getTask(far_node->first);

        double remaining_loaded_data = farNodeTask->getEdgeWeight();

        double local_data_to_unload;
        if (divisible) {
            local_data_to_unload = min(remaining_loaded_data, data_to_unload);
        } else {
            local_data_to_unload = remaining_loaded_data;
        }
        /*if it "fits", that is if the amount of data is lower than what we need to unload*/
        if (local_data_to_unload >= data_to_unload) {
            //cout<<"-------First fit success: ";
            //cout<<"break edge "<<far_node->first<<endl;
            unloaded_data += local_data_to_unload;
            unloaded_nodes.push_back(far_node->first);
            farNodeTask->breakEdge(); //break this edge
            loaded_nodes.erase(far_node);
            return unloaded_data;
        }
        far_node++;
    }

    /*if we did not unloaded enough data, call furthest node*/
    ////cout<<"-------First fit failed, go to LSNF: "<<endl;
    unloaded_data = unload_furthest_nodes(tree, unloaded_nodes, loaded_nodes, data_to_unload, divisible);

    return unloaded_data;
}

double
unload_furthest_best_fit(io_map &unloaded_nodes, schedule_traversal &loaded_nodes, const double data_to_unload,
                         double *ewghts,
                         bool divisible) {
    double unloaded_data = 0.0;
    /*unload furthest non unloaded node which is NOT in current_node children first*/
    while (unloaded_data < data_to_unload) {
        bool is_best_already_unloaded = false;
        double best_remaining_loaded_data = 0.0;
        unsigned int best_candidate = -1;
        double best_candi_score = -1.0;

        schedule_traversal::reverse_iterator far_node_id = loaded_nodes.rbegin();
        while ((far_node_id != loaded_nodes.rend())) {
            /*try to unload this node*/
            bool is_already_unloaded = false;
            double remaining_loaded_data = ewghts[*far_node_id];
            if (unloaded_nodes.find(*far_node_id) != unloaded_nodes.end()) {
                is_already_unloaded = true;
                remaining_loaded_data = max(remaining_loaded_data - unloaded_nodes[*far_node_id], 0.0);
            }
            double local_data_to_unload;
            if (divisible) {
                local_data_to_unload = min(remaining_loaded_data, data_to_unload);
            } else {
                local_data_to_unload = remaining_loaded_data;
            }
            /*if it "fits", that is if the amount of data is lower than what we need to unload*/
            if (local_data_to_unload <= data_to_unload - unloaded_data) {
                if (local_data_to_unload > best_candi_score) {
                    best_candi_score = local_data_to_unload;
                    best_candidate = *far_node_id;
                    is_best_already_unloaded = is_already_unloaded;
                    best_remaining_loaded_data = remaining_loaded_data;
                }
            }

            far_node_id++;
        }

        if (best_candi_score != -1) {
            //cerr<<"best found :"<<best_candidate<<endl;
            unloaded_data += best_candi_score;
            if (is_best_already_unloaded) {
                unloaded_nodes[best_candidate] += best_candi_score;
            } else {
                unloaded_nodes[best_candidate] = best_candi_score;
            }
            if (best_remaining_loaded_data == best_candi_score) {
                loaded_nodes.remove(best_candidate);
            }
        } else {
            break;
        }
    }
    return unloaded_data;
}

double
unload_furthest_first_fit_abs(io_map &unloaded_nodes, schedule_traversal &loaded_nodes, const double data_to_unload,
                              double *ewghts, bool divisible) {
    double unloaded_data = 0.0;
    /*unload furthest non unloaded node which is NOT in current_node children first*/
    schedule_traversal::reverse_iterator far_node_id = loaded_nodes.rbegin();
    while ((far_node_id != loaded_nodes.rend()) && (unloaded_data < data_to_unload)) {
        /*try to unload this node*/
        bool is_already_unloaded = false;
        double remaining_loaded_data = ewghts[*far_node_id];
        if (unloaded_nodes.find(*far_node_id) != unloaded_nodes.end()) {
            is_already_unloaded = true;
            remaining_loaded_data = max(remaining_loaded_data - unloaded_nodes[*far_node_id], 0.0);
        }
        double local_data_to_unload;
        if (divisible) {
            local_data_to_unload = min(remaining_loaded_data, data_to_unload);
        } else {
            local_data_to_unload = remaining_loaded_data;
        }
        /*if it "fits", that is if the amount of data is lower than what we need to unload*/
        if (local_data_to_unload >= data_to_unload - unloaded_data) {
            unloaded_data += local_data_to_unload;
            if (is_already_unloaded) {
                unloaded_nodes[*far_node_id] += local_data_to_unload;
            } else {
                unloaded_nodes[*far_node_id] = local_data_to_unload;
            }
            if (remaining_loaded_data == local_data_to_unload) {
                loaded_nodes.remove(*far_node_id);
            }
        }

        far_node_id++;
    }

    return unloaded_data;
}

double
unload_furthest_best_fit_abs(io_map &unloaded_nodes, schedule_traversal &loaded_nodes, const double data_to_unload,
                             double *ewghts, bool divisible) {
    double unloaded_data = 0.0;
    /*unload furthest non unloaded node which is NOT in current_node children first*/
    while (unloaded_data < data_to_unload) {
        bool is_best_already_unloaded = false;
        double best_remaining_loaded_data = 0.0;
        unsigned int best_candidate = -1;
        double best_candi_score = -1.0;

        schedule_traversal::reverse_iterator far_node_id = loaded_nodes.rbegin();
        while ((far_node_id != loaded_nodes.rend())) {
            /*try to unload this node*/
            bool is_already_unloaded = false;
            double remaining_loaded_data = ewghts[*far_node_id];
            if (unloaded_nodes.find(*far_node_id) != unloaded_nodes.end()) {
                is_already_unloaded = true;
                remaining_loaded_data = max(remaining_loaded_data - unloaded_nodes[*far_node_id], 0.0);
            }
            double local_data_to_unload;
            if (divisible) {
                local_data_to_unload = min(remaining_loaded_data, data_to_unload);
            } else {
                local_data_to_unload = remaining_loaded_data;
            }
            /*if it "fits", that is if the amount of data is lower than what we need to unload*/
            if (abs(data_to_unload - unloaded_data - local_data_to_unload) <
                abs(data_to_unload - unloaded_data - best_candi_score)) {
                best_candi_score = local_data_to_unload;
                best_candidate = *far_node_id;
                is_best_already_unloaded = is_already_unloaded;
                best_remaining_loaded_data = remaining_loaded_data;
            }

            far_node_id++;
        }

        if (best_candi_score != -1) {
            //cerr<<"best found :"<<best_candidate<<endl;
            unloaded_data += best_candi_score;
            if (is_best_already_unloaded) {
                unloaded_nodes[best_candidate] += best_candi_score;
            } else {
                unloaded_nodes[best_candidate] = best_candi_score;
            }
            if (best_remaining_loaded_data == best_candi_score) {
                loaded_nodes.remove(best_candidate);
            }
        } else {
            break;
        }
    }

    return unloaded_data;
}

int next_comb(int comb[], int k, int n) {
    int i = k - 1;
    ++comb[i];

    while ((i > 0) && (comb[i] >= n - k + 1 + i)) {
        --i;
        ++comb[i];
    }

    if (comb[0] > n - k) /* Combination (n-k, n-k+1, ..., n) reached */
        return 0;        /* No more combinations can be generated */

    /* comb now looks like (..., x, n, n, n, ..., n).
     *  Turn it into (..., x, x + 1, x + 2, ...) */
    for (i = i + 1; i < k; ++i)
        comb[i] = comb[i - 1] + 1;

    return 1;
}

double
unload_best_increasing_combi(io_map &unloaded_nodes, schedule_traversal &loaded_nodes, const double data_to_unload,
                             double *ewghts, bool divisible, unsigned int init_combi_size, bool quiet) {
    assert(!divisible);
    vector<unsigned int> candidates;
    double unloaded_data = 0.0;
    vector<unsigned int> best_combi;
    double best_combi_score = numeric_limits<double>::infinity();
    while (unloaded_data < data_to_unload) {
        /*unload furthest non unloaded node which is NOT in current_node children first*/
        candidates.clear();
        schedule_traversal::reverse_iterator far_node_id = loaded_nodes.rbegin();
        while ((far_node_id != loaded_nodes.rend()) && (candidates.size() <= init_combi_size)) {
            /*try to unload this node*/
            if (unloaded_nodes.find(*far_node_id) == unloaded_nodes.end()) {
                candidates.push_back(*far_node_id);
            }

            far_node_id++;
        }
        /*now we got our max_candidates candidates, compute every permutations*/

        vector<unsigned int> cur_combi;

        int n = candidates.size(); /* The size of the set; for {1, 2, 3, 4} it's 4 */
        for (int k = 1; k <= n; k++) {
            /* k is the size of the subsets; for {1, 2}, {1, 3}, ... it's 2 */
            int comb[k]; /* comb[i] is the index of the i-th element in the
                          combination */

            /* Setup comb for the initial combination */
            for (int i = 0; i < k; ++i)
                comb[i] = i;
            /*compute score*/
            cur_combi.clear();
            double cur_unloaded_data = 0.0;
            if (!quiet) {
                cerr << "cur_combi is [";
            }
            for (int i = 0; i < k; ++i) {
                unsigned int cur_node_id = candidates[comb[i]];
                cur_unloaded_data += ewghts[cur_node_id];
                cur_combi.push_back(cur_node_id);
                if (!quiet) {
                    cerr << " " << cur_node_id;
                }
            }
            if (!quiet) {
                cerr << "], its score is " << cur_unloaded_data << endl;
            }

            if (abs(data_to_unload - unloaded_data - cur_unloaded_data) <
                abs(data_to_unload - unloaded_data - best_combi_score)) {
                if (!quiet) {
                    cerr << "cur_combi is better than prev best " << cur_unloaded_data << " vs " << best_combi_score
                         << endl;
                }
                best_combi.assign(cur_combi.begin(), cur_combi.end());
                best_combi_score = cur_unloaded_data;
            }

            /* Generate and print all the other combinations */
            while (next_comb(comb, k, n)) {
                /*compute score*/
                cur_combi.clear();
                double cur_unloaded_data = 0.0;
                if (!quiet) {
                    cerr << "cur_combi is [";
                }
                for (int i = 0; i < k; ++i) {
                    unsigned int cur_node_id = candidates[comb[i]];
                    cur_unloaded_data += ewghts[cur_node_id];
                    cur_combi.push_back(cur_node_id);
                    if (!quiet) {
                        cerr << " " << cur_node_id;
                    }
                }
                if (!quiet) {
                    cerr << "], its score is " << cur_unloaded_data << endl;
                }

                if (abs(data_to_unload - unloaded_data - cur_unloaded_data) <
                    abs(data_to_unload - unloaded_data - best_combi_score)) {
                    if (!quiet) {
                        cerr << "cur_combi is better than prev best " << cur_unloaded_data << " vs " << best_combi_score
                             << endl;
                    }
                    best_combi.assign(cur_combi.begin(), cur_combi.end());
                    best_combi_score = cur_unloaded_data;
                }
            }
        }

        unloaded_data = best_combi_score;

        if (unloaded_data < data_to_unload) {
            init_combi_size *= 2;
        }
    }
    if (!quiet) {
        cerr << "best_combi is [";
    }
    for (unsigned int i = 0; i < best_combi.size(); i++) {
        unsigned int cur_id = best_combi[i];
        unloaded_nodes[cur_id] = ewghts[cur_id];
        loaded_nodes.remove(cur_id);
        if (!quiet) {
            cerr << " " << cur_id;
        }
    }
    if (!quiet) {
        cerr << "], its score is " << best_combi_score << endl;
    }

    if (!quiet) {
        cerr << "We've unloaded " << unloaded_data << " out of " << data_to_unload << " so far" << endl;
    }
    unloaded_data = best_combi_score;
    return unloaded_data;
}

double unload_best_furthest_nodes(io_map &unloaded_nodes, schedule_traversal &loaded_nodes, const double data_to_unload,
                                  double *ewghts, bool divisible, unsigned int max_candidates, bool quiet) {
    assert(!divisible);
    vector<unsigned int> candidates;
    double unloaded_data = 0.0;
    while (unloaded_data < data_to_unload) {
        /*unload furthest non unloaded node which is NOT in current_node children first*/
        candidates.clear();
        schedule_traversal::reverse_iterator far_node_id = loaded_nodes.rbegin();
        while ((far_node_id != loaded_nodes.rend()) && (candidates.size() <= max_candidates)) {
            /*try to unload this node*/
            if (unloaded_nodes.find(*far_node_id) == unloaded_nodes.end()) {
                candidates.push_back(*far_node_id);
            }

            far_node_id++;
        }
        /*now we got our max_candidates candidates, compute every permutations*/

        vector<unsigned int> best_combi;
        vector<unsigned int> cur_combi;
        double best_combi_score = numeric_limits<double>::infinity();

        int n = candidates.size(); /* The size of the set; for {1, 2, 3, 4} it's 4 */
        for (int k = 1; k <= n; k++) {
            /* k is the size of the subsets; for {1, 2}, {1, 3}, ... it's 2 */
            int comb[k]; /* comb[i] is the index of the i-th element in the
                          combination */

            /* Setup comb for the initial combination */
            for (int i = 0; i < k; ++i)
                comb[i] = i;
            /*compute score*/
            cur_combi.clear();
            double cur_unloaded_data = 0.0;
            if (!quiet) {
                cerr << "cur_combi is [";
            }
            for (int i = 0; i < k; ++i) {
                unsigned int cur_node_id = candidates[comb[i]];
                cur_unloaded_data += ewghts[cur_node_id];
                cur_combi.push_back(cur_node_id);
                if (!quiet) {
                    cerr << " " << cur_node_id;
                }
            }
            if (!quiet) {
                cerr << "], its score is " << cur_unloaded_data << endl;
            }

            if ((double) abs(data_to_unload - unloaded_data - cur_unloaded_data) <
                (double) abs(data_to_unload - unloaded_data - best_combi_score)) {
                if (!quiet) {
                    cerr << "cur_combi is better than prev best " << cur_unloaded_data << " vs " << best_combi_score
                         << endl;
                }
                best_combi.assign(cur_combi.begin(), cur_combi.end());
                best_combi_score = cur_unloaded_data;
            }

            /* Generate and print all the other combinations */
            while (next_comb(comb, k, n)) {
                /*compute score*/
                cur_combi.clear();
                double cur_unloaded_data = 0.0;
                if (!quiet) {
                    cerr << "cur_combi is [";
                }
                for (int i = 0; i < k; ++i) {
                    unsigned int cur_node_id = candidates[comb[i]];
                    cur_unloaded_data += ewghts[cur_node_id];
                    cur_combi.push_back(cur_node_id);
                    if (!quiet) {
                        cerr << " " << cur_node_id;
                    }
                }
                if (!quiet) {
                    cerr << "], its score is " << cur_unloaded_data << endl;
                }

                if (abs(data_to_unload - unloaded_data - cur_unloaded_data) <
                    abs(data_to_unload - unloaded_data - best_combi_score)) {
                    if (!quiet) {
                        cerr << "cur_combi is better than prev best " << cur_unloaded_data << " vs " << best_combi_score
                             << endl;
                    }
                    best_combi.assign(cur_combi.begin(), cur_combi.end());
                    best_combi_score = cur_unloaded_data;
                }
            }
        }

        if (!quiet) {
            cerr << "best_combi is [";
        }
        for (unsigned int i = 0; i < best_combi.size(); i++) {
            unsigned int cur_id = best_combi[i];
            unloaded_nodes[cur_id] = ewghts[cur_id];
            loaded_nodes.remove(cur_id);
            if (!quiet) {
                cerr << " " << cur_id;
            }
        }
        if (!quiet) {
            cerr << "], its score is " << best_combi_score << endl;
        }

        unloaded_data += best_combi_score;
        if (!quiet) {
            cerr << "We've unloaded " << unloaded_data << " out of " << data_to_unload << " so far" << endl;
        }
    }

    return unloaded_data;
}

list<unsigned int> getAllNodesUnderCurrent(const Tree *tree, int cur_task_id) {
    list<unsigned int> temp;
    list<unsigned int> queue;
    queue.push_back(cur_task_id);
    //cout<<"children "<<endl;
    do {
        Task *firstInQueue = tree->getTask(queue.front());
        temp.push_back(queue.front());
        queue.pop_front();

        for (Task *child: *firstInQueue->getChildren()) {
            //cout<<*(children+i)<<" ";
            queue.push_back(child->getId());
        }
    } while (!queue.empty());
    return temp;
}

//TODO check for N
double IOCounter(Tree *tree, int *schedule, bool divisible, int quiet, unsigned int &com_freq,
                 vector<unsigned int> *brokenEdges, io_method_t method) {
    double memory_occupation = tree->getTasks()->at(schedule[tree->getSize()] - 1)->getEdgeWeight();
    double io_volume = 0;
    vector<unsigned int> unloaded_nodes;
    list<node_sche> loaded_nodes;
    list<node_ew> loaded_nodes_ew;
    vector<int> schedule_vec(schedule, schedule + tree->getSize() + 1);

    int cur_task_id;
    vector<unsigned int>::iterator unloaded;
    list<unsigned int> temp;

    unsigned int subtree_size;
    list<unsigned int>::iterator iter;
    double maxoutD, memory_required, IO_sub = 0;
    schedule_traversal *schedule_f = new schedule_traversal();
    list<int>::iterator ite_sche;
    vector<unsigned int> subtreeBrokenEdges;

    /*iterates through the given permutation (schedule)*/
    for (int rank = tree->getSize(); rank >= 1; rank--) {
        cur_task_id = schedule[rank];

        //cout<<"current task id: "<<cur_task_id<<endl;
        if (cur_task_id != 0) { //0 means this node is on other subtrees
            /*if the node was unloaded*/
            Task *currentTask = tree->getTask(cur_task_id);
            unloaded = find(unloaded_nodes.begin(), unloaded_nodes.end(), cur_task_id);
            if (unloaded != unloaded_nodes.end()) { //find node cur_task_id unloaded
                //cout<<", (break) "<<endl;
                if (Cluster::getFixedCluster()->hasFreeProcessor()) {
                    Cluster::getFixedCluster()->getFirstFreeProcessor()->assignTaskId(
                            tree->getTask(cur_task_id)->getOtherSideId());
                }

                brokenEdges->push_back(tree->getTask(cur_task_id)->getOtherSideId());
//                cout << "iocounter get free proc "
                //                   << Cluster::getFixedCluster()->getFirstFreeProcessor()->getMemorySize() << " for task "
                //                  << subtree->getTask(cur_task_id)->getOtherSideId() << endl;
                ++com_freq;
                unloaded_nodes.erase(unloaded);
                temp.clear();
                temp = getAllNodesUnderCurrent(tree, cur_task_id);

                subtree_size = temp.size(); //just an approximation
                for (long i = rank - 1; i >= 0; i--) {
                    iter = find(temp.begin(), temp.end(), schedule[i]);
                    if (iter != temp.end()) {
                        schedule[i] = 0; //IO counter will pass 0;
                        temp.erase(iter);
                    }
                    if (temp.size() == 1) {
                        break;
                    }
                }

                Tree *subtree = BuildSubtree(tree, tree->getTask(cur_task_id));
                subtree_size = subtree->getSize();
                //  cout << "subtree size " << subtree_size << endl;

                int *schedule_copy = new int[subtree_size + 1];
                maxoutD = MaxOutDegree(subtree, true);
                schedule_f->clear();

                MinMem(subtree, maxoutD, memory_required, *schedule_f, true);
                schedule_copy = copyScheduleBackwards(schedule_f);
                if (Cluster::getFixedCluster()->hasFreeProcessor()) {
                    Cluster::getFixedCluster()->getFirstFreeProcessor()->setOccupiedMemorySize(
                            memory_required);
                }

                if (memory_required > Cluster::getFixedCluster()->getFirstFreeProcessorOrSmallest()->getMemorySize()) {
                    //   cout << "memory required " << memory_required << ", is larger than what is available " << available_memory << endl;
                    // cout << "----------------------Processing subtree!" << endl;
                    IO_sub = IOCounter(subtree, schedule_copy, divisible, quiet, com_freq,
                                       &subtreeBrokenEdges,
                                       method);

                    for (vector<unsigned int>::iterator iter = subtreeBrokenEdges.begin();
                         iter != subtreeBrokenEdges.end(); ++iter) {
                        //  cout<<"broken edge from subtree "<<subtree->getTask(*iter)->getId()
                        // << " otherside id "<<subtree->getTask(*iter)->getOtherSideId()
                        //   << "in tree " <<tree->getTask(*iter)->getOtherSideId()<<endl;
                        brokenEdges->push_back(tree->getTask(*iter)->getOtherSideId());
                    }
                    //   cout << "----------------------Out of Processing subtree!" << endl;
                }

                delete[] schedule_copy;
                delete subtree;

                io_volume += IO_sub;
            } else {

                double node_cost = currentTask->getEdgeWeight() + currentTask->getNodeWeight();
                for (Task *child: *currentTask->getChildren()) {
                    node_cost += child->getEdgeWeight();
                }

                double data_to_unload = memory_occupation + node_cost - currentTask->getEdgeWeight() -
                                        Cluster::getFixedCluster()->getFirstFreeProcessorOrSmallest()->getMemorySize();

                if (!quiet) {
                    cerr << "min data to unload " << data_to_unload << endl;
                }

                if (data_to_unload > 0) {
                    //cerr<<"We must commit I/O in order to process node "<<cur_task_id<<" which requires "<< memory_occupation + node_cost - ewghts[cur_task_id]<< " but has "<<available_memory<<"available"<<endl;
                    /*if we dont have enough room, unload files and update both io and occupation*/

                    switch (method) {
                        case FIRST_FIT:
                            loaded_nodes.remove(make_pair(cur_task_id, schedule_vec.end() -
                                                                       find(schedule_vec.begin(), schedule_vec.end(),
                                                                            cur_task_id)));
                            loaded_nodes.sort(sort_sche); //descending schedule order
                            break;
                        case LARGEST_FIT:
                            loaded_nodes_ew.remove(make_pair(cur_task_id, currentTask->getEdgeWeight()));
                            loaded_nodes_ew.sort(sort_ew);
                            break;

                        default:
                            break;
                    }

                    double unloaded_data = 0.0;
                    switch (method) {
                        case FIRST_FIT:
                            unloaded_data = unload_furthest_first_fit(tree, unloaded_nodes, loaded_nodes,
                                                                      data_to_unload, divisible);
                            break;
                        case LARGEST_FIT:
                            unloaded_data = unload_largest_first_fit(tree, unloaded_nodes, loaded_nodes_ew,
                                                                     data_to_unload);
                            break;
                        default:
                            unloaded_data = unload_furthest_first_fit(tree, unloaded_nodes, loaded_nodes,
                                                                      data_to_unload, divisible);
                            break;
                    }
                    io_volume += unloaded_data;
                    memory_occupation -= unloaded_data;
                }

                if (!quiet) {
                    cerr << "occupation before processing " << memory_occupation << endl;
                }
                /*if we have enough memory to process the node, update occupation*/
                memory_occupation += node_cost - 2 * currentTask->getEdgeWeight() - currentTask->getNodeWeight();
                memory_occupation = max(0.0, memory_occupation);

                if (!quiet) {
                    cerr << "processing " << cur_task_id << endl;
                    cerr << "unloading " << cur_task_id << endl;
                    cerr << "loading ";
                }

                switch (method) {
                    case FIRST_FIT:
                        loaded_nodes.remove(make_pair(cur_task_id, schedule_vec.end() -
                                                                   find(schedule_vec.begin(), schedule_vec.end(),
                                                                        cur_task_id)));
                        break;
                    case LARGEST_FIT:
                        loaded_nodes_ew.remove(make_pair(cur_task_id, currentTask->getEdgeWeight()));
                        break;

                    default:
                        break;
                }

                for (Task *child: *currentTask->getChildren()) {
                    if (!quiet) {
                        cerr << child->getId() << " ";
                    }
                    pair<unsigned int, long> idAndPosition;
                    long positionOfTaskInMinimalTraversal;
                    //vector<int>::iterator edgeWeightOnSpot;

                    switch (method) {
                        case FIRST_FIT:
                            positionOfTaskInMinimalTraversal = schedule_vec.end() -
                                                               find(schedule_vec.begin(),
                                                                    schedule_vec.end(), child->getId());
                            idAndPosition = make_pair(
                                    child->getId(), positionOfTaskInMinimalTraversal);
                            loaded_nodes.push_back(idAndPosition);
                            break;
                        case LARGEST_FIT: {
                            idAndPosition = make_pair(
                                    child->getId(), child->getEdgeWeight());
                            loaded_nodes_ew.push_back(idAndPosition);
                        }
                            break;

                        default:
                            break;
                    }
                }


                if (!quiet) {
                    cerr << endl;
                }
                if (!quiet) {
                    cerr << "New occupation after processing " << memory_occupation << endl;
                }
            }
        }
    }
    delete schedule_f;

    return io_volume;
    //    cerr<<"IO Volume "<<io_volume<<endl;
}

double
IOCounterWithVariableMem(Tree *tree, int *schedule,
                         Cluster *cluster, bool divisible, int quiet, unsigned int &com_freq,
                         vector<unsigned int> *brokenEdges, io_method_t method) {
    throw new runtime_error("not implemented");
}


double MaxOutDegree(Tree *tree, int quiet) {
    double max_out = 0;
    double max_j = 0;
    for (unsigned int j = 1; j <= tree->getTasks()->size(); j++) {
        ////cout<<j<<endl;
        double cur_out = tree->getTask(j)->getCost();
        if (cur_out >= max_out) {
            max_out = cur_out;
            max_j = tree->getTask(j)->getId();
        }
    }
    if (!quiet) {
        cerr << "Max out degree " << max_out << " at " << max_j << endl;
    }
    return max_out;
}

Tree *SubtreeRooted(Task *node) {
    Tree *subtree = new Tree();

    subtree->setRootId(1);
    subtree->setTreeId(node->getId());
    subtree->addTask(node);

    vector<Task *> visit_next;
    Task *end_node;
    if (node->isLeaf()) {
        return subtree;
    } else {
        visit_next = *(node->getChildren());
        while (!visit_next.empty()) {
            if (!visit_next.back()->isBroken()) { // this child has not been cut
                subtree->addTask(visit_next.back());
                end_node = visit_next.back();
                visit_next.pop_back();
                if (!end_node->isLeaf()) {
                    visit_next.insert(visit_next.end(), end_node->getChildren()->begin(),
                                      end_node->getChildren()->end());
                }
            } else {
                visit_next.pop_back();
            }
        }
        return subtree;
    }
}

// Subtree starting from subtreeRoot; goes down until it meets broken edges or leaves
//Subtree is a separate Tree entity with copies of original Tasks
//Each Task in Subtree knows the Task in the big Tree that it had been copied from
//The relation between a Task in the Subtree and its cunterpart  in the big Tree is OtherSideId
Tree *BuildSubtree(Tree *tree, Task *subtreeRoot) {

    pair<Task *, Task *> currentToBeExplored;
    list<pair<Task *, Task * >> toBeExplored;
    auto *tasksInNewSubtree = new vector<Task *>();
    Task *copy;
    Task *rootCopy;
    //Copy of subtreeRoot will be root in the Subtree and will have no parent
    copy = new Task(*subtreeRoot, 1, nullptr);
    rootCopy = copy;
    tasksInNewSubtree->push_back(copy);
    //original subtreeRoot knows that copy, whose id is 1
    subtreeRoot->setOtherSideId(1);
    //copy knows the id of the original subtreeRoot
    copy->setOtherSideId(subtreeRoot->getId());
    //the children of subtreeRoot need to be explored
    toBeExplored.emplace_back(subtreeRoot, copy);

    unsigned int idInSubtree = 1;
    while (!toBeExplored.empty()) {
        currentToBeExplored = toBeExplored.front();
        toBeExplored.pop_front();
        //for all children of the node in the big Tree
        for (Task *child: *currentToBeExplored.first->getChildren()) {
            if (!child->isBroken()) {
                // broken edge means node is on another subtree
                idInSubtree++;
                copy = new Task(*child, idInSubtree, currentToBeExplored.second);
                child->setOtherSideId(idInSubtree);
                copy->setOtherSideId(child->getId());

                tasksInNewSubtree->push_back(copy);

                //add as child to the copied parent, whose counterpart we are currently exploring
                currentToBeExplored.second->addChild(copy);
                //need to explore their children
                toBeExplored.emplace_back(child, copy);

            }
        }
    }

    Tree *treeobj = new Tree(tasksInNewSubtree, rootCopy, tree);
    return treeobj;
}

#endif