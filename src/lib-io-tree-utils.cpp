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

#include "../include/lib-io-tree-utils.h"
#include "../include/lib-io-tree.h"
#include "../include/lib-io-tree-minmem.h"
#include "../include/lib-io-tree-free-methods.h"
#include "../include/heuristics.h"
#include <sys/time.h>
#include <algorithm>

#ifndef CLUSTER_H
#define CLUSTER_H

using namespace std;

Tree *  Tree::originalTree = NULL;
bool Tree::originalTreeInitialized = false;

double BANDWIDTH = 1;

bool sort_sche(node_sche a, node_sche b)
{
    return (a.second > b.second);
}

bool sort_ew(node_ew a, node_ew b)
{
    return (a.second > b.second);
}

double u_wseconds(void)
{
    struct timeval tp;

    gettimeofday(&tp, NULL);

    return (double)tp.tv_sec + (double)tp.tv_usec / 1000000.0;
};


///Qtree corresponds to a whole original tree
Tree * Tree::BuildQtree()
{ //Qtree is for makespan side, so do not use it for space side
    Task *root = this->GetRoot();
    root->BreakEdge();
    this->GetRoot()->GetMSCost(true, true); //update
    size_t tree_size = this->GetNodes()->size();
    unsigned long num_subtrees = this->HowmanySubtrees(true);

    int *prnts = new int[num_subtrees + 1];
    double *ewghts = new double[num_subtrees + 1];
    double *timewghts = new double[num_subtrees + 1];
    int *brokenEdges = new int[num_subtrees + 1];

    //creat Quotient tree
    brokenEdges[1] = 1; //root node
    prnts[1] = 0;
    ewghts[1] = 0;
    timewghts[1] = root->GetSequentialPart();
    unsigned int j = 2;
    root->SetothersideID(1);

    Task *currentNode;
    for (unsigned int i = 2; i <= tree_size; ++i)
    {
        currentNode = this->GetNode(i);
        if (currentNode->IsBroken())
        {
            currentNode->SetothersideID(j); //corresponding node's ID on Qtree
            brokenEdges[j] = i;
            timewghts[j] = currentNode->GetSequentialPart();
            ewghts[j] = currentNode->GetEW();
            ++j;
        }
    }

    for (unsigned int i = 2; i <= num_subtrees; ++i)
    {
        currentNode = this->GetNode(brokenEdges[i])->GetParent();
        while (!currentNode->IsBroken())
        {
            currentNode = currentNode->GetParent();
        }
        prnts[i] = currentNode->GetothersideID();
    }

    Tree *Qtreeobj = new Tree(num_subtrees, prnts, timewghts, ewghts, timewghts); //Qtree only reprents makespan, not memory consumption

    for (unsigned int i = 1; i <= num_subtrees; i++)
    {
        Qtreeobj->GetNode(i)->BreakEdge();                    //break edge
        Qtreeobj->GetNode(i)->SetothersideID(brokenEdges[i]); //corresponding node's ID on tree
    }

    delete[] prnts;
    delete[] ewghts;
    delete[] timewghts;
    delete[] brokenEdges;

    return Qtreeobj;
}

unsigned int Tree::HowmanySubtrees(bool quiet)
{
    unsigned int number_subtrees = 0;
    this->GetRoot()->BreakEdge();
    const vector<Task *> *Nodes = this->GetNodes();
    if (quiet == false)
    {
        cout << "Broken Edges { ";
    }
    for (auto it = Nodes->begin(); it != Nodes->end(); ++it)
    {
        if ((*it)->IsBroken())
        {
            number_subtrees++;
            if (quiet == false)
            {
                cout << (*it)->GetId() << " ";
            }
        }
    }
    if (quiet == false)
    {
        cout << "}" << endl;
    }
    return number_subtrees;
}


bool Tree::increaseMS(Tree *Qtree, Task *&smallestNode, int *chstart, int *childrenID, double memory_size, bool CheckMemory)
{

    Task *currentNode;
    double diff, increase, temp;
    bool memoryEnough;
    bool feasible = false;
    Task *LargestNode;
    Task *secondLargest;
    Task *parent;
    double smallestIncrease = this->GetRoot()->GetMSCost(true, false);
    bool leaf = false;
    const vector<Task *> *subtrees = Qtree->GetNodes();
    vector<Task *> *children;

    if (subtrees->front()->GetId() != 1)
    {
        cout << "error in function increaseMs" << endl;
        return false;
    }

    vector<Task *>::const_iterator iter = subtrees->begin();
    ++iter;
    for (; iter != subtrees->end(); ++iter)
    {
        currentNode = (*iter);

        if (this->GetNode(currentNode->GetothersideID())->IsBroken() == true)
        { //this subtree has not been merged yet
            children = currentNode->GetChildren();

            if (children->empty())
            {
                leaf = true;
            }

            //check the memory cost, if merge itself to its parent subtree
            if (CheckMemory == true)
            {
                memoryEnough = this->MemoryEnough(currentNode->GetParent(), currentNode, leaf, memory_size, chstart, childrenID);
            }
            else
            {
                memoryEnough = true;
            }

            //cout<<"   subtree "<<currentNode->GetothersideID()<<" ";//print subtree's root id
            if (memoryEnough == true)
            {
                feasible = true;
                //cout<<"memory fit."<<endl;
                if (children->empty())
                {
                    children = currentNode->GetParent()->GetChildren();
                    if (children->size() == 2)
                    {
                        increase = children->front()->GetMSCost(false, false) + children->back()->GetMSCost(false, false) - currentNode->GetParent()->GetParallelPart();
                    }
                    else if (children->size() == 1)
                    {
                        increase = -currentNode->GetEW() / BANDWIDTH;
                    }
                    else
                    {
                        //cout<<"   current subtree "<<currentNode->GetothersideID()<<", parent id "<<currentNode->GetParent()->GetothersideID()<<", number of siblings "<<children->size()-1<<endl;

                        GetTwoLargestElementTypetwo(children, LargestNode, secondLargest); //no-increasing, communication counted
                        if (currentNode->GetMSCost(true, false) == LargestNode->GetMSCost(true, false))
                        {
                            if (currentNode->GetMSCost(true, false) == secondLargest->GetMSCost(true, false))
                            {
                                increase = currentNode->GetMSW();
                            }
                            else
                            {
                                increase = -currentNode->GetEW() / BANDWIDTH + secondLargest->GetMSCost(true, false);
                            }
                        }
                        else
                        {
                            increase = currentNode->GetMSW();
                        }
                    }
                }
                else
                {
                    children = currentNode->GetParent()->GetChildren();
                    //cout<<"   current subtree "<<currentNode->GetothersideID()<<", parent id "<<currentNode->GetParent()->GetothersideID()<<", number of siblings "<<children->size()-1<<endl;
                    diff = currentNode->GetMSCost(true, false) - currentNode->GetParent()->GetParallelPart();
                    if (diff < 0)
                    {
                        increase = currentNode->GetMSW();
                    }
                    else
                    { //diff=0
                        if (children->size() == 1)
                        {
                            increase = -currentNode->GetEW() / BANDWIDTH;
                        }
                        else
                        {                                                                      //children's size larger than 1
                            GetTwoLargestElementTypetwo(children, LargestNode, secondLargest); //no-increasing, communication counted
                            temp = currentNode->GetParallelPart() - secondLargest->GetMSCost(true, false);
                            if (temp >= 0)
                            {
                                increase = -currentNode->GetEW() / BANDWIDTH;
                            }
                            else
                            {
                                increase = -temp - currentNode->GetEW() / BANDWIDTH;
                            }
                        }
                    }
                }

                parent = currentNode->GetParent();
                while (parent->GetId() != 1)
                { //not the root node
                    temp = parent->GetParent()->GetParallelPart() - (parent->GetMSCost(true, false) + increase);
                    if (temp >= 0)
                    {
                        increase = 0;
                        break;
                    }
                    else
                    {
                        increase = -temp;
                        parent = parent->GetParent();
                    }
                }

                //cout<<"   merge, increase in MS(r) "<<increase<<endl;
                if (increase < smallestIncrease)
                {
                    smallestIncrease = increase;
                    smallestNode = currentNode;
                }
            }
            else
            {
                //cout<<"memory does not fit!!!"<<endl;
            }
        }
    }

    //cout<<"   ---end compute the minimum combination"<<endl;
    return feasible;
}


bool Tree::MemoryEnough(Task *Qrootone, Task *Qroottwo, bool leaf, double memory_size, int *chstart, int *children)
{
    bool enough = false;
    unsigned long new_tree_size = this->GetNodes()->size();

    Task *SubtreeRoot = this->GetNode(Qrootone->GetothersideID());

    vector<Task *> *childrenvector = Qrootone->GetChildren();
    if ((leaf == true) & (childrenvector->size() == 2))
    {
        this->GetNode(childrenvector->front()->GetothersideID())->RestoreEdge();
        this->GetNode(childrenvector->back()->GetothersideID())->RestoreEdge();
    }
    else
    {
        this->GetNode(Qroottwo->GetothersideID())->RestoreEdge(); //restore edge temporarilly
    }

    double *ewghts, *timewghts, *spacewghts;
    int *prnts;
    Tree *subtree = BuildSubtree(this, SubtreeRoot);
    delete[] ewghts;
    delete[] timewghts;
    delete[] spacewghts;
    delete[] prnts;
    double maxout, requiredMemory;
    uint64_t count = 0;
    schedule_t *schedule_f = new schedule_t();
    maxout = MaxOutDegree(subtree, true);
    MinMem(subtree, maxout, requiredMemory, *schedule_f, true, count);

    if (requiredMemory <= memory_size)
    {
        enough = true;
    }

    if ((leaf == true) & (childrenvector->size() == 2))
    {
        this->GetNode(childrenvector->front()->GetothersideID())->BreakEdge();
        this->GetNode(childrenvector->back()->GetothersideID())->BreakEdge();
    }
    else
    {
        this->GetNode(Qroottwo->GetothersideID())->BreakEdge();
    }

    delete subtree;
    delete schedule_f;

    return enough;
}



double Task::Sequence()
{
    return this->GetMSCost();
}


void parse_tree(const char *filename, int *N, int **prnts, double **nwghts, double **ewghts, double **mswghts)
{   
    ifstream OpenFile(filename);
    char begin;
    char cur_char;
    int line_index = 0;
    bool nodes_cnt_read = false;
    unsigned int nb_of_nodes = 0;
    string line;

    do
    {
        /*skip commentary lines*/
        begin = OpenFile.peek();
        if (OpenFile.good())
        {

            if (begin == '%')
            {
                do
                {
                    cur_char = OpenFile.get();
                } while (cur_char != '\n' && OpenFile.good());
            }
            else
            {
                if (!nodes_cnt_read)
                {
                    /* get the number of nodes and skip last trailing character*/
                    while (getline(OpenFile, line))
                    {
                        ++nb_of_nodes;
                    }
                    OpenFile.clear();
                    OpenFile.seekg(0, ios::beg);
                    *N = nb_of_nodes;
                    nodes_cnt_read = true;
                    /*allocate space for nodes*/
                    *prnts = new int[nb_of_nodes + 1];
                    *nwghts = new double[nb_of_nodes + 1];
                    *ewghts = new double[nb_of_nodes + 1];
                    *mswghts = new double[nb_of_nodes + 1];
                }
                else
                {
                    /*parse actual nodes*/
                    unsigned int id;
                    unsigned int parent;
                    double ew, nw, msw;

                    OpenFile >> id >> parent >> nw >> msw >> ew;
                    do
                    {
                        cur_char = OpenFile.get();
                    } while (cur_char != '\n' && OpenFile.good());
                    parent = nb_of_nodes - parent + 1; //root has the largest id in the txt file
                   // cout<<"nbnodes "<<nb_of_nodes <<"parent "<< parent<<" line index "<<line_index<<" resulting index "<<nb_of_nodes - line_index<<endl;
                    (*prnts)[nb_of_nodes - line_index] = parent;
                    (*nwghts)[nb_of_nodes - line_index] = nw;
                    (*ewghts)[nb_of_nodes - line_index] = ew;
                    (*mswghts)[nb_of_nodes - line_index] = msw;

                    line_index++;
                }
            }
        }
    } while (OpenFile.good());
    (*prnts)[1] = 0;

    OpenFile.close();
}

void poaux(const int *chstart, const int *children, int N, int r, int *por, int *label)
{
    int *stack = new int[N + 2];
    memset(stack, 0, (N + 2) * sizeof(*stack));

    int push = 1;
    stack[push] = r;
    push++;
    int *chtmp = (int *)malloc((N + 2) * sizeof(int));
    memcpy(chtmp, chstart, (N + 2) * sizeof(*chstart));

    while (push > 1)
    {
        int top = stack[push - 1];
        if (chstart[top + 1] - chtmp[top] > 0)
        {
            int ch = children[chtmp[top]];
            chtmp[top]++;
            stack[push] = ch;
            push++;
        }
        else
        {
            por[*label] = top;
            (*label)++;
            push--;
        }
    }

    delete[] chtmp;
    delete[] stack;
}

void po_construct(const int treeSize, const int *prnts, int **chstart, int **chend, int **children, int *root)
{   
    *chstart = new int[treeSize + 2];
    *children = new int[treeSize + 1];
     *chend = new int[treeSize + 2];
    memset((void *)*chstart, 0, (treeSize + 2) * sizeof(**chstart));
    memset((void *)*children, 0, (treeSize + 1) * sizeof(**children));

    *root = -1;

    for (int ii = 1; ii < treeSize + 1; ii++)
    {
        if (prnts[ii] > 0)
        {
            (*chstart)[prnts[ii]]++;
        }
        else
        {
            *root = ii;
        }
    }

    /*compute cumsum*/
    int cum_val = 1;
    for (int ii = 1; ii < treeSize + 2; ii++)
    {
        int val = cum_val;
        cum_val += (*chstart)[ii];
        (*chstart)[ii] = val;
    }

    memcpy(*chend, *chstart, (treeSize + 2) * sizeof(**chstart));

    for (int ii = 1; ii < treeSize + 1; ii++)
    {
        if (prnts[ii] > 0)
        {
            (*children)[(*chend)[prnts[ii]]] = ii;
            (*chend)[prnts[ii]]++;
        }
    }
}

double IOCounter(Tree &tree, schedule_t &sub_schedule, double available_memory, bool divisible, int quiet)
{
    double memory_occupation = 0;
    double io_volume = 0;
    io_map unloaded_nodes;
    schedule_t loaded_nodes;

    /*iterates through the given permutation (schedule)*/
    for (schedule_t::iterator cur_task_id = sub_schedule.begin(); cur_task_id != sub_schedule.end(); cur_task_id++)
    {
        Task *cur_node = tree.GetNode(*cur_task_id);

        /*if the node was unloaded*/
        if (unloaded_nodes.find(*cur_task_id) != unloaded_nodes.end())
        {
            if (!quiet)
            {
                cerr << "Loading " << unloaded_nodes[*cur_task_id] << "of " << *cur_task_id << " (IO)" << endl;
            }

            memory_occupation += unloaded_nodes[*cur_task_id];
            unloaded_nodes.erase(*cur_task_id);
        }

        double data_to_unload = memory_occupation + cur_node->GetCost() - cur_node->GetEW() - available_memory;
        if (!quiet)
        {
            cerr << "min data to unload " << data_to_unload << endl;
        }
        if (data_to_unload > 0)
        {
            /*if we dont have enough room, unload files and update both io and occupation*/
            double unloaded_data = 0;
            /*unload furthest non unloaded node which is NOT in current_node children first*/
            schedule_t::reverse_iterator far_node_id = loaded_nodes.rbegin();
            while ((far_node_id != loaded_nodes.rend()) && (unloaded_data < data_to_unload))
            {
                /*try to unload this node*/
                bool is_already_unloaded = false;
                double remaining_loaded_data = tree.GetNode(*far_node_id)->GetEW();
                if (unloaded_nodes.find(*far_node_id) != unloaded_nodes.end())
                {
                    is_already_unloaded = true;
                    remaining_loaded_data = max(remaining_loaded_data - unloaded_nodes[*far_node_id], 0.0);
                }

                double local_data_to_unload;
                if (divisible)
                {
                    local_data_to_unload = min(remaining_loaded_data, data_to_unload);
                }
                else
                {
                    local_data_to_unload = remaining_loaded_data;
                }

                if (!quiet)
                {
                    cerr << "unloading (IO) " << local_data_to_unload << " of " << *far_node_id << endl;
                }
                unloaded_data += local_data_to_unload;
                if (is_already_unloaded)
                {
                    unloaded_nodes[*far_node_id] += local_data_to_unload;
                }
                else
                {
                    unloaded_nodes[*far_node_id] = local_data_to_unload;
                }

                if (remaining_loaded_data == local_data_to_unload)
                {
                    loaded_nodes.remove(*far_node_id);
                }

                far_node_id++;
            }

            io_volume += unloaded_data;
            memory_occupation -= unloaded_data;
        }

        if (!quiet)
        {
            cerr << "occupation before processing " << memory_occupation << endl;
        }
        /*if we have enough memory to process the node, update occupation*/
        memory_occupation += cur_node->GetCost() - 2 * cur_node->GetEW() - cur_node->GetNW();
        if (!quiet)
        {
            cerr << "processing " << *cur_task_id << endl;
        }
        if (!quiet)
        {
            cerr << "unloading " << *cur_task_id << endl;
        }
        loaded_nodes.remove(*cur_task_id);

        if (!quiet)
        {
            cerr << "loading ";
        }
        for (vector<Task *>::iterator child = cur_node->GetChildren()->begin(); child != cur_node->GetChildren()->end(); child++)
        {
            if (!quiet)
            {
                cerr << (*child)->GetId() << " ";
            }
            loaded_nodes.push_back((*child)->GetId());
        }
        if (!quiet)
        {
            cerr << endl;
        }
        if (!quiet)
        {
            cerr << "New occupation after processing " << memory_occupation << endl;
        }
    }

    return io_volume;
    //    cerr<<"IO Volume "<<io_volume<<endl;
}

double unload_largest_first_fit(Tree *tree, vector<unsigned int> &unloaded_nodes, list<node_ew> &loaded_nodes, const double data_to_unload, double *ewghts)
{
    double unloaded_data = 0.0;

    /*loaded_nodes already sorted, unload largest nodes till there is enough space*/
    list<node_ew>::iterator largest_node = loaded_nodes.begin();
    while ((largest_node != loaded_nodes.end()) && (unloaded_data < data_to_unload))
    {
        largest_node = loaded_nodes.begin();
        tree->GetNode(largest_node->first)->BreakEdge(); //break this edge;
        //cout<<"******Largest first: break edge "<<largest_node->first<<endl;
        double edge_weight = largest_node->second;
        unloaded_data += edge_weight;
        unloaded_nodes.push_back(largest_node->first);
        loaded_nodes.pop_front();
    }

    return unloaded_data;
}

double unload_furthest_nodes(Tree *tree, vector<unsigned int> &unloaded_nodes, list<node_sche> &loaded_nodes, const double data_to_unload, double *ewghts, bool divisible)
{
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
    while ((far_node != loaded_nodes.end()) && (unloaded_data < data_to_unload))
    {
        // cout << far_node->first << " " << far_node->second << endl;
        /*try to unload this node*/
        far_node = loaded_nodes.begin();
        double remaining_loaded_data = ewghts[(*far_node).first];
        double local_data_to_unload;
        if (divisible)
        {
            local_data_to_unload = min(remaining_loaded_data, data_to_unload);
        }
        else
        {
            local_data_to_unload = remaining_loaded_data;
        }
#if VERBOSE
        //cerr<<"unloading (IO) "<<local_data_to_unload<<" of "<<*far_node<<endl;
#endif
        unloaded_data += local_data_to_unload;

        unloaded_nodes.push_back((*far_node).first);

        ////cout<<"-------LSNF remove "<<local_data_to_unload<<endl;
        ////cout<<"break edge "<<far_node->first<<endl;
        if (far_node->first == 0)
        {
            cout << "Problem! loaded nodes " << endl;
            list<node_sche>::iterator loaded_iterator = loaded_nodes.begin();
            if (loaded_nodes.size() < 20)
            {
                while (loaded_iterator != loaded_nodes.end())
                {
                    cout << loaded_iterator->first << endl;
                }
            }
            else
            {
                cout << "big loaded nodes";
            }
            cout << "unloaded nodes: ";
            for (int i = 0; i < unloaded_nodes.size(); i++)
            {
                cout << unloaded_nodes[i] << endl;
            }

            cout << "old loaded nodes size " << old_loaded_nodes.size() << endl;
            cout << "old loaded nodes first " << old_loaded_nodes.begin()->first << endl;
            tree->GetNode(far_node->first)->BreakEdge(); //break this edge
            loaded_nodes.pop_front();
        }

        // cout << "after pop " << far_node->first << " " << far_node->second << endl;
    }

    return unloaded_data;
}

double unload_furthest_first_fit(Tree *tree, vector<unsigned int> &unloaded_nodes, list<node_sche> &loaded_nodes, const double data_to_unload, double *ewghts, bool divisible)
{
    double unloaded_data = 0.0;
    /*unload furthest non unloaded node which is NOT in current_node children first*/
    list<node_sche>::iterator far_node = loaded_nodes.begin();
    while ((far_node != loaded_nodes.end()) && (unloaded_data < data_to_unload))
    {
        /*try to unload this node*/

        double remaining_loaded_data = ewghts[far_node->first];

        double local_data_to_unload;
        if (divisible)
        {
            local_data_to_unload = min(remaining_loaded_data, data_to_unload);
        }
        else
        {
            local_data_to_unload = remaining_loaded_data;
        }
#if VERBOSE
        //cerr<<"unloading (IO) "<<local_data_to_unload<<" of "<<*far_node<<endl;
#endif

        /*if it "fits", that is if the amount of data is lower than what we need to unload*/
        if (local_data_to_unload >= data_to_unload)
        {
            //cout<<"-------First fit success: ";
            //cout<<"break edge "<<far_node->first<<endl;
            unloaded_data += local_data_to_unload;
            unloaded_nodes.push_back(far_node->first);
            tree->GetNode(far_node->first)->BreakEdge(); //break this edge
            loaded_nodes.erase(far_node);
            return unloaded_data;
        }
        far_node++;
    }

    /*if we did not unloaded enough data, call furthest node*/
    ////cout<<"-------First fit failed, go to LSNF: "<<endl;
    unloaded_data = unload_furthest_nodes(tree, unloaded_nodes, loaded_nodes, data_to_unload, ewghts, divisible);

    return unloaded_data;
}

double unload_furthest_best_fit(io_map &unloaded_nodes, schedule_t &loaded_nodes, const double data_to_unload, double *ewghts, bool divisible)
{
    double unloaded_data = 0.0;
    /*unload furthest non unloaded node which is NOT in current_node children first*/
    while (unloaded_data < data_to_unload)
    {
        bool is_best_already_unloaded = false;
        double best_remaining_loaded_data = 0.0;
        unsigned int best_candidate = -1;
        double best_candi_score = -1.0;

        schedule_t::reverse_iterator far_node_id = loaded_nodes.rbegin();
        while ((far_node_id != loaded_nodes.rend()))
        {
            /*try to unload this node*/
            bool is_already_unloaded = false;
            double remaining_loaded_data = ewghts[*far_node_id];
            if (unloaded_nodes.find(*far_node_id) != unloaded_nodes.end())
            {
                is_already_unloaded = true;
                remaining_loaded_data = max(remaining_loaded_data - unloaded_nodes[*far_node_id], 0.0);
            }
            double local_data_to_unload;
            if (divisible)
            {
                local_data_to_unload = min(remaining_loaded_data, data_to_unload);
            }
            else
            {
                local_data_to_unload = remaining_loaded_data;
            }
            /*if it "fits", that is if the amount of data is lower than what we need to unload*/
            if (local_data_to_unload <= data_to_unload - unloaded_data)
            {
                if (local_data_to_unload > best_candi_score)
                {
                    best_candi_score = local_data_to_unload;
                    best_candidate = *far_node_id;
                    is_best_already_unloaded = is_already_unloaded;
                    best_remaining_loaded_data = remaining_loaded_data;
                }
            }

            far_node_id++;
        }

        if (best_candi_score != -1)
        {
            //cerr<<"best found :"<<best_candidate<<endl;
            unloaded_data += best_candi_score;
            if (is_best_already_unloaded)
            {
                unloaded_nodes[best_candidate] += best_candi_score;
            }
            else
            {
                unloaded_nodes[best_candidate] = best_candi_score;
            }
            if (best_remaining_loaded_data == best_candi_score)
            {
                loaded_nodes.remove(best_candidate);
            }
        }
        else
        {
            break;
        }
    }
    return unloaded_data;
}

double unload_furthest_first_fit_abs(io_map &unloaded_nodes, schedule_t &loaded_nodes, const double data_to_unload, double *ewghts, bool divisible)
{
    double unloaded_data = 0.0;
    /*unload furthest non unloaded node which is NOT in current_node children first*/
    schedule_t::reverse_iterator far_node_id = loaded_nodes.rbegin();
    while ((far_node_id != loaded_nodes.rend()) && (unloaded_data < data_to_unload))
    {
        /*try to unload this node*/
        bool is_already_unloaded = false;
        double remaining_loaded_data = ewghts[*far_node_id];
        if (unloaded_nodes.find(*far_node_id) != unloaded_nodes.end())
        {
            is_already_unloaded = true;
            remaining_loaded_data = max(remaining_loaded_data - unloaded_nodes[*far_node_id], 0.0);
        }
        double local_data_to_unload;
        if (divisible)
        {
            local_data_to_unload = min(remaining_loaded_data, data_to_unload);
        }
        else
        {
            local_data_to_unload = remaining_loaded_data;
        }
        /*if it "fits", that is if the amount of data is lower than what we need to unload*/
        if (local_data_to_unload >= data_to_unload - unloaded_data)
        {
            unloaded_data += local_data_to_unload;
            if (is_already_unloaded)
            {
                unloaded_nodes[*far_node_id] += local_data_to_unload;
            }
            else
            {
                unloaded_nodes[*far_node_id] = local_data_to_unload;
            }
            if (remaining_loaded_data == local_data_to_unload)
            {
                loaded_nodes.remove(*far_node_id);
            }
        }

        far_node_id++;
    }

    return unloaded_data;
}

double unload_furthest_best_fit_abs(io_map &unloaded_nodes, schedule_t &loaded_nodes, const double data_to_unload, double *ewghts, bool divisible)
{
    double unloaded_data = 0.0;
    /*unload furthest non unloaded node which is NOT in current_node children first*/
    while (unloaded_data < data_to_unload)
    {
        bool is_best_already_unloaded = false;
        double best_remaining_loaded_data = 0.0;
        unsigned int best_candidate = -1;
        double best_candi_score = -1.0;

        schedule_t::reverse_iterator far_node_id = loaded_nodes.rbegin();
        while ((far_node_id != loaded_nodes.rend()))
        {
            /*try to unload this node*/
            bool is_already_unloaded = false;
            double remaining_loaded_data = ewghts[*far_node_id];
            if (unloaded_nodes.find(*far_node_id) != unloaded_nodes.end())
            {
                is_already_unloaded = true;
                remaining_loaded_data = max(remaining_loaded_data - unloaded_nodes[*far_node_id], 0.0);
            }
            double local_data_to_unload;
            if (divisible)
            {
                local_data_to_unload = min(remaining_loaded_data, data_to_unload);
            }
            else
            {
                local_data_to_unload = remaining_loaded_data;
            }
            /*if it "fits", that is if the amount of data is lower than what we need to unload*/
            if (abs(data_to_unload - unloaded_data - local_data_to_unload) < abs(data_to_unload - unloaded_data - best_candi_score))
            {
                best_candi_score = local_data_to_unload;
                best_candidate = *far_node_id;
                is_best_already_unloaded = is_already_unloaded;
                best_remaining_loaded_data = remaining_loaded_data;
            }

            far_node_id++;
        }

        if (best_candi_score != -1)
        {
            //cerr<<"best found :"<<best_candidate<<endl;
            unloaded_data += best_candi_score;
            if (is_best_already_unloaded)
            {
                unloaded_nodes[best_candidate] += best_candi_score;
            }
            else
            {
                unloaded_nodes[best_candidate] = best_candi_score;
            }
            if (best_remaining_loaded_data == best_candi_score)
            {
                loaded_nodes.remove(best_candidate);
            }
        }
        else
        {
            break;
        }
    }

    return unloaded_data;
}

int next_comb(int comb[], int k, int n)
{
    int i = k - 1;
    ++comb[i];

    while ((i > 0) && (comb[i] >= n - k + 1 + i))
    {
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

double unload_best_increasing_combi(io_map &unloaded_nodes, schedule_t &loaded_nodes, const double data_to_unload, double *ewghts, bool divisible, unsigned int init_combi_size, bool quiet)
{
    assert(!divisible);
    vector<unsigned int> candidates;
    double unloaded_data = 0.0;
    vector<unsigned int> best_combi;
    double best_combi_score = numeric_limits<double>::infinity();
    while (unloaded_data < data_to_unload)
    {
        /*unload furthest non unloaded node which is NOT in current_node children first*/
        candidates.clear();
        schedule_t::reverse_iterator far_node_id = loaded_nodes.rbegin();
        while ((far_node_id != loaded_nodes.rend()) && (candidates.size() <= init_combi_size))
        {
            /*try to unload this node*/
            if (unloaded_nodes.find(*far_node_id) == unloaded_nodes.end())
            {
                candidates.push_back(*far_node_id);
            }

            far_node_id++;
        }
        /*now we got our max_candidates candidates, compute every permutations*/

        vector<unsigned int> cur_combi;

        int n = candidates.size(); /* The size of the set; for {1, 2, 3, 4} it's 4 */
        for (int k = 1; k <= n; k++)
        {
            /* k is the size of the subsets; for {1, 2}, {1, 3}, ... it's 2 */
            int comb[k]; /* comb[i] is the index of the i-th element in the
                          combination */

            /* Setup comb for the initial combination */
            for (int i = 0; i < k; ++i)
                comb[i] = i;
            /*compute score*/
            cur_combi.clear();
            double cur_unloaded_data = 0.0;
            if (!quiet)
            {
                cerr << "cur_combi is [";
            }
            for (int i = 0; i < k; ++i)
            {
                unsigned int cur_node_id = candidates[comb[i]];
                cur_unloaded_data += ewghts[cur_node_id];
                cur_combi.push_back(cur_node_id);
                if (!quiet)
                {
                    cerr << " " << cur_node_id;
                }
            }
            if (!quiet)
            {
                cerr << "], its score is " << cur_unloaded_data << endl;
            }

            if (abs(data_to_unload - unloaded_data - cur_unloaded_data) < abs(data_to_unload - unloaded_data - best_combi_score))
            {
                if (!quiet)
                {
                    cerr << "cur_combi is better than prev best " << cur_unloaded_data << " vs " << best_combi_score << endl;
                }
                best_combi.assign(cur_combi.begin(), cur_combi.end());
                best_combi_score = cur_unloaded_data;
            }

            /* Generate and print all the other combinations */
            while (next_comb(comb, k, n))
            {
                /*compute score*/
                cur_combi.clear();
                double cur_unloaded_data = 0.0;
                if (!quiet)
                {
                    cerr << "cur_combi is [";
                }
                for (int i = 0; i < k; ++i)
                {
                    unsigned int cur_node_id = candidates[comb[i]];
                    cur_unloaded_data += ewghts[cur_node_id];
                    cur_combi.push_back(cur_node_id);
                    if (!quiet)
                    {
                        cerr << " " << cur_node_id;
                    }
                }
                if (!quiet)
                {
                    cerr << "], its score is " << cur_unloaded_data << endl;
                }

                if (abs(data_to_unload - unloaded_data - cur_unloaded_data) < abs(data_to_unload - unloaded_data - best_combi_score))
                {
                    if (!quiet)
                    {
                        cerr << "cur_combi is better than prev best " << cur_unloaded_data << " vs " << best_combi_score << endl;
                    }
                    best_combi.assign(cur_combi.begin(), cur_combi.end());
                    best_combi_score = cur_unloaded_data;
                }
            }
        }

        unloaded_data = best_combi_score;

        if (unloaded_data < data_to_unload)
        {
            init_combi_size *= 2;
        }
    }
    if (!quiet)
    {
        cerr << "best_combi is [";
    }
    for (unsigned int i = 0; i < best_combi.size(); i++)
    {
        unsigned int cur_id = best_combi[i];
        unloaded_nodes[cur_id] = ewghts[cur_id];
        loaded_nodes.remove(cur_id);
        if (!quiet)
        {
            cerr << " " << cur_id;
        }
    }
    if (!quiet)
    {
        cerr << "], its score is " << best_combi_score << endl;
    }

    if (!quiet)
    {
        cerr << "We've unloaded " << unloaded_data << " out of " << data_to_unload << " so far" << endl;
    }
    unloaded_data = best_combi_score;
    return unloaded_data;
}

double unload_best_furthest_nodes(io_map &unloaded_nodes, schedule_t &loaded_nodes, const double data_to_unload, double *ewghts, bool divisible, unsigned int max_candidates, bool quiet)
{
    assert(!divisible);
    vector<unsigned int> candidates;
    double unloaded_data = 0.0;
    while (unloaded_data < data_to_unload)
    {
        /*unload furthest non unloaded node which is NOT in current_node children first*/
        candidates.clear();
        schedule_t::reverse_iterator far_node_id = loaded_nodes.rbegin();
        while ((far_node_id != loaded_nodes.rend()) && (candidates.size() <= max_candidates))
        {
            /*try to unload this node*/
            if (unloaded_nodes.find(*far_node_id) == unloaded_nodes.end())
            {
                candidates.push_back(*far_node_id);
            }

            far_node_id++;
        }
        /*now we got our max_candidates candidates, compute every permutations*/

        vector<unsigned int> best_combi;
        vector<unsigned int> cur_combi;
        double best_combi_score = numeric_limits<double>::infinity();

        int n = candidates.size(); /* The size of the set; for {1, 2, 3, 4} it's 4 */
        for (int k = 1; k <= n; k++)
        {
            /* k is the size of the subsets; for {1, 2}, {1, 3}, ... it's 2 */
            int comb[k]; /* comb[i] is the index of the i-th element in the
                          combination */

            /* Setup comb for the initial combination */
            for (int i = 0; i < k; ++i)
                comb[i] = i;
            /*compute score*/
            cur_combi.clear();
            double cur_unloaded_data = 0.0;
            if (!quiet)
            {
                cerr << "cur_combi is [";
            }
            for (int i = 0; i < k; ++i)
            {
                unsigned int cur_node_id = candidates[comb[i]];
                cur_unloaded_data += ewghts[cur_node_id];
                cur_combi.push_back(cur_node_id);
                if (!quiet)
                {
                    cerr << " " << cur_node_id;
                }
            }
            if (!quiet)
            {
                cerr << "], its score is " << cur_unloaded_data << endl;
            }

            if ((double)abs(data_to_unload - unloaded_data - cur_unloaded_data) < (double)abs(data_to_unload - unloaded_data - best_combi_score))
            {
                if (!quiet)
                {
                    cerr << "cur_combi is better than prev best " << cur_unloaded_data << " vs " << best_combi_score << endl;
                }
                best_combi.assign(cur_combi.begin(), cur_combi.end());
                best_combi_score = cur_unloaded_data;
            }

            /* Generate and print all the other combinations */
            while (next_comb(comb, k, n))
            {
                /*compute score*/
                cur_combi.clear();
                double cur_unloaded_data = 0.0;
                if (!quiet)
                {
                    cerr << "cur_combi is [";
                }
                for (int i = 0; i < k; ++i)
                {
                    unsigned int cur_node_id = candidates[comb[i]];
                    cur_unloaded_data += ewghts[cur_node_id];
                    cur_combi.push_back(cur_node_id);
                    if (!quiet)
                    {
                        cerr << " " << cur_node_id;
                    }
                }
                if (!quiet)
                {
                    cerr << "], its score is " << cur_unloaded_data << endl;
                }

                if (abs(data_to_unload - unloaded_data - cur_unloaded_data) < abs(data_to_unload - unloaded_data - best_combi_score))
                {
                    if (!quiet)
                    {
                        cerr << "cur_combi is better than prev best " << cur_unloaded_data << " vs " << best_combi_score << endl;
                    }
                    best_combi.assign(cur_combi.begin(), cur_combi.end());
                    best_combi_score = cur_unloaded_data;
                }
            }
        }

        if (!quiet)
        {
            cerr << "best_combi is [";
        }
        for (unsigned int i = 0; i < best_combi.size(); i++)
        {
            unsigned int cur_id = best_combi[i];
            unloaded_nodes[cur_id] = ewghts[cur_id];
            loaded_nodes.remove(cur_id);
            if (!quiet)
            {
                cerr << " " << cur_id;
            }
        }
        if (!quiet)
        {
            cerr << "], its score is " << best_combi_score << endl;
        }

        unloaded_data += best_combi_score;
        if (!quiet)
        {
            cerr << "We've unloaded " << unloaded_data << " out of " << data_to_unload << " so far" << endl;
        }
    }

    return unloaded_data;
}

double IOCounter(Tree *tree, int N, double *nwghts, double *ewghts, int *chstart, int *children, int *schedule, double available_memory,
                 bool divisible, int quiet, unsigned int &com_freq, vector<unsigned int> *brokenEdges, io_method_t method)
{
    double memory_occupation = ewghts[schedule[N - 1]];
    double io_volume = 0;
    vector<unsigned int> unloaded_nodes;
    list<node_sche> loaded_nodes;
    list<node_ew> loaded_nodes_ew;
    vector<int> schedule_vec(schedule, schedule + N);
    int cur_task_id;
    vector<unsigned int>::iterator unloaded;
    list<unsigned int> temp;
    list<unsigned int> queue;
    unsigned int child_start, child_end;
    unsigned int subtree_size;
    list<unsigned int>::iterator iter;
    double maxoutD, memory_required, IO_sub = 0;
    uint64_t count;
    schedule_t *schedule_f = new schedule_t();
    list<int>::iterator ite_sche;
    int rootid;
    vector<unsigned int> subtreeBrokenEdges;

    /*iterates through the given permutation (schedule)*/
    for (int rank = N - 1; rank >= 1; rank--)
    {
        cur_task_id = schedule[rank];

        //cout<<"current task id: "<<cur_task_id<<endl;
        if (cur_task_id != 0)
        { //0 means this node is on other subtrees
            /*if the node was unloaded*/
            unloaded = find(unloaded_nodes.begin(), unloaded_nodes.end(), cur_task_id);
            if (unloaded != unloaded_nodes.end())
            { //find node cur_task_id unloaded
                //cout<<", (break) "<<endl;
                brokenEdges->push_back(tree->GetNode(cur_task_id)->GetothersideID());
                ++com_freq;
                unloaded_nodes.erase(unloaded);
                temp.clear();
                queue.push_back(cur_task_id);
                //cout<<"children "<<endl;
                do
                {
                    child_start = *(chstart + queue.front());
                    child_end = *(chstart + queue.front() + 1);
                    temp.push_back(queue.front());
                    queue.pop_front();
                    for (unsigned int i = child_start; i < child_end; ++i)
                    {
                        //cout<<*(children+i)<<" ";
                        queue.push_back(*(children + i));
                    }
                } while (!queue.empty());
                //cout<<endl;

                subtree_size = temp.size(); //just an approximation
                for (long i = rank - 1; i >= 0; i--)
                {
                    iter = find(temp.begin(), temp.end(), schedule[i]);
                    if (iter != temp.end())
                    {
                        schedule[i] = 0; //IO counter will pass 0;
                        temp.erase(iter);
                    }
                    if (temp.size() == 1)
                    {
                        break;
                    }
                }

                double *ewghtssub, *timewghtssub, *spacewghtssub;
                int *prntssub;
                Tree *subtree = BuildSubtree(tree, tree->GetNode(cur_task_id));

                subtree_size = subtree->GetNodes()->size();
                cout << "subtree size " << subtree_size << endl;

                int *schedule_copy = new int[subtree_size + 1];
                maxoutD = MaxOutDegree(subtree, true);
                schedule_f->clear();
                count = 0;
                MinMem(subtree, maxoutD, memory_required, *schedule_f, true, count);
                ite_sche = schedule_f->begin();
                for (unsigned int i = subtree_size; i >= 1; --i)
                {
                    schedule_copy[i] = *ite_sche;
                    advance(ite_sche, 1);
                }
                schedule_copy[0] = subtree_size + 1;
                int *chstartsub, *childrensub, *chendsub;
                po_construct(subtree_size, prntssub, &chstartsub, &chendsub, &childrensub, &rootid);

                if (memory_required > available_memory)
                {
                    // cout << "memory required " << memory_required << ", is larger than what is available BLABLA " << available_memory << endl;
                    // cout << "----------------------Processing subtree!" << endl;
                    IO_sub = IOCounter(subtree, subtree_size + 1, spacewghtssub, ewghtssub, chstartsub, childrensub, schedule_copy, available_memory, divisible, quiet, com_freq, &subtreeBrokenEdges, method);

                    for (vector<unsigned int>::iterator iter = subtreeBrokenEdges.begin(); iter != subtreeBrokenEdges.end(); ++iter)
                    {
                        brokenEdges->push_back(tree->GetNode(*iter)->GetothersideID());
                    }
                    //   cout << "----------------------Out of Processing subtree!" << endl;
                }

                delete[] ewghtssub;
                delete[] timewghtssub;
                delete[] spacewghtssub;
                delete[] prntssub;
                delete[] chstartsub;
               
                delete[] childrensub;
                delete[] schedule_copy;
                delete subtree;

                io_volume += IO_sub;
            }
            else
            {

                double node_cost = ewghts[cur_task_id] + nwghts[cur_task_id];
                for (int j = chstart[cur_task_id]; j < chstart[cur_task_id + 1]; j++)
                {
                    node_cost += ewghts[children[j]];
                }

                double data_to_unload = memory_occupation + node_cost - ewghts[cur_task_id] - available_memory;

                if (!quiet)
                {
                    cerr << "min data to unload " << data_to_unload << endl;
                }

                if (data_to_unload > 0)
                {
                    //cerr<<"We must commit I/O in order to process node "<<cur_task_id<<" which requires "<< memory_occupation + node_cost - ewghts[cur_task_id]<< " but has "<<available_memory<<"available"<<endl;
                    /*if we dont have enough room, unload files and update both io and occupation*/

                    switch (method)
                    {
                    case FIRST_FIT:
                        loaded_nodes.remove(make_pair(cur_task_id, schedule_vec.end() - find(schedule_vec.begin(), schedule_vec.end(), cur_task_id)));
                        loaded_nodes.sort(sort_sche); //descending schedule order
                        break;
                    case LARGEST_FIT:
                        loaded_nodes_ew.remove(make_pair(cur_task_id, ewghts[cur_task_id]));
                        loaded_nodes_ew.sort(sort_ew);
                        break;

                    default:
                        break;
                    }

                    double unloaded_data = 0.0;
                    switch (method)
                    {
                    case FIRST_FIT:
                        unloaded_data = unload_furthest_first_fit(tree, unloaded_nodes, loaded_nodes, data_to_unload, ewghts, divisible);
                        break;
                    case LARGEST_FIT:
                        unloaded_data = unload_largest_first_fit(tree, unloaded_nodes, loaded_nodes_ew, data_to_unload, ewghts);
                        break;
                    default:
                        unloaded_data = unload_furthest_first_fit(tree, unloaded_nodes, loaded_nodes, data_to_unload, ewghts, divisible);
                        break;
                    }
                    io_volume += unloaded_data;
                    memory_occupation -= unloaded_data;
                }

                if (!quiet)
                {
                    cerr << "occupation before processing " << memory_occupation << endl;
                }
                /*if we have enough memory to process the node, update occupation*/
                memory_occupation += node_cost - 2 * ewghts[cur_task_id] - nwghts[cur_task_id];
                memory_occupation = max(0.0, memory_occupation);

                if (!quiet)
                {
                    cerr << "processing " << cur_task_id << endl;
                    cerr << "unloading " << cur_task_id << endl;
                    cerr << "loading ";
                }

                switch (method)
                {
                case FIRST_FIT:
                    loaded_nodes.remove(make_pair(cur_task_id, schedule_vec.end() - find(schedule_vec.begin(), schedule_vec.end(), cur_task_id)));
                    break;
                case LARGEST_FIT:
                    loaded_nodes_ew.remove(make_pair(cur_task_id, ewghts[cur_task_id]));
                    break;

                default:
                    break;
                }

                for (int j = chstart[cur_task_id]; j < chstart[cur_task_id + 1]; j++)
                {
                    int ch = children[j];
                    if (!quiet)
                    {
                        cerr << ch << " ";
                    }
                    switch (method)
                    {
                    case FIRST_FIT:
                        loaded_nodes.push_back(make_pair(ch, schedule_vec.end() - find(schedule_vec.begin(), schedule_vec.end(), ch)));
                        break;
                    case LARGEST_FIT:
                        loaded_nodes_ew.push_back(make_pair(ch, ewghts[ch]));
                        break;

                    default:
                        break;
                    }
                }

                if (!quiet)
                {
                    cerr << endl;
                }
                if (!quiet)
                {
                    cerr << "New occupation after processing " << memory_occupation << endl;
                }
            }
        }
    }
    delete schedule_f;

    return io_volume;
    //    cerr<<"IO Volume "<<io_volume<<endl;
}

double IOCounterWithVariableMem(Tree *tree, int N, double *nwghts, double *ewghts, int *chstart, int *children, int *schedule, 
                                Cluster *cluster, bool divisible, int quiet, unsigned int &com_freq, vector<unsigned int> *brokenEdges, io_method_t method)

{
    double memory_occupation = ewghts[schedule[N - 1]];
    double io_volume = 0;
    vector<unsigned int> unloaded_nodes;
    list<node_sche> loaded_nodes;
    list<node_ew> loaded_nodes_ew;
    vector<int> schedule_vec(schedule, schedule + N);
    int cur_task_id;
    vector<unsigned int>::iterator unloaded;
    list<unsigned int> temp;
    list<unsigned int> queue;
    unsigned int child_start, child_end;
    unsigned int subtree_size;
    list<unsigned int>::iterator iter;
    double maxoutD, memory_required, IO_sub = 0;
    uint64_t count;
    schedule_t *schedule_f = new schedule_t();
    list<int>::iterator ite_sche;
    int rootid;
    vector<unsigned int> subtreeBrokenEdges;

    /*iterates through the given permutation (schedule)*/
    for (int rank = N - 1; rank >= 1; rank--)
    {
        cur_task_id = schedule[rank];

        // cout<<"current task id: "<<cur_task_id<<endl;
        if (cur_task_id != 0)
        { //0 means this node is on other subtrees
            /*if the node was unloaded*/
            unloaded = find(unloaded_nodes.begin(), unloaded_nodes.end(), cur_task_id);
            if (unloaded != unloaded_nodes.end())
            { //find node cur_task_id unloaded
                //cout<<", (break) "<<endl;
                brokenEdges->push_back(tree->GetNode(cur_task_id)->GetothersideID());
                ++com_freq;
                unloaded_nodes.erase(unloaded);
                temp.clear();
                queue.push_back(cur_task_id);
                //cout<<"children "<<endl;
                do
                {
                    child_start = *(chstart + queue.front());
                    child_end = *(chstart + queue.front() + 1);
                    temp.push_back(queue.front());
                    queue.pop_front();
                    for (unsigned int i = child_start; i < child_end; ++i)
                    {
                        //cout<<*(children+i)<<" ";
                        queue.push_back(*(children + i));
                    }
                } while (!queue.empty());
                //cout<<endl;

                subtree_size = temp.size(); //just an approximation
                for (long i = rank - 1; i >= 0; i--)
                {
                    iter = find(temp.begin(), temp.end(), schedule[i]);
                    if (iter != temp.end())
                    {
                        schedule[i] = 0; //IO counter will pass 0;
                        temp.erase(iter);
                    }
                    if (temp.size() == 1)
                    {
                        break;
                    }
                }

                double *ewghtssub, *timewghtssub, *spacewghtssub;
                int *prntssub;
                Tree *subtree = BuildSubtree(tree, tree->GetNode(cur_task_id));

                subtree_size = subtree->GetNodes()->size();

                int *schedule_copy = new int[subtree_size + 1];
                maxoutD = MaxOutDegree(subtree, true);
                schedule_f->clear();
                count = 0;
                MinMem(subtree, maxoutD, memory_required, *schedule_f, true, count);
                ite_sche = schedule_f->begin();
                for (unsigned int i = subtree_size; i >= 1; --i)
                {
                    schedule_copy[i] = *ite_sche;
                    advance(ite_sche, 1);
                }
                schedule_copy[0] = subtree_size + 1;
                int *chstartsub, *childrensub, *chendsub;
                po_construct(subtree_size, prntssub, &chstartsub, &chendsub, &childrensub, &rootid);

                if (memory_required > cluster->getFirstFreeProcessor()->getMemorySize())
                {
                    //   cout << "memory required " << memory_required << ", is larger than what is available " << availableMemorySizesA2[currentProcessor] << " on proc " << currentProcessor << endl;
                    //  cout << "----------------------Processing subtree! " << cur_task_id << endl;
                    // currentProcessor++;
                    //INcrease processor??
                    IO_sub = IOCounterWithVariableMem(subtree, subtree_size + 1, spacewghtssub, ewghtssub, chstartsub, childrensub, schedule_copy, cluster, divisible, quiet, com_freq, &subtreeBrokenEdges, method);

                    //    cout << "subtree broken edges " << subtreeBrokenEdges.size() << endl;

                    for (vector<unsigned int>::iterator iter = subtreeBrokenEdges.begin(); iter != subtreeBrokenEdges.end(); ++iter)
                    {
                        brokenEdges->push_back(tree->GetNode(*iter)->GetothersideID());
                    }
                    //    cout << "----------------------Out of Processing subtree!" << endl;
                }
                else
                {
                    cluster->getFirstFreeProcessor()->assignTask(tree->GetNode(cur_task_id));
                    //   cout << "just increase proc to " << currentProcessor << endl;
                }

                //   cout << "broken edges " << brokenEdges->size() << endl;

                delete[] ewghtssub;
                delete[] timewghtssub;
                delete[] spacewghtssub;
                delete[] prntssub;
                delete[] chstartsub;                
                delete[] childrensub;
                delete[] schedule_copy;
                delete subtree;

                io_volume += IO_sub;
            }
            else
            {
                double node_cost = ewghts[cur_task_id] + nwghts[cur_task_id];
                for (int j = chstart[cur_task_id]; j < chstart[cur_task_id + 1]; j++)
                {
                    node_cost += ewghts[children[j]];
                }

                double data_to_unload = memory_occupation + node_cost - ewghts[cur_task_id] - cluster->getFirstFreeProcessor()->getMemorySize();

                if (!quiet)
                {
                    cerr << "min data to unload " << data_to_unload << endl;
                }
                if (data_to_unload > 0)
                {
                    //cerr<<"We must commit I/O in order to process node "<<cur_task_id<<" which requires "<< memory_occupation + node_cost - ewghts[cur_task_id]<< " but has "<<available_memory<<"available"<<endl;
                    /*if we dont have enough room, unload files and update both io and occupation*/

                    switch (method)
                    {
                    case FIRST_FIT:
                        loaded_nodes.remove(make_pair(cur_task_id, schedule_vec.end() - find(schedule_vec.begin(), schedule_vec.end(), cur_task_id)));
                        loaded_nodes.sort(sort_sche); //descending schedule order
                        break;
                    case LARGEST_FIT:
                        loaded_nodes_ew.remove(make_pair(cur_task_id, ewghts[cur_task_id]));
                        loaded_nodes_ew.sort(sort_ew);
                        break;

                    default:
                        break;
                    }

                    double unloaded_data = 0.0;
                    switch (method)
                    {
                    case FIRST_FIT:
                        unloaded_data = unload_furthest_first_fit(tree, unloaded_nodes, loaded_nodes, data_to_unload, ewghts, divisible);
                        break;
                    case LARGEST_FIT:
                        unloaded_data = unload_largest_first_fit(tree, unloaded_nodes, loaded_nodes_ew, data_to_unload, ewghts);
                        break;
                    default:
                        unloaded_data = unload_furthest_first_fit(tree, unloaded_nodes, loaded_nodes, data_to_unload, ewghts, divisible);
                        break;
                    }
                    io_volume += unloaded_data;
                    memory_occupation -= unloaded_data;
                }

                if (!quiet)
                {
                    cerr << "occupation before processing " << memory_occupation << endl;
                }
                /*if we have enough memory to process the node, update occupation*/
                memory_occupation += node_cost - 2 * ewghts[cur_task_id] - nwghts[cur_task_id];
                memory_occupation = max(0.0, memory_occupation);

                if (!quiet)
                {
                    cerr << "processing " << cur_task_id << endl;
                    cerr << "unloading " << cur_task_id << endl;
                    cerr << "loading ";
                }

                switch (method)
                {
                case FIRST_FIT:
                    loaded_nodes.remove(make_pair(cur_task_id, schedule_vec.end() - find(schedule_vec.begin(), schedule_vec.end(), cur_task_id)));
                    break;
                case LARGEST_FIT:
                    loaded_nodes_ew.remove(make_pair(cur_task_id, ewghts[cur_task_id]));
                    break;

                default:
                    break;
                }

                for (int j = chstart[cur_task_id]; j < chstart[cur_task_id + 1]; j++)
                {
                    int ch = children[j];
                    if (!quiet)
                    {
                        cerr << ch << " ";
                    }
                    switch (method)
                    {
                    case FIRST_FIT:
                        loaded_nodes.push_back(make_pair(ch, schedule_vec.end() - find(schedule_vec.begin(), schedule_vec.end(), ch)));
                        break;
                    case LARGEST_FIT:
                        loaded_nodes_ew.push_back(make_pair(ch, ewghts[ch]));
                        break;

                    default:
                        break;
                    }
                }

                if (!quiet)
                {
                    cerr << endl;
                }
                if (!quiet)
                {
                    cerr << "New occupation after processing " << memory_occupation << endl;
                }
            }
            cluster->getFirstFreeProcessor()->assignTask(tree->GetNode(cur_task_id));
        }
    }
    delete schedule_f;

    return io_volume;
    //    cerr<<"IO Volume "<<io_volume<<endl;
}



double MaxOutDegree(Tree *tree, int quiet)
{
    double max_out = 0;
    double max_j = 0;
    for (unsigned int j = 1; j <= tree->GetNodes()->size(); j++)
    {
        ////cout<<j<<endl;
        double cur_out = tree->GetNode(j)->GetCost();
        if (cur_out >= max_out)
        {
            max_out = cur_out;
            max_j = tree->GetNode(j)->GetId();
        }
    }
    if (!quiet)
    {
        cerr << "Max out degree " << max_out << " at " << max_j << endl;
    }
    return max_out;
}

double MaxOutDegree(int N, double *nwghts, double *ewghts, int *chstart, int *children)
{
    double max_out = 0;
    for (int j = 1; j < N + 1; j++)
    {
        double cur_out = nwghts[j] + ewghts[j];

        for (int ch = chstart[j]; ch < chstart[j + 1]; ch++)
        {
            cur_out += ewghts[children[ch]];
        }

        if (cur_out >= max_out)
        {
            max_out = cur_out;
        }
    }
    ////cout<<leaves<<" leaves"<<endl;
    ////cout<<max_j<<" is the max node"<<endl;
    return max_out;
}

double MaxOutDegree(int N, int *prnts, double *nwghts, double *ewghts)
{
    int *chstart, *children, *chend;
    int root;

    po_construct(N, prnts, &chstart, &chend, &children, &root);

    double max_out = MaxOutDegree(N, nwghts, ewghts, chstart, children);

    delete[] chstart;   
    delete[] children;

    return max_out;
}

Tree *SubtreeRooted(Task *node)
{
    Tree *subtree = new Tree();

    subtree->SetRootId(1);
    subtree->SetTreeId(node->GetId());
    subtree->AddNode(node);

    vector<Task *> visit_next;
    vector<Task *>::iterator first_node;
    Task *end_node;
    if (node->IsLeaf())
    {
        return subtree;
    }
    else
    {
        visit_next = *(node->GetChildren());
        while (!visit_next.empty())
        {
            if (!visit_next.back()->IsBroken())
            { // this child has not been cut
                subtree->AddNode(visit_next.back());
                end_node = visit_next.back();
                visit_next.pop_back();
                if (!end_node->IsLeaf())
                {
                    visit_next.insert(visit_next.end(), end_node->GetChildren()->begin(), end_node->GetChildren()->end());
                }
            }
            else
            {
                visit_next.pop_back();
            }
        }
        return subtree;
    }
}

Tree *BuildSubtreeOld(Tree *tree, Task *SubtreeRoot, unsigned int new_tree_size, int **prnts, double **ewghts, double **timewghts, double **spacewghts, int *chstart, int *children)
{
    *prnts = new int[new_tree_size + 1];
    *ewghts = new double[new_tree_size + 1];
    *timewghts = new double[new_tree_size + 1];
    *spacewghts = new double[new_tree_size + 1];
    unsigned int *originalIDs = new unsigned int[new_tree_size + 1];

    unsigned int real_tree_size = 0;
    unsigned int nodeID = 1;
    (*prnts)[1] = 0;
    originalIDs[1] = SubtreeRoot->GetId();
    (*ewghts)[1] = SubtreeRoot->GetEW();
    (*timewghts)[1] = SubtreeRoot->GetMSW();
    (*spacewghts)[1] = SubtreeRoot->GetNW();

    Task *currentNode;
    list<unsigned int> que;
    unsigned int originalID = SubtreeRoot->GetId();
    que.push_back(originalID);
    unsigned int parentID;
    unsigned int tempid = SubtreeRoot->GetothersideID();
    SubtreeRoot->SetothersideID(1);

    while (!que.empty())
    {
        originalID = que.front();
        que.pop_front();
        parentID = tree->GetNode(originalID)->GetothersideID();
        for (int j = chstart[originalID]; j < chstart[originalID + 1]; j++)
        {
            currentNode = tree->GetNode(children[j]);
            if (!currentNode->IsBroken())
            { // broken edge means node is on aother subtree
                nodeID++;
                originalIDs[nodeID] = children[j];
                (*prnts)[nodeID] = parentID;
                (*ewghts)[nodeID] = currentNode->GetEW();
                (*timewghts)[nodeID] = currentNode->GetMSW();
                (*spacewghts)[nodeID] = currentNode->GetNW();
                que.push_back(children[j]);
                tree->GetNode(children[j])->SetothersideID(nodeID);
            }
        }
    }
    real_tree_size = nodeID;

    SubtreeRoot->SetothersideID(tempid);

    Tree *treeobj = new Tree(real_tree_size, *prnts, *spacewghts, *ewghts, *timewghts);

    for (unsigned int i = 1; i <= real_tree_size; i++)
    {
        treeobj->GetNode(i)->SetothersideID(originalIDs[i]); //corresponding to the original tree's id
    }

    delete[] originalIDs;

    return treeobj;
}

Tree *BuildSubtree(Tree *tree, Task *SubtreeRoot)
{
    
    
    Task *currentNode;
    list<Task*> que;
    vector<Task*> nodesOfSubtree;
    que.push_back(SubtreeRoot);
    Task * parent;
    Task *temp = SubtreeRoot;
    SubtreeRoot->SetothersideID(1);

    while (!que.empty())
    {
        Task * currentTask = que.front();
        que.pop_front();
        parent = currentTask->GetParent();
        nodesOfSubtree.push_back(currentTask);
        for(Task * child: *(currentTask->GetChildren())){
            if (!child->IsBroken())
            {                 
                que.push_back(child); 
            }
        }      
    }
    
    Tree *treeobj = new Tree(nodesOfSubtree, tree->getOriginalTree());

    //TODO: set other side ids? 
    /*  for (unsigned int i = 1; i <= real_tree_size; i++)
    {
        treeobj->GetNode(i)->SetothersideID(originalIDs[i]); //corresponding to the original tree's id
    } */

    return treeobj;
}

#endif