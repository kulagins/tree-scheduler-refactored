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
#include "heuristics.h"
#include "lib-io-tree-utils.h"
#include "lib-io-tree-free-methods.h"

//#include <omp.h>

extern double BANDWIDTH;

bool cmp_noincreasing(Task *a, Task *b) { return (a->GetMSCost(true, false) >= b->GetMSCost(true, false)); };
bool cmp_nodecreasing(Task *a, Task *b) { return (a->GetMSCost(true, false) < b->GetMSCost(true, false)); };
bool cmp_noIn_noCommu(Task *a, Task *b) { return (a->GetMSCost(false, false) >= b->GetMSCost(false, false)); };

struct CompareMapEntries
{
    int val;
    CompareMapEntries(const int &val) : val(val) {}
};

bool operator==(const std::pair<int, int> &p, const CompareMapEntries &c)
{
    return c.val == p.second;
}
bool operator==(const CompareMapEntries &c, const std::pair<int, int> &p)
{
    return c.val == p.second;
}

template <class T, class U>
void GetTwoLargestElementTypetwo(T container, U &Largest, U &secondLargest)
{
    if (container->front()->GetMSCost(true, false) > container->at(1)->GetMSCost(true, false))
    {
        Largest = container->front();
        secondLargest = container->at(1);
    }
    else
    {
        Largest = container->at(1);
        secondLargest = container->front();
    }

    if (container->size() > 2)
    {
        vector<Task *>::iterator iter = container->begin();
        iter = iter + 2;
        for (; iter != container->end(); ++iter)
        {
            if ((*iter)->GetMSCost(true, false) > Largest->GetMSCost(true, false))
            {
                secondLargest = Largest;
                Largest = *iter;
            }
            else if ((*iter)->GetMSCost(true, false) > secondLargest->GetMSCost(true, false))
            {
                secondLargest = *iter;
            }
        }
    }
}

void GetTwoLargestElementTypethree(vector<Task *> *container, vector<Task *>::iterator &Largest, vector<Task *>::iterator &secondLargest)
{
    if (container->front()->GetMSCost(false, false) >= container->back()->GetMSCost(false, false))
    {
        Largest = container->begin();
        secondLargest = Largest;
        advance(secondLargest, 1);
    }
    else
    {
        secondLargest = container->begin();
        Largest = secondLargest;
        advance(Largest, 1);
    }

    if (container->size() > 2)
    {
        vector<Task *>::iterator iter = container->begin();
        advance(iter, 2);
        for (; iter != container->end(); ++iter)
        {
            if ((*iter)->GetMSCost(false, false) > (*Largest)->GetMSCost(false, false))
            {
                secondLargest = Largest;
                Largest = iter;
            }
            else if ((*iter)->GetMSCost(false, false) > (*secondLargest)->GetMSCost(false, false))
            {
                secondLargest = iter;
            }
        }
    }
}

void GetTwoSmallestElement(list<Task *> *container, list<Task *>::iterator &Smallest, list<Task *>::iterator &secondSmallest)
{
    if (container->front()->GetMSW() <= container->back()->GetMSW())
    {
        Smallest = container->begin();
        secondSmallest = Smallest;
        advance(secondSmallest, 1);
    }
    else
    {
        secondSmallest = container->begin();
        Smallest = secondSmallest;
        advance(Smallest, 1);
    }

    if (container->size() > 2)
    {
        list<Task *>::iterator iter = container->begin();
        advance(iter, 2);
        for (; iter != container->end(); ++iter)
        {
            if ((*iter)->GetMSW() < (*Smallest)->GetMSW())
            {
                secondSmallest = Smallest;
                Smallest = iter;
            }
            else if ((*iter)->GetMSW() < (*secondSmallest)->GetMSW())
            {
                secondSmallest = iter;
            }
        }
    }
}

double SplitSubtrees(Task *root, unsigned long num_processor, double twolevel, list<Task *> &parallelRoots, unsigned long &sequentialLength)
{
    parallelRoots.clear();
    parallelRoots.emplace_front(root);
    //cout<<"   insert root"<<endl;
    vector<double> MS(1, root->GetMSCost(true, true)); // take communication cost into account
    double MS_sequential = root->GetEW() / BANDWIDTH, Weight_more, Weight_PQ;
    unsigned long amountSubtrees;
    vector<Task *> *children;

    Task *currentNode = root;
    double temp;
    unsigned int mergetime;
    while (!currentNode->IsLeaf())
    {
        MS_sequential = MS_sequential + currentNode->GetMSW();

        Weight_PQ = 0;
        parallelRoots.remove(currentNode);
        //cout<<"pop up "<<currentNode->GetId()<<endl;

        children = currentNode->GetChildren();
        for (vector<Task *>::iterator iter = children->begin(); iter != children->end(); iter++)
        {
            if ((*iter)->IsBroken())
            {
                temp = (*iter)->GetMSCost(true, false);
                if (temp > Weight_PQ)
                {
                    Weight_PQ = temp;
                }
            }
            else
            {
                parallelRoots.push_back(*iter);
                //cout<<"   insert "<<(*iter)->GetId()<<endl;
            }
        }

        if (parallelRoots.empty())
        {
            break;
        }
        else
        {
            currentNode = *max_element(parallelRoots.begin(), parallelRoots.end(), cmp_nodecreasing); //non-decreasing
        }

        temp = currentNode->GetMSCost(true, false);
        if (temp > Weight_PQ)
        {
            Weight_PQ = temp;
        }

        Weight_more = 0;
        amountSubtrees = parallelRoots.size() + 1;
        if (amountSubtrees > num_processor)
        {
            parallelRoots.sort(cmp_noIn_noCommu); //non-increasing sort, computation weight, no communication
            list<Task *>::reverse_iterator iter = parallelRoots.rbegin();
            mergetime = amountSubtrees - num_processor;
            for (unsigned int i = 0; i < mergetime; ++i, ++iter)
            {
                Weight_more += (*iter)->GetMSCost(false, false); // no comunication cost, ImprovedSplit never goes to here.
            }
        }
        //cout<<"makespan "<<MS_sequential+Weight_more+Weight_PQ<<endl;
        MS.push_back(MS_sequential + Weight_more + Weight_PQ);
    }

    if (twolevel == true)
    {
        double makespan;
        makespan = *std::min_element(MS.begin(), MS.end());
        return makespan;
    }

    //return broken edges, i.e., root of subtrees
    vector<double>::iterator smallestMS_iter = min_element(MS.begin(), MS.end());
    unsigned long minMS_step = smallestMS_iter - MS.begin();
    //cout<<"minMS_step "<<minMS_step<<endl;
    sequentialLength = minMS_step;
    unsigned int i = 0;
    parallelRoots.clear();
    parallelRoots.push_back(root);
    currentNode = root;
    while (i < minMS_step)
    {
        parallelRoots.remove(currentNode);

        children = currentNode->GetChildren();
        for (vector<Task *>::iterator iter = children->begin(); iter != children->end(); iter++)
        {
            if (!(*iter)->IsBroken())
            {
                parallelRoots.push_back(*iter);
            }
        }

        currentNode = *max_element(parallelRoots.begin(), parallelRoots.end(), cmp_nodecreasing); //non-decreasing
        i++;
    }

    if (parallelRoots.size() > 1)
    {
        amountSubtrees = parallelRoots.size() + 1;
    }
    else
    {
        amountSubtrees = 1;
    }

    if (amountSubtrees > num_processor)
    {
        parallelRoots.sort(cmp_noIn_noCommu); //non-increasing sort, computation weight, no communication cost
        mergetime = amountSubtrees - num_processor;
        for (unsigned int i = 0; i < mergetime; ++i)
        {
            parallelRoots.pop_back();
        }
    }

    root->BreakEdge(); //root should always be broken
    for (list<Task *>::iterator iter = parallelRoots.begin(); iter != parallelRoots.end(); ++iter)
    {
        (*iter)->BreakEdge();
    }
    //    cout<<"makespan from the tree root "<<root->GetMSCost(true,true)<<endl;

    return *smallestMS_iter;
}

void ISCore(Task *root, unsigned long num_processors, bool sequentialPart)
{ //number of processors here assumed to the same as tree'size
    list<Task *> parallelRoots;
    double MS_before;
    double MS_now;
    unsigned long SF_now; //avoid dead lock
    double makespan;

    if (root->IsLeaf())
    {
        //cout<<"root is leaf, return."<<endl;
        return;
    }

    makespan = SplitSubtrees(root, num_processors, false, parallelRoots, SF_now); //SF_now will be modified in SplitSubtrees, it represents the length of sequential part, 0 means the subtree no need to partition

    if (sequentialPart == true)
    {
        if (SF_now == 0)
        {
            //cout<<"this subtree has already been fully checked, return."<<endl;
            return;
        }
    }

    parallelRoots.sort(cmp_noincreasing); //non-increasing sort, communication counted

    Task *frontNode;
    if (parallelRoots.size() > 1)
    { //==1 means there is no parallel part
        while (true)
        {
            frontNode = parallelRoots.front();
            parallelRoots.pop_front();
            MS_before = frontNode->GetMSCost(true, false);

            //cout<<"---ISCore works on Parallel root "<<frontNode->GetId()<<endl;
            ISCore(frontNode, num_processors, false);

            MS_now = frontNode->GetMSCost(true, true); //Makespan updated enforced

            if (MS_now >= MS_before)
            {
                break;
            }

            if (parallelRoots.empty())
            {
                break;
            }

            if (parallelRoots.front()->GetMSCost(true, false) <= MS_now)
            {
                break;
            }
        }

        //cout<<"---ISCore works on Sequential root "<<root->GetId()<<endl;
        ISCore(root, num_processors, true);
    }

    return;
}

double ImprovedSplit(Tree *tree, unsigned int number_processor, int *chstart, int *childrenID)
{
    unsigned long tree_size = tree->GetNodes()->size();
    Task *root = tree->GetRoot();
    //cout<<"---ISCore works on the root"<<endl;
    ISCore(root, tree_size, false);

    unsigned int numberSubtrees = tree->HowmanySubtrees(true);
    double makespan = Merge(tree, numberSubtrees, number_processor, 0, chstart, childrenID, false);
    return makespan;
}

bool cmp_merge_smallest(const pair<double, Task *> &a, const pair<double, Task *> &b) { return a.first < b.first; };

bool estimateMS(Tree *tree, Tree *Qtree, Task *&smallestNode, int *chstart, int *childrenID, double memory_size, bool CheckMemory)
{
    //cout<<"   ---start compute the minimum combination"<<endl;

    Task *currentQNode;
    double increase;
    bool memoryEnough;
    Task *LargestNode;
    Task *secondLargest;
    bool leaf = false;
    const vector<Task *> *subtrees = Qtree->GetNodes();
    vector<Task *> *children;

    if (subtrees->front()->GetId() != 1)
    { //the root is supposed to be the first element in vector nodes
        cout << "error in function estimateMS" << endl;
        return false;
    }

    vector<Task *> tempQue;
    currentQNode = Qtree->GetRoot();
    currentQNode->SetMSDiff(0);
    children = currentQNode->GetChildren();
    tempQue.insert(tempQue.end(), children->begin(), children->end());
    //cout<<"   ---compute makespan difference---"<<endl;
    while (!tempQue.empty())
    {
        currentQNode = tempQue.back();
        tempQue.pop_back();
        currentQNode->SetMSDiff(currentQNode->GetParent()->GetMSDiff() + currentQNode->GetParent()->GetParallelPart() - currentQNode->GetMSCost(true, false));
        //cout<<"   subtree "<<currentQNode->GetothersideID()<<", makespan difference: "<<currentQNode->GetMSDiff()<<endl;
        children = currentQNode->GetChildren();
        tempQue.insert(tempQue.end(), children->begin(), children->end());
    }
    //cout<<"   ----------------------------------"<<endl;

    list<pair<double, Task *>> list_increase_id;
    vector<Task *>::const_iterator iter = subtrees->begin();
    ++iter;
    unsigned long size = subtrees->size() - 1;
    //  #pragma omp parallel for
    for (unsigned int step = 0; step < size; ++step)
    {
        currentQNode = *(iter + step);

        if (tree->GetNode(currentQNode->GetothersideID())->IsBroken() == true)
        { //this subtree has not been merged yet
            children = currentQNode->GetChildren();
            if (children->empty())
            { //this is a leaf node
                children = currentQNode->GetParent()->GetChildren();
                if (children->size() == 2)
                {
                    increase = children->front()->GetMSW() + children->back()->GetMSW() - currentQNode->GetParent()->GetParallelPart();
                }
                else if (children->size() == 1)
                {
                    increase = -currentQNode->GetEW() / BANDWIDTH;
                }
                else
                {

                    GetTwoLargestElementTypetwo(children, LargestNode, secondLargest); //no-increasing, communication counted

                    if (currentQNode->GetId() == LargestNode->GetId())
                    {
                        increase = currentQNode->GetMSW() + secondLargest->GetMSCost(true, false) - currentQNode->GetMSCost(true, false);
                    }
                    else
                    {
                        increase = currentQNode->GetMSW();
                    }
                }
            }
            else
            { //not a leaf node
                children = currentQNode->GetParent()->GetChildren();

                if (children->size() == 1)
                {
                    increase = -currentQNode->GetEW() / BANDWIDTH;
                }
                else
                {
                    GetTwoLargestElementTypetwo(children, LargestNode, secondLargest); //no-increasing, communication counted

                    if (currentQNode->GetId() == LargestNode->GetId())
                    {
                        increase = currentQNode->GetMSW() + max(secondLargest->GetMSCost(true, false), currentQNode->GetParallelPart()) - currentQNode->GetMSCost(true, false);
                    }
                    else
                    {
                        increase = currentQNode->GetMSW();
                    }
                }
            }

            //now consider the "slack" between the parent of current node with siblings of its parent
            increase = increase - currentQNode->GetParent()->GetMSDiff();

            //cout<<"merge, increase in MS(r) "<<increase<<endl;
            list_increase_id.push_back(pair<double, Task *>(increase, currentQNode));
        }
    }

    bool feasible = false;
    list<pair<double, Task *>>::iterator smallest_iter;
    while (feasible == false && list_increase_id.empty() == false)
    {
        smallest_iter = min_element(list_increase_id.begin(), list_increase_id.end(), cmp_merge_smallest);
        currentQNode = (*smallest_iter).second;
        //cout<<"   increase in MS(r) estimated: "<<(*smallest_iter).first<<endl;

        children = currentQNode->GetChildren();
        if (children->empty())
        {
            leaf = true;
        }

        if (CheckMemory == true)
        {
            memoryEnough = tree->MemoryEnough( currentQNode->GetParent(), currentQNode, leaf, memory_size, chstart, childrenID);
        }
        else
        {
            memoryEnough = true;
        }

        if (memoryEnough == true)
        {
            feasible = true;
            smallestNode = currentQNode;
        }
        else
        {
            list_increase_id.erase(smallest_iter);
        }
    }

    //cout<<"   ---end compute the minimum combination"<<endl;
    return feasible;
}

double Merge(Tree *tree, unsigned int num_subtrees, unsigned int processor_number, double const memory_size, int *chstart, int *childrenID, bool CheckMemory)
{
    Task *root = tree->GetRoot();

    if (processor_number >= num_subtrees)
    {
        return root->GetMSCost(true, true);
    }

    Tree *Qtreeobj = tree->BuildQtree();

    Task *node_smallest_increase;
    Task *parent;
    int shortage = num_subtrees - processor_number;
    double temp;
    Task *nodeone;
    Task *nodetwo;
    bool memoryEnough;

    while (shortage > 0)
    { //merge subtree
        //cout<<"shortage "<<shortage<<endl;
        temp = Qtreeobj->GetRoot()->GetMSCost(true, true); //initilize ms
        temp = tree->GetRoot()->GetMSCost(true, true);     //update ms

        //memoryEnough=increaseMS(tree, Qtreeobj, node_smallest_increase, chstart, childrenID, memory_size, CheckMemory);
        memoryEnough = estimateMS(tree, Qtreeobj, node_smallest_increase, chstart, childrenID, memory_size, CheckMemory);

        //when parameter checkMemory is false, memoryEnough will always be true;
        if (memoryEnough == true)
        {
            //merge currentNode (or and its sibling) to its parent
            if (node_smallest_increase->IsLeaf())
            {
                parent = node_smallest_increase->GetParent();
                if (parent->GetChildren()->size() == 2)
                {
                    nodeone = parent->GetChildren()->front();
                    nodetwo = parent->GetChildren()->back();
                    //cout<<"Merge node "<<nodeone->GetothersideID()<<" and its sibling "<<nodetwo->GetothersideID()<<endl;
                    nodeone->MergetoParent();
                    nodetwo->MergetoParent();
                    shortage = shortage - 2;
                    tree->GetNode(nodeone->GetothersideID())->RestoreEdge();
                    tree->GetNode(nodetwo->GetothersideID())->RestoreEdge();
                }
                else
                {
                    //cout<<"Merge node "<<node_smallest_increase->GetothersideID()<<endl;
                    node_smallest_increase->MergetoParent();
                    shortage--;
                    tree->GetNode(node_smallest_increase->GetothersideID())->RestoreEdge();
                }
            }
            else
            {
                //cout<<"Merge node "<<node_smallest_increase->GetothersideID()<<endl;
                node_smallest_increase->MergetoParent();
                shortage--;
                tree->GetNode(node_smallest_increase->GetothersideID())->RestoreEdge();
            }
            //cout<<"------------------------"<<endl;
        }
        else
        {
            break;
        }
    }

    if (shortage > 0)
    { //failure
        temp = -1;
    }
    else
    {
        temp = root->GetMSCost(true, true);
    }
    delete Qtreeobj;

    return temp;
}

double MergeV2(Tree *tree, unsigned int num_subtrees, unsigned int processor_number, double const memory_size, int *chstart, int *childrenID, bool CheckMemory)
{
    if (processor_number >= num_subtrees)
    {
        return tree->GetRoot()->GetMSCost(true, true);
    }

    tree->GetRoot()->GetMSCost(true, true); //update makespan

    Tree *Qtreeobj = tree->BuildQtree();

    Task *currentNode;
    Task *Qroot = Qtreeobj->GetRoot();
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

    while (shortage > 0)
    {
        DeadBreak = true;
        firstTime = true;
        Llist.clear();
        CriticalPath.clear();
        temp = Qroot->GetMSCost(true, true);           //update ms
        temp = tree->GetRoot()->GetMSCost(true, true); //update ms

        CriticalPath.push_back(1);
        largestNode = Qroot;
        Children = largestNode->GetChildren();

        while (!Children->empty())
        { //initialize critical path
            temp = largestNode->GetParallelPart();
            for (vector<Task *>::iterator iter = Children->begin(); iter != Children->end(); ++iter)
            {
                if ((*iter)->GetMSCost(true, false) == temp)
                {
                    largestNode = (*iter);
                    break;
                }
            }
            CriticalPath.push_back(largestNode->GetId());
            Children = largestNode->GetChildren();
        }

        Children = Qroot->GetChildren();
        pathlength = CriticalPath.size();

        for (unsigned int i = 1; i < pathlength; ++i)
        { //initialize vector L
            for (vector<Task *>::iterator iter = Children->begin(); iter != Children->end(); ++iter)
            {
                if ((*iter)->GetId() != CriticalPath[i])
                {
                    queue.push_back(*iter);
                }
            }
            Children = Qtreeobj->GetNode(CriticalPath[i])->GetChildren();
        }

        while (!queue.empty())
        { //initialize vector L
            currentNode = queue.back();
            queue.pop_back();
            Children = currentNode->GetChildren();
            for (vector<Task *>::iterator iter = Children->begin(); iter != Children->end(); ++iter)
            {
                Llist.push_back(*iter);
                queue.push_back(*iter);
            }
        }

    CheckOnCritical:
        if (Llist.empty())
        {
            queue.push_back(Qroot);
            while (!queue.empty())
            { //initialize vector L
                currentNode = queue.back();
                queue.pop_back();
                Children = currentNode->GetChildren();
                for (vector<Task *>::iterator iter = Children->begin(); iter != Children->end(); ++iter)
                {
                    Llist.push_back(*iter);
                    queue.push_back(*iter);
                }
            }
        }
        do
        {
            if (Llist.size() == 1)
            {
                secondSmallest = Llist.begin();
                smallest = Llist.begin();
            }
            else
            {
                GetTwoSmallestElement(&Llist, smallest, secondSmallest);
            }

            if ((*smallest)->IsLeaf())
            {
                leaf = true;
            }
            else
            {
                leaf = false;
            }

            if (CheckMemory == true)
            {
                memoryCheckPass = tree->MemoryEnough( (*smallest)->GetParent(), (*smallest), leaf, memory_size, chstart, childrenID);
            }
            else
            {
                memoryCheckPass = true;
            }

            if (memoryCheckPass == false)
            {
                if ((*secondSmallest)->IsLeaf())
                {
                    leaf = true;
                }
                else
                {
                    leaf = false;
                }
                memoryCheckPass = tree->MemoryEnough( (*secondSmallest)->GetParent(), *secondSmallest, leaf, memory_size, chstart, childrenID);
                if (memoryCheckPass == true)
                {
                    currentNode = *secondSmallest;
                    DeadBreak = false;
                }
                else
                {
                    Llist.erase(secondSmallest);
                }
            }
            else
            {
                currentNode = *smallest;
                DeadBreak = false;
            }
        } while ((memoryCheckPass == false) && (!Llist.empty()));

        if (DeadBreak == true && firstTime == true)
        {
            Llist.clear();
            firstTime = false;
            goto CheckOnCritical;
        }

        if (DeadBreak == true)
        {
            delete Qtreeobj;
            return -1; //which means failure
        }

        //merge currentNode (or and its sibling) to its parent
        if (currentNode->IsLeaf())
        {
            if (currentNode->GetParent()->GetChildren()->size() == 2)
            {
                nodeone = currentNode->GetParent()->GetChildren()->front();
                nodetwo = currentNode->GetParent()->GetChildren()->back();
                //cout<<"Merge node "<<nodeone->GetId()<<"-"<<" and its sibling "<<nodetwo->GetId()<<"-"<<endl;
                nodeone->MergetoParent();
                nodetwo->MergetoParent();
                shortage = shortage - 2;
                tree->GetNode(nodeone->GetothersideID())->RestoreEdge();
                tree->GetNode(nodetwo->GetothersideID())->RestoreEdge();
            }
            else
            {
                //cout<<"Merge node "<<currentNode->GetId()<<"-"<<endl;
                currentNode->MergetoParent();
                shortage--;
                tree->GetNode(currentNode->GetothersideID())->RestoreEdge();
            }
        }
        else
        {
            //cout<<"Merge node "<<currentNode->GetId()<<"-"<<endl;
            currentNode->MergetoParent();
            shortage--;
            tree->GetNode(currentNode->GetothersideID())->RestoreEdge();
        }
        //cout<<"------------------------"<<endl;
    }

    temp = tree->GetRoot()->GetMSCost(true, true);
    delete Qtreeobj;

    return temp;
}

bool cmp_asapc(Task *a, Task *b) { return (a->GetMSminusComu() < b->GetMSminusComu()); };

bool cmp_asap(Task *a, Task *b) { return (a->GetMSCost(false, false) < b->GetMSCost(false, false)); };

double ASAP(Tree *tree, unsigned int num_processors)
{
    list<Task *> PriorityQue;
    vector<Task *> BrokenEdges;
    unsigned long step_minimumMS = 0;
    double minimumMS = tree->GetRoot()->GetMSCost(true, true);
    //cout<<"Excuting sequentially, makespan "<<minimumMS<<endl;
    double temp;
    Task *LargestNode;
    list<Task *>::iterator node_position;

    vector<Task *> *children = tree->GetRoot()->GetChildren();
    while (children->size() == 1)
    { //avoid the linear chain
        children = children->front()->GetChildren();
    }

    PriorityQue.insert(PriorityQue.end(), children->begin(), children->end());

    while (num_processors > 1)
    { //Breaking an edge a time
        if (PriorityQue.empty())
        { //when having more processors than nodes
            break;
        }

        if (PriorityQue.size() == 1)
        {
            LargestNode = PriorityQue.front();
        }
        else
        {
            LargestNode = *max_element(PriorityQue.begin(), PriorityQue.end(), cmp_asap); //computation weight, no communication
        }

        node_position = find(PriorityQue.begin(), PriorityQue.end(), LargestNode);
        PriorityQue.erase(node_position);

        if (LargestNode->GetParent()->GetChildren()->size() > 1)
        {
            LargestNode->BreakEdge(); //break edge
            BrokenEdges.push_back(LargestNode);
            temp = tree->GetRoot()->GetMSCost(true, true);
            //cout<<"Break edge "<<LargestNode->GetId()<<", makespan now: "<<temp;
            num_processors--;
            if (temp < minimumMS)
            {
                minimumMS = temp;
                step_minimumMS = BrokenEdges.size();
                //cout<<", makespan decreased";
            }
            //cout<<endl;
        }
        //cout<<"   pop up node "<<LargestNode->GetId()<<endl;
        children = LargestNode->GetChildren();
        PriorityQue.insert(PriorityQue.end(), children->begin(), children->end());
    }

    //cout<<"resotre edge ";
    unsigned long restore_index = BrokenEdges.size();
    while (restore_index > step_minimumMS)
    {
        BrokenEdges[restore_index - 1]->RestoreEdge();
        //cout<<BrokenEdges[restore_index-1]->GetId()<<" ";
        restore_index--;
    }
    //cout<<endl;

    return minimumMS;
}

bool EstimateDecrase(int idleP, Tree *tree, vector<Task *> *criticalPath, bool *lastsubtree, Task **node_i, Task **node_j)
{
    //cout<<"   --------------estimate decrease in makespan-----------------"<<endl;
    *lastsubtree = false;
    bool MSdecreased = false;
    vector<Task *> *children;
    vector<double> decreaseSequence;
    double temp, decrease = -1;
    vector<Task *> tempQue;
    Task *lastSubtreeRoot = tree->GetNode(criticalPath->back()->GetothersideID());

    //cout<<"   Last subtree root "<<lastSubtreeRoot->GetId()<<endl;
    //nodes on the last subtree of critical path
    if (idleP > 1)
    { //has at least 2 idle processor
        tempQue.push_back(lastSubtreeRoot);
        vector<Task *>::iterator largestNode, secondLargest;
        //cout<<"   work on the last subtree "<<criticalPath->back()->GetothersideID()<<endl;
        while (!tempQue.empty())
        {
            children = tempQue.back()->GetChildren();
            tempQue.pop_back();
            if (children->size() > 1)
            {                                                                        //has at least 2 children
                GetTwoLargestElementTypethree(children, largestNode, secondLargest); //node_i is the largest, in terms of W
                temp = min((*largestNode)->GetSequentialPart() - (*secondLargest)->GetEW() / BANDWIDTH, (*secondLargest)->GetSequentialPart() - (*largestNode)->GetEW() / BANDWIDTH);
                if (temp > decrease)
                {
                    decrease = temp;
                    *node_i = *largestNode;
                    *node_j = *secondLargest;
                }

                for (vector<Task *>::iterator it = children->begin(); it != children->end(); ++it)
                {
                    tempQue.push_back(*it);
                    if (it != largestNode)
                    {
                        temp = min((*largestNode)->GetSequentialPart() - (*secondLargest)->GetEW() / BANDWIDTH, (*secondLargest)->GetSequentialPart() - (*largestNode)->GetEW() / BANDWIDTH);
                        if (temp > decrease)
                        {
                            decrease = temp;
                            *node_i = *largestNode;
                            *node_j = *it;
                        }
                    }
                }
            }
            else
            {
                for (vector<Task *>::iterator it = children->begin(); it != children->end(); ++it)
                {
                    tempQue.push_back(*it);
                }
            }
        }
    }

    if (criticalPath->size() == 1)
    { //only one subtree on the critical path
        if (decrease >= 0)
        {
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
    do
    {
        currentNode = tree->GetNode(SubtreeT->GetothersideID());
        SubtreeT = SubtreeT->GetParent();
        //cout<<"   "<<SubtreeT->GetothersideID()<<"{ "<<endl;
        subtreeRoot = tree->GetNode(SubtreeT->GetothersideID());
        MS_t = SubtreeT->GetMSCost(true, false);
        W_t = SubtreeT->GetMSW();

        do
        {
            nodeOnPath = currentNode;
            currentNode = currentNode->GetParent();
            tempQue.push_back(currentNode);
            while (!tempQue.empty())
            {
                children = tempQue.back()->GetChildren();
                tempQue.pop_back();
                for (vector<Task *>::iterator it = children->begin(); it != children->end(); ++it)
                {
                    if ((*it)->GetId() != nodeOnPath->GetId() && (!(*it)->IsBroken()))
                    {
                        //cout<<"    "<<(*it)->GetId()<<" W_i "<<(*it)->GetSequentialPart()<<", MS(t) "<<MS_t<<", W_t "<<W_t<<", MS_tj "<<(*it)->GetParallelPart()<<endl;
                        tempQue.push_back((*it));
                        temp = min((*it)->GetSequentialPart(), MS_t - W_t - (*it)->GetEW() / BANDWIDTH - (*it)->GetParallelPart());
                        if (temp > decrease_othersubtrees)
                        {
                            decrease_othersubtrees = temp;
                            output_node = (*it);
                        }
                    }
                }
            }
        } while (!currentNode->IsBroken());
        //cout<<"   }"<<endl;
    } while (subtreeRoot->GetId() != tree->GetRootId());
    //cout<<endl;

    if (decrease_othersubtrees >= 0)
    {
        MSdecreased = true;
        if (decrease_othersubtrees < decrease)
        {
            *lastsubtree = true;
        }
        else
        {
            *node_i = output_node;
        }
    }
    else
    {
        if (decrease >= 0)
        {
            *lastsubtree = true;
            MSdecreased = true;
        }
    }

    return MSdecreased;
}


double SplitAgainV2(Tree *tree, unsigned int processor_number, unsigned int num_subtrees,  std::map<int, int>  &taskToPrc, std::map<int, bool>  &isProcBusy)
{
    double MS_now;
    Task *root = tree->GetRoot();
    Tree *Qtreeobj = tree->BuildQtree();

    vector<Task *> CriticalPath; //Q nodes on Critical Path

    Task *Qroot = Qtreeobj->GetRoot();
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
    while (idleProcessors > 0)
    {

        //cout<<"******** root id "<<tree->GetRootId()<<" ********"<<endl;
        CriticalPath.clear();
        CriticalPath.push_back(Qroot);
        MS_now = root->GetMSCost(true, true); //update makespan
        Qroot->GetMSCost(true, true);         //update critical path
        largestNode = Qroot;
        Children = Qroot->GetChildren();
        //cout<<"critical path (subtres' roots){1 ";
        while (!Children->empty())
        { //initialize critical path
            temp = largestNode->GetParallelPart();
            for (vector<Task *>::iterator iter = Children->begin(); iter != Children->end(); ++iter)
            {
                if ((*iter)->GetMSCost(true, false) == temp)
                {
                    largestNode = (*iter);
                    break;
                }
            }
            //cout<<largestNode->GetothersideID()<<" ";
            CriticalPath.push_back(largestNode);
            Children = largestNode->GetChildren();
        }
        //cout<<"}"<<endl;

        //cout<<"Idle processor now: "<<idleProcessors<<endl;
        MSReduced = EstimateDecrase(idleProcessors, tree, &CriticalPath, &onLastSubtree, &node_i, &node_j);

        if (MSReduced == true)
        {
            if (onLastSubtree == false)
            {
                // cout<<"split again cut edge "<<node_i->GetId()<<endl;
                node_i->BreakEdge(); //C<-C\cup C_k
                idleProcessors--;
                taskToPrc.at(node_i->GetId()) = currentIdleProcessor;
                isProcBusy.at(currentIdleProcessor) = true;
                //   cout<<"is busy? "<< (isProcBusy.at(currentIdleProcessor)? "true": "false")<<endl;
                currentIdleProcessor++;

                node_i->SetothersideID(Qtreeobj->GetNodes()->size() + 1);
                parent = node_i->GetParent();
                while (!parent->IsBroken())
                {
                    parent = parent->GetParent();
                }
                Task *Qparent = Qtreeobj->GetNode(parent->GetothersideID());
                Task *Qchild;
                Task *newNode = new Task(parent->GetothersideID(), 0, node_i->GetEW(), node_i->GetSequentialPart());
                newNode->SetId(Qtreeobj->GetNodes()->size() + 1);
                newNode->SetParent(Qparent);
                newNode->BreakEdge();
                newNode->SetothersideID(node_i->GetId());
                Qparent->AddChild(newNode);
                Qtreeobj->addNode(newNode);
                temp = Qparent->GetMSW();
                Qparent->SetMSW(temp - newNode->GetMSW());
                //cout<<"create new Q node "<<newNode->GetId()<<", msw "<<newNode->GetMSW()<<", its parent "<<Qparent->GetId()<<", msw "<<Qparent->GetMSW()<<endl;

                newNode->GetChildren()->clear();
                if (node_i->GetParallelPart() > 0)
                {
                    //cout<<"went to here1."<<endl;
                    tempVector.push_back(node_i);
                    while (!tempVector.empty())
                    {
                        Children = tempVector.back()->GetChildren();
                        tempVector.pop_back();
                        for (vector<Task *>::iterator iter = Children->begin(); iter != Children->end(); ++iter)
                        {
                            if ((*iter)->IsBroken())
                            {
                                //cout<<"went to here2."<<endl;
                                Qchild = Qtreeobj->GetNode((*iter)->GetothersideID());
                                newNode->AddChild(Qchild);
                                Qchild->SetParent(newNode);
                                Qchild->SetParentId(newNode->GetId());
                                Qparent->RemoveChild((*iter)->GetothersideID());
                            }
                            else
                            {
                                tempVector.push_back((*iter));
                            }
                        }
                    }
                }
            }
            else
            {
                //  cout<<"split again cut edge "<<node_i->GetId()<<" and edge "<<node_j->GetId()<<endl;
                node_i->BreakEdge(); //C<-C\cup C_k
                node_j->BreakEdge(); //C<-C\cup C_k
                idleProcessors = idleProcessors - 2;

                taskToPrc.at(node_i->GetId()) = currentIdleProcessor;
                isProcBusy.at(currentIdleProcessor) = true;
                currentIdleProcessor++;
                //    cout<<"is busy? "<<  (isProcBusy.at(currentIdleProcessor)? "true": "false")<<endl;

                taskToPrc.at(node_j->GetId()) = currentIdleProcessor;
                isProcBusy.at(currentIdleProcessor) = true;
                currentIdleProcessor++;
                //      cout<<"is busy? "<<  (isProcBusy.at(currentIdleProcessor)? "true": "false")<<endl;

                node_i->SetothersideID(Qtreeobj->GetNodes()->size() + 1);
                node_j->SetothersideID(Qtreeobj->GetNodes()->size() + 2);

                Task *newNodeone = new Task(CriticalPath.back()->GetId(), 0, node_i->GetEW(), node_i->GetSequentialPart());
                newNodeone->SetId(Qtreeobj->GetNodes()->size() + 1);
                newNodeone->GetChildren()->clear();
                newNodeone->SetParent(CriticalPath.back());
                newNodeone->BreakEdge();
                newNodeone->SetothersideID(node_i->GetId());
                CriticalPath.back()->AddChild(newNodeone);
                Qtreeobj->addNode(newNodeone);
                temp = CriticalPath.back()->GetMSW();
                temp = temp - newNodeone->GetMSW();

                Task *newNodetwo = new Task(CriticalPath.back()->GetId(), 0, node_j->GetEW(), node_j->GetSequentialPart());
                newNodetwo->SetId(Qtreeobj->GetNodes()->size() + 1);
                newNodetwo->GetChildren()->clear();
                newNodetwo->SetParent(CriticalPath.back());
                newNodetwo->BreakEdge();
                newNodetwo->SetothersideID(node_j->GetId());
                CriticalPath.back()->AddChild(newNodetwo);
                Qtreeobj->addNode(newNodetwo);
                temp = temp - newNodetwo->GetMSW();
                CriticalPath.back()->SetMSW(temp);
                //cout<<"create new Q node "<<newNodetwo->GetId()<<", msw "<<newNodetwo->GetMSW()<<" and new node "<<newNodeone->GetId()<<", msw "<<newNodeone->GetMSW()<<", their parent "<<CriticalPath.back()->GetId()<<", msw "<<CriticalPath.back()->GetMSW()<<endl;
            }
        }
        else
        {
            break;
        }
    }

    delete Qtreeobj;

    MS_now = root->GetMSCost(true, true);
    return MS_now;
}

double SplitAgain(Tree* tree, unsigned int processor_number, unsigned int num_subtrees){
    double MS_now;
    Task* root=tree->GetRoot();
    Tree* Qtreeobj =tree->BuildQtree();
    
    vector<Task*> CriticalPath;//Q nodes on Critical Path
    
    Task* Qroot=Qtreeobj->GetRoot();
    Task* largestNode;
    Task* node_i;
    Task* node_j;
    Task* parent;
    double temp;
    vector<Task*>* Children;
    bool MSReduced, onLastSubtree;
    vector<Task*> tempVector;
    
    int idleProcessors=processor_number-num_subtrees;
    while (idleProcessors>0) {
        //cout<<"******** root id "<<tree->GetRootId()<<" ********"<<endl;
        CriticalPath.clear();
        CriticalPath.push_back(Qroot);
        MS_now = root->GetMSCost(true, true); //update makespan
        Qroot->GetMSCost(true, true);         //update critical path
        largestNode = Qroot;
        Children = Qroot->GetChildren();
        //cout<<"critical path (subtres' roots){1 ";
        while (!Children->empty()) {//initialize critical path
            temp=largestNode->GetParallelPart();
            for (vector<Task*>::iterator iter=Children->begin(); iter!=Children->end(); ++iter) {
                if ((*iter)->GetMSCost(true, false)==temp) {
                    largestNode=(*iter);
                    break;
                }
            }
            //cout<<largestNode->GetothersideID()<<" ";
            CriticalPath.push_back(largestNode);
            Children = largestNode->GetChildren();
        }
        //cout<<"}"<<endl;

        //cout<<"Idle processor now: "<<idleProcessors<<endl;
        MSReduced = EstimateDecrase(idleProcessors, tree, &CriticalPath, &onLastSubtree, &node_i, &node_j);

        if (MSReduced == true)
        {
            if (onLastSubtree == false)
            {
                //cout<<"cut edge "<<node_i->GetId()<<endl;
                node_i->BreakEdge(); //C<-C\cup C_k
                idleProcessors--;

                node_i->SetothersideID(Qtreeobj->GetNodes()->size() + 1);
                parent = node_i->GetParent();
                while (!parent->IsBroken())
                {
                    parent = parent->GetParent();
                }
                Task* Qparent = Qtreeobj->GetNode(parent->GetothersideID());
                Task* Qchild;
                Task* newNode = new Task(parent->GetothersideID(), 0, node_i->GetEW(), node_i->GetSequentialPart());
                newNode->SetId(Qtreeobj->GetNodes()->size()+1);
                newNode->SetParent(Qparent);
                newNode->BreakEdge();
                newNode->SetothersideID(node_i->GetId());
                Qparent->AddChild(newNode);
                Qtreeobj->addNode(newNode);
                temp = Qparent->GetMSW();
                Qparent->SetMSW(temp - newNode->GetMSW());
                //cout<<"create new Q node "<<newNode->GetId()<<", msw "<<newNode->GetMSW()<<", its parent "<<Qparent->GetId()<<", msw "<<Qparent->GetMSW()<<endl;

                newNode->GetChildren()->clear();
                if (node_i->GetParallelPart() > 0)
                {
                    //cout<<"went to here1."<<endl;
                    tempVector.push_back(node_i);
                    while (!tempVector.empty())
                    {
                        Children = tempVector.back()->GetChildren();
                        tempVector.pop_back();
                        for (vector<Task*>::iterator iter=Children->begin(); iter!=Children->end(); ++iter){
                            if ((*iter)->IsBroken()) {
                                //cout<<"went to here2."<<endl;
                                Qchild = Qtreeobj->GetNode((*iter)->GetothersideID());
                                newNode->AddChild(Qchild);
                                Qchild->SetParent(newNode);
                                Qchild->SetParentId(newNode->GetId());
                                Qparent->RemoveChild((*iter)->GetothersideID());
                            }
                            else
                            {
                                tempVector.push_back((*iter));
                            }
                        }
                    }
                }
            }
            else
            {
                //cout<<"cut edge "<<node_i->GetId()<<" and edge "<<node_j->GetId()<<endl;
                node_i->BreakEdge();//C<-C\cup C_k
                node_j->BreakEdge();//C<-C\cup C_k
                idleProcessors=idleProcessors-2;
                
                node_i->SetothersideID(Qtreeobj->GetNodes()->size()+1);
                node_j->SetothersideID(Qtreeobj->GetNodes()->size()+2);
                
                Task* newNodeone = new Task(CriticalPath.back()->GetId(), 0, node_i->GetEW(), node_i->GetSequentialPart());
                newNodeone->SetId(Qtreeobj->GetNodes()->size()+1);
                newNodeone->GetChildren()->clear();
                newNodeone->SetParent(CriticalPath.back());
                newNodeone->BreakEdge();
                newNodeone->SetothersideID(node_i->GetId());
                CriticalPath.back()->AddChild(newNodeone);
                Qtreeobj->addNode(newNodeone);
                temp=CriticalPath.back()->GetMSW();
                temp=temp-newNodeone->GetMSW();
                
                Task* newNodetwo = new Task(CriticalPath.back()->GetId(), 0, node_j->GetEW(), node_j->GetSequentialPart());
                newNodetwo->SetId(Qtreeobj->GetNodes()->size()+1);
                newNodetwo->GetChildren()->clear();
                newNodetwo->SetParent(CriticalPath.back());
                newNodetwo->BreakEdge();
                newNodetwo->SetothersideID(node_j->GetId());
                CriticalPath.back()->AddChild(newNodetwo);
                Qtreeobj->addNode(newNodetwo);
                temp = temp - newNodetwo->GetMSW();
                CriticalPath.back()->SetMSW(temp);
                //cout<<"create new Q node "<<newNodetwo->GetId()<<", msw "<<newNodetwo->GetMSW()<<" and new node "<<newNodeone->GetId()<<", msw "<<newNodeone->GetMSW()<<", their parent "<<CriticalPath.back()->GetId()<<", msw "<<CriticalPath.back()->GetMSW()<<endl;
            }
        }
        else
        {
            break;
        }
    }

    delete Qtreeobj;

    MS_now = root->GetMSCost(true, true);
    return MS_now;
}


void Immediately(Tree *tree, unsigned long N, double *nwghts, double *ewghts, int *chstart, int *children, int *schedule, double m_availble, unsigned int &num_para_subtrees, vector<unsigned int> *brokenEdges)
{
    double memory_occupied = ewghts[schedule[N - 1]];
    list<unsigned int> allNodes;
    list<unsigned int> queue;
    unsigned long subtree_size;
    list<unsigned int>::iterator iter;
    unsigned int com_freq;
    schedule_t *schedule_f = new schedule_t();
    list<int>::iterator ite_sche;
    double maxoutD, memory_required, node_cost, data_to_unload;
    uint64_t count;
    int rootid, cur_task_id;
    unsigned int child_start, child_end;
    vector<unsigned int> subtreeBrokenEdges;

    //cout<<"current task:";
    for (unsigned long rank = N - 1; rank >= 1; rank--)
    {
        cur_task_id = schedule[rank];
        if (cur_task_id != 0)
        { //=0 means this node has already been moved to another processor
            //cout<<" "<<cur_task_id;
            node_cost = ewghts[cur_task_id] + nwghts[cur_task_id];
            for (int j = chstart[cur_task_id]; j < chstart[cur_task_id + 1]; j++)
            {
                node_cost += ewghts[children[j]];
            }

            data_to_unload = memory_occupied + node_cost - ewghts[cur_task_id] - m_availble;
            if (data_to_unload > 0)
            { // schedule the subtree that is rooted at this node onto another processor
                //cout<<"(break)"<<endl;
                tree->GetNode(cur_task_id)->BreakEdge(); // set it cut, used for building a quotient tree later
                brokenEdges->push_back(tree->GetNode(cur_task_id)->GetothersideID());
                ++num_para_subtrees;
                allNodes.clear();
                queue.clear();
                queue.push_back(cur_task_id);

                do
                {
                    child_start = *(chstart + queue.front());
                    child_end = *(chstart + queue.front() + 1);
                    allNodes.push_back(queue.front());
                    queue.pop_front();
                    for (unsigned int i = child_start; i < child_end; ++i)
                    {
                        queue.push_back(*(children + i));
                    }
                } while (!queue.empty());

                subtree_size = allNodes.size();
                for (long i = rank - 1; i >= 0; i--)
                {
                    iter = find(allNodes.begin(), allNodes.end(), schedule[i]);
                    if (iter != allNodes.end())
                    {
                        schedule[i] = 0; //IO counter will pass 0;
                        allNodes.erase(iter);
                    }
                    if (allNodes.size() == 1)
                    {
                        break;
                    }
                }

                double *ewghtssub, *timewghtssub, *spacewghtssub;
                int *prntssub, *chstartsub, *childrensub;
                Tree *subtree = BuildSubtree(tree, tree->GetNode(cur_task_id), subtree_size, &prntssub, &ewghtssub, &timewghtssub, &spacewghtssub, chstart, children);

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
                po_construct(subtree_size, prntssub, &chstartsub, &childrensub, &rootid);

                if (memory_required > m_availble)
                {
                    Immediately(subtree, subtree_size, spacewghtssub, ewghtssub, chstartsub, childrensub, schedule_copy, m_availble, com_freq, &subtreeBrokenEdges);

                    for (vector<unsigned int>::iterator iter = subtreeBrokenEdges.begin(); iter != subtreeBrokenEdges.end(); ++iter)
                    {
                        brokenEdges->push_back(tree->GetNode(*iter)->GetothersideID());
                    }
                }

                delete[] ewghtssub;
                delete[] timewghtssub;
                delete[] spacewghtssub;
                delete[] prntssub;
                delete[] chstartsub;               
                delete[] childrensub;
                delete[] schedule_copy;
                delete subtree;

                memory_occupied -= ewghts[cur_task_id];
                memory_occupied = max(0.0, memory_occupied);
            }
            else
            { //memory is enough for executing this node
                memory_occupied += node_cost - 2 * ewghts[cur_task_id] - nwghts[cur_task_id];
                memory_occupied = max(0.0, memory_occupied);
            }
        }
    }
    delete schedule_f;
    //cout<<endl;
}

void MemoryCheck(Tree *tree, int *chstart, int *children, Cluster *cluster,  io_method_t method)
{ //chstart, children are not modified
    vector<Task *> subtreeRoots;
    Task *currentnode;
    Task *subtreeRoot;
    int rootid;
    tree->GetRoot()->BreakEdge();

    //cout<<"Subtrees' roots: ";
    unsigned long treeSize = tree->GetNodes()->size();
    for (unsigned int i = treeSize; i >= 1; --i)
    {
        currentnode = tree->GetNode(i);
        if (currentnode->IsBroken())
        {
            //cout<<i<<" ";
            subtreeRoots.push_back(currentnode);
        }
    }
    //cout<<endl;

    double maxoutD, memory_required;
    schedule_t *schedule_f = new schedule_t();
    uint64_t count;
    unsigned int com_freq;
    unsigned long subtreeSize;
    list<int>::iterator ite_sche;
    vector<unsigned int> BrokenEdgesID;
    double IO_volume;
    while (!subtreeRoots.empty())
    {
        subtreeRoot = subtreeRoots.back();
        subtreeRoots.pop_back();

        double *ewghts, *timewghts, *spacewghts;
        int *prnts;
        Tree *subtree = BuildSubtree(tree, subtreeRoot, treeSize, &prnts, &ewghts, &timewghts, &spacewghts, chstart, children);

        subtreeSize = subtree->GetNodes()->size();
        int *schedule_copy = new int[subtreeSize + 1];
        maxoutD = MaxOutDegree(subtree, true);
        schedule_f->clear();
        count = 0;
        MinMem(subtree, maxoutD, memory_required, *schedule_f, true, count);

        int *chstartsub, *childrensub;
        po_construct(subtreeSize, prnts, &chstartsub, &childrensub, &rootid);

        //cout<<"Subtree "<<subtreeRoot->GetId()<<" needs memory "<<memory_required;
        if (memory_required > cluster->getProcessors().at(0)->getMemorySize())
        {
            //cout<<", larger than what is available: "<<memory_size<<endl;

            ite_sche = schedule_f->begin();
            for (unsigned int i = subtreeSize; i >= 1; --i)
            {
                schedule_copy[i] = *ite_sche;
                advance(ite_sche, 1);
            }
            schedule_copy[0] = subtreeSize + 1;
            double memory_size = cluster->getProcessors().at(0)->getMemorySize();
            switch (method)
            {
            case FIRST_FIT:
                IO_volume = IOCounter(subtree, subtreeSize + 1, spacewghts, ewghts, chstartsub, childrensub, schedule_copy, memory_size, false, true, com_freq, &BrokenEdgesID, FIRST_FIT);
                break;
            case LARGEST_FIT:
                IO_volume = IOCounter(subtree, subtreeSize + 1, spacewghts, ewghts, chstartsub, childrensub, schedule_copy, memory_size, false, true, com_freq, &BrokenEdgesID, LARGEST_FIT);
                break;
            case IMMEDIATELY:
                Immediately(subtree, subtreeSize + 1, spacewghts, ewghts, chstartsub, childrensub, schedule_copy, memory_size, com_freq, &BrokenEdgesID);
                break;

            default:
                break;
            }
        }
        //cout<<endl;

        delete[] ewghts;
        delete[] timewghts;
        delete[] spacewghts;
        delete[] prnts;
        delete[] schedule_copy;
        delete[] chstartsub;
        delete[] childrensub;
        delete subtree;
    }
    delete schedule_f;

    for (vector<unsigned int>::iterator iter = BrokenEdgesID.begin(); iter != BrokenEdgesID.end(); ++iter)
    {
        tree->GetNode(*iter)->BreakEdge();
    }
}
std::map<int, int> MemoryCheckA2(Tree *tree, Cluster *cluster, io_method_t method, bool skipBig)

{ //chstart, children are not modified
    vector<Task *> subtreeRoots;
    vector<Task *> subtreeRootsSkipped;
    Task *currentnode;
    Task *subtreeRoot;
    int rootid;
    tree->GetRoot()->BreakEdge();

    //cout<<"Subtrees' roots: ";
    unsigned long treeSize = tree->GetNodes()->size();
    for (unsigned int i = treeSize; i >= 1; --i)
    {
        currentnode = tree->GetNode(i);
        if (currentnode->IsBroken())
        {
            //cout<<i<<" ";
            subtreeRoots.push_back(currentnode);
            //cout << "root " << currentnode->GetMSCost() << endl;
        }
    }
    //cout<<endl;
    sort(subtreeRoots.begin(), subtreeRoots.end(), [](Task *lhs, Task *rhs)
         { return lhs->GetMSCost() < rhs->GetMSCost(); });

    double maxoutD, memory_required;
    schedule_t *schedule_f = new schedule_t();
    uint64_t count;
    unsigned int com_freq;
    unsigned long subtreeSize;
    list<int>::iterator ite_sche;
    vector<unsigned int> BrokenEdgesID;
    double IO_volume;
    Processor *currentProcessor = cluster->getFirstFreeProcessor();

    while (!subtreeRoots.empty())
    {
        double currentMem = currentProcessor->getMemorySize();
        subtreeRoot = subtreeRoots.back();
        subtreeRoots.pop_back();

        double *ewghts, *timewghts, *spacewghts;
        int *prnts;
        Tree *subtree = BuildSubtree(tree, subtreeRoot, treeSize, &prnts, &ewghts, &timewghts, &spacewghts, chstart, children);

        subtreeSize = subtree->GetNodes()->size();
        int *schedule_copy = new int[subtreeSize + 1];
        maxoutD = MaxOutDegree(subtree, true);
        schedule_f->clear();
        count = 0;
        MinMem(subtree, maxoutD, memory_required, *schedule_f, true, count);

        int *chstartsub, *childrensub;
        po_construct(subtreeSize, prnts, &chstartsub, &childrensub, &rootid);

        //  cout << "Subtree " << subtreeRoot->GetId() << " needs memory " << memory_required;
        if (memory_required > currentMem)
        {
            if (!skipBig)
            {
                //   cout << ", larger than what is available: " << currentMem << " on proc " << currentProcessor << endl;

                ite_sche = schedule_f->begin();
                for (unsigned int i = subtreeSize; i >= 1; --i)
                {
                    schedule_copy[i] = *ite_sche;
                    advance(ite_sche, 1);
                }

                schedule_copy[0] = subtreeSize + 1;
    
                switch (method)
                {
                case FIRST_FIT:
                    IO_volume = IOCounterWithVariableMem(subtree, subtreeSize + 1, spacewghts, ewghts, chstartsub, childrensub, schedule_copy, cluster, false, true, com_freq, &BrokenEdgesID, FIRST_FIT);
                    break;
                case LARGEST_FIT:
                    IO_volume = IOCounterWithVariableMem(subtree, subtreeSize + 1, spacewghts, ewghts, chstartsub, childrensub, schedule_copy, cluster, false, true, com_freq, &BrokenEdgesID, LARGEST_FIT);
                    break;
                case IMMEDIATELY:
                    Immediately(subtree, subtreeSize + 1, spacewghts, ewghts, chstartsub, childrensub, schedule_copy, currentMem, com_freq, &BrokenEdgesID);
                    break;

                default:
                    break;
                }
            }
            else
            {
                subtreeRootsSkipped.push_back(subtreeRoot);
            }
        }
        else
        {
            currentProcessor->assignTask(subtreeRoot);
            currentProcessor = cluster->getFirstFreeProcessor();
        }
        //cout<<endl;
        //
        //   cout << "Now big trees" << endl;
        while (!subtreeRootsSkipped.empty())
        {
            double currentMem = currentProcessor->getMemorySize();
            subtreeRoot = subtreeRootsSkipped.back();
            subtreeRootsSkipped.pop_back();
            //     cout << ", larger than what is available: " << currentMem << " on proc " << currentProcessor << endl;

            ite_sche = schedule_f->begin();
            for (unsigned int i = subtreeSize; i >= 1; --i)
            {
                schedule_copy[i] = *ite_sche;
                advance(ite_sche, 1);
            }

            schedule_copy[0] = subtreeSize + 1;

            switch (method)
            {
            case FIRST_FIT:
                IO_volume = IOCounterWithVariableMem(subtree, subtreeSize + 1, spacewghts, ewghts, chstartsub, childrensub, schedule_copy, cluster, false, true, com_freq, &BrokenEdgesID, FIRST_FIT);
                break;
            case LARGEST_FIT:
                IO_volume = IOCounterWithVariableMem(subtree, subtreeSize + 1, spacewghts, ewghts, chstartsub, childrensub, schedule_copy, cluster, false, true, com_freq, &BrokenEdgesID, LARGEST_FIT);
                break;
            case IMMEDIATELY:
                Immediately(subtree, subtreeSize + 1, spacewghts, ewghts, chstartsub, childrensub, schedule_copy, currentMem, com_freq, &BrokenEdgesID);
                break;

            default:
                break;
            }
        }

        delete[] ewghts;
        delete[] timewghts;
        delete[] spacewghts;
        delete[] prnts;
        delete[] schedule_copy;
        delete[] chstartsub;
        delete[] childrensub;
        delete subtree;
    }
    delete schedule_f;

    for (vector<unsigned int>::iterator iter = BrokenEdgesID.begin(); iter != BrokenEdgesID.end(); ++iter)
    {
        tree->GetNode(*iter)->BreakEdge();
    }
}

void SetBandwidth(double CCR, unsigned long tree_size, double *ewghts, double *timewghts)
{
    double sum_edges = 0;
    double sum_weights = 0;
    for (unsigned int i = 1; i <= tree_size; ++i)
    {
        sum_edges = sum_edges + ewghts[i];
        sum_weights = sum_weights + timewghts[i];
    }
    BANDWIDTH = sum_edges / (sum_weights * CCR);
}
