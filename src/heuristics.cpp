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
//#include <omp.h>

extern double BANDWIDTH;

bool cmp_noincreasing(Cnode *a, Cnode *b) { return (a->GetMSCost(true, false) >= b->GetMSCost(true, false)); };
bool cmp_nodecreasing(Cnode *a, Cnode *b) { return (a->GetMSCost(true, false) < b->GetMSCost(true, false)); };
//bool cmp_noIn_minusCommu (Cnode* a, Cnode* b){return (a->GetMSminusComu()>=b->GetMSminusComu());};
bool cmp_noIn_noCommu(Cnode *a, Cnode *b) { return (a->GetMSCost(false, false) >= b->GetMSCost(false, false)); };
bool cmp_c_nodecrea(Cnode *a, Cnode *b) { return (a->GetMSW() < b->GetMSW()); };

//void GetTwoLargestElementTypeone(vector<Cnode*>* container, Cnode* & Largest, Cnode* & secondLargest) {
//    if (container->front()->GetMSminusComu()>=container->at(1)->GetMSminusComu()) {
//        Largest=container->front();
//        secondLargest=container->at(1);
//    }else{
//        Largest=container->at(1);
//        secondLargest=container->front();
//    }
//
//    if (container->size()>2) {
//        vector<Cnode*>::iterator iter=container->begin();
//        iter=iter+2;
//        for (;iter!=container->end(); ++iter) {
//            if ((*iter)->GetMSminusComu()>Largest->GetMSminusComu()) {
//                secondLargest=Largest;
//                Largest=*iter;
//            } else if((*iter)->GetMSminusComu()>secondLargest->GetMSminusComu()){
//                secondLargest=*iter;
//            }
//        }
//    }
//}

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
        vector<Cnode *>::iterator iter = container->begin();
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

void GetTwoLargestElementTypethree(vector<Cnode *> *container, vector<Cnode *>::iterator &Largest, vector<Cnode *>::iterator &secondLargest)
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
        vector<Cnode *>::iterator iter = container->begin();
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

void GetTwoSmallestElement(list<Cnode *> *container, list<Cnode *>::iterator &Smallest, list<Cnode *>::iterator &secondSmallest)
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
        list<Cnode *>::iterator iter = container->begin();
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

double SplitSubtrees(Cnode *root, unsigned long num_processor, double twolevel, list<Cnode *> &parallelRoots, unsigned long &sequentialLength)
{
    parallelRoots.clear();
    parallelRoots.emplace_front(root);
    //cout<<"   insert root"<<endl;
    vector<double> MS(1, root->GetMSCost(true, true)); // take communication cost into account
    double MS_sequential = root->GetEW() / BANDWIDTH, Weight_more, Weight_PQ;
    unsigned long amountSubtrees;
    vector<Cnode *> *children;

    Cnode *currentNode = root;
    double temp;
    unsigned int mergetime;
    //list<Cnode*>::iterator iter;
    //unsigned int target;
    //unsigned int round=0;
    while (!currentNode->IsLeaf())
    {
        MS_sequential = MS_sequential + currentNode->GetMSW();

        Weight_PQ = 0;
        parallelRoots.remove(currentNode);
        //cout<<"pop up "<<currentNode->GetId()<<endl;

        children = currentNode->GetChildren();
        for (vector<Cnode *>::iterator iter = children->begin(); iter != children->end(); iter++)
        {
            if ((*iter)->IsBorken())
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
            list<Cnode *>::reverse_iterator iter = parallelRoots.rbegin();
            mergetime = amountSubtrees - num_processor;
            for (unsigned int i = 0; i < mergetime; ++i, ++iter)
            {
                Weight_more += (*iter)->GetMSCost(false, false); // no comunication cost, ImprovedSplit never goes to here.
            }
        }

        //round++;
        //cout<<"round "<<round<<"---MS Now "<<MS_sequential+Weight_more+Weight_PQ<<", edges broken{ ";
        //        iter=parallelRoots.begin();
        //        if (parallelRoots.size()>(num_processor-1)) {
        //            target = num_processor-1;
        //        }else{
        //            target = parallelRoots.size();
        //        }
        //        if (!parallelRoots.empty()) {
        //            for (unsigned int count=0; count<target; ++iter, ++count) {
        //                cout<<(*iter)->GetId()<<" ";
        //            }
        //        }
        //        cout<<"}"<<endl;

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
        for (vector<Cnode *>::iterator iter = children->begin(); iter != children->end(); iter++)
        {
            if (!(*iter)->IsBorken())
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
    for (list<Cnode *>::iterator iter = parallelRoots.begin(); iter != parallelRoots.end(); ++iter)
    {
        (*iter)->BreakEdge();
    }

    //    cout<<"   broken edges: ";
    //    for (list<Cnode*>::iterator iter=parallelRoots.begin(); iter!=parallelRoots.end(); ++iter) {
    //        cout<<(*iter)->GetId()<<" ";
    //        if (!(*iter)->IsBorken()) {
    //            cout<<"(error) ";
    //        }
    //    }
    //    cout<<endl;

    //    cout<<"makespan from the tree root "<<root->GetMSCost(true,true)<<endl;

    return *smallestMS_iter;
}

void ISCore(Cnode *root, unsigned long num_processors, bool sequentialPart)
{ //number of processors here assumed to the same as tree'size
    list<Cnode *> parallelRoots;
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

    Cnode *frontNode;
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

double ImprovedSplit(Ctree *tree, unsigned int number_processor, int *chstart, int *childrenID)
{
    //double ImprovedSplit(Ctree* tree, unsigned int number_processor){
    unsigned long tree_size = tree->GetNodes()->size();
    Cnode *root = tree->GetRoot();
    //cout<<"---ISCore works on the root"<<endl;
    ISCore(root, tree_size, false);

    //    Ctree* Qtreeobj = BuildQtree(tree);
    //    long index=Qtreeobj->GetNodes()->size()-number_processor;
    //    if (index>0) {
    //        list<Cnode*> C;
    //        C.assign(Qtreeobj->GetNodes()->begin(),Qtreeobj->GetNodes()->end());
    //        C.pop_front();//pop up the root
    //        C.sort(cmp_c_nodecrea);
    //        list<Cnode*>::iterator iter=C.begin();
    //        while (index>0) {
    //            tree->GetNode((*iter)->GetothersideID())->RestoreEdge();
    //            advance(iter, 1);
    //            index--;
    //        }
    //    }

    unsigned int numberSubtrees = HowmanySubtrees(tree, true);
    double makespan = Merge(tree, numberSubtrees, number_processor, 0, chstart, childrenID, false);

    //double makespan = root->GetMSCost(true,true);
    //delete Qtreeobj;
    return makespan;
}

//double ImprovedSplit(Ctree* tree, unsigned int processor_number){//first implementation
//    unsigned long tree_size=tree->GetNodes()->size();
//    Cnode* root=tree->GetRoot();
//    ISCore(root, tree_size, false, 0);
//
////    cout<<"Broken Edges: ";
////    for (unsigned int i=1; i<=tree_size; ++i) {
////        if (tree->GetNode(i)->IsBorken()) {
////            cout<<i<<" ";
////        }
////    }
////    cout<<endl;
//
//    unsigned int num_subtrees=0;
//    num_subtrees=HowmanySubtrees(tree,true);
//
//    if (processor_number>=num_subtrees) {
//        return root->GetMSCost(true, true);
//    }
//
//    Ctree* Qtreeobj = BuildQtree(tree);//makespan will also be updated in BuildQtree //root->GetMSCost(true, true);
//    Cnode* currentNode;
//    double temp;
//    unsigned int shortage=num_subtrees-processor_number;
//    unsigned int subtree_root_id;
//    double smallestParaPart;
//    vector<Cnode*>* Children;
//    while (shortage>0) {
//        currentNode=Qtreeobj->GetRoot();
//        smallestParaPart=currentNode->GetMSCost(true, true);
//        while (!currentNode->IsLeaf()) {
//            Children=currentNode->GetChildren();
//            for (vector<Cnode*>::iterator iter=Children->begin(); iter!=Children->end(); iter++) {
//                temp=(*iter)->GetMSCost(true, false);
//                if (temp<smallestParaPart) {
//                    smallestParaPart=temp;
//                    currentNode=(*iter);
//                }
//            }
//        }//now current node is the leaf node and is the smallest one from the smallest one iteratively
//
//        subtree_root_id = currentNode->GetothersideID();
//        tree->GetNode(subtree_root_id)->RestoreEdge();
//        currentNode->MergetoParent();
//        shortage--;
//
////        for (unsigned int i=1; i<=Qtreeobj->GetNodes()->size(); ++i) {
////            cout<<i<<" "<<Qtreeobj->GetNode(i)->GetParentId()<<" "<<Qtreeobj->GetNode(i)->GetMSW()<<" "<<Qtreeobj->GetNode(i)->GetEW()<<endl;
////        }
////        cout<<"-----------------------------"<<endl;
////        cout<<"Restore edge "<<subtree_root_id<<", shortage "<<shortage<<endl;
//    }
//
//    double Makespan=root->GetMSCost(true, true);
//
//    delete Qtreeobj;
//
//    return Makespan;
//}

bool MemoryEnough(Ctree *tree, Cnode *Qrootone, Cnode *Qroottwo, bool leaf, double memory_size, int *chstart, int *children)
{
    bool enough = false;
    unsigned long new_tree_size = tree->GetNodes()->size();

    Cnode *SubtreeRoot = tree->GetNode(Qrootone->GetothersideID());

    vector<Cnode *> *childrenvector = Qrootone->GetChildren();
    if ((leaf == true) & (childrenvector->size() == 2))
    {
        tree->GetNode(childrenvector->front()->GetothersideID())->RestoreEdge();
        tree->GetNode(childrenvector->back()->GetothersideID())->RestoreEdge();
    }
    else
    {
        tree->GetNode(Qroottwo->GetothersideID())->RestoreEdge(); //restore edge temporarilly
    }

    //clock_t time;
    //time=clock();

    double *ewghts, *timewghts, *spacewghts;
    int *prnts;
    Ctree *subtree = BuildSubtree(tree, SubtreeRoot, new_tree_size, &prnts, &ewghts, &timewghts, &spacewghts, chstart, children);
    delete[] ewghts;
    delete[] timewghts;
    delete[] spacewghts;
    delete[] prnts;

    //time=clock()-time;
    //time=time/10000;
    //printf("  building subtree took me %d *10^4 clicks. \n",time);

    //time=clock();
    double maxout, requiredMemory;
    uint64_t count = 0;
    schedule_t *schedule_f = new schedule_t();
    maxout = MaxOutDegree(subtree, true);
    MinMem(subtree, maxout, requiredMemory, *schedule_f, true, count);

    //time=clock()-time;
    //time=time/10000;
    //printf("  memory checking took me %d *10^4 clicks. \n",time);

    if (requiredMemory <= memory_size)
    {
        enough = true;
    }

    if ((leaf == true) & (childrenvector->size() == 2))
    {
        tree->GetNode(childrenvector->front()->GetothersideID())->BreakEdge();
        tree->GetNode(childrenvector->back()->GetothersideID())->BreakEdge();
    }
    else
    {
        tree->GetNode(Qroottwo->GetothersideID())->BreakEdge();
    }

    delete subtree;
    delete schedule_f;

    return enough;
}

///Qtree corresponds to a whole original tree
Ctree *BuildQtree(Ctree *tree)
{ //Qtree is for makespan side, so do not use it for space side
    Cnode *root = tree->GetRoot();
    root->BreakEdge();
    tree->GetRoot()->GetMSCost(true, true); //update
    size_t tree_size = tree->GetNodes()->size();
    unsigned long num_subtrees = HowmanySubtrees(tree, true);

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

    Cnode *currentNode;
    for (unsigned int i = 2; i <= tree_size; ++i)
    {
        currentNode = tree->GetNode(i);
        if (currentNode->IsBorken())
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
        currentNode = tree->GetNode(brokenEdges[i])->GetParent();
        while (!currentNode->IsBorken())
        {
            currentNode = currentNode->GetParent();
        }
        prnts[i] = currentNode->GetothersideID();
    }

    Ctree *Qtreeobj = new Ctree(num_subtrees, prnts, timewghts, ewghts, timewghts); //Qtree only reprents makespan, not memory consumption

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

bool increaseMS(Ctree *tree, Ctree *Qtree, Cnode *&smallestNode, int *chstart, int *childrenID, double memory_size, bool CheckMemory)
{
    //cout<<"   ---start compute the minimum combination"<<endl;
    //vector<Cnode*> que;
    //vector<Cnode*>* children=Qrooot->GetChildren();
    //que.insert(que.end(), children->begin(),children->end());

    Cnode *currentNode;
    double diff, increase, temp;
    bool memoryEnough;
    bool feasible = false;
    Cnode *LargestNode;
    Cnode *secondLargest;
    Cnode *parent;
    double smallestIncrease = tree->GetRoot()->GetMSCost(true, false);
    bool leaf = false;
    const vector<Cnode *> *subtrees = Qtree->GetNodes();
    vector<Cnode *> *children;

    if (subtrees->front()->GetId() != 1)
    {
        cout << "error in function increaseMs" << endl;
        return false;
    }

    vector<Cnode *>::const_iterator iter = subtrees->begin();
    ++iter;
    for (; iter != subtrees->end(); ++iter)
    {
        currentNode = (*iter);

        if (tree->GetNode(currentNode->GetothersideID())->IsBorken() == true)
        { //this subtree has not been merged yet
            children = currentNode->GetChildren();

            if (children->empty())
            {
                leaf = true;
            }

            //check the memory cost, if merge itself to its parent subtree
            if (CheckMemory == true)
            {
                memoryEnough = MemoryEnough(tree, currentNode->GetParent(), currentNode, leaf, memory_size, chstart, childrenID);
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

bool cmp_merge_smallest(const pair<double, Cnode *> &a, const pair<double, Cnode *> &b) { return a.first < b.first; };

bool estimateMS(Ctree *tree, Ctree *Qtree, Cnode *&smallestNode, int *chstart, int *childrenID, double memory_size, bool CheckMemory)
{
    //cout<<"   ---start compute the minimum combination"<<endl;

    Cnode *currentQNode;
    double increase;
    bool memoryEnough;
    Cnode *LargestNode;
    Cnode *secondLargest;
    bool leaf = false;
    const vector<Cnode *> *subtrees = Qtree->GetNodes();
    vector<Cnode *> *children;

    if (subtrees->front()->GetId() != 1)
    { //the root is supposed to be the first element in vector nodes
        cout << "error in function estimateMS" << endl;
        return false;
    }

    vector<Cnode *> tempQue;
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

    list<pair<double, Cnode *>> list_increase_id;
    vector<Cnode *>::const_iterator iter = subtrees->begin();
    ++iter;
    unsigned long size = subtrees->size() - 1;
    //  #pragma omp parallel for
    for (unsigned int step = 0; step < size; ++step)
    {
        currentQNode = *(iter + step);

        if (tree->GetNode(currentQNode->GetothersideID())->IsBorken() == true)
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
            list_increase_id.push_back(pair<double, Cnode *>(increase, currentQNode));
        }
    }

    bool feasible = false;
    list<pair<double, Cnode *>>::iterator smallest_iter;
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
            memoryEnough = MemoryEnough(tree, currentQNode->GetParent(), currentQNode, leaf, memory_size, chstart, childrenID);
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

double Merge(Ctree *tree, unsigned int num_subtrees, unsigned int processor_number, double const memory_size, int *chstart, int *childrenID, bool CheckMemory)
{
    Cnode *root = tree->GetRoot();

    if (processor_number >= num_subtrees)
    {
        return root->GetMSCost(true, true);
    }

    Ctree *Qtreeobj = BuildQtree(tree);

    Cnode *node_smallest_increase;
    Cnode *parent;
    int shortage = num_subtrees - processor_number;
    double temp;
    Cnode *nodeone;
    Cnode *nodetwo;
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

double MergeV2(Ctree *tree, unsigned int num_subtrees, unsigned int processor_number, double const memory_size, int *chstart, int *childrenID, bool CheckMemory)
{
    if (processor_number >= num_subtrees)
    {
        return tree->GetRoot()->GetMSCost(true, true);
    }

    tree->GetRoot()->GetMSCost(true, true); //update makespan

    Ctree *Qtreeobj = BuildQtree(tree);

    Cnode *currentNode;
    Cnode *Qroot = Qtreeobj->GetRoot();
    int shortage = num_subtrees - processor_number;
    list<Cnode *> Llist;
    vector<unsigned int> CriticalPath;
    double temp;
    Cnode *largestNode;
    //list<Cnode*>::iterator largest;
    //list<Cnode*>::iterator secondLargest;
    list<Cnode *>::iterator smallest;
    list<Cnode *>::iterator secondSmallest;
    vector<Cnode *> *Children;
    long pathlength;
    vector<Cnode *> queue;
    bool memoryCheckPass = false, leaf = false;
    Cnode *nodeone;
    Cnode *nodetwo;
    bool DeadBreak, firstTime;

    //clock_t time;
    //time = clock();
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

        //time=clock();
        while (!Children->empty())
        { //initialize critical path
            temp = largestNode->GetParallelPart();
            for (vector<Cnode *>::iterator iter = Children->begin(); iter != Children->end(); ++iter)
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
        //time=clock()-time;
        //time=time/10000;
        //printf("  initializing critical path took me %d *10^4 clicks. \n",time);

        //        cout<<"Critical Path: ";
        //        for (vector<unsigned int>::iterator iter=CriticalPath.begin(); iter!=CriticalPath.end(); ++iter) {
        //            cout<<(*iter)<<" ";
        //        }
        //        cout<<endl;
        //        cout<<endl;

        Children = Qroot->GetChildren();
        pathlength = CriticalPath.size();

        //time=clock();
        for (unsigned int i = 1; i < pathlength; ++i)
        { //initialize vector L
            for (vector<Cnode *>::iterator iter = Children->begin(); iter != Children->end(); ++iter)
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
            for (vector<Cnode *>::iterator iter = Children->begin(); iter != Children->end(); ++iter)
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
                for (vector<Cnode *>::iterator iter = Children->begin(); iter != Children->end(); ++iter)
                {
                    Llist.push_back(*iter);
                    queue.push_back(*iter);
                }
            }
        }
        //time=clock()-time;
        //time=time/10000;
        //printf("  initializing L took me %d *10^4 clicks. \n",time);

        //        cout<<"List L: ";
        //        for (list<Cnode*>::iterator iter=Llist.begin(); iter!=Llist.end(); ++iter) {
        //            cout<<(*iter)->GetId()<<" ";
        //        }
        //        cout<<endl;
        //        cout<<endl;

        //cout<<"List L size: "<<Llist.size()<<", shortage: "<<shortage<<endl;
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
                memoryCheckPass = MemoryEnough(tree, (*smallest)->GetParent(), (*smallest), leaf, memory_size, chstart, childrenID);
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
                memoryCheckPass = MemoryEnough(tree, (*secondSmallest)->GetParent(), *secondSmallest, leaf, memory_size, chstart, childrenID);
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
                //time=clock();
                currentNode->MergetoParent();
                //time=clock()-time;
                //time=time/10000;
                //printf("  merge node to its parent took me %d *10^4 clicks. \n",time);
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

bool cmp_asapc(Cnode *a, Cnode *b) { return (a->GetMSminusComu() < b->GetMSminusComu()); };

double ASAP(Ctree *tree, unsigned int num_processors, unsigned int depth)
{ //depth should be at least 1
    list<Cnode *> PriorityQue;
    vector<Cnode *> BrokenEdges;
    unsigned long step_minimumMS = 0;
    double minimumMS = tree->GetRoot()->GetMSCost(true, true);
    //cout<<"Excuting sequentially, makespan "<<minimumMS<<endl;
    double temp;
    list<Cnode *> children_buffer;
    Cnode *currentNode;
    unsigned int add_child_deepth;
    list<Cnode *>::iterator node_position;
    unsigned int i;
    Cnode *LargestNode;

    vector<Cnode *> *children = tree->GetRoot()->GetChildren();

    while (children->size() == 1)
    { //avoid the linear chain
        children = children->front()->GetChildren();
    }

    for (vector<Cnode *>::iterator iter = children->begin(); iter != children->end(); ++iter)
    {
        PriorityQue.push_back(*iter);
        //cout<<"   add node "<<(*iter)->GetId()<<" into PQ."<<endl;
        children_buffer.push_back(*iter);
    }

    unsigned long buffer_size;
    unsigned int generation = 1;
    while (generation < depth)
    {
        buffer_size = children_buffer.size();
        while (buffer_size > 0)
        {
            children = children_buffer.front()->GetChildren();
            children_buffer.pop_front();
            for (vector<Cnode *>::iterator child = children->begin(); child != children->end(); ++child)
            {
                children_buffer.push_back(*child);
                PriorityQue.push_back(*child);
                //cout<<"   add node "<<(*child)->GetId()<<" into PQ."<<endl;
            }
            --buffer_size;
        }
        ++generation;
    }

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
            LargestNode = *max_element(PriorityQue.begin(), PriorityQue.end(), cmp_asapc); //computation weight minus communication cost
        }

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

        currentNode = LargestNode;
        node_position = find(PriorityQue.begin(), PriorityQue.end(), currentNode);
        PriorityQue.erase(node_position);
        //cout<<"   pop up node "<<currentNode->GetId()<<endl;
        children_buffer.clear();
        children_buffer.push_back(currentNode);
        add_child_deepth = 1;

        node_position = find(PriorityQue.begin(), PriorityQue.end(), currentNode->GetParent());
        while (node_position != PriorityQue.end())
        {
            PriorityQue.erase(node_position);
            ++add_child_deepth;
            currentNode = currentNode->GetParent();
            node_position = find(PriorityQue.begin(), PriorityQue.end(), currentNode->GetParent());
        }

        i = 0;
        while (i < depth - add_child_deepth)
        { //avoid insert nodes that are already in PQ, when depth is larger than 1
            ++i;
            buffer_size = children_buffer.size();
            while (buffer_size > 0)
            {
                children = children_buffer.front()->GetChildren();
                children_buffer.pop_front();
                for (vector<Cnode *>::iterator child = children->begin(); child != children->end(); ++child)
                {
                    children_buffer.push_back(*child);
                }
                --buffer_size;
            }
        }

        while (add_child_deepth >= 1)
        {
            buffer_size = children_buffer.size();
            while (buffer_size > 0)
            {
                children = children_buffer.front()->GetChildren();
                children_buffer.pop_front();
                for (vector<Cnode *>::iterator child = children->begin(); child != children->end(); ++child)
                {
                    children_buffer.push_back(*child);
                    PriorityQue.push_back(*child);
                    //cout<<"   add node "<<(*child)->GetId()<<" into PQ."<<endl;
                }
                --buffer_size;
            }
            --add_child_deepth;
        }
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

bool cmp_asap(Cnode *a, Cnode *b) { return (a->GetMSCost(false, false) < b->GetMSCost(false, false)); };

double ASAP(Ctree *tree, unsigned int num_processors)
{
    list<Cnode *> PriorityQue;
    vector<Cnode *> BrokenEdges;
    unsigned long step_minimumMS = 0;
    double minimumMS = tree->GetRoot()->GetMSCost(true, true);
    //cout<<"Excuting sequentially, makespan "<<minimumMS<<endl;
    double temp;
    Cnode *LargestNode;
    list<Cnode *>::iterator node_position;

    vector<Cnode *> *children = tree->GetRoot()->GetChildren();
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

//unsigned long AvoidChain(Ctree* tree){//it works, the first implementation
//    Cnode* root=tree->GetRoot();
//
//    root->BreakEdge();
//    unsigned long num_subtrees=0;
//    num_subtrees = HowmanySubtrees(tree, true);
//
//    Ctree* Qtreeobj = BuildQtree(tree);
//
//    //find the chain
//    Cnode* currentNode;
//    forward_list<Cnode*> Que;
//    Que.push_front(Qtreeobj->GetRoot());
//    vector<Cnode*>* children;
//    while (!Que.empty()) {
//        currentNode=Que.front();
//        Que.pop_front();
//
//        while (currentNode->GetChildren()->size()==1) {
//            currentNode=currentNode->GetChildren()->front();
//            tree->GetNode(currentNode->GetothersideID())->RestoreEdge();// Restore Edge
//            //cout<<"Restore edge "<<currentNode->GetothersideID()<<endl;
//            num_subtrees--;
//        }
//
//        children=currentNode->GetChildren();
//        for (vector<Cnode*>::iterator iter=children->begin();iter!=children->end();iter++){
//            Que.push_front((*iter));
//        }
//    }
//
//    delete Qtreeobj;
//
//    return num_subtrees;
//}

unsigned long AvoidChain(Ctree *tree)
{
    Cnode *root = tree->GetRoot();
    root->BreakEdge();
    unsigned long num_subtrees = 0;
    num_subtrees = HowmanySubtrees(tree, true);

    Ctree *Qtreeobj = BuildQtree(tree);
    const vector<Cnode *> *AllNodes = Qtreeobj->GetNodes();
    vector<Cnode *> *children;
    for (vector<Cnode *>::const_iterator iter = AllNodes->begin(); iter != AllNodes->end(); iter++)
    {
        children = (*iter)->GetChildren();
        if (children->size() == 1)
        {
            tree->GetNode(children->front()->GetothersideID())->RestoreEdge(); //restore edge
            //cout<<"resotre edge "<<children->front()->GetothersideID()<<endl;
            num_subtrees--;
        }
    }

    delete Qtreeobj;

    return num_subtrees;
}

bool cmp_larSav(Cnode *a, Cnode *b) { return (a->GetMSminusComu() > b->GetMSminusComu()); };

//double LarSav(Ctree* tree, unsigned int processor_number, unsigned int num_subtrees){
//    Cnode* root=tree->GetRoot();
//    if (processor_number<=num_subtrees) {
//        return root->GetMSCost(true, true);
//    }
//
//    root->GetMSCost(true, true);//update makespan
//
//    Ctree* Qtreeobj = BuildQtree(tree);
//
////    cout<<"id parentId ew msw nw"<<endl;
////    for (int i=1; i<=100; ++i) {
////        cout<<Qtreeobj->GetNode(i)->GetId()<<" "<<Qtreeobj->GetNode(i)->GetParentId()<<" "<<Qtreeobj->GetNode(i)->GetEW()<<" "<<Qtreeobj->GetNode(i)->GetMSW()<<" "<<Qtreeobj->GetNode(i)->GetNW()<<endl;
////        Qtreeobj->GetNode(i)->BreakEdge();
////    }
//
//    vector<Cnode*> CriticalPath;
//    list<Cnode*> listL;
//    list<Cnode*> que;
//    int idleProcessors=processor_number-num_subtrees;
//    Cnode* largestNode=root;
//    vector<Cnode*>* Children;
//    double temp;
//    bool childSubtreeEnd;
//    Cnode* Largest;
//    Cnode* secondLargest;
//    unsigned int tempid;
//    Cnode* currentNode;
//    Cnode* Qroot=Qtreeobj->GetRoot();
//    long QtreeSize;
////    clock_t time;
//    vector<Cnode*> tempStore;
//    vector<Cnode*> nodesLastSubtree;
//    bool DeadLock=true;
//    while (idleProcessors>0) {
//        DeadLock=true;
//        childSubtreeEnd=false;
//        CriticalPath.clear();
//        CriticalPath.push_back(root);
//        Children=Qroot->GetChildren();
//        Qroot->GetMSCost(true,true);//update makespan
//        largestNode=Qroot;
////        time = clock();
//        while (!Children->empty()) {//initialize critical path
//            temp=largestNode->GetParallelPart();
//            for (vector<Cnode*>::iterator iter=Children->begin(); iter!=Children->end(); ++iter) {
//                //cout<<(*iter)->GetId()<<" "<<endl;
//                if ((*iter)->GetMSCost(true, false)==temp) {
//                    largestNode=(*iter);
//                    break;
//                }
//            }
//            CriticalPath.push_back(tree->GetNode(largestNode->GetothersideID()));
//            Children=largestNode->GetChildren();
//        }
//
////        time=clock()-time;
////        time=time/10000;
////        printf("Build CriPat took me %d *10^4 clicks. \n",time);
//
////        cout<<endl;
////        cout<<"Critical path: ";
////        for (vector<Cnode*>::iterator iter=CriticalPath.begin(); iter!=CriticalPath.end(); ++iter) {
////            cout<<(*iter)->GetId()<<" ";
////        }
////        cout<<endl;
////        cout<<"Critical path size: "<<CriticalPath.size()<<endl;
////
////        cout<<"Idle processor number: "<<idleProcessors<<endl;
//
//        listL.clear();
////        time=clock();
//
//        QtreeSize=Qtreeobj->GetNodes()->size();
//        for (unsigned int i=2; i<=QtreeSize; ++i) {//remove ancestors of roots
//            currentNode=tree->GetNode(Qtreeobj->GetNode(i)->GetothersideID())->GetParent();
//            while (!currentNode->IsBorken()) {
//                currentNode->BreakEdge();
//                tempStore.push_back(currentNode);
//                que.push_back(currentNode);
//                currentNode=currentNode->GetParent();
//            }
//        }
//
//        for (unsigned int i=0; i<CriticalPath.size(); ++i) {//initialize set L
//            que.push_back(CriticalPath[i]);
//            while (!que.empty()) {
//                Children=que.front()->GetChildren();
//                que.pop_front();
//                for (vector<Cnode*>::iterator iter=Children->begin(); iter!=Children->end(); ++iter) {
//                    if (!(*iter)->IsBorken()) {
//                        listL.push_back((*iter));
//                        que.push_back((*iter));
//                    }
//                }
//            }
//        }
//
//        while (!tempStore.empty()) {
//            currentNode=tempStore.back();
//            currentNode->RestoreEdge();
//            tempStore.pop_back();
//        }
//
////        cout<<endl;
////        cout<<"List L: ";
////        for (list<Cnode*>::iterator iter=listL.begin(); iter!=listL.end(); ++iter) {
////            cout<<(*iter)->GetId()<<" ";
////        }
////        cout<<endl;
////        cout<<"List L size: "<<listL.size()<<endl;
//
//        if (listL.empty()) {
//            break;
//        }
//
////        time=clock()-time;
////        time=time/10000;
////        printf("Build listL took me %d *10^4 clicks. \n",time);
//
////        time=clock();
////        largestNode = *max_element(listL.begin(), listL.end(), cmp_larSav);//computation weight minus communication cost
//        listL.sort(cmp_larSav);//computation weight minus communication cost, non-increasing
////        cout<<endl;
////        for (list<Cnode*>::iterator iter=listL.begin(); iter!=listL.end(); advance(iter, 1)) {
////            cout<<(*iter)->GetId()<<"_"<<(*iter)->GetMSminusComu()<<" ";
////        }
////        cout<<endl;
////        time=clock()-time;
////        time=time/10000;
////        printf("Get max element took me %d *10^4 clicks. \n",time);
//
//        //cout<<"Nodes of the last subtree on Critical path: ";
//        Children=CriticalPath.back()->GetChildren();
//        que.clear();
//        nodesLastSubtree.clear();
//        que.assign(Children->begin(),Children->end());
//        while (!que.empty()) {
//            for (vector<Cnode*>::iterator iter=Children->begin(); iter!=Children->end(); ++iter) {
//                nodesLastSubtree.push_back((*iter));
//                que.push_back((*iter));
//                //cout<<(*iter)->GetId()<<" ";
//            }
//            Children=que.front()->GetChildren();
//            que.pop_front();
//        }
//        //cout<<endl;
//
//
//        while (!listL.empty()) {
//            tempid = listL.front()->GetId();
//            childSubtreeEnd=false;
//            for (vector<Cnode*>::iterator iter=nodesLastSubtree.begin(); iter!=nodesLastSubtree.end(); ++iter) {
//                if ((*iter)->GetId()==tempid) {
//                    childSubtreeEnd=true;
//                    break;
//                }
//            }
//
//            if (childSubtreeEnd==true) {//on the last subtree of critical path
//                currentNode=listL.front();
//                if ((idleProcessors>1)&&(currentNode->GetParent()->GetChildren()->size()>1)) {
//                    GetTwoLargestElementTypeone(currentNode->GetParent()->GetChildren(), Largest, secondLargest);
//
//                    temp=Qtreeobj->GetNode(CriticalPath.back()->GetothersideID())->GetMSW();
//                    Largest->BreakEdge();
//                    Largest->SetothersideID(Qtreeobj->GetNodes()->size()+1);
//                    secondLargest->BreakEdge();
//                    secondLargest->SetothersideID(Qtreeobj->GetNodes()->size()+2);
//
//                    //cout<<"Break edge "<<Largest->GetId()<<", "<<secondLargest->GetId()<<" of subtree "<<CriticalPath.back()->GetothersideID()<<", create Q node ";
//
//                    Cnode* newNodeone = new Cnode(CriticalPath.back()->GetothersideID(), 0, Largest->GetEW(), Largest->GetMSCost(false, false));//bug here
//                    newNodeone->SetId(Qtreeobj->GetNodes()->size()+1);
//                    //cout<<Qtreeobj->GetNodes()->size()+1<<", ";
//                    newNodeone->GetChildren()->clear();
//                    newNodeone->SetParent(Qtreeobj->GetNode(CriticalPath.back()->GetothersideID()));
//                    newNodeone->BreakEdge();
//                    newNodeone->SetothersideID(Largest->GetId());
//                    Qtreeobj->GetNode(CriticalPath.back()->GetothersideID())->AddChild(newNodeone);
//                    Qtreeobj->addNode(newNodeone);
//                    temp=temp-newNodeone->GetMSW();
//
//                    Cnode* newNodetwo = new Cnode(CriticalPath.back()->GetothersideID(), 0, secondLargest->GetEW(), secondLargest->GetMSCost(false, false));//bug here
//                    newNodetwo->SetId(Qtreeobj->GetNodes()->size()+1);
//                    //cout<<Qtreeobj->GetNodes()->size()+1<<endl;
//                    newNodetwo->GetChildren()->clear();
//                    newNodetwo->SetParent(Qtreeobj->GetNode(CriticalPath.back()->GetothersideID()));
//                    newNodetwo->BreakEdge();
//                    newNodetwo->SetothersideID(secondLargest->GetId());
//                    Qtreeobj->GetNode(CriticalPath.back()->GetothersideID())->AddChild(newNodetwo);
//                    Qtreeobj->addNode(newNodetwo);
//                    temp=temp-newNodetwo->GetMSW();
//                    Qtreeobj->GetNode(CriticalPath.back()->GetothersideID())->SetMSW(temp);
//
//                    //cout<<"Node "<<newNodeone->GetId()<<", msweight "<<newNodeone->GetMSW()<<", "<<"Node "<<newNodetwo->GetId()<<", msweight "<<newNodetwo->GetMSW()<<endl;
//
//                    idleProcessors=idleProcessors-2;
//                    DeadLock=false;
//                    break;
//                }else{
//                    //cout<<"list L pop node "<<listL.front()->GetId()<<endl;
//                    listL.pop_front();
//                }
//            }else{//not on the last subtree
//                largestNode=listL.front();
//                largestNode->BreakEdge();
//                largestNode->SetothersideID(Qtreeobj->GetNodes()->size()+1);
//                currentNode=largestNode->GetParent();
//                while (!currentNode->IsBorken()) {
//                    currentNode=currentNode->GetParent();
//                }
//                Cnode* newNode = new Cnode(currentNode->GetothersideID(), 0, largestNode->GetEW(), largestNode->GetMSCost(false, false));//bug here
//                //cout<<"Break edge "<<listL.front()->GetId()<<" of subtree "<<currentNode->GetothersideID()<<", ";
//                newNode->SetId(Qtreeobj->GetNodes()->size()+1);
//                //cout<<"create Q node "<<Qtreeobj->GetNodes()->size()+1<<endl;
//                newNode->GetChildren()->clear();
//                newNode->SetParent(Qtreeobj->GetNode(currentNode->GetothersideID()));
//                newNode->BreakEdge();
//                newNode->SetothersideID(largestNode->GetId());
//                Qtreeobj->GetNode(currentNode->GetothersideID())->AddChild(newNode);
//                Qtreeobj->addNode(newNode);
//                temp=Qtreeobj->GetNode(currentNode->GetothersideID())->GetMSW();
//                Qtreeobj->GetNode(currentNode->GetothersideID())->SetMSW(temp-newNode->GetMSW());
//
//                //cout<<"Node "<<newNode->GetId()<<", msweight "<<newNode->GetMSW()<<endl;
//
//                idleProcessors--;
//                DeadLock=false;
//                break;
//            }
//        }
//
//        if (DeadLock==true) {
//            break;
//        }
//
////        cout<<"---------------------------------------------"<<endl;
//
////        Qroot->GetMSCost(true,true);//update makespan
////        for (int i=1; i<=Qtreeobj->GetNodes()->size(); ++i) {
////            cout<<"node "<<i<<", MScost "<<Qtreeobj->GetNode(i)->GetMSCost(true, false)<<", SequentialPart "<<Qtreeobj->GetNode(i)->GetSequentialPart()<<", ParallelPart "<<Qtreeobj->GetNode(i)->GetParallelPart()<<endl;
////            for (vector<Cnode*>::iterator iter=Qtreeobj->GetNode(i)->GetChildren()->begin(); iter!=Qtreeobj->GetNode(i)->GetChildren()->end(); ++iter) {
////                cout<<"   child "<<(*iter)->GetId()<<", MScost "<<(*iter)->GetMSCost(true, false)<<endl;
////            }
////        }
//
//    }
//
//    temp=tree->GetRoot()->GetMSCost(true, true);
//    delete Qtreeobj;
//
//    return temp;
//}

bool EstimateDecrase(int idleP, Ctree *tree, vector<Cnode *> *criticalPath, bool *lastsubtree, Cnode **node_i, Cnode **node_j)
{
    //cout<<"   --------------estimate decrease in makespan-----------------"<<endl;
    *lastsubtree = false;
    bool MSdecreased = false;
    vector<Cnode *> *children;
    vector<double> decreaseSequence;
    double temp, decrease = -1;
    vector<Cnode *> tempQue;
    Cnode *lastSubtreeRoot = tree->GetNode(criticalPath->back()->GetothersideID());

    //cout<<"   Last subtree root "<<lastSubtreeRoot->GetId()<<endl;
    //nodes on the last subtree of critical path
    if (idleP > 1)
    { //has at least 2 idle processor
        tempQue.push_back(lastSubtreeRoot);
        vector<Cnode *>::iterator largestNode, secondLargest;
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

                for (vector<Cnode *>::iterator it = children->begin(); it != children->end(); ++it)
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
                for (vector<Cnode *>::iterator it = children->begin(); it != children->end(); ++it)
                {
                    tempQue.push_back(*it);
                }
            }
        }
    }

    //    if (decrease>=0) {
    //        cout<<"   if cut edge "<<(*node_i)->GetId()<<" and "<<(*node_j)->GetId()<<" of last subtree, decrease is "<<decrease<<endl;
    //    }else{
    //        cout<<"   fail on the last subtree."<<endl;
    //    }
    //
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
    Cnode *output_node;
    Cnode *subtreeRoot;
    Cnode *currentNode; //current node is on the path composed of critial path nodes
    Cnode *nodeOnPath;
    Cnode *SubtreeT = criticalPath->back();
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
                for (vector<Cnode *>::iterator it = children->begin(); it != children->end(); ++it)
                {
                    if ((*it)->GetId() != nodeOnPath->GetId() && (!(*it)->IsBorken()))
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
        } while (!currentNode->IsBorken());
        //cout<<"   }"<<endl;
    } while (subtreeRoot->GetId() != tree->GetRootId());
    //cout<<endl;

    //    if (decrease_othersubtrees>=0) {
    //        cout<<"   if cut edge "<<output_node->GetId()<<", decrease is "<<decrease_othersubtrees<<endl;
    //    }else{
    //        cout<<"   cut an edge, fail"<<endl;
    //    }

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

//Cnode* GetLargestSibling(Cnode* node){
//    unsigned int id=node->GetId();
//    Cnode* node_return;
//    vector<Cnode*>* Children=node->GetParent()->GetChildren();
//    double temp = 0;
//    double cost_sibling;
//    for (vector<Cnode*>::iterator iter=Children->begin(); iter!=Children->end(); ++iter) {
//        cost_sibling = (*iter)->GetMSCost(false,false);
//        if ((cost_sibling>=temp)&&((*iter)->GetId()!=id)) {
//            node_return = (*iter);
//            temp = cost_sibling;
//        }
//    }
//
//    return node_return;
//}

double SplitAgain(Ctree *tree, unsigned int processor_number, unsigned int num_subtrees)
{
    double MS_now;
    Cnode *root = tree->GetRoot();
    Ctree *Qtreeobj = BuildQtree(tree);

    vector<Cnode *> CriticalPath; //Q nodes on Critical Path

    Cnode *Qroot = Qtreeobj->GetRoot();
    Cnode *largestNode;
    Cnode *node_i;
    Cnode *node_j;
    Cnode *parent;
    double temp;
    vector<Cnode *> *Children;
    bool MSReduced, onLastSubtree;
    vector<Cnode *> tempVector;

    int idleProcessors = processor_number - num_subtrees;
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
            for (vector<Cnode *>::iterator iter = Children->begin(); iter != Children->end(); ++iter)
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
                //cout<<"cut edge "<<node_i->GetId()<<endl;
                node_i->BreakEdge(); //C<-C\cup C_k
                idleProcessors--;

                node_i->SetothersideID(Qtreeobj->GetNodes()->size() + 1);
                parent = node_i->GetParent();
                while (!parent->IsBorken())
                {
                    parent = parent->GetParent();
                }
                Cnode *Qparent = Qtreeobj->GetNode(parent->GetothersideID());
                Cnode *Qchild;
                Cnode *newNode = new Cnode(parent->GetothersideID(), 0, node_i->GetEW(), node_i->GetSequentialPart());
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
                        for (vector<Cnode *>::iterator iter = Children->begin(); iter != Children->end(); ++iter)
                        {
                            if ((*iter)->IsBorken())
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
                //cout<<"cut edge "<<node_i->GetId()<<" and edge "<<node_j->GetId()<<endl;
                node_i->BreakEdge(); //C<-C\cup C_k
                node_j->BreakEdge(); //C<-C\cup C_k
                idleProcessors = idleProcessors - 2;

                node_i->SetothersideID(Qtreeobj->GetNodes()->size() + 1);
                node_j->SetothersideID(Qtreeobj->GetNodes()->size() + 2);

                Cnode *newNodeone = new Cnode(CriticalPath.back()->GetId(), 0, node_i->GetEW(), node_i->GetSequentialPart());
                newNodeone->SetId(Qtreeobj->GetNodes()->size() + 1);
                newNodeone->GetChildren()->clear();
                newNodeone->SetParent(CriticalPath.back());
                newNodeone->BreakEdge();
                newNodeone->SetothersideID(node_i->GetId());
                CriticalPath.back()->AddChild(newNodeone);
                Qtreeobj->addNode(newNodeone);
                temp = CriticalPath.back()->GetMSW();
                temp = temp - newNodeone->GetMSW();

                Cnode *newNodetwo = new Cnode(CriticalPath.back()->GetId(), 0, node_j->GetEW(), node_j->GetSequentialPart());
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

void Immediately(Ctree *tree, unsigned long N, double *nwghts, double *ewghts, int *chstart, int *children, int *schedule, double m_availble, unsigned int &num_para_subtrees, vector<unsigned int> *brokenEdges)
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
                int *prntssub, *chstartsub, *chendsub, *childrensub;
                Ctree *subtree = BuildSubtree(tree, tree->GetNode(cur_task_id), subtree_size, &prntssub, &ewghtssub, &timewghtssub, &spacewghtssub, chstart, children);

                subtree_size = subtree->GetNodes()->size();
                //                        for (unsigned int index=1; index<=subtree_size; ++index) {
                //                            cout<<index<<" "<<subtree->GetNode(index)->GetParentId()<<" "<<subtree->GetNode(index)->GetNW()<<" "<<subtree->GetNode(index)->GetMSW()<<" "<<subtree->GetNode(index)->GetEW()<<endl;
                //                        }

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
                po_construct(subtree_size, prntssub, &chstartsub, &chendsub, &childrensub, &rootid);

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
                delete[] chendsub;
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

void MemoryCheck(Ctree *tree, int *chstart, int *children, double const memory_size, io_method_t method)
{ //chstart, children are not modified
    vector<Cnode *> subtreeRoots;
    Cnode *currentnode;
    Cnode *subtreeRoot;
    int rootid;
    tree->GetRoot()->BreakEdge();

    //cout<<"Subtrees' roots: ";
    unsigned long treeSize = tree->GetNodes()->size();
    for (unsigned int i = treeSize; i >= 1; --i)
    {
        currentnode = tree->GetNode(i);
        if (currentnode->IsBorken())
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
        Ctree *subtree = BuildSubtree(tree, subtreeRoot, treeSize, &prnts, &ewghts, &timewghts, &spacewghts, chstart, children);

        subtreeSize = subtree->GetNodes()->size();
        int *schedule_copy = new int[subtreeSize + 1];
        maxoutD = MaxOutDegree(subtree, true);
        schedule_f->clear();
        count = 0;
        MinMem(subtree, maxoutD, memory_required, *schedule_f, true, count);

        int *chstartsub, *chendsub, *childrensub;
        po_construct(subtreeSize, prnts, &chstartsub, &chendsub, &childrensub, &rootid);

        //cout<<"Subtree "<<subtreeRoot->GetId()<<" needs memory "<<memory_required;
        if (memory_required > memory_size)
        {
            //cout<<", larger than what is available: "<<memory_size<<endl;

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
        delete[] chendsub;
        delete[] childrensub;
        delete subtree;
    }
    delete schedule_f;

    for (vector<unsigned int>::iterator iter = BrokenEdgesID.begin(); iter != BrokenEdgesID.end(); ++iter)
    {
        tree->GetNode(*iter)->BreakEdge();
    }
}
void breakSubtreeFurther(int *schedule_copy, int subtreeSize)
{
}

std::map<int, int> MemoryCheckA2(Ctree *tree, int *chstart, int *children, vector<double> memory_sizes, io_method_t method, bool skipBig)
{ //chstart, children are not modified
    vector<Cnode *> subtreeRoots;
    vector<Cnode *> subtreeRootsSkipped;
    Cnode *currentnode;
    Cnode *subtreeRoot;
    int rootid;
    tree->GetRoot()->BreakEdge();

    //cout<<"Subtrees' roots: ";
    unsigned long treeSize = tree->GetNodes()->size();
    for (unsigned int i = treeSize; i >= 1; --i)
    {
        currentnode = tree->GetNode(i);
        if (currentnode->IsBorken())
        {
            //cout<<i<<" ";
            subtreeRoots.push_back(currentnode);
            cout << "root " << currentnode->GetMSCost() << endl;
        }
    }
    //cout<<endl;
    sort(subtreeRoots.begin(), subtreeRoots.end(), [](Cnode *lhs, Cnode *rhs)
         { return lhs->GetMSCost() < rhs->GetMSCost(); });
    sort(memory_sizes.begin(), memory_sizes.end());
    //std::sort(vec.begin(), vec.end());

    std::map<int, int> taskToPrc;

    for (int i = 0; i < treeSize; i++)
    {
        taskToPrc.insert(pair<int, int>(i, -1));
    }

    double maxoutD, memory_required;
    schedule_t *schedule_f = new schedule_t();
    uint64_t count;
    unsigned int com_freq;
    unsigned long subtreeSize;
    list<int>::iterator ite_sche;
    vector<unsigned int> BrokenEdgesID;
    double IO_volume;
    int currentProcessor;
    currentProcessor = 0;
    cout<<endl<<"first small trees"<<endl;
    while (!subtreeRoots.empty())
    {
        double currentMem = memory_sizes[currentProcessor];
        subtreeRoot = subtreeRoots.back();
        subtreeRoots.pop_back();

        double *ewghts, *timewghts, *spacewghts;
        int *prnts;
        Ctree *subtree = BuildSubtree(tree, subtreeRoot, treeSize, &prnts, &ewghts, &timewghts, &spacewghts, chstart, children);

        subtreeSize = subtree->GetNodes()->size();
        int *schedule_copy = new int[subtreeSize + 1];
        maxoutD = MaxOutDegree(subtree, true);
        schedule_f->clear();
        count = 0;
        MinMem(subtree, maxoutD, memory_required, *schedule_f, true, count);

        int *chstartsub, *chendsub, *childrensub;
        po_construct(subtreeSize, prnts, &chstartsub, &chendsub, &childrensub, &rootid);

        cout << "Subtree " << subtreeRoot->GetId() << " needs memory " << memory_required;
        if (memory_required > currentMem)
        {
            if (!skipBig)
            {
                cout << ", larger than what is available: " << currentMem << " on proc " << currentProcessor << endl;

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
                    IO_volume = IOCounterWithVariableMem(subtree, subtreeSize + 1, spacewghts, ewghts, chstartsub, childrensub, schedule_copy, memory_sizes, currentProcessor, taskToPrc, false, true, com_freq, &BrokenEdgesID, FIRST_FIT);
                    break;
                case LARGEST_FIT:
                    IO_volume = IOCounterWithVariableMem(subtree, subtreeSize + 1, spacewghts, ewghts, chstartsub, childrensub, schedule_copy, memory_sizes, currentProcessor, taskToPrc, false, true, com_freq, &BrokenEdgesID, LARGEST_FIT);
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
            taskToPrc.at(subtreeRoot->GetId()) = currentProcessor;
            currentProcessor++;
        }
        //cout<<endl;
        //
        cout << "Now big trees" << endl;
        while (!subtreeRootsSkipped.empty())
        {
            double currentMem = memory_sizes[currentProcessor];
            subtreeRoot = subtreeRootsSkipped.back();
            subtreeRootsSkipped.pop_back();
            cout << ", larger than what is available: " << currentMem << " on proc " << currentProcessor << endl;

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
                IO_volume = IOCounterWithVariableMem(subtree, subtreeSize + 1, spacewghts, ewghts, chstartsub, childrensub, schedule_copy, memory_sizes, currentProcessor, taskToPrc, false, true, com_freq, &BrokenEdgesID, FIRST_FIT);
                break;
            case LARGEST_FIT:
                IO_volume = IOCounterWithVariableMem(subtree, subtreeSize + 1, spacewghts, ewghts, chstartsub, childrensub, schedule_copy, memory_sizes, currentProcessor, taskToPrc, false, true, com_freq, &BrokenEdgesID, LARGEST_FIT);
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
        delete[] chendsub;
        delete[] childrensub;
        delete subtree;
    }
    delete schedule_f;

    for (vector<unsigned int>::iterator iter = BrokenEdgesID.begin(); iter != BrokenEdgesID.end(); ++iter)
    {
        tree->GetNode(*iter)->BreakEdge();
    }

    // cout <<"task to proc "<<endl;
    //for(int i=0; i<taskToPrc.size(); i++){
    //  if(taskToPrc.at(i)!=-1)
    // cout<< i<<" " <<taskToPrc.at(i)<< " \t ";
    //}

    int count1 = std::count(taskToPrc.begin(), taskToPrc.end(), CompareMapEntries(1));
    int count2 = std::count(taskToPrc.begin(), taskToPrc.end(), CompareMapEntries(2));
    int count3 = std::count(taskToPrc.begin(), taskToPrc.end(), CompareMapEntries(3));
    int count4 = std::count(taskToPrc.begin(), taskToPrc.end(), CompareMapEntries(4));
    int count5 = std::count(taskToPrc.begin(), taskToPrc.end(), CompareMapEntries(5));
    int count6 = std::count(taskToPrc.begin(), taskToPrc.end(), CompareMapEntries(6));
    int count7 = std::count(taskToPrc.begin(), taskToPrc.end(), CompareMapEntries(7));
    cout << "counts " << count1 << " " << count2 << " " << count3 << " " << count4 << " " << count5 << " " << count6 << " " << count7 << endl;
    return taskToPrc;
}

unsigned int HowmanySubtrees(const Ctree *tree, bool quiet)
{
    unsigned int number_subtrees = 0;
    tree->GetRoot()->BreakEdge();
    const vector<Cnode *> *Nodes = tree->GetNodes();
    if (quiet == false)
    {
        cout << "Broken Edges { ";
    }
    for (auto it = Nodes->begin(); it != Nodes->end(); ++it)
    {
        if ((*it)->IsBorken())
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

double Sequence(Cnode *root)
{
    return root->GetMSCost();
}
