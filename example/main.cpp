// simulation sequence+stage2+splitAgain/Merge
//  Created by changjiang GOU on 10/05/2018.
//  Copyright © 2018 Changjiang GOU. All rights reserved.

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <algorithm>
#include <stdlib.h>
#include "../include/lib-io-tree.h"
#include "../include/heuristics.h"

void RunWithClusterConfig(bool skipBigTrees, int *chstart, int *children, Tree *treeobj,
                          io_method_t method)

{
    if (Cluster::getFixedCluster()->isHomogeneous())
        MemoryCheck(treeobj, chstart, children, Cluster::getFixedCluster(), method);
    else
        MemoryCheckA2(treeobj, chstart, children,  Cluster::getFixedCluster(), method, skipBigTrees);
}

void firstStep(double CCR, unsigned int num_processors, double *ewghts, double *spacewghts, double *timewghts,
               int *prnts, int tree_size, int *&chstart, int *&chend, int *&children, string &stage2heuristic,
               vector<double> &memorySizes) {
    Cluster *cluster;
    cluster= new Cluster(num_processors, true);
    clock_t time;
    unsigned int number_subtrees;
    int root;

    double minMem;
    uint64_t count;
    list<Task *> parallelSubtrees;
    map<int, int> processor_speeds = Cluster::buildProcessorSpeeds(num_processors);
    map<int, int> taskToPrc;
    map<int, bool> isProcBusy;

    double makespan;
    double maxoutd;
    cluster->SetBandwidth(CCR, tree_size, ewghts, timewghts);
    Cluster::setFixedCluster(cluster);

    Tree *treeobj = new Tree(tree_size, prnts, spacewghts, ewghts, timewghts);
    Tree::setOriginalTree(treeobj);

    maxoutd = MaxOutDegree(treeobj, true);

    po_construct(tree_size, prnts, &chstart, &chend, &children, &root);
    time = clock();
    unsigned long sequentialLen;
    makespan = treeobj->GetRoot()->SplitSubtrees(false, parallelSubtrees, sequentialLen);
    time = clock() - time;
    cout<<"makespan "<<makespan<<endl;

    schedule_t *schedule_f = new schedule_t();
    count = 0;
    MinMem(treeobj, maxoutd, minMem, *schedule_f, true, count);
    delete schedule_f;
    delete treeobj;

    memorySizes = Cluster::buildMemorySizes(maxoutd, minMem, num_processors);
    cluster->setMemorySizes(memorySizes);
}

void thirdStep(unsigned int num_processors, clock_t time, unsigned int number_subtrees, int *chstart, int *children,
               const vector<double> &memorySizes, double makespan, Tree *treeobj) {
    makespan = treeobj->GetRoot()->GetMSCost(true, true);
    number_subtrees = treeobj->HowmanySubtrees(true);
    // std::cout << "after 2nd step "
//       << number_subtrees << " " << num_processors << " " << makespan << " " << stage2heuristic << "+Nothing " << 0 << endl;

    if (number_subtrees > num_processors)
    {
        time = clock();
        makespan = treeobj->MergeV2(number_subtrees, num_processors, memorySizes[0], chstart, children, true);
        time = clock() - time;
        number_subtrees = treeobj->HowmanySubtrees(true);
        cout << "w merge "
                  << "#subtrees: " << number_subtrees << ", #numberProcessors; " << num_processors << " makespan: " << makespan << endl;
    }
    else if (number_subtrees == num_processors)
    {
        cout << "w equal "
                  << "#subtrees: " << number_subtrees << ", #numberProcessors; " << num_processors << " makespan: " << makespan << endl;
    }
    else
    {
        time = clock();
        makespan = treeobj->SplitAgain();
        time = clock() - time;
        number_subtrees = treeobj->HowmanySubtrees(true);
        cout << "w split "
                  << "#subtrees: " << number_subtrees << ", #numberProcessors; " << num_processors << " makespan: " << makespan << endl;
    }

    treeobj->printBrokenEdges();
    delete treeobj;
}

void secondStep(bool skipBigTrees, int *chstart, int *children, string &stage2heuristic, int stage2Method, Tree *treeobj) {
    switch (stage2Method)
    {
        case 0:
            stage2heuristic = "FIRST_FIT";
            RunWithClusterConfig(skipBigTrees, chstart, children, treeobj, FIRST_FIT);
            break;
        case 1:
            stage2heuristic = "LARGEST_FIT";
            RunWithClusterConfig(skipBigTrees, chstart, children, treeobj,  LARGEST_FIT);
            break;
        case 2:
            stage2heuristic = "IMMEDIATELY";
            RunWithClusterConfig(skipBigTrees, chstart, children, treeobj,  IMMEDIATELY);
            break;

        default:
            stage2heuristic = "FIRST_FIT";
            RunWithClusterConfig(skipBigTrees, chstart, children, treeobj,  IMMEDIATELY);
            break;
    }
}

int main(int argc, const char *argv[])
{
    int tree_size = 0;
    int *prnts;
    clock_t time;
    unsigned int number_subtrees;
    unsigned int num_processors;
    string dir = argv[1];
    string treename;
    double *ewghts, *spacewghts, *timewghts;
    double CCR = atof(argv[3]);
    double NPR = atof(argv[4]);

    bool skipBigTrees = (atoi(argv[6]) == 1);
    cout.precision(0);
    cout.setf(ios::fixed);

    ifstream OpenFile(dir + argv[2]);
    do
    {
        OpenFile >> treename;
        cout << treename << endl;
        for (int clusterConfig = 1; clusterConfig <= 2; clusterConfig++)
        {
            cout << "clusterConfig: " << clusterConfig << endl;

            parse_tree((dir + treename).c_str(), &tree_size, &prnts, &spacewghts, &ewghts, &timewghts);

            num_processors = ceil(tree_size / NPR);
            if (num_processors < 3)
            {
                num_processors = 3;
            }
            int *chstart;
            int *chend;
            int *children;
            string stage2heuristic;
            vector<double> memorySizes;
            double makespan;

            firstStep(CCR, num_processors, ewghts, spacewghts, timewghts, prnts, tree_size, chstart, chend, children,
                      stage2heuristic,
                      memorySizes);

            for (int stage2Method = 0; stage2Method < 1; ++stage2Method)
            {
                Tree *treeobj = new Tree(tree_size, prnts, spacewghts, ewghts, timewghts);
                time = clock();
                secondStep(skipBigTrees, chstart, children, stage2heuristic, stage2Method, treeobj);
                time = clock() - time;
                thirdStep(num_processors, time, number_subtrees, chstart, children, memorySizes, makespan, treeobj);

                delete[] chstart;
                delete[] chend;
                delete[] children;
            }
            delete[] prnts;
            delete[] ewghts;
            delete[] spacewghts;
            delete[] timewghts;
        }
    } while (OpenFile.good());
    OpenFile.close();

    return 0;
}
