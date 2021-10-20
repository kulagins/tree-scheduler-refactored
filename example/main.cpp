// simulation sequence+stage2+splitAgain/Merge
//  Created by changjiang GOU on 10/05/2018.
//  Copyright © 2018 Changjiang GOU. All rights reserved.

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <algorithm>
#include <stdlib.h>
#include "lib-io-tree.h"
#include "heuristics.h"




void RunWithClusterConfig(bool skipBigTrees, int *chstart, int *children, Tree *treeobj,
                          Cluster *cluster, io_method_t method)

{
    if (cluster->isHomogeneous())
        MemoryCheck(treeobj, chstart, children, cluster, method);
    else
        MemoryCheckA2(treeobj, chstart, children,  cluster, method, skipBigTrees);
}



void actualActions(double CCR, unsigned int num_processors, double *ewghts, double *spacewghts, double *timewghts, int *prnts, int tree_size, bool skipBigTrees)
{
    clock_t time;
    unsigned int number_subtrees;

    int *chstart, *chend, *children;
    int root;

    double minMem;
    uint64_t count;
    string stage2heuristic;
    vector<double> memorySizes;
    list<Task *> parallelSubtrees;
    std::map<int, int> processor_speeds = Cluster::buildProcessorSpeeds(num_processors);
    std::map<int, int> taskToPrc;
    std::map<int, bool> isProcBusy;

    double makespan;
    double maxoutd;

    SetBandwidth(CCR, tree_size, ewghts, timewghts);

    Tree *treeobj = new Tree(tree_size, prnts, spacewghts, ewghts, timewghts);
    treeobj->setOriginalTree(treeobj);

    maxoutd = MaxOutDegree(treeobj, true);

    po_construct(tree_size, prnts, &chstart, &chend, &children, &root);
    time = clock();
    makespan = treeobj->GetRoot()->GetMSCost();
    number_subtrees = 1;
    time = clock() - time;

    //<< " " << NPR << " " << CCR << " NA " << number_subtrees << " " << num_processors << " " << makespan << " Sequence " << time << endl;

    schedule_t *schedule_f = new schedule_t();
    count = 0;
    MinMem(treeobj, maxoutd, minMem, *schedule_f, true, count);
    delete schedule_f;
    delete treeobj;

    memorySizes = Cluster::buildMemorySizes(maxoutd, minMem, num_processors);
    Cluster *cluster = new Cluster(memorySizes);
    for (int stage2Method = 0; stage2Method < 1; ++stage2Method)
    {

        Tree *treeobj = new Tree(tree_size, prnts, spacewghts, ewghts, timewghts);
        time = clock();
        switch (stage2Method)
        {
        case 0:
            stage2heuristic = "FIRST_FIT";
            RunWithClusterConfig(skipBigTrees, chstart, children, treeobj, cluster, FIRST_FIT);
            break;
        case 1:
            stage2heuristic = "LARGEST_FIT";
            RunWithClusterConfig(skipBigTrees, chstart, children, treeobj, cluster, LARGEST_FIT);
            break;
        case 2:
            stage2heuristic = "IMMEDIATELY";
            RunWithClusterConfig(skipBigTrees, chstart, children, treeobj, cluster, IMMEDIATELY);
            break;

        default:
            stage2heuristic = "FIRST_FIT";
            RunWithClusterConfig(skipBigTrees, chstart, children, treeobj, cluster, IMMEDIATELY);
            break;
        }

        time = clock() - time;

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
            std::cout << "w merge "
                      << "#subtrees: " << number_subtrees << ", #numberProcessors; " << num_processors << " makespan: " << makespan << endl;
        }
        else if (number_subtrees == num_processors)
        {
            std::cout << "w equal "
                      << "#subtrees: " << number_subtrees << ", #numberProcessors; " << num_processors << " makespan: " << makespan << endl;
        }
        else
        {
            time = clock();
            makespan = treeobj->SplitAgain(num_processors, number_subtrees);
            time = clock() - time;
            number_subtrees = treeobj->HowmanySubtrees(true);
            std::cout << "w split "
                      << "#subtrees: " << number_subtrees << ", #numberProcessors; " << num_processors << " makespan: " << makespan << endl;
        }

        treeobj->printBrokenEdges();
        delete treeobj;

        delete[] chstart;
        delete[] chend;
        delete[] children;
    }
}

int main(int argc, const char *argv[])
{
    int tree_size = 0;
    int *prnts;
    unsigned int num_processors;
    string dir = argv[1];
    string treename;
    double *ewghts, *spacewghts, *timewghts;
    double CCR = atof(argv[3]);
    double NPR = atof(argv[4]);

    bool skipBigTrees = (atoi(argv[6]) == 1);
    cout.precision(0);
    cout.setf(ios::fixed);

    //  std::cout << " AmountSubtrees AmountProcessors Makespan Heuristic TimeConsuming" << std::endl;

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

            actualActions(CCR, num_processors, ewghts, spacewghts, timewghts, prnts, tree_size, skipBigTrees);
            delete[] prnts;
            delete[] ewghts;
            delete[] spacewghts;
            delete[] timewghts;
        }
    } while (OpenFile.good());
    OpenFile.close();

    return 0;
}
