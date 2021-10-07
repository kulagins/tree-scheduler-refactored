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

vector<double> buildMemorySizes(double maxoutd, double minMem, int num_processors)
{
    cout << "max deg " << maxoutd << ", MinMem " << minMem << endl;
    double cumulativeMem = 0;
    vector<double> memSizes(num_processors);
    memSizes.resize(num_processors);
    maxoutd = maxoutd / 4;
    //cout << "minProc " << maxoutd << " " << (maxoutd + minMem) / 2 << " " << minMem << endl;
    for (int k = 0; k < num_processors / 3; k++)
    {
        memSizes[k] = maxoutd; 
        cumulativeMem += memSizes[k];
    }
    for (int k = num_processors / 3; k < 2 * num_processors / 3; k++)
    {
        memSizes[k] = (maxoutd + minMem) / 2;
        cumulativeMem += memSizes[k];
    }
    for (int k = 2 * num_processors / 3; k < num_processors; k++)
    {
        memSizes[k] = minMem;
        cumulativeMem += memSizes[k];
    }
    cout << "cumulative mem in system: " << cumulativeMem << endl;
    return memSizes;
}

std::map<int, int> buildProcessorSpeeds(int num_processors)
{
    std::map<int, int> procSpeeds;
    for (int k = 0; k < num_processors / 3; k++)
    {
        procSpeeds.insert(pair<int, int>(k, 1));
    }
    for (int k = num_processors / 3; k < 2 * num_processors / 3; k++)
    {
        procSpeeds.insert(pair<int, int>(k, 2));
    }
    for (int k = 2 * num_processors / 3 + 1; k < num_processors; k++)
    {
        procSpeeds.insert(pair<int, int>(k, 3));
    }

    return procSpeeds;
}

void RunWithClusterConfig(int clusterConfig, bool skipBigTrees, int *chstart, int *children, Ctree *treeobj, vector<double> memorySizesA2, std::map<int, int> &taskToPrc, std::map<int, bool> &isProcBusy, io_method_t method)
{
    switch (clusterConfig)
    {
    case 1:
        MemoryCheck(treeobj, chstart, children, memorySizesA2[0], method);
        break;
    case 2:
        MemoryCheckA2(treeobj, chstart, children, memorySizesA2, method, skipBigTrees, taskToPrc, isProcBusy);
        break;
    case 3:
    default:
        throw std::invalid_argument("not implemented");
    }
}

void printBrokenEdges(Ctree *tree)
{
    cout << "Print broken edges" << endl;
    unsigned long treeSize = tree->GetNodes()->size();
    for (unsigned int i = treeSize; i >= 1; --i)
    {
        Cnode *currentnode = tree->GetNode(i);
        if (currentnode->IsBroken())
        {
            cout << i << " ";

            //cout << "root " << currentnode->GetMSCost() << endl;
        }
    }
    cout << "End" << endl;
}

void actualActions(double CCR, double NPR, unsigned int num_processors, double *ewghts, double *spacewghts, double *timewghts, int *prnts, int tree_size, bool skipBigTrees, int clusterConfig)
{
    clock_t time;
    unsigned int number_subtrees;

    int *chstart, *chend, *children;
    int root;

    double minMem;
    uint64_t count;
    string stage2heuristic;
    vector<double> memorySizes;
    list<Cnode *> parallelSubtrees;
    unsigned long sequentialLen;
    std::map<int, int> processor_speeds = buildProcessorSpeeds(num_processors);
    std::map<int, int> taskToPrc;
    std::map<int, bool> isProcBusy;

    double makespan;
    double maxoutd;

    SetBandwidth(CCR, tree_size, ewghts, timewghts);

    Ctree *treeobj = new Ctree(tree_size, prnts, spacewghts, ewghts, timewghts);
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

    memorySizes = buildMemorySizes(maxoutd, minMem, num_processors);
    for (int stage2Method = 0; stage2Method < 1; ++stage2Method)
    {

        Ctree *treeobj = new Ctree(tree_size, prnts, spacewghts, ewghts, timewghts);

        time = clock();
        switch (stage2Method)
        {
        case 0:
            stage2heuristic = "FIRST_FIT";
            RunWithClusterConfig(clusterConfig, skipBigTrees, chstart, children, treeobj, memorySizes, taskToPrc, isProcBusy, FIRST_FIT);
            break;
        case 1:
            stage2heuristic = "LARGEST_FIT";
            RunWithClusterConfig(clusterConfig, skipBigTrees, chstart, children, treeobj, memorySizes, taskToPrc, isProcBusy, LARGEST_FIT);
            break;
        case 2:
            stage2heuristic = "IMMEDIATELY";
            RunWithClusterConfig(clusterConfig, skipBigTrees, chstart, children, treeobj, memorySizes, taskToPrc, isProcBusy, IMMEDIATELY);
            break;

        default:
            stage2heuristic = "FIRST_FIT";
            RunWithClusterConfig(clusterConfig, skipBigTrees, chstart, children, treeobj, memorySizes, taskToPrc, isProcBusy, IMMEDIATELY);
            break;
        }

        time = clock() - time;

        makespan = treeobj->GetRoot()->GetMSCost(true, true);
        number_subtrees = HowmanySubtrees(treeobj, true);
        // std::cout << "after 2nd step "
        //       << number_subtrees << " " << num_processors << " " << makespan << " " << stage2heuristic << "+Nothing " << 0 << endl;

        if (number_subtrees > num_processors)
        {
            time = clock();
            makespan = MergeV2(treeobj, number_subtrees, num_processors, memorySizes[0], chstart, children, true);
            time = clock() - time;
            number_subtrees = HowmanySubtrees(treeobj, true);
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
            makespan = SplitAgain(treeobj, num_processors, number_subtrees);
            time = clock() - time;
            number_subtrees = HowmanySubtrees(treeobj, true);
            std::cout << "w split "
                      << "#subtrees: " << number_subtrees << ", #numberProcessors; " << num_processors << " makespan: " << makespan << endl;
        }

        printBrokenEdges(treeobj);
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

            actualActions(CCR, NPR, num_processors, ewghts, spacewghts, timewghts, prnts, tree_size, skipBigTrees, clusterConfig);
            delete[] prnts;
            delete[] ewghts;
            delete[] spacewghts;
            delete[] timewghts;
        }
    } while (OpenFile.good());
    OpenFile.close();

    return 0;
}