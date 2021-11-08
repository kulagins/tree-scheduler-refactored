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

void RunWithClusterConfig(bool skipBigTrees, Tree *treeobj,
                          io_method_t method) {
    if (Cluster::getFixedCluster()->isHomogeneous())
        MemoryCheck(treeobj, method);
    else
        MemoryCheckA2(treeobj, Cluster::getFixedCluster(), method, skipBigTrees);
}

void buildFixedClusterWithSpeedsAndMemory(double CCR, unsigned int num_processors, Tree *treeobj) {
    Cluster *cluster = new Cluster(num_processors, true);
    double minMem;
    map<int, int> processor_speeds = Cluster::buildProcessorSpeeds(num_processors);
    cluster->SetBandwidth(CCR, treeobj);
    Cluster::setFixedCluster(cluster);

    Cluster::getFixedCluster()->SetBandwidth(CCR, treeobj);
    double maxoutd = MaxOutDegree(treeobj, true);
    schedule_t *temp_schedule = new schedule_t();
    MinMem(treeobj, maxoutd, minMem, *temp_schedule, true);

    vector<double> memorySizes = Cluster::buildMemorySizes(maxoutd, minMem, num_processors);
    //Fix, for now we consider the non-homog cluster homogeneuos
    Cluster::getFixedCluster()->setMemorySizes(memorySizes);
    delete temp_schedule;
}

int main(int argc, char **argv) {
    string stage1, stage2 = "FirstFit", stage3;

    string dir = argv[1];
    double CCR = atof(argv[3]);
    double NPR = atof(argv[4]);
   // bool skipBigTrees = (atoi(argv[6]) == 1);
    //  cout.precision(0);
    // cout.setf(ios::fixed);

    ifstream OpenFile(dir + argv[2]);
    bool quiet = false;

    int tree_size = 0;

    list<Task *> parallelSubtrees;
    string treename;
    unsigned int num_processors;
    double makespan, maxoutd, minMem;
    clock_t time;

    unsigned int number_subtrees;

    cout.precision(20);

    std::cout << "TreeName " << "NPR " << "CCR " << "MemoryConstraint " << "AmountSubtrees " << "AmountProcessors "
              << "Makespan " << "Stage1 " << "Stage2 " << "Stage3 " << "TimeConsuming" << std::endl;

    std::vector<int> brokenEdges;
    do {
        OpenFile >> treename;

        Tree *tree = read_tree((dir + treename).c_str());
        Tree *untouchedTree = read_tree((dir + treename).c_str());
        Tree::setOriginalTree(untouchedTree);

        num_processors = ceil(tree_size / NPR);
        if (num_processors < 3) {
            num_processors = 3;
        }

        buildFixedClusterWithSpeedsAndMemory(CCR, num_processors, tree);

        time = clock();
        //stage1==SplitSubtrees
        unsigned long sequentialLen;
        makespan = tree->getRoot()->SplitSubtrees(false, parallelSubtrees,
                                                  sequentialLen);// for counting how many subtrees produced, twolevel is set as false
        // stage1 == "ImprovedSplit") {
        //po_construct(tree_size, prnts, &chstart, &chend, &children, &root);
        makespan = tree->ImprovedSplit();
        //makespan = ImprovedSplit(tree, num_processors);
        //stage1 == "AvoidChain") {
        makespan = tree->ASAP();
        //   number_subtrees = HowmanySubtrees(tree, true);
        //  time = clock() - time;
        //  std::cout << treename << " " << NPR << " " << CCR << " NA " << number_subtrees << " " << num_processors
        //             << " " << makespan << " " << "ASAP NA NA " << time << std::endl;

        //   time = clock();
        //number_subtrees = AvoidChain(tree);
        //  makespan = tree->GetRoot()->GetMSCost(true, true);
        // }
        time = clock() - time;

        number_subtrees = tree->HowmanySubtrees(quiet);
        std::cout << treename << " " << NPR << " " << CCR << " NA " << number_subtrees << " " << num_processors << " "
                  << makespan << " " << stage1 << " NA NA " << time << std::endl;


        maxoutd = MaxOutDegree(tree, true);
        schedule_t *schedule_f = new schedule_t();

        MinMem(tree, maxoutd, minMem, *schedule_f, true);
        delete schedule_f;
        vector<Tree *> trees;
        Tree *treeobj2 = read_tree((dir + treename).c_str());;
        Tree *treeobj3 = read_tree((dir + treename).c_str());;
        trees.push_back(treeobj3);
        trees.push_back(treeobj2);
        trees.push_back(tree);

        for (int i = 0; i < 3; ++i) {
            tree = trees.back();
            trees.pop_back();

            time = clock();
            if (stage2 == "LargestFirst") {
                MemoryCheck(tree, LARGEST_FIT);
            } else if (stage2 == "FirstFit") {
                MemoryCheck(tree, FIRST_FIT);
            } else if (stage2 == "Immediately") {
                MemoryCheck(tree, IMMEDIATELY);
            }

            time = clock() - time;
            number_subtrees = tree->HowmanySubtrees(true);
            makespan = tree->getRoot()->getMakespanCost(true, true);
            std::cout << treename << " " << NPR << " " << CCR << " " << " memConstraint " << " " << number_subtrees
                      << " " << num_processors << " " << makespan << " " << stage1 << " " << stage2 << " NA "
                      << time
                      << std::endl;

            time = clock();
            if (number_subtrees > num_processors) {
                stage3 = "Merge";
                makespan = tree->MergeV2(number_subtrees, num_processors,
                                         Cluster::getFixedCluster()->getFirstFreeProcessor()->getMemorySize(),
                                         true);
                //Merge(tree, number_subtrees, num_processors, memorySize, chstart, children, true);
            } else if (number_subtrees == num_processors) {
                stage3 = "Nothing";
            } else {
                stage3 = "SplitAgain";
                makespan = tree->SplitAgain();
            }
            time = clock() - time;

            number_subtrees = tree->HowmanySubtrees(true);
            std::cout << treename << " " << NPR << " " << CCR << " " << "mem constraint" << " " << number_subtrees
                      << " " << num_processors << " " << makespan << " " << stage1 << " " << stage2 << " " << stage3
                      << " " << time << std::endl;
            delete tree;
            delete untouchedTree;
        }
    } while (OpenFile.good());
    OpenFile.close();
    exit(EXIT_SUCCESS);
}