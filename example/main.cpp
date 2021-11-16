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
#include "../include/inputParser.h"

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

    vector<double> memorySizes = Cluster::buildHomogeneousMemorySizes(maxoutd, num_processors);
    //Fix, for now we consider the non-homog cluster homogeneuos
    Cluster::getFixedCluster()->setMemorySizes(memorySizes);
    delete temp_schedule;
}

int main(int argc, char **argv) {
    InputParser * input = new InputParser(argc, argv);
    string stage1, stage2 = "FirstFit", stage3;

    ifstream OpenFile(input->getPathToTreeList());

    
    list<Task *> parallelSubtrees;
    string treename;
    unsigned int num_processors;
    double makespan, maxoutd, minMem;
    clock_t time;
    unsigned int number_subtrees;
    float CCR=0;
    float NPR=0;
    switch (input->getClusteringMode())
    {
    case staticClustering:
    {
        float numberOfProcessors = input->getNumberOfProcessors();
        float processorMemory = input->getProcessorMemory();
        throw "Not implemented yet";
        break;
    }
    case treeDependent:
    {
        float CCR = input->getCCR();
        float NPR = input->getNPR();
        break;
    }
    default:
        break;
    }

    cout.precision(20);

    std::cout << "AmountSubtrees " << "AmountProcessors "
              << "Makespan " << "Stage1 " << "Stage2 " << "Stage3 " << "TimeConsuming" << std::endl;
    std::vector<int> brokenEdges;
    do {
        OpenFile >> treename;
        cout << treename << endl;
        Tree *tree = read_tree((input->getWorkingDirectory() + treename).c_str());
        Tree *untouchedTree = read_tree((input->getWorkingDirectory() + treename).c_str());
        Tree::setOriginalTree(untouchedTree);
        num_processors = ceil(tree->getSize() / NPR);
        if (num_processors < 3) {
            num_processors = 3;
        }

        buildFixedClusterWithSpeedsAndMemory(CCR, num_processors, tree);

        time = clock();

        unsigned long sequentialLen;
        makespan = tree->getRoot()->SplitSubtrees(false, parallelSubtrees,
                                                  sequentialLen);// for counting how many subtrees produced, twolevel is set as false
        // stage1 == "ImprovedSplit") {
        //po_construct(tree_size, prnts, &chstart, &chend, &children, &root);
        //   makespan = tree->ImprovedSplit();
        //makespan = ImprovedSplit(tree, num_processors);
        //stage1 == "AvoidChain") {
        //   makespan = tree->ASAP();
        //   number_subtrees = HowmanySubtrees(tree, true);
        //  time = clock() - time;
        //  std::cout << treename << " " << NPR << " " << CCR << " NA " << number_subtrees << " " << num_processors
        //             << " " << makespan << " " << "ASAP NA NA " << time << std::endl;

        //   time = clock();
        //number_subtrees = AvoidChain(tree);
        //  makespan = tree->GetRoot()->GetMSCost(true, true);
        // }
        time = clock() - time;
        cout << "1 step ready " << endl;
        number_subtrees = tree->HowmanySubtrees(false);

        std::cout << number_subtrees << " " << num_processors << " "
                  << makespan << " " << stage1 << " NA NA " << time << std::endl;


        maxoutd = MaxOutDegree(tree, true);
        schedule_t *schedule_f = new schedule_t();

        MinMem(tree, maxoutd, minMem, *schedule_f, true);
        delete schedule_f;

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
        std::cout << number_subtrees << " " << num_processors << " " << makespan << " " << stage1 << " " << stage2
                  << " NA " << time << std::endl;

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

        number_subtrees = tree->HowmanySubtrees(false);
        std::cout << number_subtrees << " " << num_processors << " " << makespan << " " << stage1 << " " << stage2
                  << " " << stage3 << " " << time << std::endl;

        delete tree;
        delete untouchedTree;
    } while (OpenFile.good());
    OpenFile.close();
    exit(EXIT_SUCCESS);
}