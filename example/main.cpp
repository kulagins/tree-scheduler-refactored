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

void buildHomogeneousBandwidths(double CCR, unsigned int num_processors, Tree *treeobj, double &minMem, double &maxoutd,
                                schedule_traversal *&temp_schedule);


void buildHomogeneousCluster(double CCR, unsigned int num_processors, Tree *treeobj, HeterogeneousAdaptationMode mode) {
    // mode 0: build a cluster that uses all nodes with smallest memory
    // mode 1: build a cluster that uses 2/3 of nodes with middle amount of memory
    // mode 2: build a cluster that uses 1/3 of nodes with big memory

    double minMem;
    double maxoutd;
    schedule_traversal *temp_schedule;
    buildHomogeneousBandwidths(CCR, num_processors, treeobj, minMem, maxoutd, temp_schedule);
    vector<double> memorySizes;
    switch (mode) {
        case noAdaptation:
            memorySizes = Cluster::buildHomogeneousMemorySizes(maxoutd, num_processors);
            break;
        case manySmall:
            memorySizes = Cluster::buildHomogeneousMemorySizes(min(maxoutd, minMem), num_processors);
            break;
        case average:
            memorySizes = Cluster::buildHomogeneousMemorySizes((maxoutd + minMem) / 2, num_processors * 2 / 3);
            break;
        case fewBig:
            memorySizes = Cluster::buildHomogeneousMemorySizes(max(maxoutd, minMem), num_processors / 3);
            break;
        default:
            throw std::runtime_error("Bad mode for creating homogeneous cluster.");

    }

    //Fix, for now we consider the non-homog cluster homogeneuos
    Cluster::getFixedCluster()->setMemorySizes(memorySizes);
    delete temp_schedule;
}

void buildMemHeterogeneousCluster(double CCR, unsigned int num_processors, Tree *treeobj) {
    double minMem;
    double maxoutd;
    schedule_traversal *temp_schedule;
    buildHomogeneousBandwidths(CCR, num_processors, treeobj, minMem, maxoutd, temp_schedule);

    vector<double> memorySizes = Cluster::build3LevelMemorySizes(min(maxoutd, minMem), max(maxoutd, minMem),
                                                                 num_processors);
    Cluster::getFixedCluster()->setMemorySizes(memorySizes);
    delete temp_schedule;
}


void buildHomogeneousBandwidths(double CCR, unsigned int num_processors, Tree *treeobj, double &minMem, double &maxoutd,
                                schedule_traversal *&temp_schedule) {
    maxoutd = MaxOutDegree(treeobj, true);
    temp_schedule = new schedule_traversal();
    Cluster *cluster = new Cluster(num_processors, true);
    map<int, int> processor_speeds = Cluster::buildProcessorSpeeds(num_processors);
    cluster->SetBandwidth(CCR, treeobj);
    Cluster::setFixedCluster(cluster);

    Cluster::getFixedCluster()->SetBandwidth(CCR, treeobj);
    MinMem(treeobj, maxoutd, minMem, *temp_schedule, true);
}

int main(int argc, char **argv) {
    InputParser *input = new InputParser(argc, argv);
    string stage1, stage2 = "FirstFit", stage3;

    ifstream OpenFile(input->getPathToTreeList());


    list<Task *> parallelSubtrees;
    string treename;
    unsigned int num_processors;
    double makespan, maxoutd, minMem;
    clock_t time;
    unsigned int number_subtrees;
    float CCR = 0;
    float NPR = 0;
    float numberOfProcessors = 0;
    float processorMemory = 0;
    switch (input->getClusteringMode()) {
        case staticClustering: {
            numberOfProcessors = input->getNumberOfProcessors();
            processorMemory = input->getProcessorMemory();
            throw "Not implemented yet";
            break;
        }
        case treeDependent: {
            CCR = input->getCCR();
            NPR = input->getNPR();
            break;
        }
        default:
            break;
    }

    //  cout.precision(20);

    //  std::cout << "AmountSubtrees " << "AmountProcessors "
    //       << "Makespan " << "Stage1 " << "Stage2 " << "Stage3 " << "TimeConsuming" << std::endl;
    std::vector<int> brokenEdges;
    do {
        OpenFile >> treename;
        cout << treename << "\t";
        Tree *tree = read_tree((input->getWorkingDirectory() + treename).c_str());
        Tree *untouchedTree = read_tree((input->getWorkingDirectory() + treename).c_str());
        Tree::setOriginalTree(untouchedTree);


        if (input->getClusteringMode() == treeDependent) {
            num_processors = ceil(tree->getSize() / NPR);
            if (num_processors < 3) {
                num_processors = 3;
            }

            if (input->getHeterogenityLevel() == homogeneus) {
                buildHomogeneousCluster(CCR, num_processors, tree,
                                        input->getAdaptationMode() ? input->getAdaptationMode() : noAdaptation);
            } else {
                buildMemHeterogeneousCluster(CCR, num_processors, tree);
            }
            //   Cluster::getFixedCluster()->printProcessors();
        } else {
            throw std::runtime_error("No static cluster sizes yet.");
        }


        time = clock();
        /// tree->Print(cout);

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
        //  cout << "1 step ready " << endl;
        //  number_subtrees = tree->HowmanySubtrees(false);
        // cout << makespan << endl;
        //    std::cout << number_subtrees << " " << num_processors << " "
        //         << makespan << " " << stage1 << " NA NA " << time << std::endl;


        maxoutd = MaxOutDegree(tree, true);
        schedule_traversal *schedule_f = new schedule_traversal();

        MinMem(tree, maxoutd, minMem, *schedule_f, true);
        delete schedule_f;

        time = clock();
        int ret = 0;
        if (stage2 == "LargestFirst") {
            ret = MemoryCheck(tree, LARGEST_FIT);
        } else if (stage2 == "FirstFit") {
            ret = MemoryCheck(tree, FIRST_FIT);
        } else if (stage2 == "Immediately") {
            ret = MemoryCheck(tree, IMMEDIATELY);
        }
        if (ret == -1) {
            cout << "unsolvable currently" << endl;
            continue;
        }
        //  cout << "2 step ready " << endl;
        //number_subtrees = tree->HowmanySubtrees(false);

        time = clock() - time;
        number_subtrees = tree->HowmanySubtrees(true);
        makespan = tree->getRoot()->getMakespanCost(true, true);
        //   std::cout << number_subtrees << " " << num_processors << " " << makespan << " " << stage1 << " " << stage2
        //             << " NA " << time << std::endl;
        // cout << makespan << endl;
        time = clock();
        if (number_subtrees > num_processors) {
            stage3 = "Merge";
            makespan = tree->Merge(true);
              //  makespan = tree->MergeV2(number_subtrees, num_processors,
                //                       Cluster::getFixedCluster()->getFirstFreeProcessorOrSmallest()->getMemorySize(),
                  //                   true);
            //Merge(tree, number_subtrees, num_processors, memorySize, chstart, children, true);
        } else if (number_subtrees == num_processors) {
            stage3 = "Nothing";
        } else {
            stage3 = "SplitAgain";
            makespan = tree->SplitAgain();
        }
        time = clock() - time;
        //  Cluster::getFixedCluster()->printProcessors();
        // number_subtrees = tree->HowmanySubtrees(false);
        cout << makespan << endl;
        //  std::cout << number_subtrees << " " << num_processors << " " << makespan << " " << stage1 << " " << stage2
        //           << " " << stage3 << " " << time << std::endl;

        delete tree;
        delete untouchedTree;
    } while (OpenFile.good());
    OpenFile.close();
    exit(EXIT_SUCCESS);
}