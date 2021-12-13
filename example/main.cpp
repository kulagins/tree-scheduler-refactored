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

const bool verbose = false;

void buildFixedClusterWithSpeedsAndMemory(double CCR, unsigned int num_processors, Tree *treeobj) {
    Cluster *cluster = new Cluster(num_processors, true);
    double minMem;
    map<int, int> processor_speeds = Cluster::buildProcessorSpeeds(num_processors);
    cluster->SetBandwidth(CCR, treeobj);
    Cluster::setFixedCluster(cluster);

    Cluster::getFixedCluster()->SetBandwidth(CCR, treeobj);
    double maxoutd = MaxOutDegree(treeobj, true);
    schedule_traversal *temp_schedule = new schedule_traversal();
    MinMem(treeobj, maxoutd, minMem, *temp_schedule, true);

    vector<double> memorySizes = Cluster::buildHomogeneousMemorySizes(maxoutd, num_processors);
    //Fix, for now we consider the non-homog cluster homogeneuos
    Cluster::getFixedCluster()->setMemorySizes(memorySizes);
    delete temp_schedule;
}
void initOutput(){
    if (!verbose){
        cout.setstate(std::ios_base::failbit);
    }
}

void quietPrint(string text){
    cout.clear();
    cout << text << endl;
    initOutput();
}

int main(int argc, char **argv) {
    initOutput();
    string stage1, stage2 = "FirstFit", stage3;

    ifstream OpenFile(input->getPathToTreeList());
    ifstream OpenFilePreliminary(input->getPathToTreeList());
    list<Task *> parallelSubtrees;
    string treename;
    unsigned int num_processors;
    double makespan, maxoutd, minMem;
    clock_t time;
    unsigned int number_subtrees;
    float CCR = 0;
    float NPR = 0;
    int clusterConfigurationNumber = 0;

    switch (input->getClusteringMode()) {
        case staticClustering: {
            clusterConfigurationNumber = (int) input->getStaticClusterConfigurationNumber();
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

    if (input->getClusteringMode() == staticClustering) {
        double maxMinMem = 0;
        double maxEdgesToMakespanWeights = 0;
        double sum_edges = 0;
        double sum_weights = 0;

        do {
            schedule_traversal *schedule_f = new schedule_traversal();
            OpenFilePreliminary >> treename;
            Tree *tree = read_tree((input->getWorkingDirectory() + treename).c_str());
            maxoutd = MaxOutDegree(tree, true);
            MinMem(tree, maxoutd, minMem, *schedule_f, true);
            if (minMem > maxMinMem) {
                maxMinMem = minMem;
                for (Task *task: *tree->getTasks()) {
                    sum_edges += task->getEdgeWeight();
                    sum_weights += task->getMakespanWeight();
                }
                maxEdgesToMakespanWeights = sum_edges / sum_weights;
            }
            delete tree;
            delete schedule_f;
        } while (OpenFilePreliminary.good());
        OpenFilePreliminary.close();
        switch (clusterConfigurationNumber) {
            case 1:
                if (input->getHeterogenityLevel() == homogeneus) {
                    Cluster::buildHomStatic2LevelCluster(maxMinMem, maxEdgesToMakespanWeights,
                                                         input->getAdaptationMode() ? input->getAdaptationMode()
                                                                                    : noAdaptation);
                } else {
                    Cluster::buildStatic2LevelCluster(maxMinMem, maxEdgesToMakespanWeights);
                }
                break;
            case 2:
                if (input->getHeterogenityLevel() == homogeneus) {
                    Cluster::buildHomStatic3LevelCluster(maxMinMem, maxEdgesToMakespanWeights,
                                                         input->getAdaptationMode() ? input->getAdaptationMode()
                                                                                    : noAdaptation);
                } else {
                    Cluster::buildStatic3LevelCluster(maxMinMem, maxEdgesToMakespanWeights);
                }
                break;
            default:
                throw "No such cluster configuration is implemented: " + to_string(clusterConfigurationNumber);
        }
        Cluster::getFixedCluster()->printInfo();
    }

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
                Cluster::buildHomogeneousCluster(CCR, num_processors, tree,
                                                 input->getAdaptationMode() ? input->getAdaptationMode()
                                                                            : noAdaptation);
            } else {
                Cluster::buildMemHetTreeDepCluster(CCR, num_processors, tree);
            }
            //   Cluster::getFixedCluster()->printProcessors();
        }

        quietPrint(Cluster::getFixedCluster()->getPrettyClusterString());

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
        
        quietPrint(treename+" "+to_string(makespan)+" "+to_string(time)+"\n\n");
        //  Cluster::getFixedCluster()->printProcessors();
        // number_subtrees = tree->HowmanySubtrees(false);
        //  std::cout << number_subtrees << " " << num_processors << " " << makespan << " " << stage1 << " " << stage2
        //           << " " << stage3 << " " << time << std::endl;

        delete tree;
        delete untouchedTree;
    } while (OpenFile.good());
    OpenFile.close();
    exit(EXIT_SUCCESS);
}