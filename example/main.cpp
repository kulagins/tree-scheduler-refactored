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

const bool verbose = true;

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

void initOutput() {
    if (!verbose) {
        cout.setstate(std::ios_base::failbit);
    }
}

void quietPrint(string text) {
    cout.clear();
    cout << text << endl;
    initOutput();
}

int main(int argc, char **argv) {
    InputParser *input = new InputParser(argc, argv);

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
            clusterConfigurationNumber = (int) input->getStaticClusterConfigurationNumber();
            break;
        }
        default:
            break;
    }

    if (input->getClusteringMode() == staticClustering) {
        double maxMinMem = 0;
        double maxMaxoutD = 0;
        double maxEdgesToMakespanWeights = 0;
        double sum_edges = 0;
        double sum_weights = 0;

        do {
            schedule_traversal *schedule_f = new schedule_traversal();
            OpenFilePreliminary >> treename;
            Tree *tree = read_tree((input->getWorkingDirectory() + treename).c_str());
            maxoutd = MaxOutDegree(tree, true);
            if (maxoutd > maxMaxoutD) maxMaxoutD = maxoutd;
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
                    Cluster::buildHomStatic2LevelCluster(maxMaxoutD, maxEdgesToMakespanWeights,
                                                         input->getAdaptationMode() ? input->getAdaptationMode()
                                                                                    : noAdaptation);
                } else {
                    Cluster::buildStatic2LevelCluster(maxMaxoutD, maxEdgesToMakespanWeights);
                }
                break;
            case 2:
                if (input->getHeterogenityLevel() == homogeneus) {
                    Cluster::buildHomStatic3LevelCluster(maxMaxoutD, maxEdgesToMakespanWeights,
                                                         input->getAdaptationMode() ? input->getAdaptationMode()
                                                                                    : noAdaptation);
                } else {
                    Cluster::buildStatic3LevelCluster(maxMaxoutD, maxEdgesToMakespanWeights);
                }
                break;
            default:
                throw "No such cluster configuration is implemented: " + to_string(clusterConfigurationNumber);
        }
        Cluster::getFixedCluster()->printInfo();
        num_processors = Cluster::getFixedCluster()->getProcessors().size();
    }

    std::vector<int> brokenEdges;
    do {
        OpenFile >> treename;
        Tree *tree = read_tree((input->getWorkingDirectory() + treename).c_str());
        Tree *untouchedTree = read_tree((input->getWorkingDirectory() + treename).c_str());
        Tree::setOriginalTree(untouchedTree);
        maxoutd = MaxOutDegree(tree, true);
        quietPrint("maxoutD" + to_string(maxoutd));
        /*    schedule_traversal *schedule_f = new schedule_traversal();

             MinMem(tree, maxoutd, minMem, *schedule_f, true);
             cout<< minMem<< " "<< maxoutd<< " "<< minMem/maxoutd<< " "<< maxoutd*100/minMem<<"%"<<endl;
   */
        if (input->getClusteringMode() == treeDependent) {
            switch (clusterConfigurationNumber) {
                case 1:
                    if (input->getHeterogenityLevel() == homogeneus) {
                        Cluster::buildHomStatic2LevelCluster(maxoutd, 610,
                                                             input->getAdaptationMode() ? input->getAdaptationMode()
                                                                                        : noAdaptation);
                    } else {
                        Cluster::buildStatic2LevelCluster(maxoutd, 610);
                    }
                    break;
                case 2:
                    if (input->getHeterogenityLevel() == homogeneus) {
                        Cluster::buildHomStatic3LevelCluster(maxoutd, 610,
                                                             input->getAdaptationMode() ? input->getAdaptationMode()
                                                                                        : noAdaptation);
                    } else {
                        Cluster::buildStatic3LevelCluster(maxoutd, 610);
                    }
                    break;
                default:
                    throw "No such cluster configuration is implemented: " + to_string(clusterConfigurationNumber);
            }
            Cluster::getFixedCluster()->printInfo();
            num_processors = Cluster::getFixedCluster()->getProcessors().size();

        }

        time = clock();
        unsigned long sequentialLen;
        // for counting how many subtrees produced, twolevel is set as false
        makespan = tree->getRoot()->SplitSubtrees(false, parallelSubtrees, sequentialLen, -1);
        // tree->ImprovedSplit();

        time = clock() - time;
        cout << "1 step ready " << endl;
        number_subtrees = tree->HowmanySubtrees(false);
        cout << "Makespan " << makespan << " #trees: " << number_subtrees << endl;


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
        cout << "2 step ready " << endl;
        makespan = tree->getRoot()->getMakespanCost(true, true);
        number_subtrees = tree->HowmanySubtrees(false);
        cout << "Makespan " << makespan << " #trees: " << number_subtrees << endl;


        if (number_subtrees > num_processors) {
            stage3 = "Merge";
            makespan = tree->Merge(true);
        } else if (number_subtrees == num_processors) {
            stage3 = "Nothing";
        } else {
            stage3 = "SplitAgain";
            makespan = tree->SplitAgain();
        }
        time = clock() - time;

        quietPrint(treename + " " + to_string(makespan) + " " + to_string(time));
        quietPrint(Cluster::getFixedCluster()->getPrettyClusterString());
        quietPrint(Cluster::getFixedCluster()->getAverageLoadAndNumberOfUsedProcessors());
        //quietPrint(Cluster::getFixedCluster()->getUsageString());
        delete tree;
        delete untouchedTree;
        Cluster::getFixedCluster()->clean();
    } while (OpenFile.good());
    OpenFile.close();
    exit(EXIT_SUCCESS);
}