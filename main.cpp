// simulation sequence+stage2+splitAgain/Merge
//  Created by changjiang GOU on 10/05/2018.
//  Copyright © 2018 Changjiang GOU. All rights reserved.

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <algorithm>
#include <stdlib.h>
#include <lib-io-tree.h>
//

//#include "../include/cluster.h"

#include <inputParser.h>
#include "OutputPrinter.h"

#include <lib-io-tree-minmem.h>
#include "../include/heuristics.h"
#include "../include/inputParser.h"
#include "../include/OutputPrinter.h"
#include "../include/tree.h"
#include "../include/lib-io-tree-minmem.h"

void
buildTreeDependentCluster(InputParser *input, Tree *tree);

OutputPrinter printer;

double a2Steps(Tree *tree, OutputPrinter *printer) {
    unsigned int number_subtrees = 0;
    double makespan;
    clock_t time;
    time = clock();

    /* double maxoutd = MaxOutDegree(tree, true);
     double minMem;
     schedule_traversal *schedule_f = new schedule_traversal();
     MinMem(tree, maxoutd, minMem, *schedule_f, true);
     delete schedule_f;
     */
    tree->getRoot()->precomputeMinMems(tree);
    number_subtrees = tree->HowmanySubtrees(true);
    //makespan = tree->getRoot()->getMakespanCostWithSpeeds(true, true);
    cout << "1 step: " << clock() - time << " ";

    // int numberUnfeasibleTasks = 0;

    /* for (Task *task: *tree->getTasks()) {
         if (task->getFeasibleProcessors()->empty()) {
             numberUnfeasibleTasks++;
         }
     }*/

// cout << "#unfeasible tasks: " << numberUnfeasibleTasks << endl;

    time = clock();

    try {
        seqSetAndFeasSets(tree);

//        number_subtrees = tree->HowmanySubtrees(true);
        //      makespan = tree->getRoot()->getMakespanCostWithSpeeds(true, true);
        //    cout << "Makespan " << makespan << " #trees: " << number_subtrees << endl;

        //  cout<<"2 step: "<<clock()-time<<endl;
        //time = clock();

        makespan = tree->getRoot()->getMakespanCostWithSpeeds(true, true);
        number_subtrees = tree->HowmanySubtrees(true);
        cout << "2&3 step: " << clock() - time <<endl; //<<" #trees: " << number_subtrees << endl;
    }
    catch (exception e) {
        printer->quietPrint("An error has occurred: ");//+ e.what());  // Not executed
    }
    catch (const char *str) {
        printer->quietPrint(str);
        printer->quietPrint("No solution"); //<< str << endl;
        makespan = -1;
    }
    return makespan;

}

double threeSteps(Tree *tree, OutputPrinter *printer) {
    string stage2 = "FirstFit";
    unsigned int number_subtrees = 0;
    unsigned long sequentialLen;
    unsigned int num_processors = Cluster::getFixedCluster()->getNumberProcessors();
    list < Task * > parallelSubtrees;
    double makespan;
    // for counting how many subtrees produced, twolevel is set as false
    makespan = tree->ASAP();
    tree->getRoot()->SplitSubtrees(false, parallelSubtrees, sequentialLen, -1);
    // tree->ImprovedSplit();


    cout << "1 step ready " << endl;
    number_subtrees = tree->HowmanySubtrees(true);
    cout << "Makespan " << makespan << " #trees: " << number_subtrees << endl;

    int ret = 0;
    double processorMemorySize = Cluster::getFixedCluster()->getProcessors().at(0)->getMemorySize();
    if (stage2 == "LargestFirst") {
        if (Cluster::getFixedCluster()->isMemoryHomogeneous()) {
            printer->quietPrint("run homp memcheck");
            ret = MemoryCheckHomp(tree, LARGEST_FIT, processorMemorySize);
        } else {
            ret = MemoryCheck(tree, LARGEST_FIT, false);
        }
    } else if (stage2 == "FirstFit") {
        if (Cluster::getFixedCluster()->isMemoryHomogeneous()) {
            printer->quietPrint("run homp memcheck");
            ret = MemoryCheckHomp(tree, FIRST_FIT, processorMemorySize);
        } else {
            ret = MemoryCheck(tree, FIRST_FIT, false);
        }
    } else if (stage2 == "Immediately") {
        ret = Cluster::getFixedCluster()->isMemoryHomogeneous() ? MemoryCheckHomp(tree, IMMEDIATELY,
                                                                                  processorMemorySize) : MemoryCheck(
                tree, IMMEDIATELY, false);
    }
    if (ret == -1) {
        cout << "unsolvable currently" << endl;
        return -2;
    }
    cout << "2 step ready " << endl;
    makespan = tree->getRoot()->getMakespanCost(true, true);
    number_subtrees = tree->HowmanySubtrees(false);
    cout << "Makespan " << makespan << " #trees: " << number_subtrees << "weight " << endl;


    if (number_subtrees > num_processors) {
        cout << "merge" << endl;
        if (Cluster::getFixedCluster()->isMemoryHomogeneous()) {
            printer->quietPrint("run homp");
            // nsigned int num_subtrees, unsigned int processor_number, double const memory_size, bool CheckMemory)
            makespan = tree->MergeOld(number_subtrees, num_processors, processorMemorySize, true);
        } else
            makespan = tree->Merge(true);
    } else if (number_subtrees == num_processors) {
        cout << "nothing" << endl;
    } else {
        cout << "splitAgain" << endl;
        //  makespan = tree->SplitAgain();
        makespan = tree->SplitAgain();
        // SplitAgainOld(tree, num_processors, tree->HowmanySubtrees(true));
    }
    return makespan;
}

int main(int argc, char **argv) {
    InputParser *input = new InputParser(argc, argv);
    OutputPrinter *printer = new OutputPrinter;
    printer->setVerbose(input->getVerbosity());
    printer->initOutput();
    string treesToRerun = "";

    ifstream OpenFile(input->getPathToTreeList());
    ifstream OpenFilePreliminary(input->getPathToTreeList());

    string treename;

    double makespan, maxoutd, minMem;
    clock_t time;


    if (input->getClusteringMode() == staticClustering) {
        double maxMinMem = 0, maxMaxoutD = 0, sum_edges = 0, sum_weights = 0, maxEdgesToMakespanWeights = 0;
        do {

            OpenFilePreliminary >> treename;
            Tree *tree = read_tree((input->getWorkingDirectory() + treename).c_str());

            maxoutd = MaxOutDegree(tree, true);
            if (maxoutd > maxMaxoutD) maxMaxoutD = maxoutd;

            schedule_traversal *schedule_f = new schedule_traversal();
            MinMem(tree, maxoutd, minMem, *schedule_f, true);
            delete schedule_f;

            if (minMem > maxMinMem) {
                maxMinMem = minMem;
                for (Task *task: *tree->getTasks()) {
                    sum_edges += task->getEdgeWeight();
                    sum_weights += task->getMakespanWeight();
                }
                maxEdgesToMakespanWeights = sum_edges / sum_weights;
            }
            delete tree;

        } while (OpenFilePreliminary.good());
        OpenFilePreliminary.close();

        input->setClusterFromFile(1);
    }

    string header_column = "treename\t";
    do {
        header_column += input->getPathToCluster() + "\t";
    } while (input->nextCluster());

    input->resetClusterIterator();
    printer->quietPrint(header_column);
    double processorUtilizationOverall = 0;
    int numTrees = 0;
    std::vector<int> brokenEdges;
    do {
        OpenFile >> treename;
        //printer->quietPrint(treename);
        // cout << treename << endl;
        string extraSlash = input->getWorkingDirectory().back() != '/' ? "/" : "";
        Tree *tree = read_tree((input->getWorkingDirectory() + extraSlash +
                                treename).c_str());
        if (tree->getSize() == 0) {
            printer->quietPrint(
                    "Read empty tree with directory " + input->getWorkingDirectory() + " and file " + treename);
            continue;
        }
        Tree *untouchedTree = read_tree((input->getWorkingDirectory() + "/" + treename).c_str());
        Tree::setOriginalTree(untouchedTree);
        tree->getRoot()->breakEdge();
        numTrees++;
        //  const vector<double> fanouts = maxAndAvgFanout(tree);
        //   cout << treename << " Fanout: Max: " << fanouts[0] << ",  Avg: " << fanouts[1] <<
        //      " Max depth " << maxDepth(tree->getRoot()) << " num leaves " << tree->numberOfLeaves()
        //       << " #tasks: " << tree->getSize() << endl;



        string tree_column = treename + "\t";
        do {
            try {

                if (input->getClusteringMode() == treeDependent) {
                    buildTreeDependentCluster(input, tree);
                }
            }
            catch (const char *str) {
                printer->quietPrint(str); //<< str << endl
                continue;
            }

            time = clock();
            makespan = input->getRunA1() ? threeSteps(tree, printer) : a2Steps(tree, printer);
            //maxoutd = MaxOutDegree(tree, true);

            //schedule_traversal *schedule_f = new schedule_traversal();
            //MinMem(tree, maxoutd, minMem, *schedule_f, true);
            //delete schedule_f;

            //makespan = (maxoutd/(minMem+maxoutd))*100;
            time = clock() - time;

            if (makespan == -1) {
                cout << "no solution" << endl;
            }
            // cout<<"makespan "<<makespan<<endl;
            int num_subtrees = tree->HowmanySubtrees(true);

            //  makespan = tree->getRoot()->getMakespanCost(true, true);
            // cout<<"makespan "<<makespan<<endl;
            //quietPrint("&& " + treename + " " + to_string(makespan) + " " + to_string(time));
            // quietPrint(Cluster::getFixedCluster()->getPrettyClusterString());
            // quietPrint(Cluster::getFixedCluster()->getAverageLoadAndNumberOfUsedProcessors());
            //quietPrint(Cluster::getFixedCluster()->getUsageString());
            //  quietPrint(Cluster::getFixedCluster()->printProcessors());
            tree_column += to_string(makespan) + "\t" + to_string(num_subtrees) + "\t";
            double processorUtilization = 0;
            for (Processor *proc: (Cluster::getFixedCluster()->getProcessors())) {
                if (proc->isBusy) {
                    processorUtilization++;
                }
            }
            processorUtilization /= Cluster::getFixedCluster()->getProcessors().size();
            processorUtilizationOverall += processorUtilization;
            tree->clearComputedValues();
        } while (input->nextCluster());
        printer->quietPrint(tree_column);
        delete tree;
        delete untouchedTree;
        Cluster::getFixedCluster()->clean();
        input->resetClusterIterator();
    } while (OpenFile.good());
    OpenFile.close();
    processorUtilizationOverall /= numTrees;
    //quietPrint("proc util "+ to_string(processorUtilizationOverall));
    cout << treesToRerun << endl;
    exit(EXIT_SUCCESS);

}

void
buildTreeDependentCluster(InputParser *input, Tree *tree) {
    double maxoutd, minMem;
    maxoutd = MaxOutDegree(tree, true);

    schedule_traversal *schedule_f = new schedule_traversal();
    MinMem(tree, maxoutd, minMem, *schedule_f, true);

    input->setClusterFromFile(maxoutd);

    if (minMem > Cluster::getFixedCluster()->getCumulativeMemory()) {
        throw "Cluster too small: cumulative memory: " + to_string(Cluster::getFixedCluster()->getCumulativeMemory()) +
              " vs required " +
              to_string(minMem);
    }
}
