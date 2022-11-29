// simulation sequence+stage2+splitAgain/Merge
//  Created by changjiang GOU on 10/05/2018.
//  Copyright © 2018 Changjiang GOU. All rights reserved.

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <algorithm>
#include <stdlib.h>


#include "../include/heuristics.h"
#include "include/inputParser.h"
#include "include/tree.h"
#include "include/lib-io-tree-minmem.h"
#include "include/OutputPrinter.h"


/*
 * ### How to run

#### ./main trees-directory trees-list clusters-list clustering-mode run-a1 verbosity choose-subtree choose-task choose-subtree-assign

The executable requires the following input parameters:

- `clustering-mode` : Static clustering (memories are fixed) vs tree-dependent (memory given to each tree is individual)
  - `clustering-mode = 0`: the clustering-mode is dependent on the tree instance.
  - `clustering-mode = 1`: the clustering-mode is static.
- `run-a1`:
  - `run-a1 = 0`: do not run A1 code (run A2)
  - `run-a1 = 1`: run A1 code
- `verbosity`: Should all of the debugging output be printed when the program runs
- `choose-subtree`: how we choose multi-level subtree to cut
  - `choose-subtree = LMW`: Largest M<sub>i</sub> * W<sub>i</sub> of the subtree
  - `choose-subtree = CP`: First (unprocessed) subtree on the critical path
- `choose-task`: how we choose the task *within the chosen subtree* to cut his children
  - `choose-task = EX`: Exhaustive, try cutting all tasks and shoose the one with best makespan
  - `choose-task = FFT`: First from top that gives an improvement
  - `choose-task = M`: Task on the middle level that gives the best makespan
- `choose-subtree-assign`: how we choose the subtree to assign best processors
  - `choose-subtree-assign = LW`: Largest W<sub>i</sub> of the subtree
  - `choose-subtree-assign = MD`: Maximum number of descendants
  - `choose-subtree-assign = CP`: First (unprocessed) subtree on the critical path
 *
 */

void
buildTreeDependentCluster(InputParser *input, Tree *tree);

void buildStaticCluster(InputParser *input, ifstream &OpenFilePreliminary, string &treename,
                        double maxoutd, double minMem);

double computeProcesorUtilization(double processorUtilizationOverall);

OutputPrinter printer;

string a2WithNewMinMem(Tree *tree, OutputPrinter *printer, double &makespan, InputParser *pParser) {

    string result = "";
    clock_t time;
    time = clock();

    tree->mergeLinearChains();
    result += "reprocessing: " + to_string((clock() - time) / CLOCKS_PER_SEC);

    try {
        result += partitionHeuristics(tree, pParser->getChooseSubtree(), pParser->getChooseNode(),
                                      pParser->getAssignChooseSubtree(), pParser->getCutWhat());

    }
    catch (const char *str) {
        printer->quietPrint(str);
        printer->quietPrint("No solution"); //<< str << endl;
        makespan = -1;
    }
    catch (...) {
        printer->quietPrint("No solution"); //<< str << endl;
        makespan = -1;
    }
    return result;

}

double threeSteps(Tree *tree, OutputPrinter *printer) {
    string stage2 = "FirstFit";
    unsigned int number_subtrees = 0;
    unsigned long sequentialLen;
    unsigned int num_processors = Cluster::getFixedCluster()->getNumberProcessors();
    list<Task *> parallelSubtrees;
    double makespan;
    // for counting how many subtrees produced, twolevel is set as false
    makespan = tree->ASAP();
    tree->getRoot()->SplitSubtrees(false, parallelSubtrees, sequentialLen, -1);
    // tree->ImprovedSplit();


    // cout << "1 step ready " << endl;
    number_subtrees = tree->HowmanySubtrees(true);
    // cout << "Makespan " << makespan << " #trees: " << number_subtrees << endl;

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
    //  cout << "2 step ready " << endl;
    makespan = tree->getRoot()->getMakespanCost(true, true);
    number_subtrees = tree->HowmanySubtrees(true);
    //  cout << "Makespan " << makespan << " #trees: " << number_subtrees << "weight " << endl;


    if (number_subtrees > num_processors) {
        //cout << "merge" << endl;
        if (Cluster::getFixedCluster()->isMemoryHomogeneous()) {
            printer->quietPrint("run homp");
            // nsigned int num_subtrees, unsigned int processor_number, double const memory_size, bool CheckMemory)
            makespan = tree->MergeOld(number_subtrees, num_processors, processorMemorySize, true);
        } else
            makespan = tree->Merge(true);
    } else if (number_subtrees == num_processors) {
        //cout << "nothing" << endl;
    } else {
        // cout << "splitAgain" << endl;
        //  makespan = tree->SplitAgain();
        makespan = tree->SplitAgain();
        // SplitAgainOld(tree, num_processors, tree->HowmanySubtrees(true));
    }

    for (const auto &item: Cluster::getFixedCluster()->getProcessors()) {
        if (item->getAssignedTaskId() != -1 && !item->getAssignedTask()->isBroken()) {
            item->isBusy = false;
            item->setAssignedTaskId(-1);
            item->setAssignedTask(NULL);
            item->setOccupiedMemorySize(0);
        }
    }

    for (const auto &item: tree->getBrokenTasks()) {
        if (item->getAssignedProcessor() == NULL) {

            double maxoutd, minMem;
            Tree *subtree = BuildSubtree(tree, item);
            maxoutd = MaxOutDegree(subtree, true);
            schedule_traversal *schedule_f = new schedule_traversal();
            MinMem(subtree, maxoutd, minMem, *schedule_f, true);

            delete schedule_f;
            delete subtree;
            for (const auto &proc: Cluster::getFixedCluster()->getProcessors()) {
                if (!proc->isBusy && proc->getMemorySize() >= minMem) {
                    proc->assignTask(item);
                    break;
                }
            }
        }
    }
    for (const auto &item: tree->getBrokenTasks()) {
        if (item->getAssignedProcessor() == NULL) {
            Cluster::getFixedCluster()->getBiggestFreeProcessor()->assignTask(item);
            cout << "assigning afterwards! " << item->getId() << " "
                 << Cluster::getFixedCluster()->getNumberFreeProcessors() << endl;
        }
        assert(item->getAssignedProcessor() != NULL);
    }
    for (Task *brokenTask: tree->getBrokenTasks()) {
        vector<Task *> allTasksInSubtree = brokenTask->getTasksInSubtreeRootedHere();
        for (Task *taskInSubtree: allTasksInSubtree) {
            taskInSubtree->setAssignedProcessor(brokenTask->getAssignedProcessor());
        }
    }
    for (const auto &item: *tree->getTasks()) {
        assert(item->getAssignedProcessor() != NULL);
    }
    return tree->getRoot()->getMakespanCost(true, true);
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
    clock_t time, time1;


    if (input->getClusteringMode() == staticClustering) {
        buildStaticCluster(input, OpenFilePreliminary, treename, maxoutd, minMem);
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
    while (OpenFile.good()) {
        OpenFile >> treename;
        if (OpenFile.eof()) break;
        //printer->quietPrint(treename);
        //cout << treename << endl;
        string extraSlash = input->getWorkingDirectory().back() != '/' ? "/" : "";
        Tree *tree = read_tree((input->getWorkingDirectory() + extraSlash +
                                treename).c_str());
        if (tree->getSize() == 0) {
            printer->quietPrint(
                    "Read empty tree with directory " + input->getWorkingDirectory() + " and file " + treename);
            continue;
        }
        Tree *untouchedTree = read_tree((input->getWorkingDirectory() + "/" + treename).c_str());
        veryOriginalTree = untouchedTree;
        tree->getRoot()->breakEdge();
        numTrees++;

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

            makespan = threeSteps(tree, printer);

            Tree *qtree = tree->BuildQtree();
            assignAllCorrespondingTreeTasks(tree, qtree);
            makespan = tree->getRoot()->getMakespanCostWithSpeeds(true, true);
            time = (clock() - time) / CLOCKS_PER_SEC;
            time1 = clock();
            double makespan1 = swapUntilBest(tree);
            if (makespan1 == numeric_limits<double>::infinity()) {
                printer->quietPrint("No MS");
                makespan1 = makespan;
            }

            string result = to_string(makespan1);

            time1 = (clock() - time1)/CLOCKS_PER_SEC;

            if (makespan == -1) {
                cout << "no solution" << endl;
            }


            //tree->HowmanySubtreesAndWeights(false);
            tree_column += " " + to_string(makespan) + "\t" + to_string(tree->HowmanySubtrees(true)) + "\t" +
                           // to_string(time )+ " " + to_string(CLOCKS_PER_SEC);
                           result + " " + to_string(time ) + " "+ to_string(time1);
            /*for (Processor *proc: (Cluster::getFixedCluster()->getProcessors())) {
                if (proc->isBusy) {
                    cout<<proc->getMemorySize()<<" "<<proc->getAssignedTaskId()<<endl;
                }
                else cout<<proc->getMemorySize()<<" free"<<endl;
            }*/

            processorUtilizationOverall = computeProcesorUtilization(processorUtilizationOverall);
            tree->clearComputedValues();
        } while (input->nextCluster());
        printer->quietPrint(tree_column);
        delete tree;
        delete veryOriginalTree;
        Cluster::getFixedCluster()->clean();
        input->resetClusterIterator();
    };
    OpenFile.close();
    processorUtilizationOverall /= numTrees;
    //quietPrint("proc util "+ to_string(processorUtilizationOverall));
    cout << treesToRerun << endl;
    exit(EXIT_SUCCESS);

}

double computeProcesorUtilization(double processorUtilizationOverall) {
    double processorUtilization = 0;
    for (Processor *proc: (Cluster::getFixedCluster()->getProcessors())) {
        if (proc->isBusy) {
            processorUtilization++;
        }
    }
    processorUtilization /= Cluster::getFixedCluster()->getProcessors().size();
    processorUtilizationOverall += processorUtilization;
    return processorUtilizationOverall;
}

void buildStaticCluster(InputParser *input, ifstream &OpenFilePreliminary, string &treename,
                        double maxoutd, double minMem) {
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
