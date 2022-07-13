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

OutputPrinter printer;

string a2MultiLevel(Tree *tree, OutputPrinter *printer, double &makespan, InputParser *pParser) {
    unsigned int number_subtrees = 0;

    clock_t time;
    time = clock();

    tree->getRoot()->precomputeMinMems(tree);
    string result = /*"0 step: " +*/ to_string(clock() - time) + " ";


    time = clock();

    try {
        result += partitionHeuristics(tree, pParser->getChooseSubtree(), pParser->getChooseNode(),
                                      pParser->getAssignChooseSubtree());
        // double makespan1 = tree->getRoot()->getMakespanCostWithSpeeds(true, true);
        // assert(makespan==makespan1);

        number_subtrees = tree->HowmanySubtrees(true);
        result +=/* " 2&3 step: " + */
               // to_string(clock() - time) + " "/*+" #trees: "*/+
            " "+   to_string(number_subtrees) + "\n";
    }
    catch (exception e) {
        printer->quietPrint("An error has occurred: ");
        printer->quietPrint(e.what());
    }
    catch (const char *str) {
        printer->quietPrint(str);
        printer->quietPrint("No solution"); //<< str << endl;
        makespan = -1;
    }
    return result;

}

string a2Steps(Tree *tree, OutputPrinter *printer, double &makespan, InputParser *pParser) {
    unsigned int number_subtrees = 0;

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
    string result = "1 step: " + to_string(clock() - time) + " ";


    time = clock();

    try {
        result += seqSetAndFeasSets(tree);
        makespan = tree->getRoot()->getMakespanCostWithSpeeds(true, true);
        number_subtrees = tree->HowmanySubtrees(true);
        result += " 2&3 step: " + to_string(clock() - time) + " "; //<<" #trees: " << number_subtrees << endl;
    }
    catch (exception e) {
        printer->quietPrint("An error has occurred: ");//+ e.what());  // Not executed
    }
    catch (const char *str) {
        printer->quietPrint(str);
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
            //  makespan = input->getRunA1() ? threeSteps(tree, printer) : a2Steps(tree, printer);
            string result = a2MultiLevel(tree, printer, makespan, input);
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

            makespan = tree->getRoot()->getMakespanCostWithSpeeds(true, true);
            tree_column += " "+to_string(makespan) + "\t" + to_string(tree->HowmanySubtrees(true)) + "\t" + result;
            /*for (Processor *proc: (Cluster::getFixedCluster()->getProcessors())) {
                if (proc->isBusy) {
                    cout<<proc->getMemorySize()<<" "<<proc->getAssignedTaskId()<<endl;
                }
                else cout<<proc->getMemorySize()<<" free"<<endl;
            }*/

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
