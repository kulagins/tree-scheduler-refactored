//
//  heuristics.h
//  memCom2
//
//  Created by changjiang GOU on 10/05/2018.
//  Copyright © 2018 ROMA. All rights reserved.
//

#ifndef heuristics_h
#define heuristics_h

#include "lib-io-tree.h"
//#include "lib-io-tree-minmem.h"
//#include "cluster.h"

template<class T, class U>
void GetTwoLargestElementTypetwo(T container, U &Largest, U &secondLargest);

static double timeForAssignment = 0;

int MemoryCheck(Tree *tree, io_method_t method, bool useMinimalAvailableProvcessor);

int MemoryCheckHomp(Tree *tree, io_method_t method, double processor_memory_size);

//void
//MemoryCheckA2(Tree *tree, Cluster *cluster, io_method_t method, bool skipBig);

void SetBandwidth(double CCR, unsigned long tree_size, double *ewghts, double *timewghts);

void distributeProcessors(Tree *qTree);

string seqSetAndFeasSets(Tree *tree);

double assignToBestProcessors(Tree *tree, string assignSubtreeChoiceCode = "N");

void removeProcessorFromAllFeasSets(Processor *processor, Tree *tree);

string partitionHeuristics(Tree *tree, string subtreeChoiceCode, string nodeChoiceCode, string assignSubtreeChoiceCode);

Task *chooseSubtree(string subtreeChoiceCode, Tree * tree, vector<Task*> candidates);

Task *chooseNode(Task *root, Tree *tree, string nodeChoiceCode, string assignSubtreeChoiceCode);

void chooseAssignSubtree(string assignSubtreeChoiceCode, Tree * tree);

Task *findBestCutAmong(Tree *tree, vector<Task *> candidates, string assignSubtreeChoiceCode, double initMS = -1);

Task * CutTaskWithMaxImprovement(Tree *tree, string assignSubtreeChoiceCode);

void SiftInfTmaxUpPreserveOrder(vector<Task *> *taskHeap);

vector<Task *> criticalPath(Tree *tree);
vector<Task *> buildCriticalPath(Task* root);

#endif /* heuristics_h */
