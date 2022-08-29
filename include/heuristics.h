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

template<class T, class U>

 extern double timeForAssignment;
extern double timeChooseNode;
extern double timeChooseTree;
extern double timeBestCutInNodeChoice;

bool cmp_noincreasing(Task *a, Task *b);

bool cmp_nodecreasing(Task *a, Task *b);

bool cmp_noIn_noCommu(Task *a, Task *b);

bool cmp_Mem_nodecreasing(Task *a, Task *b);

bool cmp_asap(Task *a, Task *b);


int MemoryCheck(Tree *tree, io_method_t method, bool useMinimalAvailableProvcessor);

int MemoryCheckHomp(Tree *tree, io_method_t method, double processor_memory_size);

//void
//MemoryCheckA2(Tree *tree, Cluster *cluster, io_method_t method, bool skipBig);

void distributeProcessors(Tree *qTree);

string seqSetAndFeasSets(Tree *tree);

double assignToBestProcessors(Tree *tree, vector<Task *> newlyBroken, string assignSubtreeChoiceCode = "N");

void removeProcessorFromAllFeasSets(Processor *processor, Tree *tree);

string partitionHeuristics(Tree *tree, string subtreeChoiceCode, string nodeChoiceCode, string assignSubtreeChoiceCode);

string
partitionHeuristicsNoPreprocessing(Tree *tree, string subtreeChoiceCode, string nodeChoiceCode,
                                   string assignSubtreeChoiceCode);

Task *chooseSubtree(string subtreeChoiceCode, Tree *tree, vector<Task *> candidates);

Task *chooseTask(Task *root, Tree *tree, string nodeChoiceCode, string assignSubtreeChoiceCode);

vector<Task *> buildCandidatesForNode(Tree *tree, const string &nodeChoiceCode, Task *root);

void chooseAssignSubtree(string assignSubtreeChoiceCode, Tree *tree);

void assignCorrespondingTreeTasks(Tree *tree, Tree *qTree);

pair<Task *, double>
findBestCutAmong(Tree *tree, vector<Task *> candidates, string assignSubtreeChoiceCode, double initMS = -1);

Task *CutTaskWithMaxImprovement(Tree *tree, string assignSubtreeChoiceCode);

double CutTaskWithMaxImprovementHeuristicChoice(Tree *tree, string assignSubtreeChoiceCode);

void SiftInfTmaxUpPreserveOrder(vector<Task *> *taskHeap);

vector<Task *> criticalPath(Tree *tree);

vector<Task *> buildCriticalPath(Task *root);

#endif /* heuristics_h */
