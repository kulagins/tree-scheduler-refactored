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
#include "lib-io-tree-minmem.h"
#include "cluster.h"

template <class T, class U>
void GetTwoLargestElementTypetwo(T container, U &Largest, U &secondLargest);


double SplitSubtrees(Task* root, unsigned long num_processor,  double twolevel, list<Task*>& parallelRoots, unsigned long & sequentialLength);
double ImprovedSplit(Tree* tree, unsigned int number_processor, int* chstart, int* childrenID);
double Merge(Tree* tree, unsigned int num_subtrees, unsigned int processor_number, double const memory_size,int * chstart,int * childrenID, bool CheckMemory);
double MergeV2(Tree* tree, unsigned int num_subtrees, unsigned int processor_number, double const memory_size,int * chstart,int * childrenID, bool CheckMemory);
double ASAP(Tree* tree, unsigned int num_processors);

double SplitAgainV2(Tree* tree, unsigned int processor_number, unsigned int num_subtrees,  std::map<int, int>  &taskToPrc, std::map<int, bool>  &isProcBusy);
double SplitAgain(Tree* tree, unsigned int processor_number, unsigned int num_subtrees);


void MemoryCheck(Tree* tree, int* chstart, int*children,  Cluster *cluster, io_method_t method);
std::map<int, int> MemoryCheckA2(Tree* tree, Cluster *cluster,io_method_t method, bool skipBig);
void SetBandwidth(double CCR, unsigned long tree_size, double * ewghts, double * timewghts);


#endif /* heuristics_h */
