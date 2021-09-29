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

typedef enum {TIME=1,SPACE} subtreeType;

double SplitSubtrees(Cnode* root, unsigned long num_processor,  double twolevel, list<Cnode*>& parallelRoots, unsigned long & sequentialLength);
double SplitSubtreesV3(Cnode* root, unsigned long num_processor,  std::map<int, int> processor_speeds, double twolevel, list<Cnode*>& parallelRoots, unsigned long & sequentialLength);
double ImprovedSplit(Ctree* tree, unsigned int number_processor, int* chstart, int* childrenID);
//double ImprovedSplit(Ctree* tree, unsigned int number_processor);
double Merge(Ctree* tree, unsigned int num_subtrees, unsigned int processor_number, double const memory_size,int * chstart,int * childrenID, bool CheckMemory);
double MergeV2(Ctree* tree, unsigned int num_subtrees, unsigned int processor_number, double const memory_size,int * chstart,int * childrenID, bool CheckMemory);
double ASAP(Ctree* tree, unsigned int num_processors, unsigned int depth);
double ASAP(Ctree* tree, unsigned int num_processors);
unsigned long AvoidChain(Ctree* tree);
//double LarSav(Ctree* tree, unsigned int processor_number, unsigned int num_subtrees);
double SplitAgainV2(Ctree* tree, unsigned int processor_number, unsigned int num_subtrees,  std::map<int, int>  &taskToPrc, std::map<int, bool>  &isProcBusy);
double SplitAgain(Ctree* tree, unsigned int processor_number, unsigned int num_subtrees);
double Sequence(Cnode* root);

Ctree* BuildQtree(Ctree* tree);
void MemoryCheck(Ctree* tree, int* chstart, int*children, double const memory_size, io_method_t method);
std::map<int, int> MemoryCheckA2(Ctree* tree, int* chstart, int*children, vector<double> const memory_sizes, io_method_t method, bool skipBig,  std::map<int, int>  &taskToPrc, std::map<int, bool>  &isProcBusy);
unsigned int HowmanySubtrees(const Ctree* tree, bool quiet);
void SetBandwidth(double CCR, unsigned long tree_size, double * ewghts, double * timewghts);


#endif /* heuristics_h */
