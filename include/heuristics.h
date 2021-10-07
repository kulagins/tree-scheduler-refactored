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

double SplitSubtrees(Node* root, unsigned long num_processor,  double twolevel, list<Node*>& parallelRoots, unsigned long & sequentialLength);
double SplitSubtreesV3(Node* root, unsigned long num_processor,  std::map<int, int> processor_speeds, double twolevel, list<Node*>& parallelRoots, unsigned long & sequentialLength);
double ImprovedSplit(Tree* tree, unsigned int number_processor, int* chstart, int* childrenID);
//double ImprovedSplit(Tree* tree, unsigned int number_processor);
double Merge(Tree* tree, unsigned int num_subtrees, unsigned int processor_number, double const memory_size,int * chstart,int * childrenID, bool CheckMemory);
double MergeV2(Tree* tree, unsigned int num_subtrees, unsigned int processor_number, double const memory_size,int * chstart,int * childrenID, bool CheckMemory);
double ASAP(Tree* tree, unsigned int num_processors, unsigned int depth);
double ASAP(Tree* tree, unsigned int num_processors);
unsigned long AvoidChain(Tree* tree);
//double LarSav(Tree* tree, unsigned int processor_number, unsigned int num_subtrees);
double SplitAgainV2(Tree* tree, unsigned int processor_number, unsigned int num_subtrees,  std::map<int, int>  &taskToPrc, std::map<int, bool>  &isProcBusy);
double SplitAgain(Tree* tree, unsigned int processor_number, unsigned int num_subtrees);
double Sequence(Node* root);

Tree* BuildQtree(Tree* tree);
void MemoryCheck(Tree* tree, int* chstart, int*children, double const memory_size, io_method_t method);
std::map<int, int> MemoryCheckA2(Tree* tree, int* chstart, int*children, vector<double> const memory_sizes, io_method_t method, bool skipBig,  std::map<int, int>  &taskToPrc, std::map<int, bool>  &isProcBusy);
unsigned int HowmanySubtrees(const Tree* tree, bool quiet);
void SetBandwidth(double CCR, unsigned long tree_size, double * ewghts, double * timewghts);


#endif /* heuristics_h */
