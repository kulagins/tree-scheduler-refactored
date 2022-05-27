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




int MemoryCheck(Tree *tree, io_method_t method, bool useMinimalAvailableProvcessor);
int MemoryCheckHomp(Tree *tree, io_method_t method, double processor_memory_size);

//void
//MemoryCheckA2(Tree *tree, Cluster *cluster, io_method_t method, bool skipBig);

void SetBandwidth(double CCR, unsigned long tree_size, double *ewghts, double *timewghts);
void distributeProcessors(Tree * qTree);

void seqSetAndFeasSets(Tree * tree);
double assignToBestProcessors(Tree *tree);
void removeProcessorFromAllFeasSets(Processor * processor, Tree * tree);

#endif /* heuristics_h */
