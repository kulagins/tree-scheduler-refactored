#include "tree.h"
#include "cluster.h"

double IOCounter(Tree *tree, int N, double *nwghts, double *ewghts, int *chstart, int *children, int *schedule, double available_memory, bool divisible, int quiet, unsigned int &com_freq, vector<unsigned int> *brokenEdges, io_method_t method);
double IOCounterWithVariableMem(Tree *tree, int N, double *nwghts, double *ewghts, int *chstart, int *children, int *schedule,
                                Cluster *cluster, bool divisible, int quiet, unsigned int &com_freq, vector<unsigned int> *brokenEdges, io_method_t method);
