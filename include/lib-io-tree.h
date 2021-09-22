#ifndef LIB_IO_TREE_H
#define LIB_IO_TREE_H



#ifdef __cplusplus


#include <assert.h>

#include "lib-io-tree-utils.h"


#define MIN_ACCURACY 0.001

struct iter_node_t;

extern "C" {
#endif

#ifndef __cplusplus
typedef enum {FURTHEST_NODE=1, BEST_K_COMBI, BEST_FIT_ABS, FIRST_FIT_ABS, BEST_FIT, FIRST_FIT,BEST_INC_COMBI, BEST_COMBI,LARGEST_FIT} io_method_t;
#endif

  
/* Bora's functions (mex interface) */
double MaxOutDegree(int N, int * prnts, double * nwghts, double * ewghts);
#ifndef __cplusplus
//double IOCounter(int N, int * prnts, double * nwghts, double * ewghts, int * schedule, double available_memory,int divisible,io_method_t method);
#else
//double IOCounter(int N, int * prnts, double * nwghts, double * ewghts, int * schedule, double available_memory,int divisible,io_method_t method=FURTHEST_NODE);
#endif

double PostOrderRecurAlgorithm(int N, int *prnts, double *nwghts, double *ewghts, int *schedule);
double PostOrderRecurAlgorithm_timed(int N, int *prnts, double *nwghts, double *ewghts, int *schedule,double * usec,int quiet);


double PostOrderIterAlgorithm(int N, int *prnts, double *nwghts, double *ewghts, int *schedule);
double PostOrderIterAlgorithm_timed(int N, int *prnts, double *nwghts, double *ewghts, int *schedule, double * usec,int quiet);

double PebbleOrderingRecurAlgorithm(int N, int *prnts, double *nwghts, double *ewghts, int *schedule);
double PebbleOrderingRecurAlgorithm_timed(int N, int *prnts, double *nwghts, double *ewghts, int *schedule,double * usec,int quiet);

double PebbleOrderingIterAlgorithm(int N, int *prnts, double *nwghts, double *ewghts, int *schedule);
double PebbleOrderingIterAlgorithm_timed(int N, int *prnts, double *nwghts, double *ewghts, int *schedule,double * usec,int quiet);

#ifndef DEBUG_USING_MINMEM
double MinMemRecurAlgorithm( int N, int *prnts, double *nwghts, double *ewghts, int *schedule);
double MinMemRecurAlgorithm_timed(int N, int *prnts, double *nwghts, double *ewghts, int *schedule,double * usec,int quiet);
#else
double MinMemRecurAlgorithm( int N, int *prnts, double *nwghts, double *ewghts, int *schedule,iter_node_t * minmem_trace);
double MinMemRecurAlgorithm_timed(int N, int *prnts, double *nwghts, double *ewghts, int *schedule,double * usec,int quiet,iter_node_t * minmem_trace);
#endif
	
double MinMemArrayRecurAlgorithm(int N, int *prnts, double *nwghts, double *ewghts, int *schedule);
double MinMemArrayRecurAlgorithm_timed(int N, int *prnts, double *nwghts, double *ewghts, int *schedule,double * usec,int quiet);
#ifdef __cplusplus
} /* closing brace for extern "C" */
#endif


#endif
