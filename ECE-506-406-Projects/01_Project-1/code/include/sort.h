#ifndef SORT_H
#define SORT_H

#include "edgelist.h"
#include "graph.h"
// Order edges by id of a source vertex,
// using the Counting Sort
// Complexity: O(E + V)

#ifdef OPENMP_HARNESS
struct Graph *radixSortEdgesBySourceOpenMP (struct Graph *graph);
#endif

#ifdef MPI_HARNESS
struct Graph *radixSortEdgesBySourceMPI (struct Graph *graph);
#endif

#ifdef HYBRID_HARNESS
struct Graph *radixSortEdgesBySourceHybrid (struct Graph *graph);
#endif
//struct Graph *countSortSerialEdgesBySource(struct Graph *graph);
//struct Graph *countSortParallelEdgesBySource(struct Graph *graph);
struct Graph *countSortEdgesBySource (struct Graph *graph,int mask,int bits);
struct Graph *radixSortEdgesBySource (struct Graph *graph);

extern int numThreads;

#endif
