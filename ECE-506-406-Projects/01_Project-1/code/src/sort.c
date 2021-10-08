#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#ifdef OPENMP_HARNESS
#include <omp.h>
#endif

#ifdef MPI_HARNESS
#include <mpi.h>
#endif

#ifdef HYBRID_HARNESS
#include <omp.h>
#include <mpi.h>
#endif

#include "sort.h"
#include "graph.h"

struct Graph *countSortEdgesBySource (struct Graph *graph,int mask,int bits)
{

    int i;
    int j;
    int key;
    int pos;
    int base=0;
    int tid;
    int end_offset;
    int start_offset;
    int P;
    int m = mask;
    int b = bits;
    int **vertex_count;
    struct Edge *sorted_edges_array = newEdgeArray(graph->num_edges);
    #pragma omp parallel private(i,j,key,pos,start_offset,end_offset,tid,base)shared(graph,sorted_edges_array)
    {
        tid = omp_get_thread_num();
        //printf("tid: %d\n",tid);
        P = omp_get_num_threads();
        //printf("P: %d\n",P);
        if(tid==0)
        {
            vertex_count = (int**)malloc(P*sizeof(int*));
            for(i=0;i<P;i++)
            {
                vertex_count[i] = (int*)malloc(graph->num_vertices*sizeof(int));
                //printf("tid: %d\n",tid);
                //printf("vc: %d\n",vertex_count[i]);
            }
        }
        //printf("E: %d\n",graph->num_edges);
        //printf("V: %d\n",graph->num_vertices);
        #pragma omp barrier
        for(i=0;i<graph->num_vertices;++i)
        {
            vertex_count[tid][i] = 0;
            //count++;
            //printf("tid: %d\n",tid);
            //printf("L1: %d\n",vertex_count[tid][i]);
            //printf("c: %d\n",count);
        }
        start_offset = tid*(graph->num_edges/P);
        //printf("S: %d\n",start_offset);
        if(tid==(P-1))
        {
            end_offset = start_offset+(graph->num_edges/P)+(graph->num_edges%P);
        }
        else
        {
            end_offset = start_offset+(graph->num_edges/P);
        }
        //printf("tid: start: end: %d %d %d\n",tid,start_offset,end_offset);
        //printf("E: %d\n",end_offset);
        for(i=start_offset;i<end_offset;++i)
        {
            key = (graph->sorted_edges_array[i].src & m)>>b;
            //printf("tid: %d\n",tid);
            //printf("k: %d\n",key);
            vertex_count[tid][key]++;
            //printf("tid: %d\n",tid);
            //printf("Vc: %d\n",vertex_count[tid][key]);
        }
        /*for(i=0;i<graph->num_vertices;i++)
        {
            printf("tid: vc: %d %d\n",tid,vertex_count[tid][i]);
        }*/
        #pragma omp barrier
        printf("id = %d, base: %d\n",tid, base);
        #pragma omp barrier
        if(tid==0)
        {
            for(i=0;i<graph->num_vertices;++i)
            {
                for(j=0;j<P;++j)
                {
                    base = base + vertex_count[j][i];
                    vertex_count[j][i] = base;
                    //printf("tid: i: j: base: vertexcount: %d %d %d %d %d\n",tid,i,j,base,vertex_count[j][i]);
                }
            }
        }
        #pragma omp barrier
        /*for(i=0;i<graph->num_vertices;i++)
        {
            for(j=0;j<P;j++)
            {
                printf("vertexxount: %d\n",vertex_count[j][i]);
            }
            //count++;
            //printf("count: %d\n",count);
            //printf("tid: %d\n",tid);
        }*/
        //printf("tid: start: end: %d %d %d\n",tid,start_offset,end_offset);
        //printf("s: %d\n",start_offset);
        //printf("e: %d\n",end_offset);
        for(i=end_offset-1;i>=start_offset;--i)
        {
            //printf("tid: start: end: %d %d %d\n",tid,start_offset,end_offset);
            key = (graph->sorted_edges_array[i].src & m)>>b;
            //printf("K: %d\n",key);
            //printf("tid: key: VC: %d %d %d\n",tid,key,vertex_count[tid][key]);
            //printf("tid: key: %d %d\n",tid,key);
            pos = vertex_count[tid][key]-1;
            //printf("pos: i: sortededgesarray: %d %d %d\n",pos,i,sorted_edges_array[i]);
            //printf("G: %d\n",graph->sorted_edges_array[i]);
            sorted_edges_array[pos] = graph->sorted_edges_array[i];
            vertex_count[tid][key]--;
        }
    }
    free(vertex_count);
    //free(graph->sorted_edges_array);
    graph->sorted_edges_array = sorted_edges_array;

    return graph;

}

struct Graph *radixSortEdgesBySource (struct Graph *graph)
{
    int i;
    int mask1 = 0x000000ff;
    int mask2 = 0x0000ff00;
    int mask3 = 0x00ff0000;
    int mask4 = 0xff000000;
    int bit1 = 0;
    int bit2 = 8;
    int bit3 = 16;
    int bit4 = 24;
    for(i=0;i<4;++i)
    {
        if(i==0)
        {
            countSortEdgesBySource(graph,mask1,bit1);
        }
        else if(i==1)
        {
            countSortEdgesBySource(graph,mask2,bit2);
        }
        else if(i==2)
        {
            countSortEdgesBySource(graph,mask3,bit3);
        }
        else if(i==3)
        {
            countSortEdgesBySource(graph,mask4,bit4);
        }

    }

}
#ifdef OPENMP_HARNESS
struct Graph *radixSortEdgesBySourceOpenMP (struct Graph *graph)
{


    printf("*** START Radix Sort Edges By Source OpenMP *** \n");
    return graph;
}
#endif

#ifdef MPI_HARNESS
struct Graph *radixSortEdgesBySourceMPI (struct Graph *graph)
{

    printf("*** START Radix Sort Edges By Source MPI*** \n");
    return graph;
}
#endif

#ifdef HYBRID_HARNESS
struct Graph *radixSortEdgesBySourceHybrid (struct Graph *graph)
{

    printf("*** START Radix Sort Edges By Source Hybrid*** \n");
    return graph;
}
#endif