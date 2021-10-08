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
struct Graph *countSortSerialEdgesBySource (struct Graph *graph)
{

    int i;
    int key;
    int pos;
    struct Edge *sorted_edges_array = newEdgeArray(graph->num_edges);

    // auxiliary arrays, allocated at the start up of the program
    int *vertex_count = (int *)malloc(graph->num_vertices * sizeof(int)); // needed for Counting Sort

    for(i = 0; i < graph->num_vertices; ++i)
    {
        vertex_count[i] = 0;
    }

    // count occurrence of key: id of a source vertex
    for(i = 0; i < graph->num_edges; ++i)
    {
        key = graph->sorted_edges_array[i].src;
        vertex_count[key]++;
    }

    // transform to cumulative sum
    for(i = 1; i < graph->num_vertices; ++i)
    {
        vertex_count[i] += vertex_count[i - 1];
    }

    // fill-in the sorted array of edges
    for(i = graph->num_edges - 1; i >= 0; --i)
    {
        key = graph->sorted_edges_array[i].src;
        pos = vertex_count[key] - 1;
        sorted_edges_array[pos] = graph->sorted_edges_array[i];
        vertex_count[key]--;
    }



    free(vertex_count);
    free(graph->sorted_edges_array);

    graph->sorted_edges_array = sorted_edges_array;

    return graph;

}
struct Graph *countSortParallelEdgesBySource (struct Graph *graph)
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
            key = graph->sorted_edges_array[i].src;
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
            key = graph->sorted_edges_array[i].src;
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
    int radixnumber = 8;
    int mask1L1 = 0x00000003;
    int mask2L1 = 0x0000000c;
    int mask3L1 = 0x00000030;
    int mask4L1 = 0x000000c0;
    int mask5L1 = 0x00000300;
    int mask6L1 = 0x00000c00;
    int mask7L1 = 0x00003000;
    int mask8L1 = 0x0000c000;
    int mask9L1 = 0x00030000;
    int mask10L1 = 0x000c0000;
    int mask11L1 = 0x00300000;
    int mask12L1 = 0x00c00000;
    int mask13L1 = 0x03000000;
    int mask14L1 = 0x0c000000;
    int mask15L1 = 0x30000000;
    int mask16L1 = 0xc0000000;
    int bit1L1 = 0;
    int bit2L1 = 2;
    int bit3L1 = 4;
    int bit4L1 = 6;
    int bit5L1 = 8;
    int bit6L1 = 10;
    int bit7L1 = 12;
    int bit8L1 = 14;
    int bit9L1 = 16;
    int bit10L1 = 18;
    int bit11L1 = 20;
    int bit12L1 = 22;
    int bit13L1 = 24;
    int bit14L1 = 26;
    int bit15L1 = 28;
    int bit16L1 = 30;
    int mask1L2 = 0x0000000f;
    int mask2L2 = 0x000000f0;
    int mask3L2 = 0x00000f00;
    int mask4L2 = 0x0000f000;
    int mask5L2 = 0x000f0000;
    int mask6L2 = 0x00f00000;
    int mask7L2 = 0x0f000000;
    int mask8L2 = 0xf0000000;
    int bit1L2 = 0;
    int bit2L2 = 4;
    int bit3L2 = 8;
    int bit4L2 = 12;
    int bit5L2 = 16;
    int bit6L2 = 20;
    int bit7L2 = 24;
    int bit8L2 = 28;
    int mask1L3 = 0x000000ff;
    int mask2L3 = 0x0000ff00;
    int mask3L3 = 0x00ff0000;
    int mask4L3 = 0xff000000;
    int bit1L3 = 0;
    int bit2L3 = 8;
    int bit3L3 = 16;
    int bit4L3 = 24;
    int mask1L4 = 0x0000ffff;
    int mask2L4 = 0xffff0000;
    int bit1L4 = 0;
    int bit2L4 = 16;
    if(radixnumber==2)
    {
        for(i=0;i<16;++i)
        {
            
            if(i==0)
            {
                countSortEdgesBySource(graph,mask1L1,bit1L1);
            }
            else if(i==1)
            {
                countSortEdgesBySource(graph,mask2L1,bit2L1);
            }
            else if(i==2)
            {
                countSortEdgesBySource(graph,mask3L1,bit3L1);
            }
            else if(i==3)
            {
                countSortEdgesBySource(graph,mask4L1,bit4L1);
            }
            else if(i==4)
            {
                countSortEdgesBySource(graph,mask5L1,bit5L1);
            }
            else if(i==5)
            {
                countSortEdgesBySource(graph,mask6L1,bit6L1);
            }
            else if(i==6)
            {
                countSortEdgesBySource(graph,mask7L1,bit7L1);
            }
            else if(i==7)
            {
                countSortEdgesBySource(graph,mask8L1,bit8L1);
            }
            else if(i==8)
            {
                countSortEdgesBySource(graph,mask9L1,bit9L1);
            }
            else if(i==9)
            {
                countSortEdgesBySource(graph,mask10L1,bit10L1);
            }
            else if(i==10)
            {
                countSortEdgesBySource(graph,mask11L1,bit11L1);
            }
            else if(i==11)
            {
                countSortEdgesBySource(graph,mask12L1,bit12L1);
            }
            else if(i==12)
            {
                countSortEdgesBySource(graph,mask13L1,bit13L1);
            }
            else if(i==13)
            {
                countSortEdgesBySource(graph,mask14L1,bit14L1);
            }
            else if(i==14)
            {
                countSortEdgesBySource(graph,mask15L1,bit15L1);
            }
            else if(i==15)
            {
                countSortEdgesBySource(graph,mask16L1,bit16L1);
            }
        }
    }
    if(radixnumber==4)
    {
         for(i=0;i<8;++i)
        {
            if(i==0)
            {
                countSortEdgesBySource(graph,mask1L2,bit1L2);
            }
            else if(i==1)
            {
                countSortEdgesBySource(graph,mask2L2,bit2L2);
            }
            else if(i==2)
            {
                countSortEdgesBySource(graph,mask3L2,bit3L2);
            }
            else if(i==3)
            {
                countSortEdgesBySource(graph,mask4L2,bit4L2);
            }
            else if(i==4)
            {
                countSortEdgesBySource(graph,mask5L2,bit5L2);
            }
            else if(i==5)
            {
                countSortEdgesBySource(graph,mask6L2,bit6L2);
            }
            else if(i==6)
            {
                countSortEdgesBySource(graph,mask7L2,bit7L2);
            }
            else if(i==7)
            {
                countSortEdgesBySource(graph,mask8L2,bit8L2);
            }
        }

    }
    if(radixnumber==8)
    {
        for(i=0;i<4;++i)
        {
            if(i==0)
            {
                countSortEdgesBySource(graph,mask1L3,bit1L3);
            }
            else if(i==1)
            {
                countSortEdgesBySource(graph,mask2L3,bit2L3);
            }
            else if(i==2)
            {
                countSortEdgesBySource(graph,mask3L3,bit3L3);
            }
            else if(i==3)
            {
                countSortEdgesBySource(graph,mask4L3,bit4L3);
            }
        }
    }
    if(radixnumber==16)
    {
        for(i=0;i<2;i++)
        {
            if(i==0)
            {
                countSortEdgesBySource(graph,mask1L4,bit1L4);
            }
            else if(i==1)
            {
                countSortEdgesBySource(graph,mask2L4,bit2L4);
            }
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
