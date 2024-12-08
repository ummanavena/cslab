#include "newgpusssp1global.h"
 #include "HGraph.h"

 #include<sys/time.h>
#include<unistd.h>
struct struct_graph { 
    int   *olddist ;//has to given size of property type
    bool  *updated ;//has to given size of property type
    int   *dist ;//has to given size of property type
};
void alloc_extra_graph(HGraph &graph,int flag,int npoints)/*symtab37*/ {
     if(flag==0)graph.extra=(struct struct_graph  *)malloc(sizeof(struct struct_graph )) ;
    ((struct struct_graph  *)graph.extra)->olddist=(int  *)malloc(sizeof(int ) * npoints) ;
    ((struct struct_graph  *)graph.extra)->updated=(bool *)malloc(sizeof(bool) * npoints) ;
    ((struct struct_graph  *)graph.extra)->dist=(int  *)malloc(sizeof(int ) * npoints) ;
}
/*void alloc_extra_graph_on_device(HGraph &graph, int npoints) {
    // Allocate struct_graph on host first
    struct struct_graph *host_struct_graph = (struct struct_graph *)malloc(sizeof(struct struct_graph));

    // Allocate memory for each member on the device
    cudaMalloc(&(host_struct_graph->olddist), sizeof(int) * npoints);
    cudaMalloc(&(host_struct_graph->updated), sizeof(bool) * npoints);
    cudaMalloc(&(host_struct_graph->dist), sizeof(int) * npoints);

    // Copy struct_graph to device
    cudaMalloc(&(graph.extra), sizeof(struct struct_graph));
    cudaMemcpy(graph.extra, host_struct_graph, sizeof(struct struct_graph), cudaMemcpyHostToDevice);

    // Free the temporary host copy
    free(host_struct_graph);
}*/

void read_and_allocate_graph(HGraph  &graph ){
printf("enter number of points and edges");
 scanf("%d%d",&(graph.npoints),&(graph.nedges));
 graph.points=(union float_int *)malloc(sizeof(union float_int)*graph.npoints);
graph.edges=(union float_int *)malloc(sizeof(union float_int)*graph.nedges*2);
alloc_extra_graph(graph,0,graph.npoints);
}
//resetgraph graph

