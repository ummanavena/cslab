#include "newgpusssp1global.h"
 #include "gpu.h"
 #include<sys/time.h>
#include<unistd.h>
struct struct_graph { 
    int   *olddist ;//has to given size of property type
    bool  *updated ;//has to given size of property type
    int   *dist ;//has to given size of property type
};
#include<sys/time.h>
#include </usr/local/cuda/include/cuda.h>
#include </usr/local/cuda/include/cuda_runtime_api.h>
#include<unistd.h>

//void alloc_extra_graph(HGraph &graph,int flag,int npoints)/*symtab37*/ {
//     if(flag==0)graph.extra=(struct struct_graph  *)malloc(sizeof(struct struct_graph )) ;
//    ((struct struct_graph  *)graph.extra)->olddist=(int  *)malloc(sizeof(int ) * npoints) ;
//    ((struct struct_graph  *)graph.extra)->updated=(bool *)malloc(sizeof(bool) * npoints) ;
//    ((struct struct_graph  *)graph.extra)->dist=(int  *)malloc(sizeof(int ) * npoints) ;
//}


void alloc_extra_hgraph(HGraph &graph,int flag,int npoints)/*symtab37*/ {
     if(flag==0)graph.extra=(struct struct_graph  *)malloc(sizeof(struct struct_graph )) ;
    ((struct struct_graph  *)graph.extra)->olddist=(int  *)malloc(sizeof(int ) * npoints) ;
    ((struct struct_graph  *)graph.extra)->updated=(bool *)malloc(sizeof(bool) * npoints) ;
    ((struct struct_graph  *)graph.extra)->dist=(int  *)malloc(sizeof(int ) * npoints) ;
}
void alloc_extra_graph(GGraph &graph,int flag) {
struct struct_graph temp;
if(flag==0)cudaMalloc((void **)&(graph.extra),sizeof(struct struct_graph ));
	cudaMemcpy(&temp,((struct struct_graph *)(graph.extra)),sizeof(struct struct_graph),cudaMemcpyDeviceToHost);
	cudaMalloc( (void **)&(temp.olddist),sizeof(int )* graph.npoints);
	cudaMalloc( (void **)&(temp.updated),sizeof(bool)* graph.npoints);
	cudaMalloc( (void **)&(temp.dist),sizeof(int )* graph.npoints);
if(cudaMemcpy(graph.extra,&temp,sizeof(struct struct_graph ),cudaMemcpyHostToDevice)!=cudaSuccess)printf("memcpyerror 0");
}
void alloc_extra_graph1(GGraph &graph1,int flag) {
cudaSetDevice(1);
struct struct_graph temp;
if(flag==0)cudaMalloc((void **)&(graph1.extra),sizeof(struct struct_graph ));
	cudaMemcpy(&temp,((struct struct_graph *)(graph1.extra)),sizeof(struct struct_graph ),cudaMemcpyDeviceToHost);
	cudaMalloc( (void **)&(temp.olddist),sizeof(int )* graph1.npoints);
        cudaMalloc( (void **)&(temp.updated),sizeof(bool)* graph1.npoints);
        cudaMalloc( (void **)&(temp.dist),sizeof(int )* graph1.npoints);
if(cudaMemcpy(graph1.extra,&temp,sizeof(struct struct_graph ),cudaMemcpyHostToDevice)!=cudaSuccess)printf("memcpyerror 1");
cudaSetDevice(0);
}
void read_and_allocate_graph(HGraph  &hgraph ){
printf("enter number of points and edges");
 scanf("%d%d",&(hgraph.npoints),&(hgraph.nedges));
 hgraph.points=(union float_int *)malloc(sizeof(union float_int)*hgraph.npoints);
hgraph.edges=(union float_int *)malloc(sizeof(union float_int)*hgraph.nedges*2);
alloc_extra_hgraph(hgraph,0,hgraph.npoints);
//resetgraph graph
}

