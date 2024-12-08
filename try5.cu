#include "new.h"
#include<sys/time.h>
#include </usr/local/cuda/include/cuda.h>
#include </usr/local/cuda/include/cuda_runtime_api.h>
#include<unistd.h>

#define t (((struct struct_graph *)(graph.extra)))
#define t2 (((struct struct_graph *)(graph1.extra)))
__device__ int   changed =0, hchanged =0,changed2=0;

__global__ void relaxgraph(GGraph graph,int x) {
        int id = blockIdx.x * blockDim.x + threadIdx.x+x;
        if (id <graph.npoints){
        if (t->updated[id] == true){
        t->updated[id]=false;
        int ind0 = graph.index[id];
        int ind1 = graph.index[id+1]-graph.index[id];
        for (int ind2 = 0; ind2 < ind1; ind2++) {
                int ut0 = 2 * (ind0 + ind2); //edge index
                int ut1 = graph.edges[ut0].ipe; //dest point
                int ut2 = graph.edges[ut0 + 1].ipe;
                GMIN(&t->dist[ut1], t->dist[id] + ut2, changed);
		//changed=1;
		//printf("changed to %d\n",t->dist[ut1]);
        }
        } }
}
__global__ void relaxgraph2(GGraph graph1,int x) {
        int id = blockIdx.x * blockDim.x + threadIdx.x+x;
        if (id <graph1.npoints){
        if (t2->updated[id] == true){
        t2->updated[id]=false;
        int ind0 = graph1.index[id];
        int ind1 = graph1.index[id+1]-graph1.index[id];
        for (int ind2 = 0; ind2 < ind1; ind2++) {
                int ut0 = 2 * (ind0 + ind2); //edge index
                int ut1 = graph1.edges[ut0].ipe; //dest point
                int ut2 = graph1.edges[ut0 + 1].ipe;
                GMIN(&t2->dist[ut1], t2->dist[id] + ut2, changed2);
		//printf("changed to %d\n",t2->dist[ut1]);
	
        }
        } }
}


__global__ void   reset ( GGraph  graph,int x )
 {
    int id = blockIdx.x * blockDim.x + threadIdx.x + x;
    if(id<graph.npoints){
    t->dist[id] = 1234567890;
    t->olddist[id] = 1234567890;
    t->updated[id] = false;

    }

 }
 __global__ void   reset2 ( GGraph  graph1,int x )
 {
    int id = blockIdx.x * blockDim.x + threadIdx.x + x;
    if(id<graph1.npoints){
    t2->dist[id] = 1234567890;
    t2->olddist[id] = 1234567890;
    t2->updated[id] = false;

    }

 }
__global__ void   reset1 (GGraph graph,int x )
 {
 int id = blockIdx.x * blockDim.x + threadIdx.x + x;
 if(id<graph.npoints){
 if (t->dist[id] < t->olddist[id]) {
       t->updated[id] = true;
 }
 t->olddist[id] = t->dist[id];
 }
}
__global__ void   reset3 (GGraph graph1,int x )
 {
 int id = blockIdx.x * blockDim.x + threadIdx.x + x;
 if(id<graph1.npoints){
 if (t2->dist[id] < t2->olddist[id]) {
       t2->updated[id] = true;
 }
 t2->olddist[id] = t2->dist[id];
 }
}
void   SSSP ( char    *  name )
 {
	HGraph hgraph ;
	GGraph graph;
	GGraph graph1;
	hgraph.read2(name);
	
	int hosthgraph=1;
	hgraph.extra=(struct struct_graph *)malloc(sizeof(struct struct_graph ));
	


	alloc_extra_hgraph(hgraph,hosthgraph,hgraph.npoints);
	hgraph.cloneGPU(graph,0 );
	int TPB=1024;
	int kb;
	if ((graph.npoints / TPB + 1) > (32 * 1024))
	    kb = (32 * 1024);
	else
	    kb = (graph.npoints / TPB + 1);

	int graphflag=0;
	cudaSetDevice(0);
	cudaMalloc((void **)(&graph.extra),sizeof(struct struct_graph ));
	struct struct_graph ftemp1;
	if(cudaMemcpy(&ftemp1,graph.extra,sizeof(struct struct_graph ),cudaMemcpyDeviceToHost)!=cudaSuccess)printf("memcpyerror 2");

	if(cudaMemcpy(graph.extra,&ftemp1,sizeof(struct struct_graph ),cudaMemcpyHostToDevice)!=cudaSuccess)printf("memcpyerror 3");
	cudaSetDevice(0);
	graphflag=1;
	alloc_extra_graph(graph,graphflag);
	int falcvt1;
	hgraph.cloneGPU(graph1,1);
	int TPB1=1024;

	int kb1;
	if ((graph.npoints / TPB + 1) > (32 * 1024))
	    kb1 = (32 * 1024);
	else
	    kb1 = (graph.npoints / TPB + 1);
	int graph1flag=0;
	cudaSetDevice(1);
	cudaMalloc((void **)(&graph1.extra),sizeof(struct struct_graph ));
	struct struct_graph ftemp4;
	if(cudaMemcpy(&ftemp4,graph1.extra,sizeof(struct struct_graph ),cudaMemcpyDeviceToHost)!=cudaSuccess)printf("memcpyerror 5");
	
	if(cudaMemcpy(graph1.extra,&ftemp4,sizeof(struct struct_graph ),cudaMemcpyHostToDevice)!=cudaSuccess)printf("memcpyerror 6");
	cudaSetDevice(0);
	graph1flag=1;
	alloc_extra_graph1(graph1,graph1flag);
	cudaSetDevice(0);
	for(int i=0;i<graph.npoints;i+=kb*TPB){
		reset<<<kb,TPB>>>(graph,i);
	}
	cudaDeviceSynchronize();
	cudaSetDevice(0);
	cudaSetDevice(1);
	for(int i=0;i<graph1.npoints;i+=kb1*TPB1){
		reset2<<<kb1,TPB1>>>(graph1,i);
	}
	cudaDeviceSynchronize();
	cudaSetDevice(0);
	int f2;
	f2=0;
	struct struct_graph cpy5;
	cudaMemcpy(&cpy5,((struct struct_graph *)(graph.extra)),sizeof(struct struct_graph ),cudaMemcpyDeviceToHost);
	if(cudaMemcpy(&(cpy5.dist[0]),&(f2),sizeof(int ),cudaMemcpyHostToDevice)!=cudaSuccess)printf("memcpyerror 7");
	int f3;
	f3=0;
	struct struct_graph cpy6;
	cudaMemcpy(&cpy6,((struct struct_graph *)(graph1.extra)),sizeof(struct struct_graph ),cudaMemcpyDeviceToHost);
	if(cudaMemcpy(&(cpy6.dist[0]),&(f3),sizeof(int ),cudaMemcpyHostToDevice)!=cudaSuccess)printf("memcpyerror 8");

	bool updt=true;
	if (cudaMemcpy(&(cpy5.updated[0]), &updt, sizeof(bool), cudaMemcpyHostToDevice) != cudaSuccess)
	    printf("memcpyerror 1");
	if (cudaMemcpy(&(cpy6.updated[0]), &updt, sizeof(bool), cudaMemcpyHostToDevice) != cudaSuccess)
	    printf("memcpyerror 1");


	#pragma omp parallel sections
	{
		#pragma omp sections
		{
		while(1) {
			int f4;
			f4=0;
			//struct struct_graph cpy7;
			
			//cudaMemcpy(&cpy7,((struct struct_graph *)(graph1.extra)),sizeof(struct struct_graph ),cudaMemcpyDeviceToHost);
			if(cudaMemcpyToSymbol(changed2,&f4,sizeof(int ),0,cudaMemcpyHostToDevice)!=cudaSuccess)printf("memcpyerror 11");
			cudaSetDevice(1);
			
			for(int i=0;i<graph1.npoints/2;i+=kb1*TPB1){
				//printf("enter 2");
				relaxgraph2<<<kb1,TPB1>>>(graph1,i);
			}
			cudaDeviceSynchronize();
			cudaSetDevice(0);
			int f5;
			if(cudaMemcpyFromSymbol(&f5,(changed2),sizeof(int ),0,cudaMemcpyDeviceToHost)!=cudaSuccess)printf("memcpyerror 13");
			if( f5==0 ){printf("yes2");break;}
			cudaSetDevice(0);
	    		for (int pointIdx = 0; pointIdx < graph1.npoints; pointIdx += kb1 * TPB1) {
				reset3<<<kb1, TPB1>>>(graph1, pointIdx);
	    		}
			}
		}
		#pragma omp sections
		{
		while(1) {
			int f6;
			f6=0;
			
			if(cudaMemcpyToSymbol(changed,&(f6),sizeof(int ),0,cudaMemcpyHostToDevice)!=cudaSuccess)printf("memcpyerror 11");
			cudaSetDevice(0);
			for(int i=0;i<graph.npoints;i+=kb*TPB){
				//printf("enter");
				relaxgraph<<<kb,TPB>>>(graph,i);
			}
			cudaDeviceSynchronize();
			cudaSetDevice(0);
			int f7;
			if(cudaMemcpyFromSymbol(&f7,changed,sizeof(int ),0,cudaMemcpyDeviceToHost)!=cudaSuccess)printf("memcpyerror 13");
			if( f7==0 ){break;}
			cudaSetDevice(0);
	    		for (int pointIdx = 0; pointIdx < graph.npoints; pointIdx += kb * TPB) {
				reset1<<<kb, TPB>>>(graph, pointIdx);
	    		}
		}
		}
	}
	struct struct_graph cpy10;
	int *finalDistanceArray = (int *)malloc(sizeof(int) * graph.npoints);
	cudaMemcpy(&cpy10,((struct struct_graph *)(graph.extra)),sizeof(struct struct_graph ),cudaMemcpyDeviceToHost);
	if(cudaMemcpy(finalDistanceArray,(cpy10.dist),sizeof(int)*graph.npoints,cudaMemcpyDeviceToHost)!=cudaSuccess)printf("memcpyerror 14");
	struct struct_graph cpy11;
	int *finalDistanceArray1 = (int *)malloc(sizeof(int) * graph1.npoints);
	cudaMemcpy(&cpy11,((struct struct_graph *)(graph1.extra)),sizeof(struct struct_graph ),cudaMemcpyDeviceToHost);
	if(cudaMemcpy(finalDistanceArray1,(cpy11.dist),sizeof(int)*graph.npoints,cudaMemcpyDeviceToHost)!=cudaSuccess)printf("memcpyerror 15");
	for (int i = 0; i < graph1.npoints; i++)
	    printf("%d \n", min(finalDistanceArray1[i],finalDistanceArray[i]));


	return ;
	}
		
int   main ( int   argc ,char    *  argv [ ] )
 {
if(argc>2)FALC_THREADS=atoi(argv[2]);

 if( argc!=4  )
{
printf("error:-exec -t threads  file");
return 1;
}
SSSP(argv[3]);

}

