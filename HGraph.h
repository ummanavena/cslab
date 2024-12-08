#ifndef HGRAPH_H
#define HGRAPH_H
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include<iostream>
#include<fstream>
#include <sys/mman.h>
#include <inttypes.h>
#include "Graph.h"

#include </usr/local/cuda/include/cuda.h>
#include </usr/local/cuda/include/cuda_runtime_api.h>
#define THRESHOLDDEGREE		10
#define INVALID               1234567890
int HEAVIEST;
struct Point {
 float x;
 float y;
 float z;
};
class HGraph: public Graph{
public:
int maxnpoints,*pnpoints,maxnedges,pnedges;
 int init();
int  readNodes(char *filename,int dims);
int  printNodes(char *filename,int dims);
int readEdges(char *filename);
int readEdges1(char *filename);
void printEdges(char *filename);
int read(char *filename);
int read2(char *filename);
int read1(char *filename);
int sortPoints();
int getPoint(int pos);
};



int HGraph::read2(char *file){//read graph from binary or .txt file
	if (!strcmp(file + strlen(file) - 3, "txt"))
	{
		readEdges(file);//if file ends with txt format, read from edgelist format using readedges
		return 1;
	}
//read graph stored in .gr format.
	std::ifstream cfile;
	cfile.open(file);
	int zerocount = 0;
	int masterFD = open(file, O_RDONLY);
	if (masterFD == -1)
	{
		printf("FileGraph::structureFromFile: unable to open %s.\n", file);
		return 1;
	}

	struct stat buf;
	int f = fstat(masterFD, &buf);
	if (f == -1)
	{
		printf("FileGraph::structureFromFile: unable to stat %s.\n", file);
                abort();
	}	

	size_t masterLength = buf.st_size;

	int _MAP_BASE = MAP_PRIVATE;

	void *m = mmap(0, masterLength, PROT_READ, _MAP_BASE, masterFD, 0);
	if (m == MAP_FAILED)
	{
		m = 0;
		printf("FileGraph::structureFromFile: mmap failed.\n");
		abort();
	}

	uint64_t *fptr = (uint64_t*) m;
	__attribute__((unused)) uint64_t version = le64toh(*fptr++);
	//assert(version == 1);
	uint64_t sizeEdgeTy = le64toh(*fptr++);
	printf("SIZEOF EDGE TYPE=%d ", sizeEdgeTy);
	uint64_t numNodes = le64toh(*fptr++);
	uint64_t numEdges = le64toh(*fptr++);
	uint64_t *outIdx = fptr;
	fptr += numNodes;
	uint32_t *fptr32 = (uint32_t*) fptr;
	uint32_t *outs = fptr32;
	fptr32 += numEdges;
	if (numEdges % 2) fptr32 += 1;
	unsigned *edgeData = (unsigned*) fptr32;
	npoints = numNodes;
	nedges = numEdges + 1;
	int max = 0, min = 99999, total = 0;
	edims = 2;
	printf("npoints=%d, nedges=%d.\n", (npoints), (nedges));
	index = (int*) malloc(sizeof(int) *(npoints + 1));
	points = (union float_int *) malloc(sizeof(union float_int) *(npoints + 1) *pdims);
	edges = (union float_int *) malloc(sizeof(union float_int) *(nedges) *edims);
	for (int i = 0; i < npoints; i++)
	{
		points[i].ipe = i;
		if (i > 0)
		{
			index[i] = le64toh(outIdx[i - 1]) + 1;
			total = le64toh(outIdx[i]) - le64toh(outIdx[i - 1]);
		}
		else
		{
			index[0] = 1;
			total = le64toh(outIdx[0]);
		}

		if (total > max) max = total;
		if (total < min) min = total;
		for (unsigned jj = 0; jj < total; ++jj)
		{
			unsigned edgeindex = index[i] + jj;
			if (edgeindex > nedges)
			{
				printf("%d %d %d\n", nedges, edgeindex, i);
				exit(0);
			}

			unsigned dst = le32toh(outs[edgeindex - 1]);
			if(i<2){ printf(" %d %d %d \n",i,dst, edgeData[edgeindex-1] );}
			edges[2 *edgeindex].ipe = dst;
			edges[2 *edgeindex + 1].ipe = edgeData[edgeindex - 1];
			//edges[2*edgeindex+1].ipe=1;
		}
	}

	index[npoints] = index[npoints - 1] + total;
	printf("ZERO COUNT=%d\n ", zerocount);
	cfile.close();
	printf("READ OVER");
	printf("MAX=%d MIN=%d\n ", max, min);
	return 0;
}






int HGraph::init(){
}
int HGraph::getPoint(int pos){
if(pos<npoints)return points[pos].ipe;
return -1;
}

int HGraph::readNodes(char *filename,int dims){
char temp[100];
char *fnode=(char *)malloc(sizeof(char)*100);
sprintf(fnode,"%s.node",filename);
FILE *fp=fopen(fnode,"r");
int index,npoints1;
float x,y;
pdims=dims;
pnpoints=(int *)malloc(sizeof(int));
fscanf(fp,"%d",&npoints1);
printf("HERE %d \n",npoints1);
points=(union float_int *)(malloc(sizeof(union float_int )*npoints1*15*pdims));
npoints=npoints1;
*pnpoints=npoints1;
maxnpoints=npoints1*3;
fgets(temp,100,fp);
if(dims==2){
for(int i=0;i<npoints1;i++){
fscanf(fp,"%d%f%f",&index,&x,&y);
if(i<10)
fprintf(stderr,"%d %f %f\n",index,x,y);
fgets(temp,100,fp);
points[2*index].fpe=x;
points[2*index+1].fpe=y;
}
}
printf("HERE %d \n",npoints1);
}

int HGraph::printNodes(char *filename,int dims){
FILE *fp=fopen(filename,"w+");
for(int i=0;i<npoints;i++){
fprintf(fp,"%d %f %f \n",i,points[2*i].fpe,points[2*i+1].fpe);
}
}
int HGraph::readEdges(char *filename){                                                                                                                       
FILE *fp=fopen(filename,"r");
fscanf(fp,"%d%d",&npoints,&nedges);
nedges++;
int src,dst,weight;
edims=2;
points=(union float_int *)malloc(sizeof(union float_int )*(npoints+1)*1);
index=(int *)malloc(sizeof(int )*(npoints+1));
      edges=(union float_int *)malloc(sizeof(union float_int)*(nedges)*edims); 
int currsrc=-1,ttindex=1,currindex=0;
for(int i=0;i<(npoints);i++){ points[i].ipe=i;}
for(int i=0;i<nedges;i++){
fscanf(fp,"%d%d%d",&src,&dst,&weight);
if(src!=currsrc){
for(int kk=currsrc+1;kk<src;kk++){ //if some points/nodes does not have any edges
index[kk]=ttindex;
}
currsrc=src;
index[src]=ttindex;
currindex=ttindex;
edges[2*ttindex].ipe=dst;
edges[2*ttindex+1].ipe=weight;
ttindex++;
}
else {
edges[2*ttindex].ipe=dst;
edges[2*ttindex+1].ipe=weight;
ttindex++;
}
}
for(int i=currsrc+1;i<=npoints;i++)index[i]=ttindex;
fclose(fp);
}
int HGraph::read1(char *file){
	std::ifstream cfile;
	cfile.open(file);

	int masterFD = open(file, O_RDONLY);
	if (masterFD == -1) {
	printf("FileGraph::structureFromFile: unable to open %s.\n", file);
	return 1;
	}

	struct stat buf;
	int f = fstat(masterFD, &buf);
	if (f == -1) {
			printf("FileGraph::structureFromFile: unable to stat %s.\n", file);
			abort();
	}
	size_t masterLength = buf.st_size;

	int _MAP_BASE = MAP_PRIVATE;

	void* m = mmap(0, masterLength, PROT_READ, _MAP_BASE, masterFD, 0);
	if (m == MAP_FAILED) {
			m = 0;
			printf("FileGraph::structureFromFile: mmap failed.\n");
			abort();
	}
	uint64_t* fptr = (uint64_t*)m;
	__attribute__((unused)) uint64_t version = le64toh(*fptr++);
	uint64_t sizeEdgeTy = le64toh(*fptr++);
printf("SIZEOF EDGE TYPE=%d ",sizeEdgeTy);
	uint64_t numNodes = le64toh(*fptr++);
	uint64_t numEdges = le64toh(*fptr++);
	uint64_t *outIdx = fptr;
	fptr += numNodes;
	uint32_t *fptr32 = (uint32_t*)fptr;
	uint32_t *outs = fptr32; 
	fptr32 += numEdges;
	if (numEdges % 2) fptr32 += 1;
	unsigned  *edgeData = (unsigned *)fptr32;
	npoints = numNodes;
	nedges = numEdges+1;
int *count,*incountar;
count=(int *)malloc(sizeof(int)*npoints);
incountar=(int *)malloc(sizeof(int)*npoints);
for(int i=0;i<npoints;i++){count[i]=incountar[i]=0;}

	printf("npoints=%d, nedges=%d.\n", (npoints), (nedges));
index=(int *)malloc(sizeof(int)*npoints);
points=(union float_int *)malloc(sizeof(union float_int )*(npoints+1)*pdims);
      edges=(union float_int *)malloc(sizeof(union float_int)*(nedges)*edims);
unsigned int total;
for(int i=0;i<(npoints);i++){
points[i].ipe=i;
if(i>0){
index[i]=le64toh(outIdx[i - 1]) + 1;
total=le64toh(outIdx[i]) - le64toh(outIdx[i - 1]);

}
else {
  index[i]=1;
 total=le64toh(outIdx[0]);
}
index[npoints]=index[npoints-1]+total;
for (unsigned jj = 0; jj <total; ++jj) {
unsigned edgeindex = index[i] + jj;
int  dst =le32toh(outs[edgeindex - 1]);
if(i<=5)printf("Edge from %d to %d with cost %d \n",i,dst,edgeData[edgeindex - 1]);
if(dst==1)printf("dst=%d src=%d \n",dst,i);
if (dst >=(npoints)) printf("\tinvalid edge\n");
edges[edgeindex*2].ipe=dst;
edges[edgeindex*2+1].ipe=edgeData[edgeindex - 1];
if(dst<npoints)(count[dst])++;
}
}
inindex=(int *)malloc(sizeof(int )*(npoints+1));
unsigned long int incount=0;
for(int i=0;i<npoints;i++){inindex[i]=incount; incount+=count[i];}
inindex[npoints]=index[npoints];
inedges=(union float_int *)malloc(sizeof(union float_int)*(nedges+1));
for(int i=0;i<npoints;i++)count[i]=0;
for(int i=1;i<nedges;i++){
int p= edges[2*i+1].ipe;
inedges[inindex[i]+count[p]].ipe=i;
inedges[inindex[i]+count[p]+1].ipe=edges[2*i+1].ipe;
(count[p])++;

}
cfile.close();
return 1;
}
int HGraph::read(char *file){
	std::ifstream cfile;
	cfile.open(file);

int zerocount=0;
	int masterFD = open(file, O_RDONLY);
  	if (masterFD == -1) {
	printf("FileGraph::structureFromFile: unable to open %s.\n", file);
	return 1;
  	}

  	struct stat buf;
	int f = fstat(masterFD, &buf);
  	if (f == -1) {
    		printf("FileGraph::structureFromFile: unable to stat %s.\n", file);
    		abort();
  	}
  	size_t masterLength = buf.st_size;

  	int _MAP_BASE = MAP_PRIVATE;

  	void* m = mmap(0, masterLength, PROT_READ, _MAP_BASE, masterFD, 0);
  	if (m == MAP_FAILED) {
    		m = 0;
    		printf("FileGraph::structureFromFile: mmap failed.\n");
    		abort();
  	}
  	uint64_t* fptr = (uint64_t*)m;
  	__attribute__((unused)) uint64_t version = le64toh(*fptr++);
  	//assert(version == 1);
  	uint64_t sizeEdgeTy = le64toh(*fptr++);
printf("SIZEOF EDGE TYPE=%d ",sizeEdgeTy);
  	uint64_t numNodes = le64toh(*fptr++);
  	uint64_t numEdges = le64toh(*fptr++);
  	uint64_t *outIdx = fptr;
  	fptr += numNodes;
  	uint32_t *fptr32 = (uint32_t*)fptr;
  	uint32_t *outs = fptr32; 
  	fptr32 += numEdges;
  	if (numEdges % 2) fptr32 += 1;
  	unsigned  *edgeData = (unsigned *)fptr32;
	
	// cuda.
	npoints = numNodes;
	nedges = numEdges+1;
int *count;
int max=0,min=99999,maxedge=0,pendcnt=0;
unsigned int total=0;
int nbrmax=0,nbrmin=9999;
count=(int *)malloc(sizeof(int)*npoints);
for(int i=0;i<npoints;i++)count[i]=0;
	printf("npoints=%d, nedges=%d.\n", (npoints), (nedges));
index=(int  *)malloc(sizeof(int)*(npoints+1));
points=(union float_int *)malloc(sizeof(union float_int )*(npoints+1)*pdims);
      edges=(union float_int *)malloc(sizeof(union float_int)*(nedges)*edims);
for(int i=0;i<(npoints);i++){
points[i].ipe=i;
if(i>0){
 index[i]=le64toh(outIdx[i - 1]) + 1;
 total=le64toh(outIdx[i]) - le64toh(outIdx[i - 1]);
}
else {
  index[0]=1;
 total=le64toh(outIdx[0]);
}
for (unsigned jj = 0; jj <total; ++jj) {
unsigned edgeindex = index[i] + jj;
int dst =le32toh(outs[edgeindex - 1]);
if (dst >=(npoints)) printf("\tinvalid edge\n");
edges[edgeindex*2].ipe=dst;
edges[edgeindex*2+1].ipe=edgeData[edgeindex - 1];
}
}
index[npoints]=index[npoints-1]+total;
cfile.close();
printf("READ OVER");
return 1;
}
struct tempforsort{
int point;
int total;
}*pointlist;
int comparefunction( int p1, int p2){
//if(p1.total==p2.total) return p1.index <p2.index;
 return p1 > p2;
}

void inline  HMIN(int *t1,int t2,int id,int &result){
if(*t1 > t2){
  if (*t1 > t2){
volatile int ch=*t1;
__sync_val_compare_and_swap(t1,ch,t2);
if(*t1<ch)result=1; 
}
}
//return result;
}
  int  HMAX(int *t1,int t2,int id,int &result){
if(*t1 <  t2){
  if (*t1 < t2){
volatile int ch=*t1;
__sync_val_compare_and_swap(t1,ch,t2);
result=1; 
}
}
return result;
}
 int HADD(int *t1,int t2){
int t4=__sync_fetch_and_add(t1,t2);
return t4;
}
 int HINCR(int *t1,int t2){
int t4=__sync_fetch_and_add(t1,t2);
return t4;
}
 int HDECR(int *t1,int t2){
int t4=__sync_fetch_and_sub(t1,t2);
return t4;
}
 int HSUB(int *t1,int t2){
int t4=__sync_fetch_and_sub(t1,t2);
return t4;
}
 int HAND(int *t1,int t2){
int t4=__sync_fetch_and_and(t1,t2);
return t4;
}

 int HOR(int *t1,int t2){
int t4=__sync_fetch_and_or(t1,t2);
return t4;
}
 int HXOR(int *t1,int t2){
int t4=__sync_fetch_and_xor(t1,t2);
return t4;
}
double rtclock()
{
    struct timezone Tzp;
    struct timeval Tp;
    int stat;
    stat = gettimeofday (&Tp, &Tzp);
    if (stat != 0) printf("Error return from gettimeofday: %d",stat);
    return(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}
#endif

