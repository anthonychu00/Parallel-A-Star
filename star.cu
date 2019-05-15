//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//You need to compile the program as "nvcc -G star.cu"
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include <assert.h>
#include <curand.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <math.h>
#include "timerc.h"
__device__ __host__ void shiftqueue(int* array_sh, int dim,int k);
__device__ __host__ int checkIfQueueEmpty(int * array_sh, int dim, int k);
__device__ __host__ int heuristic(int dim, int current, int end);
__device__ int findNextInsertionPoint(int * sizes, int k);
__device__ int duplicateAdjacents(int * adjacents, int a, int k, int currentPos);
__device__ void organizeQueue(int * queue, int targetLocation, int * heuristics, int h, int which, int dim);
__device__ __host__ void print_board(int * game, int dim);

__global__ void traverse(int * grid, int start, int end, int dim, int k, int * result){//k is number of queues
									//6 and 5 in this example
	__shared__ int prevNode [36];//previous Node to print out route
	__shared__ int lowestCost [36];//current lowest cost to get to that Node
	__shared__ int array [30];//simulates priority queues, k marker at each step (dim*k)=(16*5)
	__shared__ int heuristics[30]; //will store the heuristics corresponding to a
	__shared__ int sizes [5];//stores current sizes of priority queues
	__shared__ int flag;
	__shared__ int expandedNodes[20];//k*4 for 4 directions
	int extracted;
	flag = 0;
	lowestCost[end] = 0;
	prevNode[end] = -5;
	
		
	//presets shared memory
	if (threadIdx.x == 0)
	{
		for(int i = 0; i < dim * k ;i++)
		{
			array[i] = -1;	
		}
		
		for(int i = 0; i < dim * dim ;i++)
		{
			lowestCost[i] = 2000;//whatever max int is
			
		}
		for(int i = 0; i < k ; i++)
		{
			sizes[i] = 0;
		}
		for(int i =0; i < 4 * k;i++)
		{
			expandedNodes[i] = -1;
		}
	}

	__syncthreads();
	
		

	array[0] = start;
	sizes[0] = 1;
	lowestCost[start] = 0;

	
	while(checkIfQueueEmpty(array, dim, k) != 0)
	{	

		if(flag == 1)
		{
			break;
		}

		__syncthreads();
		
		
		//stuff expandedNodes array
		if(array[dim * threadIdx.x] != -1)//check if corresponding priority queue is empty
		{
			extracted = array[dim * threadIdx.x];//set the extracted Node to the front 										of this P.Queue
			printf("Extracted %d in %d\n", extracted, threadIdx.x);
			shiftqueue(array, dim, threadIdx.x);
			sizes[threadIdx.x]--;
			//extraction and other stuff here
			

			if(extracted == end || flag == 1)
			{
				flag = 1;
				break;
			}
			
			
			//add to expanded nodes here
			int top = extracted-dim;
			int bottom = extracted+dim;
			int left = extracted-1;
			int right = extracted+1;
			if(extracted % dim == 0)
				left = -1;
			if(extracted % dim == dim-1)
				right = -1;
			
			
			//dumps adjacent squares into array
			expandedNodes[threadIdx.x * 4 + 0] = top;
			expandedNodes[threadIdx.x* 4 + 1] = bottom;
			expandedNodes[threadIdx.x* 4 + 2] = left;
			expandedNodes[threadIdx.x* 4 + 3] = right;		
			
			//checks adjacent squares
			//deduplicates list
			for(int i = 0; i < 4 ;i++)
			{
				int curNum = expandedNodes[threadIdx.x * 4 + i];
				
			      if(curNum<0 || curNum>dim*dim|| grid[curNum]==0||lowestCost[extracted]+1>lowestCost[curNum])
				{					
					expandedNodes[threadIdx.x * 4 + i] =- 1;
					continue;
				}//checks for invalid indices
			

				if(expandedNodes[threadIdx.x * 4 + i] != -1)
				{
					
					//route is shorter, therefore update cost and previous node
					lowestCost[curNum] = lowestCost[extracted] + 1;
					prevNode[curNum] = extracted;
				}
				

				
			}
			//printf("After dedup: %d %d %d %d in thread %d\n ",expandedNodes[threadIdx.x*4+0],expandedNodes[threadIdx.x*4+1],expandedNodes[threadIdx.x*4+2],expandedNodes[threadIdx.x*4+3],threadIdx.x);
			
			

			//start heuristic of Nodes not -1 in expandedNodes 
			
			for(int i = 0;i < 4 ; i++)
			{
				int r = duplicateAdjacents(expandedNodes, expandedNodes[threadIdx.x * 4 + i],k,threadIdx.x * 4 + i);
				if(r != -1)
				{
					atomicExch(&expandedNodes[r], -1);
				}
				
				if(expandedNodes[threadIdx.x*4+i] != -1)
				{
					
					int h = lowestCost[expandedNodes[threadIdx.x*4+i]]+				 							heuristic(dim,expandedNodes[threadIdx.x*4+i],end);
					int check = 0; 
					while(check == 0)
					{
						int targetLocation = findNextInsertionPoint(sizes, k);
						if(atomicCAS(&array[targetLocation * dim + sizes[targetLocation]], -1, expandedNodes[threadIdx.x * 4 + i]) == -1)
						{
							
							heuristics[targetLocation * dim + sizes[targetLocation]] = h;
							sizes[targetLocation]++;	
							organizeQueue(array, targetLocation, heuristics, h, sizes[targetLocation], dim);
							check = 1;
							
						}
					}
				
				}
								
			}			
			
				
		}//end of the larger if statement, coded this way to prevent warp divergence
		
		__syncthreads();

	}//end of while loop, shoould have found route or not by here
	if(threadIdx.x == 0){
		printf("prev : \n%d\n", prevNode[end]);
		int current = end;	
		while(current != start)
		{
			if(prevNode[current] == -5)
				break;
			grid[current] = -1;
			current = prevNode[current];
		}
		if(current == start)
		{
			grid[current] = -1;
			printf("The path is marked by -1s\n");
			print_board(grid, dim);
		}
		else
			printf("No route found!\n");
	
	}

}

__device__ void organizeQueue(int * queue, int targetLocation, int * heuristics, int h, int which, int dim)
{
	int temp1;
	int temp2;
	for(int i = which-1; i >= 0 ;i--)
	{	
		if(i < 0)
			break;
		 if(heuristics[dim * targetLocation + i] > h)
		{
			temp1 = heuristics[dim * targetLocation + i];
			temp2 = queue[dim * targetLocation + i];
			heuristics[dim * targetLocation + i] = h;
			queue[dim * targetLocation + i] = queue[dim * targetLocation + which];
			queue[dim * targetLocation + which] = temp2;
			heuristics[dim * targetLocation + which] = temp1;
			which--;
		}

	}
}

__device__ int duplicateAdjacents(int * adjacents, int a, int k, int currentPos)
{
	for(int i = 1; i < 4 * k ;i++)
	{
		if(adjacents[i] == a && i != currentPos)
		{
			return i;
		}
	}
	return -1;
}

__device__ int findNextInsertionPoint(int * sizes, int k){
	
	int smallestQueue = -100;
	int curSmallest = 5000;
	for(int i = 0; i < k ; i++)
	{
		if(sizes[i] < curSmallest)
		{
			smallestQueue = i;
			curSmallest = sizes[i];
		}
	} 
	
	return smallestQueue;

}

__device__ __host__ int heuristic(int dim, int current, int end){//heurstics (Manhattan distance)
	
	int curX = current % dim;
	int curY = current / dim;
	int endX = end % dim;
	int endY = end / dim;
	
	return (int)(fabsf(curX - endX) + fabsf(curY - endY));
	
}

 __device__  void shiftqueue(int* array_sh, int dim, int k)
{
	
	for(int i = 0; i < dim; i++)
	{
		if(i == dim - 1)
			array_sh[dim * k + i] = -1;
		else
			array_sh[dim * k + i] = array_sh[dim * k + i + 1];
		
	}
}

__device__  int checkIfQueueEmpty(int* array_sh, int dim, int k)
{
	int check = 0;//0 if all queues empty, 1 otherwise
	for(int i = 0 ; i < k ; i++)
	{
		if(array_sh[dim * i] != -1){
			check = 1;
			break;
		}
	}
	return check;
}

__device__ __host__ void print_board(int *game, int dim){

	for (int y = 0; y < dim; y++){
                for (int x = 0; x < dim; x++){
			printf("%d ", game[x + dim*y]);
		}
		printf("\n");
	}
}

int main () {
	time_t t;
  	//int num_iter = 10;
	int dim = 6;
	int *grid = (int *) malloc(dim * dim * sizeof(int));
	int *results = (int *) malloc(dim * dim * sizeof(int));
        srand((unsigned) time(&t));
	double p = 0.5;

	for (int i = 0; i < dim*dim; i++){
		grid[i] = p < (double)rand()/(double)(RAND_MAX);
		
	}
	int * dev_grid;
	cudaMalloc((void**)&dev_grid, dim*dim*sizeof(int));
	cudaMemcpy(dev_grid, grid, dim*dim*sizeof(int), cudaMemcpyHostToDevice);

	print_board(grid, dim);	
	//1s are spaces you can walk on
	int startPoint = -1;
	while (startPoint == - 1){
		int rando = rand() % (dim*dim);//dim*dim - 1 is highest number( 0 to dim*dim-1)
		if (grid[rando] == 1)
			startPoint = rando;
	}
	
	int endPoint = -1;
	while (endPoint == - 1){
		int rando = rand() % (dim*dim);
		if (grid[rando] == 1 && rando!=startPoint)
			endPoint = rando;
	}
	
	printf("%d %d\n", startPoint, endPoint);
	float gpu_time;
	gstart();
	traverse<<<1,5>>>(dev_grid, startPoint, endPoint, dim, 5, results);//last value is # of p-queues
	gend(&gpu_time);
	printf("GPU time = %f\n", gpu_time);
	cudaDeviceSynchronize();
	return 0;


}












