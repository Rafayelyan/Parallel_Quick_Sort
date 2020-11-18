#include <iostream>
#include <cmath>
#include <mpi.h>

int check(int number, int i);

int main(int argc, char *argv[])
{
	int n;
	n = 8;
	
	int x[8] = {1,4,2,5,8,6,9,3};
	int *block = new int[n];

	int tid,nthreads;

	MPI_Init(&argc,&argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &tid);

	MPI_Comm_size(MPI_COMM_WORLD, &nthreads);
	
	int N = log2(nthreads - 1) + 1;
	int delta = nthreads / 2;

	MPI_Scatter(&x[0], n / nthreads, MPI_INT, &block[0], n / nthreads, MPI_INT, 0, MPI_COMM_WORLD);
	int recive_pivot;
	int pivot;
	if(tid == 0)
	{
		//MPI_Gather(&x, n/nthreads, MPI_INT, &block, n/nthreads, MPI_INT, 0, MPI_COMM_WORLD);
		//MPI_Recv(&block, n / nthreads, MPI_INT,0,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		pivot = 0;
		for(int i=0;i < n / nthreads; i++)
		{
			pivot += x[i];
		}		
		pivot = pivot / (n / nthreads);
		int pivotArr[nthreads];
		for(int i=0; i<nthreads; i++)	
		{
			pivotArr[i] = pivot;
		}
		MPI_Scatter(&pivot, 1, MPI_INT, &recive_pivot, 1, MPI_INT, tid, MPI_COMM_WORLD);
	}


	int checkNumber = 0b10000000;
	//for(int i =  log2( nthreads - 1 ) + 1; i > 0; i--)

	MPI_Allgather(&pivot, 1, MPI_INT, &recive_pivot, 1, MPI_INT, MPI_COMM_WORLD);

	int leftArraySizes[nthreads];
	int rightArraySizes[nthreads];
	for(int i=0; i < nthreads; i++)
	{
		leftArraySizes[i] = -1;
		rightArraySizes[i] = -1;
	}
	int leftSize, rightSize;

	if(check(tid, N) &&  (tid + delta) <= (nthreads - 1) )
	{			
		for(int i=0;i<n/nthreads;i++)
		{
			std::cout<<block[i]<<'\t';
		}
		std::cout<<std::endl;
		int leftCounter = -1, rightCounter = -1;
		int leftArray[n];
		int leftArray1[n];
		int rightArray[n];
		int rightArray1[n];
		
		for(int i=0;i<n/nthreads;i++)
		{
			if( block[i] <= pivot )
			{
				leftArray[++leftCounter] = block[i];
			}
			else
				rightArray[++rightCounter] = block[i];
		}

		rightArraySizes[tid] = rightCounter;

		MPI_Send(&rightArray, rightCounter, MPI_INT, tid + delta, MPI_ANY_TAG,  MPI_COMM_WORLD);
		
		MPI_Scatter(&rightArraySizes, nthreads, MPI_INT, &rightSize, 1, MPI_INT, tid, MPI_COMM_WORLD);
		MPI_Status* status;

		MPI_Recv(&leftArray1, leftSize, MPI_INT, tid + delta, MPI_ANY_TAG, MPI_COMM_WORLD, status);
		
		for( int i=0; i < leftCounter; i++ )
		{
			block[i] = leftArray[i];
		}
		
		for( int i=leftCounter; i < leftSize; i++ )
		{
			block[i] = leftArray1[i];
		}
	}


	if(!( check(tid, N) ) &&  (tid - delta) >= 0 )
	{			
		int leftCounter = -1, rightCounter = -1;
		int leftArray[n];
		int leftArray1[n];
		int rightArray[n];
		int rightArray1[n];
		
		for(int i=0;i<n/nthreads;i++)
		{
			if( block[i] <= pivot )
			{
				leftArray[++leftCounter] = block[i];
			}
			else
		int rightArray[n];
				rightArray[++rightCounter] = block[i];
		}
		
		leftArraySizes[tid] = leftCounter;

		MPI_Send(&leftArray, leftCounter, MPI_INT, tid - delta, MPI_ANY_TAG, MPI_COMM_WORLD);

		MPI_Scatter(&leftArraySizes, nthreads, MPI_INT, &leftSize, 1, MPI_INT, tid, MPI_COMM_WORLD);

		MPI_Status status;
		MPI_Recv(&rightArray1, rightSize, MPI_INT, tid - delta, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		
		for( int i=0; i < rightCounter; i++ )
		{
			block[i] = rightArray[i];
		}
		
		for( int i=rightCounter; i < leftSize; i++ )
		{
			block[i] = rightArray1[i];
		}
	}



	/*
	 * 0 00
	 * 1 01
	 * 2 10 
	 * 3
	 */



	for(int i=0;i<n;i++)
	{
		std::cout<<block[i]<<'\t';
	}
	std::cout<<'\n';








	//MPI_Gather(&x1[0], n/nthreads, MPI_INT, &x[0], n/nthreads, MPI_INT, 0, MPI_COMM_WORLD);

	
	MPI_Finalize();
	return 0;
}

int check(int number, int i)
{
        int index = 0b10000000;
        //std::cout<<"number : "<<number<<'\n';
        //std::cout<<"index : "<<index<<'\n';
        index = index >> (8 - i);
        if(! ( number & index ) )
                return 1;
	else
		return 0;
}

