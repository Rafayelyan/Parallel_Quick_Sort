#include <iostream>
#include <cmath>
#include <mpi.h>

int check(int number, int i);

int main(int argc, char *argv[])
{
	int tid,nthreads;

	MPI_Init(&argc,&argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &tid);

	MPI_Comm_size(MPI_COMM_WORLD, &nthreads);
	
	int 		n;
	int 		*x;
	int 		stride;  
	int 		*send_counts;
	int		*displs;
	int		*block;
	int 		N;
	int 		delta;
	int 		pivot;
	int		*left_buffer;
	int		*right_buffer;
	int		left_counter;
	int		right_counter;
	MPI_Status	status;

	n = 8;
	N = log2(nthreads - 1) + 1;
	delta = nthreads / 2;
	left_counter = 0;
	right_counter = 0;
	stride = n / nthreads;
	x = new int[n];
	send_counts = new int[nthreads];
	displs = new int[nthreads];
	block = new int[n / nthreads];
	left_buffer = new int[n];
	right_buffer = new int[n];

	for(int i=0; i < nthreads; ++i)
	{
		displs[i] = i * stride;
		send_counts[i] = n / nthreads;
	}

	if (tid == 0)
	{
		x[0] = 1;
		x[1] = 4;
		x[2] = 2;
		x[3] = 5;
		x[4] = 8;
		x[5] = 6;
		x[6] = 9;
		x[7] = 3;
	}
	MPI_Scatterv(&x[0], send_counts, displs, MPI_INT, &block[0], n / nthreads, MPI_INT, 0, MPI_COMM_WORLD);

	if(tid == 0)
	{
		pivot = 0;
		for(int i=0;i < n / nthreads; i++)
		{
			pivot += block[i];
		}		
		pivot = pivot / (n / nthreads);
	}
	MPI_Bcast(&pivot, 1, MPI_INT, 0, MPI_COMM_WORLD);

	for(int count =  log2( nthreads - 1 ) + 1; count > 0; count--)
        {
		if(check(tid, N) &&  (tid + delta) <= (nthreads - 1) )
		{
			std::cout<<"tid : "<< tid <<" pivot : "<<pivot<<" N : "<<N<<'\n';
			for (int i=0; i<(n / nthreads); ++i)
			{
				std::cout<<block[i]<<'\n';
			}
			right_counter = 0;
			int *l_left_buffer = new int[n];
			int l_left_counter = 0;

			for (int i=0; i<(n / nthreads); ++i)
			{
				if (block[i] > pivot)
				{
					right_buffer[right_counter++] = block[i];	
				}
				else
				{
					l_left_buffer[l_left_counter++] = block[i];
				}
			}
			MPI_Send(&right_counter, 1, MPI_INT, tid + delta, tid + delta, MPI_COMM_WORLD);
			if (right_counter > 0)
			{
				MPI_Send(right_buffer, right_counter, MPI_INT, tid + delta, tid + delta, MPI_COMM_WORLD);
			}
			
			MPI_Recv(&left_counter, 1, MPI_INT, tid + delta, tid, MPI_COMM_WORLD, &status);
			
			int index = 0;
			if (l_left_counter > 0)
			{
				for(index=0; index<l_left_counter; ++index)
				{
					block[index] = l_left_buffer[index];
				}
			}

			else if(!( check(tid, N) ) &&  (tid - delta) >= 0 )
			{
				MPI_Recv(left_buffer, left_counter, MPI_INT, tid + delta, tid, MPI_COMM_WORLD, &status);

				for(int i=0; i<left_counter; ++i)
				{
					block[index++] = left_buffer[i];
				}
			}
			n = left_counter + l_left_counter; // kaskacner kan
			
			std::cout<<"tid\t:\t"<<tid<<std::endl;
			for (int i=0; i<n; ++i)
			{
				std::cout<<'\t'<<block[i]<<'\t';
			}
			std::cout<<'\n';
		}
		else if (!( check(tid, N) ) &&  (tid - delta) >= 0 )
		{
			std::cout<<"tid : "<< tid <<" pivot : "<<pivot<<" N : "<<N<<'\n';
			for (int i=0; i<(n / nthreads); ++i)
			{
				std::cout<<block[i]<<'\n';
			}
			left_counter = 0;
			int *l_right_buffer = new int[n];
			int l_right_counter = 0;

			for (int i=0; i<(n / nthreads); ++i)
			{
				if (block[i] <= pivot)
				{
					left_buffer[left_counter++] = block[i];	
				}
				else
				{
					l_right_buffer[l_right_counter++] = block[i];
				}
			}


			int index = 0;
			if (l_right_counter > 0)
			{
				for (index=0; index<l_right_counter; ++index)
				{
					block[index] = l_right_buffer[index];
				}
			}
			MPI_Recv(&right_counter, 1, MPI_INT, tid - delta, tid, MPI_COMM_WORLD, &status);
			MPI_Recv(right_buffer, right_counter, MPI_INT, tid - delta, tid, MPI_COMM_WORLD, &status);
			std::cout<<"tid : "<<tid<<" recived right_counter : "<<right_counter<<std::endl;
			if(right_counter > 0)
			{
				for(int i=0; i<right_counter; ++i)
				{
					block[index++] = right_buffer[i];
				}
			}

			n = right_counter + l_right_counter; // kaskacner kan

			MPI_Send(&left_counter, 1, MPI_INT, tid - delta, tid - delta, MPI_COMM_WORLD);
			if (left_counter > 0)
			{
				MPI_Send(left_buffer, left_counter, MPI_INT, tid - delta, tid - delta, MPI_COMM_WORLD);
			}

			std::cout<<"tid\t:\t"<<tid<<std::endl;
			for (int i=0; i<n; ++i)
			{
				std::cout<<'\t'<<block[i]<<'\t';
			}
			std::cout<<'\n';
		}

		MPI_Barrier(MPI_COMM_WORLD);
                --delta;
                --N;
        }


	
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
