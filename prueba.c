#include <mpi.h>
#include <stdio.h>
#include <omp.h>
int main(int argc, char* argv[]){
    int rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    int N = 4;
    int T = 2;
    int i;
    int k;
    double temp2;   
    #pragma omp parallel
    {
        #pragma omp for nowait schedule(dynamic,T)
        for(i=0;i < N/T;i++){
            for(k=0;k<i+5+1;k++){
                temp2++;
                printf("Hello MPI rank %d thread %d\n",rank,omp_get_thread_num());
            }     
        }
    }
    MPI_Finalize();
return(0);
}