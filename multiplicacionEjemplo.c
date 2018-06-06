#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>

//MAIN
int main(int argc,char*argv[]){
    double *A,*B,*C,*D;
    double *parcialAB,*parcialAB_SUB;
    //DECLARACIONES MPI
    int world_size;
    int ID;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;


    MPI_Status status;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &ID);
    MPI_Get_processor_name(processor_name, &name_len);

    //DECLARACION DE FUNCIONES UTILIZADAS EN EL PROGRAMA
    int i,j,k;
    int N = 4; //TAMANIO DE LA MATRIZ
    int T = 2; //NRO PROCESADORES
    A=(double*)malloc(sizeof(double)*N*N);
    B=(double*)malloc(sizeof(double)*N*N);
    C=(double*)malloc(sizeof(double)*N*N);
    D=(double*)malloc(sizeof(double)*N*N);
    parcialAB=(double*)malloc(sizeof(double)*N*N);
    parcialAB_SUB=(double*)malloc(sizeof(double)*N*N);

    printf("Esta funcion es la %d \n",ID);
    for(i=0;i<N;i++){
            for(j=0;j<N;j++){
                parcialAB[i*N+j]=0;
                parcialAB_SUB[i*N+j]=0;
                A[i*N+j] = 1;
                B[j*N+i] = 1;//ORDENXCOLUMNAS
                C[j*N+i]=1; //ORDENXCOLUMNAS
                D[i*N+j]=1; //ORDENXFILAS
            }
    }

    MPI_Scatter(A,(N*N)/T, MPI_INT, C, (N*N)/T, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(parcialAB,(N*N)/T, MPI_INT, parcialAB_SUB, (N*N)/T, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(B,N,MPI_DOUBLE,0,MPI_COMM_WORLD); // Comunicador utilizado (En este caso, el global)

    int inicio = ((N)/T) * (ID);
    int parcial = ((N)/T) * (ID+1);

        for(int i=inicio;i<parcial;i++){
            for(int j=0;j<N;j++){
                for(int k = 0;k<N;k++){
                    parcialAB_SUB[i*N+j] += C[i*N+k]*B[k*N+j];
                }
            }
        }

    MPI_Gather(parcialAB_SUB,(N*N)/T,MPI_DOUBLE,parcialAB,(N*N)/T,MPI_DOUBLE,0,MPI_COMM_WORLD); 
    
    MPI_Scatter(A,(N*N)/T, MPI_INT, C, (N*N)/T, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(parcialAB,(N*N)/T, MPI_INT, parcialAB_SUB, (N*N)/T, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(B,N,MPI_DOUBLE,0,MPI_COMM_WORLD); // Comunicador utilizado (En este caso, el global)


        for(int i=inicio;i<parcial;i++){
            for(int j=0;j<N;j++){
                for(int k = 0;k<N;k++){
                    parcialAB_SUB[i*N+j] += C[i*N+k]*B[k*N+j];
                }
            }
        }

    MPI_Gather(parcialAB_SUB,(N*N)/T,MPI_DOUBLE,parcialAB,(N*N)/T,MPI_DOUBLE,0,MPI_COMM_WORLD); 
    
MPI_Scatter(A,(N*N)/T, MPI_INT, C, (N*N)/T, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(parcialAB,(N*N)/T, MPI_INT, parcialAB_SUB, (N*N)/T, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(B,N,MPI_DOUBLE,0,MPI_COMM_WORLD); // Comunicador utilizado (En este caso, el global)

  

        for(int i=inicio;i<parcial;i++){
            for(int j=0;j<N;j++){
                for(int k = 0;k<N;k++){
                    parcialAB_SUB[i*N+j] += C[i*N+k]*B[k*N+j];
                }
            }
        }

    MPI_Gather(parcialAB_SUB,(N*N)/T,MPI_DOUBLE,parcialAB,(N*N)/T,MPI_DOUBLE,0,MPI_COMM_WORLD); 
        
printf(" \n");
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            printf("%f  ",parcialAB_SUB[i*N+j]);
        }
        printf(" \n");
    }
    

    free(A);
    free(B);
    free(C);
    free(D);
    free(parcialAB);
    free(parcialAB_SUB);
    MPI_Finalize();
}
