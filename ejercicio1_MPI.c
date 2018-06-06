#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>

//MAIN
int main(int argc,char*argv[]){
  //int N; // tam de la matriz
  double *A,*B,*C,*D;
  double *L,*U; // matriz triangular L superior U inferior
  double *parcialAB,*parcialLC,*parcialDU; //matrices parciales
  double *parcialAB_SUB,*parcialLC_SUB,*parcialDU_SUB;
  double *parcialA,*parcialL,*parcialD,*parcialM;
  double *M; //resultado final
  double ul,u,l,b; // promedios de las matrices U L y B
  
//DECLARACION DE FUNCIONES UTILIZADAS EN EL PROGRAMA
  int i,j,k;

  int N = 4;
  int T = 2;
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
  
   //ALOCACION DE MEMORIA PARA LAS MATRICES
    A=(double*)malloc(sizeof(double)*N*N);
    B=(double*)malloc(sizeof(double)*N*N);
    C=(double*)malloc(sizeof(double)*N*N);
    D=(double*)malloc(sizeof(double)*N*N);
    U=(double*)malloc(sizeof(double)*((N*(N+1))/2)); //superior
    L=(double*)malloc(sizeof(double)*((N*(N+1))/2)); //inferior
    M=(double*)malloc(sizeof(double)*N*N);

    parcialAB =(double*)malloc(sizeof(double)*N*N);
    parcialLC =(double*)malloc(sizeof(double)*N*N);
    parcialDU =(double*)malloc(sizeof(double)*N*N);


    parcialM =(double*)malloc(sizeof(double)*N*N);
    parcialA =(double*)malloc(sizeof(double)*N*N);
    parcialL =(double*)malloc(sizeof(double)*N*N);
    parcialD =(double*)malloc(sizeof(double)*N*N);
    
    parcialAB_SUB =(double*)malloc(sizeof(double)*N*N);
    parcialLC_SUB =(double*)malloc(sizeof(double)*N*N);
    parcialDU_SUB =(double*)malloc(sizeof(double)*N*N);
    
    //double num = 0;
    //INICIALIZACION DE LAS MATRICES
  //  gettimeofday(&tv,NULL);
    //sec = tv.tv_sec + tv.tv_usec/1000000.0;
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            parcialAB[i*N+j]=0;
            parcialLC[i*N+j]=0;
            parcialDU[i*N+j]=0;
            parcialA[i*N+j]=0;
            parcialL[i*N+j]=0;
            parcialM[i*N+j]=0;
            parcialD[i*N+j]=0;
            A[i*N+j] = 1;
            B[j*N+i] = 1;//ORDENXCOLUMNAS
            C[j*N+i]=1; //ORDENXCOLUMNAS
            D[i*N+j]=1; //ORDENXFILAS
            //Matrices triangulares, PENSAR EN ARMAR LA TRANSPUESTA
            if (i > j){
                L[j+((i*(i+1))/2)] = 1; //Almaceno L por filas
            }
            else if (i < j){
                U[i+((j*(j+1))/2)] = 1; //Almaceno U por columnas
            }
            else{
                L[j+((i*(i+1))/2)] = 1;
                U[i+((j*(j+1))/2)] = 1;
            }
        }
    }

    
    // SACO PROMEDIOS QUE NECESITO
    // PROMEDIO b
    double temp,temp1,temp2;
    temp=0;
    temp1=0;
    temp2=0;

    //HAY Q recorer todo asi que no importa la forma
    
    for(int i=0;i<N;i++){
          for(int j=0;j<N;j++){
              temp+=B[i*N+j];
            }
    }
    for(int i=0;i<N;i++){
        for(int j=i;j<N;j++){
            temp1+=U[i+((j*(j+1))/2)];
        }
    }

	for(int i=0;i<N;i++){
		for(int j=0;j<(i+1);j++)
		{
			temp2+=L[j+((i*(i+1)/2))];
		}
	}

    b=(temp/(N*N));
    u=(temp1/(N*N));
    l=(temp2/(N*N));
    ul=u*l;
    printf("u = %f  l = %f \n", u,l);


    int inicio = ((N)/T) * (ID);
    int parcial = ((N)/T) * (ID+1);

    //--------------------------------------------
    MPI_Scatter(A,(N*N)/T, MPI_DOUBLE, parcialA, (N*N)/T, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(parcialAB,(N*N)/T, MPI_DOUBLE, parcialAB_SUB, (N*N)/T, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(B,N,MPI_DOUBLE,0,MPI_COMM_WORLD); // Comunicador utilizado (En este caso, el global)

    for(int i=inicio;i<parcial;i++){
        for(int j=0;j<N;j++){
            for(k = 0;k<N;k++){
                parcialAB_SUB[i*N+j] += parcialA[i*N+k]*B[k*N+j];
            }
        }
    }

    MPI_Gather(parcialAB_SUB,(N*N)/T,MPI_DOUBLE,parcialAB,(N*N)/T,MPI_DOUBLE,0,MPI_COMM_WORLD); 
    //--------------------------------------------

    //---------------------------------------------
    MPI_Scatter(L,(((N*(N+1))/2))/T, MPI_DOUBLE, parcialL, (((N*(N+1))/2))/T, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(parcialLC,(N*N)/T, MPI_DOUBLE, parcialLC_SUB, (N*N)/T, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(C,N,MPI_DOUBLE,0,MPI_COMM_WORLD); // Comunicador utilizado (En este caso, el global)

    //L es inferior, recorrido parcial
    for(int i=inicio;i<parcial;i++){
        for(int j=0;j<N;j++){
            for(k = 0;k<i+1;k++){
                parcialLC_SUB[i*N+j] +=parcialL[k+((i*(i+1)/2))]*C[j*N+k];
            }
        }
    }

    MPI_Gather(parcialLC_SUB,(N*N)/T,MPI_DOUBLE,parcialLC,(N*N)/T,MPI_DOUBLE,0,MPI_COMM_WORLD);
    //-----------------------------------

    MPI_Scatter(D,(N*N)/T, MPI_DOUBLE, parcialD, (N*N)/T, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(parcialDU,(N*N)/T, MPI_DOUBLE, parcialDU_SUB, (N*N)/T, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(U,N,MPI_DOUBLE,0,MPI_COMM_WORLD); // Comunicador utilizado (En este caso, el global)

    //U es superior
    for(int i=inicio;i<parcial;i++){
        for(int j=0;j<N;j++){
            for(k = 0;k
            
            <j+1;k++){
                parcialDU_SUB[i*N+j] +=parcialD[i*N+k]*U[k+((j*(j+1))/2)];
            }
        }
    }
    MPI_Gather(parcialDU_SUB,(N*N)/T,MPI_DOUBLE,parcialDU,(N*N)/T,MPI_DOUBLE,0,MPI_COMM_WORLD); 
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Scatter(parcialDU,(N*N)/T, MPI_DOUBLE, parcialDU_SUB, (N*N)/T, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(parcialLC,(N*N)/T, MPI_DOUBLE, parcialLC_SUB, (N*N)/T, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(parcialAB,(N*N)/T, MPI_DOUBLE, parcialAB_SUB, (N*N)/T, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(M,(N*N)/T, MPI_DOUBLE, parcialM, (N*N)/T, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            for(k = 0;k<N;k++){
                parcialM[i*N+j] = ul*(parcialAB_SUB[i*N+j]+parcialLC_SUB[i*N+j]+parcialDU_SUB[i*N+j]);
            }
        }
    }

    MPI_Gather(M,(N*N)/T,MPI_DOUBLE,parcialM,(N*N)/T,MPI_DOUBLE,0,MPI_COMM_WORLD);


  /*gettimeofday(&tv,NULL);
    timetick = tv.tv_sec + tv.tv_usec/1000000.0;
    printf("Tiempo en segundos %f\n", timetick - sec);
*/
    /*
    for (i = 0; i < N; i++){
        for (j = 0; j < i+1; j++)
        {
            printf("U[%d] = %f\n", j+((i*(i+1))/2), U[j+((i*(i+1))/2)]);
        }
    }
        printf("\n");
    for (i = 0; i < N; i++){
        for (j = 0; j < i+1; j++)
        {
            printf("L[%d] = %f\n", j+((i*(i+1))/2), L[j+((i*(i+1))/2)]);
        }
    }
    
    printf("\n");
    imprimeMatriz(parcialAB,1);
    imprimeMatriz(parcialLC,1);
    imprimeMatriz(parcialDU,1);
    */
    free(A);
    free(B);
    free(C);
    free(D);
    free(U);
    free(L);
    free(parcialAB);
    free(parcialLC);
    free(parcialDU);
    MPI_Barrier(MPI_COMM_WORLD);

    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            printf("%f  ",M[i*N+j]);
        }
        printf(" \n");
    }
    printf(" \n");

    
    free(M);
    MPI_Finalize();
}
