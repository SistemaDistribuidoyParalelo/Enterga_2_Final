#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>

double dwalltime(){
        double sec;
        struct timeval tv;
        int tiempo = gettimeofday(&tv,NULL);
        sec = tv.tv_sec + tv.tv_usec/1000000.0;
        return sec;
}

//MAIN
int main(int argc,char*argv[]){
double *A,*B,*C,*D,*L,*U; // Matrices iniciales
double *parcialAB,*parcialLC,*parcialDU,*parcialM; //matrices resultados parciales
double *parcialAB_SUB,*parcialLC_SUB,*parcialDU_SUB; //matrices distribuidos parciales
double *pruebaA,*pruebaB,*pruebaL,*pruebaD;// temporales partidas de cada matriz
double *M; //resultado final
double sumTemp,sumTemp2,temp,temp1,temp2,u,l,b,ul;
int i,j,k;
int N; //TAMANIO DE LA MATRIZ
int T; //NROPROCESADORES

// VARIABLES DE TIEMPO
double startComunication, cuentas;
//DECLARACION DE FUNCIONES UTILIZADAS EN EL PROGRAMA

if ((argc != 2))
   {
     printf("\nUsar: %s n\n  n: Dimension de la matriz (nxn X nxn)\n", argv[0]);
     exit(1);
   }else{
       N=atoi(argv[1]);
   }

//DECLARACIONES MPI
int world_size;
int ID;
char processor_name[MPI_MAX_PROCESSOR_NAME];
int name_len;

//INICIO VARIABLES
temp=0;
temp1=0;
temp2=0;

MPI_Status status;
MPI_Init(NULL, NULL);
MPI_Comm_size(MPI_COMM_WORLD, &world_size);
MPI_Comm_rank(MPI_COMM_WORLD, &ID);
MPI_Get_processor_name(processor_name, &name_len);

T = world_size;
int inicio = ((N)/T) * (ID);


A=(double*)malloc(sizeof(double)*N*N);
B=(double*)malloc(sizeof(double)*N*N);
C=(double*)malloc(sizeof(double)*N*N);
D=(double*)malloc(sizeof(double)*N*N);
U=(double*)malloc(sizeof(double)*((N*(N+1))/2)); //superior
L=(double*)malloc(sizeof(double)*N*N); //inferior
M=(double*)malloc(sizeof(double)*N*N);

pruebaA=(double*)malloc(sizeof(double)*N*N);
pruebaB=(double*)malloc(sizeof(double)*N*N);
pruebaD=(double*)malloc(sizeof(double)*N*N);
pruebaL=(double*)malloc(sizeof(double)*N*N); //inferior


parcialAB=(double*)malloc(sizeof(double)*N*N);
parcialLC=(double*)malloc(sizeof(double)*N*N);
parcialDU=(double*)malloc(sizeof(double)*N*N);
parcialM=(double*)malloc(sizeof(double)*N*N);

parcialAB_SUB=(double*)malloc(sizeof(double)*N*N);
parcialLC_SUB=(double*)malloc(sizeof(double)*N*N);
parcialDU_SUB=(double*)malloc(sizeof(double)*N*N);

//printf("Esta funcion es la %d ",ID);
printf(" \n");
for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            A[i*N+j] = 1;
            B[j*N+i] = 1;//ORDENXCOLUMNAS
            C[j*N+i] = 1;
            D[i*N+j]=1; //ORDENXFILAS
            M[i*N+j]=0; //ORDENXFILAS

            parcialLC[i*N+j]=0;
            parcialAB[i*N+j]=0;
            parcialDU[i*N+j]=0;
            parcialM[i*N+j]=0;

            parcialAB_SUB[i*N+j]=0;
            parcialLC_SUB[i*N+j]=0;
            parcialDU_SUB[i*N+j]=0;


            pruebaA[j*N+i]=1; //ORDENXCOLUMNAS
            pruebaB[i*N+j]=1; //ORDENXFILAS
            pruebaL[j*N+i]=1; //ORDENXCOLUMNAS
            pruebaD[j*N+i]=1; //ORDENXCOLUMNAS

            if (i > j){
                L[i*N+j] = 1; //Almaceno L por filas
            }
            else if (i < j){
                U[i+((j*(j+1))/2)] = 1; //Almaceno U por columnas
                L[i*N+j] = 0;
            }
            else{
                L[i*N+j] = 1;
                U[i+((j*(j+1))/2)] = 1;
            }
        }
}

 
//printf("u = %f  l = %f \n", u,l);
startComunication = dwalltime();

cuentas = dwalltime();
//COMUNICACION
MPI_Scatter(A,(N*N)/T, MPI_DOUBLE, pruebaA, (N*N)/T, MPI_DOUBLE, 0, MPI_COMM_WORLD);
MPI_Scatter(L,(N*N)/T, MPI_DOUBLE, pruebaL, (N*N)/T, MPI_DOUBLE, 0, MPI_COMM_WORLD);
MPI_Scatter(D,(N*N)/T, MPI_DOUBLE, pruebaD, (N*N)/T, MPI_DOUBLE, 0, MPI_COMM_WORLD);
MPI_Bcast(B,N*N,MPI_DOUBLE,0,MPI_COMM_WORLD); // Comunicador utilizado (En este caso, el global)
MPI_Bcast(C,N*N,MPI_DOUBLE,0,MPI_COMM_WORLD); // Comunicador utilizado (En este caso, el global)
MPI_Bcast(U,((N*(N+1))/2),MPI_DOUBLE,0,MPI_COMM_WORLD);
//TOMO EL TIEMPO DE INICIO
// SACO PROMEDIOS QUE NECESITO
    // PROMEDIO b

//HAY Q recorer todo asi que no importa la forma

//printf("%f \n", sumTemp);
for(i=0;i<N/T;i++){
    for (j=i;j<N;j++){
        temp1+=U[i+((j*(j+1))/2)];
    }
}
//printf("%f\n", temp1);

for(i=0;i < N/T;i++){
    for(k=0;k<i+((N/T)*(ID))+1;k++){
        temp2+=pruebaL[i*N+k];
	}
}
MPI_Allreduce(&temp1,&sumTemp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
MPI_Allreduce(&temp2,&sumTemp2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
//printf("%f\n", sumTemp2);

u=(sumTemp/(N*N));
l=(sumTemp2/(N*N));
ul=u*l;
//MULTIPLICACION  A*B

    for(i=0;i<N/T;i++){
        for(j=0;j<N;j++){
            for(k = 0;k<N;k++){
                parcialAB[i*N+j] += pruebaA[i*N+k]*B[k*N+j];
            }
        }
    }

/*if(ID==0){
for( i=0;i<N;i++){
      for(int j=0;j<N;j++){
          printf("%f  ",parcialAB[i*N+j]);
      }
      printf(" \n");
}
}
      printf(" \n");
*/
//MULTIPLICACION L*C



//L es inferior, recorrido parcial
for(i=0;i<N/T;i++){
    for(j=0;j<N;j++){
        for(k = 0;k<i+inicio+1;k++){
            parcialLC[i*N+j] +=pruebaL[i*N+k]*C[j*N+k];
        }
    }
}

/*
if(ID==0){
for( i=0;i<N;i++){
      for(int j=0;j<N;j++){
          printf("%f  ",parcialLC[i*N+j]);
      }
      printf(" \n");
}
}*/
/*for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
          printf("%f  ",parcialAB_SUB[i*N+j]);
      }
      printf(" \n");
}*/


//MULTIPLICACION D*U


//U es superior, recorrido parcial
for(i=0;i<N/T;i++){
    for(j=0;j<N;j++){
        for(k = 0;k<j+1;k++){
            parcialDU[i*N+j] +=pruebaD[i*N+k]*U[k+((j*(j+1))/2)];
        }
    }
}


/*if (ID==0) {
  for( i=0;i<N;i++){
        for(int j=0;j<N;j++){
            printf("%f  ",parcialDU[i*N+j]);
        }
        printf(" \n");
  }}*/
  /*for( i=0;i<N;i++){
      for(int j=0;j<N;j++){
          printf("%f  ",parcialDU_SUB[i*N+j]);
      }
      printf(" \n");
}*/





for(i=0;i<N;i++){
    for(j=0;j<N;j++){
        for(k=0;k<N;k++){
            parcialM[i*N+j] = ul*(parcialAB[i*N+j]+parcialLC[i*N+j]+parcialDU[i*N+j]);
        }
    }
}

MPI_Gather(parcialM,(N*N)/T,MPI_DOUBLE,M,(N*N)/T,MPI_DOUBLE,0,MPI_COMM_WORLD);

if (ID==0){
   printf("Tiempo en segundos %f\n", cuentas - dwalltime());
}



if (ID==0){
  for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            printf("%f  ",M[i*N+j]);
        }
        printf(" \n");
  }}

free(A);
free(B);
free(C);
free(D);
free(U);
free(L);
free(M);
free(parcialAB);
free(parcialLC);
free(parcialDU);
free(parcialM);
free(parcialAB_SUB);
free(parcialDU_SUB);
free(parcialLC_SUB);
free(pruebaA);
free(pruebaB);
free(pruebaL);
free(pruebaD);

MPI_Finalize();
return 0;
}