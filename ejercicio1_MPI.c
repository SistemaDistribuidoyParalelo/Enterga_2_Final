#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>


//declaracion de variables

void imprimeMatriz(double *S, int tipo_fc) {
// Imprime la matriz pasada por parametro
// N es la cantidad de bloques, r dimension de cada bloque
  printf(" \n");
  printf(" \n");
  if(tipo_fc == 0){
    for(int i=0;i<N;i++){
            for(int j=0;j<N;j++){
                printf("%f  ",S[i+j*N]);
            }
            printf(" \n");
    }
  }else{
      for(int i=0;i<N;i++){
            for(int j=0;j<N;j++){
                printf("%f  ",S[i*N+j]);
            }
            printf(" \n");
    }
  }
}


//MAIN
int main(int argc,char*argv[]){
  int N; // tam de la matriz
  double *A,*B,*C,*D;
  double *L,*U; // matriz triangular L superior U inferior
  double *parcialAB,*parcialLC,*parcialDU; //matrices parciales
  double *M; //resultado final
  double ul,u,l,b; // promedios de las matrices U L y B
  double timetick;
  double sec;
  struct timeval tv;

//DECLARACION DE FUNCIONES UTILIZADAS EN EL PROGRAMA
  int i,j,k;
  if ((argc != 2) || ((N = atoi(argv[1])) <= 0) )
   {
     printf("\nUsar: %s n\n  n: Dimension de la matriz (nxn X nxn)\n", argv[0]);
     exit(1);
   }

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
    //double num = 0;
    //INICIALIZACION DE LAS MATRICES
    gettimeofday(&tv,NULL);
    sec = tv.tv_sec + tv.tv_usec/1000000.0;
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            parcialAB[i*N+j]=0;
            parcialLC[i*N+j]=0;
            parcialDU[i*N+j]=0;
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
    double temp,temp1,temp2,u,l,b;
    temp=0;
    temp1=0;
    temp2=0;
    //HAY Q recorer todo asi que no importa la forma
    for(int i=0;i<N;i++){
          for(int j=0;j<N;j++){
              temp+=B[i*N+j];
            }
    }

    for(i=0;i<N;i++){
        for (j=i;j<N;j++){
            temp1+=U[i+((j*(j+1))/2)];
        }
    }

	for(i=0;i<N;i++){
		for(j=0;j<(i+1);j++)
		{
			temp2+=L[j+((i*(i+1)/2))];
		}
	}

    b=(temp/(N*N));
    u=(temp1/(N*N));
    l=(temp2/(N*N));
    ul=u*l;
    printf("u = %f  l = %f", u,l);
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            for(k = 0;k<N;k++){
                parcialAB[i*N+j] +=A[i*N+k]*B[k*N+j];
            }
        }
    }

    //L es inferior, recorrido parcial
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            for(k = 0;k<i+1;k++){
                parcialLC[i*N+j] +=L[k+((i*(i+1)/2))]*C[j*N+k];
            }
        }
    }

    //U es superior
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            for(k = 0;k<j+1;k++){
                parcialDU[i*N+j] +=D[i*N+k]*U[k+((j*(j+1))/2)];
            }
        }
    }

    //Resultado Final
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            for(k = 0;k<N;k++){
                M[i*N+j] = ul*(parcialAB[i*N+j]+parcialLC[i*N+j]+parcialDU[i*N+j]);
            }
        }
    }

    gettimeofday(&tv,NULL);
    timetick = tv.tv_sec + tv.tv_usec/1000000.0;
    printf("Tiempo en segundos %f\n", timetick - sec);
    imprimeMatriz(M,1);

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
    free(M);
    return(0);
}
