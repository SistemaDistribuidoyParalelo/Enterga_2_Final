#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "mpi.h"
#include <omp.h>

//Para calcular tiempo
double dwalltime(){
	double sec;
	struct timeval tv;

	gettimeofday(&tv,NULL);
	sec = tv.tv_sec + tv.tv_usec/1000000.0;
	return sec;
}

int main(int argc, char** argv){
	int miID; int cantidadDeProcesos;
	int N; // Dimension de la matriz
	int sizeMatrix; // Cantidad total de datos matriz 	
	int sizePart; // Cantidad de elementos por proceso
	int sizeTrian;
	double *A_buf, *L_buf, *D_buf,*ab_temp,*lc_temp,*du_temp;
	double *A,*B,*C,*D,*L,*U,*M;
	double u=0.0,l=0.0;
	double tiempo_paral, tiempo_balance;


	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&miID);
	MPI_Comm_size(MPI_COMM_WORLD,&cantidadDeProcesos);

	/***** PROGRAMA *****/
	if (argc < 2){
		printf("\n Falta un parametro ");
		printf("\n 1. Dimension de la matriz ");
		return 0;
	}

	N = atoi(argv[1]);
 	sizeMatrix=N*N; // Cantidad total de datos matriz
 	sizePart=sizeMatrix/cantidadDeProcesos; // Cantidad de elementos por bloque x Cantidad de bloques por proceso 	
	sizeTrian = (N+1)*N/2;
	int posicion = (N/cantidadDeProcesos)*miID;

	int CANT_THREADS = atoi(argv[2]);
    omp_set_num_threads(CANT_THREADS);

 	A_buf=(double*)malloc(sizeof(double)*sizePart);
	L_buf=(double*)malloc(sizeof(double)*sizePart);
	D_buf=(double*)malloc(sizeof(double)*sizePart);

	ab_temp=(double*)malloc(sizeof(double)*sizePart);
	lc_temp=(double*)malloc(sizeof(double)*sizePart);
	du_temp=(double*)malloc(sizeof(double)*sizePart);
	
	B=(double*)malloc(sizeof(double)*sizeMatrix);
	C=(double*)malloc(sizeof(double)*sizeMatrix);
	U=(double*)malloc(sizeof(double)*sizeTrian);      

	if(miID==0) { // El proceso con ID=0 inicializa y distribuye los datos
	
		A=(double*)malloc(sizeof(double)*sizeMatrix);
		L=(double*)malloc(sizeof(double)*sizeMatrix);
		D=(double*)malloc(sizeof(double)*sizeMatrix);
		M=(double*)malloc(sizeof(double)*sizeMatrix);

	  	for(int i=0;i<N;i++){
     		for(int j=0;j<N;j++){
     			//Ordenadas por "fila"
       			A[i*N+j]=1;
       			D[i*N+j]=1;
       			//inicializa una matriz triangular inferior
	   			if(i>=j){
	   				L[i*N+j]=1;
				}else{
					L[i*N+j]=0;	
				}
       			//Ordenadas por "columna"
       			B[i+N*j]=1;
       			C[i+N*j]=1;
	   			//inicializa una matriz triangular superior
		   		U[i+j*(j+1)/2]=1;
     		}
   		}
	}

 	tiempo_paral = dwalltime();

	MPI_Scatter(A,sizePart,MPI_DOUBLE,A_buf,sizePart,MPI_DOUBLE,0,MPI_COMM_WORLD);

	MPI_Bcast(B, sizeMatrix, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Scatter(L,sizePart,MPI_DOUBLE,L_buf,sizePart,MPI_DOUBLE,0,MPI_COMM_WORLD);

	MPI_Bcast(C, sizeMatrix, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Scatter(D,sizePart,MPI_DOUBLE,D_buf,sizePart,MPI_DOUBLE,0,MPI_COMM_WORLD);

	MPI_Bcast(U, sizeTrian, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	tiempo_balance = dwalltime();

	// Promedio de U
	double temp=0;
	
	#pragma omp parallel for ordered reduction(+ : temp) schedule(static)
	for(int i=0;i<sizeTrian;i++)
		temp+=U[i];

	temp/=sizeMatrix;	

	MPI_Allreduce(&temp, &u, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	//printf("Proceso %d, Promedio u = %lf \n",miID,u);

	// Promedio de L
	temp=0;
	#pragma omp parallel for ordered reduction(+ : temp) schedule(static)
	for(int i=0;i<sizePart/N;i++){
		for (int j = 0; j < posicion+i+1; j++){
			temp+=L_buf[i*N+j];
		}
	}

	temp/=sizeMatrix;

	MPI_Allreduce(&temp, &l, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	//printf("Proceso %d, Promedio l = %lf \n",miID,l);

	#pragma omp parallel
	{
	#pragma omp for
	for(int i=0;i<sizePart/N;i++){
		for(int j=0;j<N;j++){
			ab_temp[i*N+j]=0;
			for(int k=0;k<N;k++){
				ab_temp[i*N+j]= ab_temp[i*N+j] + A_buf[i*N+k]*B[k+j*N];
			}
		}
	}

	#pragma omp for
	for (int i=0;i<sizePart;i++){
		ab_temp[i]*=u*l;
	}

	#pragma omp for
	for(int i=0;i<sizePart/N;i++){
		for(int j=0;j<N;j++){
			lc_temp[i*N+j]=0;
			for(int k=0;k<posicion+i+1;k++){
				lc_temp[i*N+j]= lc_temp[i*N+j] + L_buf[i*N+k]*C[k+j*N];
			}
		}
	}

	#pragma omp for
	for (int i=0;i<sizePart;i++){
		lc_temp[i]*=u*l;
	}

	#pragma omp for
	for(int i=0;i<sizePart/N;i++){
		for(int j=0;j<N;j++){
			du_temp[i*N+j]=0;
			for(int k=0;k<j+1;k++){
				du_temp[i*N+j]= du_temp[i*N+j] + D_buf[i*N+k]*U[k+j*(j+1)/2];
			}
		}
	}

	#pragma omp for
	for (int i=0;i<sizePart;i++){
		du_temp[i]*=u*l;
	}

	#pragma omp for
	for (int i=0;i<sizePart;i++){
		ab_temp[i]+=lc_temp[i]+du_temp[i];
	}
	}
	MPI_Gather(ab_temp, sizePart, MPI_DOUBLE, M, sizePart, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	tiempo_paral = dwalltime() - tiempo_paral;
	tiempo_balance = dwalltime() - tiempo_balance;
	printf("Tiempo en segundos paralelo %f \n", tiempo_paral);
	printf("Tiempo en segundos balance %f \n", tiempo_balance);

	if(miID==0){
		for (int j=0; j<N; j++){
			for (int i=0; i<N; i++){
				printf ("%lf \t",M[j*N+i]);
			}
			printf ("\n");
		}
	}

	//***** FIN PROGRAMA ****
	
	free(B);
	free(C);
	free(U);
	free(ab_temp);
	free(lc_temp);
	free(du_temp);
	free(A_buf);
	free(L_buf);
	free(D_buf);
	if(miID==0){
		free(A);
		free(L);
		free(D);
		free(M);		
	}
	MPI_Finalize();
	return(0);
}
