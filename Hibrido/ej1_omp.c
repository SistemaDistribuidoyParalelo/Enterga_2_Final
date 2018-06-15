#include<omp.h>
#include<stdio.h>
#include <stdlib.h>
#include <string.h>
void omp_iniciar(int T){
    omp_set_num_threads(T);
}
//Suma Parcial para realizar el promedio de U
double omp_sumaTemp1(double *U,int N, int T){
    int i;
    int j;
    double temp1;
    
    #pragma omp parallel for reduction(+:temp1) private(i,j)
    for(i=0;i<N;i++){
        for (j=i;j<N;j++){
            temp1+=U[i+((j*(j+1))/2)];
        }
    }
    
    return temp1;
}

//Suma parcial para realizar el promedio de L
double omp_sumaTemp2(double *pruebaL,int N, int T,int ID){
    int i;
    int k;
    double temp2;   
    
    #pragma omp parallel for reduction(+:temp2) private(i,k)
    for(i=0;i < N/T;i++){
        for(k=0;k<i+((N/T)*(ID))+1;k++){
            temp2+=pruebaL[i*N+k];
	    }     
    }
    
    
    return temp2;
}

//Primera multiplicacion Parcial A*B
void omp_parcialAB(double *parcialAB_SUB, double *pruebaA, double *B ,int N, int T){
    int i;
    int k;
    int j;
   //#pragma omp parallel for shared(parcialAB_SUB,pruebaA,B) private(i,j,k)
    for(i=0;i<N/T;i++){
        #pragma omp parallel for firstprivate(i)     
        for(j=0;j<N;j++){
            for(k = 0;k<N;k++){
                parcialAB_SUB[i*N+j] += pruebaA[i*N+k]*B[k*N+j];
            }
        }
    }
    //return parcialAB_SUB
}

//Segunda multiplicacion Parcial LC
void omp_parcialLC(double *parcialLC_SUB, double *pruebaL, double *C ,int N, int T,int ID){
    int i;
    int j;
    int k;
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,T)
        for(i=0;i<N/T;i++){
            for(j=0;j<N;j++){
                for(k = 0;k<i+(N/T*ID)+1;k++){
                    parcialLC_SUB[i*N+j] +=pruebaL[i*N+k]*C[j*N+k];
                }
            }
        }
    }
}

//Tercera multiplicacion Parcial DU
void omp_parcialDU(double *parcialDU_SUB, double *pruebaD, double *U ,int N, int T){
    int i;
    int j;
    int k;
    
    //U es superior, recorrido parcial
    for(i=0;i<N/T;i++){
        #pragma omp parallel for firstprivate(i)
        for(j=0;j<N;j++){
            for(k = 0;k<j+1;k++){
                parcialDU_SUB[i*N+j] +=pruebaD[i*N+k]*U[k+((j*(j+1))/2)];
            }
        }
    }


}

//Suma final, incluyendo el valor ul como parametro, promedio de cada una
void omp_parcialM(double *parcialM,double *parcialAB_SUB,double *parcialLC_SUB,double *parcialDU_SUB, double ul,int N, int T){
    int i;
    int j;
    int k;
    //EN ESTE PRAGMA, SE LE PUEDE COLOCAR EL NOWAIT, como es la suma total, cuando terminen buenisimo.
    #pragma omp parallel 
    {
        #pragma omp for schedule(dynamic,T)
        for(i=0;i<N;i++){
            for(j=0;j<N;j++){
                for(k=0;k<N;k++){
                    parcialM[i*N+j] = ul*(parcialAB_SUB[i*N+j]+parcialLC_SUB[i*N+j]+parcialDU_SUB[i*N+j]);
                }
            }
        }
    }
}
