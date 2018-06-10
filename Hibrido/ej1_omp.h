#ifndef __OMP_FUNCTIONS__
#define __OMP_FUNCTIONS__

extern void omp_function(int rank);
double omp_sumaTemp(double *pruebaB, int N, int T);
double omp_sumaTemp1(double *U,int N, int T);
double omp_sumaTemp2(double *pruebaL,int N, int T,int ID);
void omp_parcialAB(double *parcialAB_SUB, double *pruebaA, double *B, int N, int T);
void omp_parcialLC(double *parcialLC_SUB, double *pruebaL, double *C ,int N, int T,int ID);
void omp_parcialDU(double *parcialDU_SUB, double *pruebaD, double *U ,int N, int T);
void omp_parcialM();

#endif