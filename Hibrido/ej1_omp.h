#ifndef __OMP_FUNCTIONS__
#define __OMP_FUNCTIONS__

extern void omp_iniciar(int T);
extern double omp_sumaTemp(double *pruebaB, int N, int T);
extern double omp_sumaTemp1(double *U,int N, int T);
extern double omp_sumaTemp2(double *pruebaL,int N, int T,int ID);
extern void omp_parcialAB(double *parcialAB_SUB, double *pruebaA, double *B, int N, int T);
extern void omp_parcialLC(double *parcialLC_SUB, double *pruebaL, double *C ,int N, int T,int ID);
extern void omp_parcialDU(double *parcialDU_SUB, double *pruebaD, double *U ,int N, int T);
extern void omp_parcialM(double *parcialM,double *parcialAB_SUB,double *parcialLC_SUB,double *parcialDU_SUB, double ul,int N, int T);

#endif