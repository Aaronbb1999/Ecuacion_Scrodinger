#include <stdio.h>
#include <math.h>
#include "complex.h"

#define PI 3.1415926535 
#define N 1000
#define nciclos 50
#define lambda 3
#define n 500

int main(void)
{
    //Definición de constantes y parametros
    int i, j, k;
    double s, k0;
    double mod, arg;
    FILE *fout, *fnorma;

    fout=fopen("schrodinger.dat","w");
    fnorma=fopen("norma_3.dat","w");

    //Definición vectores de la simulacion
    double norma;
    double V[N+1], modphi[N+1];
    fcomplex alfa[N], beta[N], gamma[N];
    fcomplex phi[N+1], xi[N+1], b[N+1];

    //Calculamos primeros parámetros y condiciones iniciales
    k0=(2.0*PI*nciclos)/N;
    s=1/(k0*k0*4.0);

    alfa[N-1]=Complex(0.0,0.0);
    beta[N-1]=Complex(0.0,0.0);

    phi[0]=Complex(0.0,0.0);
    phi[N]=Complex(0.0,0.0);
    xi[0]=Complex(0.0,0.0);
    xi[N]=Complex(0.0,0.0);

    for(j=1;j<N;j++)
    {
        mod=exp(-8.0*((4.0*j-N)*(4.0*j-N))/(N*N));
        arg=j*k0;
        phi[j]=Cgauss(arg,mod);

        phi[j]=RCmul(0.5,phi[j]); //Renormalicización para mostrar mejor la animación
    }

    //Potencial V
    for(j=0;j<N+1;j++)
    {
        if((j<2.0*N/5)||(j>3.0*N/5))
        {
            V[j]=0;
        }
        else
        {
            V[j]=lambda*k0*k0;
        }
    }

    //Calculamos alfa y gamma

    fcomplex num, den, aux;

    for(j=N-2;j>0;j--)
    {
        aux.r=-2-V[j+1];
        aux.i=2/s;

        den=Cadd(aux,alfa[j+1]);
        num=Complex(1.0,0.0);

        gamma[j]=Cdiv(num,den);
        alfa[j]=RCmul(-1.0,gamma[j]);
    }

    //Rellenamos el fichero con las condiciones iniciales del modulo
    for(j=0;j<=N;j++)
    {
        modphi[j]=Cabs(phi[j])*Cabs(phi[j]);

        fprintf(fout,"%i, %lf, %lf\n",j,modphi[j],V[j]);
    }
    fprintf(fout, "\n");

    //Rellenamos el fichero con la norma inicial
    for(j=0;j<=N;j++)
    {
        norma+=modphi[j];
    }
    fprintf(fnorma,"%lf\n", norma);


    //Bucle principal

    fcomplex cb, cbeta;

    cb=Complex(0.0,4/s);

    for (k=0; k<n; k++)
    {
        //Calculamos b
        for (j=0; j<=N; j++)
        {
            b[j]=Cmul(cb,phi[j]);
        }

        //Calculamos Beta
        for(j=N-2; j>=0; j--)
        {
            cbeta=Csub(b[j+1],beta[j+1]);
            beta[j]=Cmul(cbeta,gamma[j]);
        }

        //Calculamos xi
        for(j=0; j<N; j++)
        {
            aux=Cmul(alfa[j],xi[j]);
            xi[j+1]=Cadd(aux,beta[j]);
        }

        //Calculamos phi y escribimos módulo en el fichero
        for(j=0; j<=N; j++)
        {
            phi[j]=Csub(xi[j],phi[j]);
            modphi[j]=Cabs(phi[j])*Cabs(phi[j]);

            fprintf(fout,"%i,%lf,%lf\n",j,modphi[j],V[j]);
        }
        fprintf(fout, "\n");

        //Calculamos la norma también y la sacamos al fichero
        norma=0;
        for (j=0; j<=N; j++)
        {
            norma+=modphi[j];
        }
        fprintf(fnorma,"%lf\n", norma);
    }

fclose(fout);
fclose(fnorma);

}