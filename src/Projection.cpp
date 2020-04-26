#include "Projection.hpp"


void projection_L1Ball(double *y, double *x, int length, double a){
	double*	aux = new double[length];
	double*  aux0 = aux;
	int		auxlength=1; 
	int		auxlengthold=-1;	
	double	tau=(*aux=*y)-a;
	int 	i=1;
	for (; i<length; i++) 
		if (y[i]>tau) {
			if ((tau+=((aux[auxlength]=y[i])-tau)/(auxlength-auxlengthold))<=y[i]-a) {
				tau=y[i]-a;
				auxlengthold=auxlength-1;
			}
			auxlength++;
		} 
	if (auxlengthold>=0) {
		auxlength-=++auxlengthold;
		aux+=auxlengthold;
		while (--auxlengthold>=0) 
			if (aux0[auxlengthold]>tau) 
				tau+=((*(--aux)=aux0[auxlengthold])-tau)/(++auxlength);
	}
	do {
		auxlengthold=auxlength-1;
		for (i=auxlength=0; i<=auxlengthold; i++)
			if (aux[i]>tau) 
				aux[auxlength++]=aux[i];	
			else 
				tau+=(tau-aux[i])/(auxlengthold-i+auxlength);
	} while (auxlength<=auxlengthold);
// printf("%f\n",tau);

	for (i=0; i<length; i++){
		x[i]=(y[i]>tau ? y[i]-tau : 0.0); }
	if (x==y) {delete[] aux0;}
} 
