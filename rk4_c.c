/*******************************************************************
*
* Mode d'emploi: 
*
* Cet algorithme vous permet de simuler le mouvement
* de particules au sein des anneaux de Saturne à l'aide de l'algorithme
* Runge-Kutta 4. Vous remarquerez les constantes physiques correspondants
* à notre modèle.
* 
* Afin d'augmenter le nombre d'astéroïdes orbitant au sein
* des anneaux, il suffit d'augmenter la variable N.
* Un fichier sera créé nous le nom de "rk4_0.19.txt"
* avec deux colonnes représentant la position en coordonnées
* cartésiennes.
*
**********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 1
#define D 5000

double f(double x, double y){
  return (-(4*pow(acos(-1),2)*x)/(pow(pow(x,2)+pow(y,2),1.5)));
}

void rk4(double T, double h){
  FILE* file;
  file = fopen("rk4_0.19.txt","w+");
  double y[4],k1[4],k2[4],k3[4],k4[4],pos[2][D][N],t=0;
  
  //Initialisation du vecteur y
  y[0]=0.5;// x
  y[1]=0;//y
  y[2]=0;//vx
  y[3]=11.5;//v y

  while(t<=T){
    
    k1[0]=y[2];
    k1[1]=y[3];
    k1[2]=f(y[0],y[1]);
    k1[3]=f(y[1],y[0]);

    k2[0]=y[2];
    k2[1]=y[3];
    k2[2]=f(y[0]+h*(k1[0]/2),y[1]+h*(k1[1]/2));
    k2[3]=f(y[1]+h*(k1[1]/2),y[0]+h*(k1[0]/2));

    k3[0]=y[2];
    k3[1]=y[3];
    k3[2]=f(y[0]+h*(k2[0]/2),y[1]+h*(k2[1]/2));
    k3[3]=f(y[1]+h*(k2[1]/2),y[0]+h*(k2[0]/2));

    k4[0]=y[2];
    k4[1]=y[3];
    k4[2]=f(y[0]+h*(k3[0]/2),y[1]+h*(k3[1]/2));
    k4[3]=f(y[1]+h*(k3[1]/2),y[0]+h*(k3[0]/2));

    y[0]+= (h/6)*(k1[0]+2*k2[0]+2*k3[0]+k4[0]);
    y[1]+= (h/6)*(k1[1]+2*k2[1]+2*k3[1]+k4[1]);
    y[2]+= (h/6)*(k1[2]+2*k2[2]+2*k3[2]+k4[2]);
    y[3]+= (h/6)*(k1[3]+2*k2[3]+2*k3[3]+k4[3]);

    fprintf(file,"%f\t%f\n",y[0],y[1]);
    
    t+=h;
  }
}

int main(){

  rk4(25.95,0.010);

  return 0;
}
