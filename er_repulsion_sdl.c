/*******************************************************************
*
* Mode d'emploi: 
*
* Cet algorithme vous permet de simuler le mouvement
* de particules au sein des anneaux de Saturne. Vous
* remarquerez les constantes physiques correspondantes
* à notre modèle.
* 
* Afin d'augmenter le nombre d'astéroïdes orbitant au sein
* des anneaux, il suffit d'augmenter la variable N.
* Un fichier sera créé nous le nom de "euler_richardson.txt"
* avec deux colonnes représentant la position en coordonnées
* cartésiennes.
*
**********************************************************************/



#include <stdio.h>
#include <stdlib.h>  /
#include <math.h>    
#include <SDL/SDL.h>


#define N 10000
#define D 100000
#define G 6.67408e-10
#define Ms 5.683e26
#define Mast 500.
#define k_raideur 50000.


void pause()
{
  int continuer = 1;
  SDL_Event event;
 
  while (continuer)
    {
      SDL_WaitEvent(&event);
      switch(event.type)
        {
	case SDL_QUIT:
	  continuer = 0;
        }
    }
}

double distance(double x, double y, double x_tab[2][D], double y_tab[2][D], int l){  
  double g(double x, double y, double x_tab[2][D], double y_tab[2][D], int i, int j){
    return (( fabs(x-x_tab[i][j])  )/( pow(pow(x-x_tab[i][j],2) + pow(y-y_tab[i][j],2),1.5) ));
  }  
  double d = g(x,y,x_tab,y_tab,0,l-1)+g(x,y,x_tab,y_tab,0,l)+g(x,y,x_tab,y_tab,0,l+1);
  return d;
}


double f(double x, double y, double x_tab[2][D], double y_tab[2][D], int l){
  return -(G*Ms*x)/(pow(pow(x,2)+pow(y,2),1.5)) +(k_raideur*Mast)*distance(x,y,x_tab,y_tab,l);
}


void euler_richardson(double T, double h, double eps_seuil){
  FILE *file;  
  double y[4], k[4], dk[4], alpha, eps, t0 = 0, h_temp = h;
  double x_tab[2][D], y_tab[2][D]; 
  int i,l,m,n,COUNT=0;

  //Configuration de la fenêtre 
  SDL_Surface *ecran = NULL, *rectangle = SDL_CreateRGBSurface(SDL_HWSURFACE, 1, 1, 32, 0, 0, 0, 0); 
  SDL_Rect posi; // coord[N][D+1] 
  SDL_Init(SDL_INIT_VIDEO);
  SDL_WM_SetCaption("Euler Richardson - Anneaux de Saturne", NULL);  
  ecran = SDL_SetVideoMode(900, 600, 32, SDL_HWSURFACE | SDL_DOUBLEBUF);


  file=fopen("euler_richardson.txt","w+");
  
  for (n=0;n<D;n++){
    x_tab[0][n]=0;
    y_tab[0][n]=0;
  }
   
  for (m=1;m<N;m++){
    y[0]=66000+0.005*m;//séparation de 5m entre astéroides
    y[1]=0;
    y[2]=0;
    y[3]=57600;
    l=0;
    t0=0;
    h = h_temp;
    while(t0<=T){
       if (l>=D-1)
	 break;
      COUNT++;      

      k[0]=y[2];
      k[1]=y[3];
      k[2]= f(y[0],y[1],x_tab,y_tab,l); 
      k[3]= f(y[1],y[0],y_tab,x_tab,l); 
    
      dk[0]=y[2]+(h/2)*k[2];
      dk[1]=y[3]+(h/2)*k[3];
      dk[2]=f(y[0]+(h/2)*k[0],y[1]+(h/2)*k[1],x_tab,y_tab,l); 
      dk[3]=f(y[1]+(h/2)*k[1],y[0]+(h/2)*k[0],y_tab,x_tab,l);

      eps  = (h/2)*sqrt(pow(dk[0]-k[0],2)+pow(dk[1]-k[1],2)+pow(dk[2]-k[2],2)+pow(dk[3]-k[3],2));
      alpha= eps/eps_seuil;

      if (alpha>1)
	h*=0.9/sqrt(alpha);
      else if(alpha<1){
	y[0]+=h*dk[0];
	y[1]+=h*dk[1];
	y[2]+=h*dk[2];
	y[3]+=h*dk[3];      
	h*=0.9/sqrt(alpha);
	t0+=h;      
      }
      x_tab[1][l]=y[0];
      y_tab[1][l]=y[1];
      fprintf(file,"%f\t%f\n",y[0],y[1]);
      l++;
      posi.x = y[0]/90+90;
      posi.y = y[1]/20+300;
      SDL_FillRect(rectangle, NULL, SDL_MapRGB(ecran->format, 255, 255, 255)); 
      SDL_BlitSurface(rectangle, NULL, ecran, &posi);
      SDL_Flip(ecran);
    }
    for (n=0;n<D;n++){
      x_tab[0][n]=x_tab[1][n];
      y_tab[0][n]=y_tab[1][n];
    }
  }  
  printf("%d\n",COUNT);
  pause();
  SDL_FreeSurface(rectangle);
  SDL_Quit();
  fclose(file);
}


int main(int argc, char** argv){  
  euler_richardson(.08,1,1000);  
  return 0;
}
