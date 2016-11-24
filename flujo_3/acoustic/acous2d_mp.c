 // http://www.reproducibility.org/RSF/book/rsf/manual/manual_html/node1907.html
#include<omp.h>
#include<rsf.h>
#include<stdlib.h>
#define Id(A,i,j) (A)[(i)*nz+(j)]
void swapptr(float *a[],float *b[]){ //swap arrays
  float *ptr=*a;
  *a = *b;
  *b = ptr;
}
int main(int argc, char* argv[]){
  
  int it,i,j;     // index variables 
  int nt,nz,nx;   // dimensions
  int Iz,Ix;      // source index coordinates   
  int skip;
  float dt,dz,dx; // cells, & steptime sizes 
  float Sx,Sz;    // source coordinates in meters
  float G;        // auxiliar variables
   
  // Initialize RSF 
  sf_init(argc,argv);

   // setup I/O files 
  sf_file Fv = sf_input ("in");   //velocity media
  sf_file Fw = sf_input ("wav" ); //wavelet input
  sf_file Fo = sf_output("out");  //P-field outputo


  //Read axes-rsf-data parameters 
  if (!sf_histint  (Fw,"n1",&nt)) sf_error("No n1= in wav");
  if (!sf_histfloat(Fw,"d1",&dt)) sf_error("No d1= in wav");

  if (!sf_histint  (Fv,"n1",&nz)) sf_error("No nz= in inp");
  if (!sf_histint  (Fv,"n2",&nx)) sf_error("No nx= in inp");
  if (!sf_histfloat(Fv,"d1",&dz)) sf_error("No dz= in inp");
  if (!sf_histfloat(Fv,"d2",&dx)) sf_error("No dx= in inp");
  
  // parameter from the command line 
  if (!sf_getfloat("Sx",&Sx)) 
    sf_error("Need Sx="); // source x-coord [m] 
   
  if (!sf_getfloat("Sz",&Sz)) 
    sf_error("Need Sz="); // source z-coord [m] 
  
  if (!sf_getint("Skip",&skip)) 
    sf_error("Need Skip="); // Interval between frames 

  //Write axes-rsf-data parameters 
  sf_putint(Fo,"n3",nt/skip);
  sf_putfloat(Fo,"d3",dt*skip);
  sf_putfloat(Fo,"o2",0.0);


//read wavelet, velocity & density
  float *wavelet=sf_floatalloc(nt); 
  sf_floatread(wavelet,nt,Fw);
   
  float *vel=sf_floatalloc(nz*nx); 
  sf_floatread(vel,nz*nx,Fv);

  // allocate fdtd arrays 
  float  *p1=sf_floatalloc(nz*nx);
  float  *p2=sf_floatalloc(nz*nx);
   
  Ix = (int)Sx/dx; Iz = (int)Sz/dz;

  // Begin time-loop 
   for (it=0; it<nt; it++) {

   #pragma omp parallel for default(shared) private(i,j,G)
     for (i=2; i<nx-2; i++){
       for (j=2; j<nz-2; j++){
         G = Id(vel,i,j)*dt/dx;
         Id(p1,i,j) = 2*Id(p2,i,j)-Id(p1,i,j) + 
                  G*G*( (4/3.)*(Id(p2,i+1,j)+Id(p2,i-1,j)) + 
                        (4/3.)*(Id(p2,i,j+1)+Id(p2,i,j-1)) +
                      (-1/12.)*(Id(p2,i+2,j)+Id(p2,i-2,j)) +
                      (-1/12.)*(Id(p2,i,j+2)+Id(p2,i,j-2)) +
                          (-5)*Id(p2,i,j) );
		     }
      } 
      // step source perturbation 
      Id(p2,Ix,Iz) += wavelet[it];

      swapptr(&p1,&p2);
      
      // write wavefield to output
      if(it%skip==0){
         sf_floatwrite(p2,nz*nx,Fo);
         sf_warning("%d\%;",it*100/nt); 
      }
   }// End time-loop 

void abc(int Nx,int Nz,int cells,float a,float *field) {
        int i,j;
        float taper;


#pragma omp parallel for default(shared) private(i,j,taper)
        for(i=0;i<cells;i++){
        	        for(j=0;j<Nz;j++){
                taper = exp(-(a*(cells-i))*(a*(cells-i))); //left side
                Id(field,i,j) *= taper;
		}
	}
	for(i=0;i<Nx;i++){
			for(j=Nz-cells;j<Nz;j++){	
                taper = exp(-(a*(Nz-cells-j))*(a*(Nz-cells-j))); //left side
                Id(field,i,j) *= taper;
	}
	}
        for(i=Nx-cells;i<Nx;i++){
                        for(j=0;j<Nz;j++){
                taper = exp(-(a*(Nx-cells-i))*(a*(Nx-cells-i))); //left side
                Id(field,i,j) *= taper;
        }
        }
}
 
   sf_close();
   return(0);
}

