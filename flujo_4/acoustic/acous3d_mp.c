 // http://www.reproducibility.org/RSF/book/rsf/manual/manual_html/node1907.html
#include<omp.h>
#include<rsf.h>
#include<stdlib.h>
#define Id(A,i,j,k) (A)[ (i)*ny*nz + (j)*nz + (k) ]
void swapptr(float *a[],float *b[]){ //swap arrays
  float *ptr=*a;
  *a = *b;
  *b = ptr;
}
int main(int argc, char* argv[]){
  
  int it,i,j,k;      // index variables 
  int nt,nz,ny,nx;   // dimensions
  int Iz,Iy,Ix;      // source index coordinates   
  int skip;
  float dt,dz,dy,dx; // cells, & steptime sizes 
  float Sz,Sy,Sx;    // source coordinates in meters
  float G;           // Couran stability parameter
   
  // Initialize RSF 
  sf_init(argc,argv);

   // setup I/O files 
  sf_file Fv = sf_input ("in");   //velocity media
  sf_file Fw = sf_input ("wav" ); //wavelet input
  sf_file Fo = sf_output("out");  //P-field outputo


  //Read axes-rsf-data parameters 
  if (!sf_histint (Fw,"n1",&nt)) sf_error("No n1= in wav");
  if (!sf_histfloat(Fw,"d1",&dt)) sf_error("No d1= in wav");
  if (!sf_histint (Fv,"n1",&nz)) sf_error("No nz= in inp");
  if (!sf_histint (Fv,"n2",&ny)) sf_error("No ny= in inp");
  if (!sf_histint (Fv,"n3",&nx)) sf_error("No nx= in inp");
  if (!sf_histfloat(Fv,"d1",&dz)) sf_error("No dz= in inp");
  if (!sf_histfloat(Fv,"d2",&dy)) sf_error("No dy= in inp");
  if (!sf_histfloat(Fv,"d3",&dx)) sf_error("No dx= in inp");

  // parameter from the command line 
  if (!sf_getfloat("Sx",&Sx)) 
    sf_error("Need Sx="); // source x-coord [m] 
  
  if (!sf_getfloat("Sy",&Sy)) 
    sf_error("Need Sy="); // source y-coord [m] 
   
  if (!sf_getfloat("Sz",&Sz)) 
    sf_error("Need Sz="); // source z-coord [m] 
  
  if (!sf_getint("Skip",&skip)) 
    sf_error("Need Skip="); // Interval between frames 

  
  //Write axes-rsf-data parameters (4D)
  sf_putint(Fo,"n4",nt/skip);
  sf_putfloat(Fo,"d4",dt*skip);
  sf_putfloat(Fo,"o4",0.0);

  //read wavelet, velocity
  float *wavelet=sf_floatalloc(nt); 
  sf_floatread(wavelet,nt,Fw);
   
  float *vel=sf_floatalloc(nz*ny*nx); 
  sf_floatread(vel,nz*ny*nx,Fv);

  // allocate fdtd arrays 
  float  *p1=sf_floatalloc(nz*ny*nx);
  float  *p2=sf_floatalloc(nz*ny*nx);
   
  Ix = (int)(Sx/dx); Iy = (int)(Sy/dy); Iz = (int)(Sz/dz);

  // Begin time-loop 
   for (it=0; it<nt; it++) {

   #pragma omp parallel for default(shared) private(i,j,k,G)
     for (i=2; i<nx-2; i++){
       for (j=2; j<ny-2; j++){
         for (k=2; k<nz-2; k++){
           G = Id(vel,i,j,k)*dt/dx;
		     Id(p1,i,j,k) = 2*Id(p2,i,j,k)-Id(p1,i,j,k) +
                  G*G*( (4/3.)*(Id(p2,i+1,j,k)+Id(p2,i-1,j,k)) +
                        (4/3.)*(Id(p2,i,j+1,k)+Id(p2,i,j-1,k)) +
                      (-1/12.)*(Id(p2,i+2,j,k)+Id(p2,i-2,j,k)) +
                      (-1/12.)*(Id(p2,i,j+2,k)+Id(p2,i,j-2,k)) +
                         (4/3.)*(Id(p2,i,j,k+1)+Id(p2,i,j,k-1)) +
		       (-1/12.)*(Id(p2,i,j,k-2)+Id(p2,i,j,k+2)) +
			  (-15/2.)*Id(p2,i,j,k)
			 );
//	sf_warning("%d\%;",Id(p1,i,j,k));	
	}
       }
     } 
     
    // step source perturbation 
    Id(p2,Ix,Iy,Iz) += wavelet[it];

    swapptr(&p1,&p2);
      
    // write wavefield to output
    if(it%skip==0){
         sf_floatwrite(p2,nz*ny*nx,Fo);
         sf_warning("%d\%;",it*100/nt); 
      }
   }// End time-loop 
   
   sf_close();
   return(0);
}

