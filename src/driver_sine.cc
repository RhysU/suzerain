/*
! This sample program illustrates the 
! use of P3DFFT library for highly scalable parallel 3D FFT. 
!
! This program initializes a 3D array with a 3D sine wave, then 
! performs 3D forward Fourier transform, then backward transform, 
! and checks that 
! the results are correct, namely the same as in the start except 
! for a normalization factor. It can be used both as a correctness
! test and for timing the library functions. 
!
! The program expects 'stdin' file in the working directory, with 
! a single line of numbers : Nx,Ny,Nz,Ndim,Nrep. Here Nx,Ny,Nz
! are box dimensions, Ndim is the dimentionality of processor grid
! (1 or 2), and Nrep is the number of repititions. Optionally
! a file named 'dims' can also be provided to guide in the choice 
! of processor geometry in case of 2D decomposition. It should contain 
! two numbers in a line, with their product equal to the total number
! of tasks. Otherwise processor grid geometry is chosen automatically.
! For better performance, experiment with this setting, varying 
! iproc and jproc. In many cases, minimizing iproc gives best results. 
! Setting it to 1 corresponds to one-dimensional decomposition.
!
! If you have questions please contact Dmitry Pekurovsky, dmitry@sdsc.edu
*/

#include "p3dfft_cxx.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void print_all(double *,long int,int,long int);
void mult_array(double *,long int,double);
double FORTNAME(t1),FORTNAME(t2),FORTNAME(t3),FORTNAME(t4),FORTNAME(tp1);
/* double t1,t2,t3,t4,tp1; */

int main(int argc,char **argv)
{
   double *A,*B,*p,*C;
   int i,j,k,x,y,z,nx,ny,nz,proc_id,nproc,dims[2],ndim,nu;
   int istart[3],isize[3],iend[3];
   int fstart[3],fsize[3],fend[3];
   int iproc,jproc,ng[3],kmax,iex,conf,m,n;
   long int Nglob,Ntot;
   double pi,twopi,sinyz;
   double *sinx,*siny,*sinz,factor;
   double rtime1,rtime2,gt1,gt2,gt3,gt4,gtp1,gtcomm,tcomm;
   double cdiff,ccdiff,ans;
   FILE *fp;

   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD,&nproc);
   MPI_Comm_rank(MPI_COMM_WORLD,&proc_id);

   pi = atan(1.0)*4.0;
   twopi = 2.0*pi;

   gt1=gt2=gt3=gt4=gtp1=0.0;

 
   if(proc_id == 0) {
     fp = fopen("stdin","r");
     ndim = 2;
     fscanf(fp,"%d,%d,%d,%d,%d\n",&nx,&ny,&nz,&ndim,&n);
     fclose(fp);
     printf("Single precision\n (%d %d %d) grid\n %d proc. dimensions\n%d repetitions\n",nx,ny,nz,ndim,n);
   }
   MPI_Bcast(&nx,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&ny,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nz,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&ndim,1,MPI_INT,0,MPI_COMM_WORLD);
   
   if(ndim == 1) {
     dims[0] = 1; dims[1] = nproc;
   }
   else if(ndim == 2) {
     fp = fopen("dims","r");
     if(fp != NULL) {
       if(proc_id == 0)
         printf("Reading proc. grid from file dims\n");
       fscanf(fp,"%d %d\n",dims,dims+1);
       fclose(fp);
       if(dims[0]*dims[1] != nproc) 
          dims[1] = nproc / dims[0];
     }
     else {
       if(proc_id == 0) 
          printf("Creating proc. grid with mpi_dims_create\n");
       dims[0]=dims[1]=0;
       MPI_Dims_create(nproc,2,dims);
       if(dims[0] > dims[1]) {
          dims[0] = dims[1];
          dims[1] = nproc/dims[0];
       }
     }
   }

   if(proc_id == 0) 
      printf("Using processor grid %d x %d\n",dims[0],dims[1]);

   /* Initialize P3DFFT */
   p3dfft_setup(dims,nx,ny,nz,1);
   /* Get dimensions for input and output arrays */
   conf = 1;
   get_dims(istart,iend,isize,conf);
   conf = 2;
   get_dims(fstart,fend,fsize,conf);

   sinx = (double *) malloc(sizeof(double)*nx);
   siny = (double *) malloc(sizeof(double)*ny);
   sinz = (double *) malloc(sizeof(double)*nz);

   for(z=0;z < isize[2];z++)
     sinz[z] = sin((z+istart[2]-1)*twopi/nz);
   for(y=0;y < isize[1];y++)
     siny[y] = sin((y+istart[1]-1)*twopi/ny);
   for(x=0;x < isize[0];x++)
     sinx[x] = sin((x+istart[0]-1)*twopi/nx);

   /* Allocate and Initialize */
   A = (double *) malloc(sizeof(double) * isize[0]*isize[1]*isize[2]);
   B = (double *) malloc(sizeof(double) * fsize[0]*fsize[1]*fsize[2]*2);
   C = (double *) malloc(sizeof(double) * isize[0]*isize[1]*isize[2]);

   p = A;
   for(z=0;z < isize[2];z++)
     for(y=0;y < isize[1];y++) {
       sinyz = siny[y]*sinz[z];
       for(x=0;x < isize[0];x++)
          *p++ = sinx[x]*sinyz;
     }

   Ntot = fsize[0]*fsize[1];
   Ntot *= fsize[2]*2;
   Nglob = nx * ny;
   Nglob *= nz;
   factor = 1.0/Nglob;

   rtime1 = 0.0;
   for(m=0;m < n;m++) {

     MPI_Barrier(MPI_COMM_WORLD);
     rtime1 = rtime1 - MPI_Wtime();
     if(proc_id == 0) 
        printf("Iteration %d\n",m);
     /* compute forward Fourier transform on A, store results in B */
     p3dfft_ftran_r2c(A,B);
     rtime1 = rtime1 + MPI_Wtime();

     if(proc_id == 0) 
        printf("Result of forward transform\n");

     print_all(B,Ntot,proc_id,Nglob);
     /* normalize */
     mult_array(B,Ntot,factor);

     /* Compute backward transform on B, store results in A */
     MPI_Barrier(MPI_COMM_WORLD);
     rtime1 = rtime1 - MPI_Wtime();
     p3dfft_btran_c2r(B,C);
     rtime1 = rtime1 + MPI_Wtime();

   } 
   /* free work space */
  p3dfft_clean();
  
  /* Check results */
  cdiff = 0.0; p = C;
  for(z=0;z < isize[2];z++)
    for(y=0;y < isize[1];y++)  {
       sinyz =siny[y]*sinz[z];
       for(x=0;x < isize[0];x++) {
          ans = sinx[x]*sinyz;
          if(cdiff < fabs(*p - ans))
           cdiff = fabs(*p - ans);
	   p++;
        }
    }

   MPI_Reduce(&cdiff,&ccdiff,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

  if(proc_id == 0)
    printf("max diff =%g\n",ccdiff);

  /* Gather timing statistics */
  MPI_Reduce(&rtime1,&rtime2,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

  /*
  MPI_Reduce(&FORTNAME(t1),&gt1,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  MPI_Reduce(&FORTNAME(t2),&gt2,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  MPI_Reduce(&FORTNAME(t3),&gt3,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  MPI_Reduce(&FORTNAME(t4),&gt4,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  MPI_Reduce(&FORTNAME(tp1),&gtp1,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  tcomm = FORTNAME(t1)+FORTNAME(t2)+FORTNAME(t3)+FORTNAME(t4);
  MPI_Reduce(&tcomm,&gtcomm,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

  gt1 = gt1 / ((double) n);
  gt2 = gt2 / ((double) n);
  gt3 = gt3 / ((double) n);
  gt4 = gt4 / ((double) n);
  gtp1 = gtp1 / ((double) n);
  gtcomm = gtcomm / ((double) n);
  */

  if(proc_id == 0) {
     printf("Time per loop=%lg\n",rtime2/((double) n));
     /*
     printf("Total comm: %g",gtcomm);
     printf("t1=%lg, t2=%lg, t3=%lg, t4=%lg, tp1=%lg\n",gt1,gt2,gt3,gt4,gtp1);
     */
  }


  MPI_Finalize();

}

void mult_array(double *A,long int nar,double f)
{
  long int i;

  for(i=0;i < nar;i++)
    A[i] *= f;
}

void print_all(double *A,long int nar,int proc_id,long int Nglob)
{
  int x,y,z,conf,Fstart[3],Fsize[3],Fend[3];
  long int i;

  conf = 2;
  get_dims(Fstart,Fend,Fsize,conf);
  Fsize[0] *= 2;
  Fstart[0] = 1 + (Fstart[0]-1)*2;
  for(i=0;i < nar;i++)
    if(fabs(A[i]) > Nglob *1.25e-4) {
      z = i/(Fsize[0]*Fsize[1]);
      y = i/Fsize[0] - z*Fsize[1];
      x = i-z*Fsize[0]*Fsize[1] - y*Fsize[0];
      printf("(%d,%d,%d) %f\n",x+Fstart[0],y+Fstart[1],z+Fstart[2],A[i]);
    }
}

  
