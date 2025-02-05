
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "Lattice.h"
#include "Substrate.h"
#include "Simulation.h"
#include <glib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_rng.h>
#define PI (3.141592653589793)


void Gen_Adsorbate(struct Adsorbate * A, struct Lattice *L, struct Parameters* P)
{
	// Generate adsorbate and set some parameters
	        
	A->num_of_samples = P->Simulation_Time/P->Adsorbate_Samplerate; // How many times we sample the adsorbate's trajectory
	A->b = P->b;
	A->mass = P->Amass;
	A->Amp = P->Amp;
	A->R0 = P->R0;
	A->Force_Amplitude = A->Amp*A->b*A->b*2; // multiply by b twice cause later we multiply R by b and need to compenstae
	A->preexp = exp(A->b*A->R0);
	A->dt_over_m = P->Timestep/(6*A->mass); // used in calculating the velocity
        A->dt2_over_m = P->Timestep*P->Timestep/(6*A->mass); // used in calculating the displacement
        if (P->Track_Adsorbate_Potential || P->GetPES2) // do we save the potential
	{
	  A->savepot = 1;
	  A->GetPES2 = P->GetPES2;  // set type of quneching. 0 means no quenching. 1 mean fixed adsorbate. 2 means adsorbate free to move only in z direction.
	}
	else A->GetPES2 = 0; // no quenching
}


void Adsorbate_Init_position(struct Adsorbate * A, size_t seed, double xmax, double ymax)
{
      // generate random initial x,y position of adsorbate
  
       const gsl_rng_type * T;
       gsl_rng * r;
          
       gsl_rng_env_setup();
     
       T = gsl_rng_default;
       r = gsl_rng_alloc (T);
       	
       gsl_rng_set(r,(unsigned long int)seed);

       double ux = gsl_rng_uniform (r);
       ux*=xmax;
       gsl_vector_set(A->Rt,0,ux);
       double uy = gsl_rng_uniform (r);
       uy*=ymax;
       gsl_vector_set(A->Rt,1,uy);         	
       gsl_rng_free (r); 
}  



void Init_Adsorbate(struct Adsorbate * A, struct Lattice *L, size_t seed)
{
    // Further initializing

    A->Adsorbate_Substrate = gsl_matrix_calloc(L->number_of_atoms ,3); // Distance between substrate atoms and adsorbates
    A->Rt = gsl_vector_calloc(3); // allocate memory for adsorbate's position
    A->RtGlobal = gsl_vector_calloc(3); // allocate memory for adsorbate's position over lattice
    A->dR = gsl_vector_calloc(3); // allocate memory for adsorbate's change in position
    A->Vt = gsl_vector_calloc(3); // allocate memory for adsorbate's velocity
    A->Ft = gsl_vector_calloc(3); // allocate memory for adsorbate's force
    
    if (A->GetPES2 == 0)
    {
      // set initial adsorbate's position
      Adsorbate_Init_position(A, seed, L->xmax, L->ymax);    
      gsl_vector_set(A->Rt,2,L->zmax+2.0); //adsorbate starts at (0,0) and 2A above substrate
      gsl_vector_memcpy(A->RtGlobal, A->Rt);
      double x = gsl_vector_get(A->RtGlobal,0);
      double y = gsl_vector_get(A->RtGlobal,1);
      double z = gsl_vector_get(A->RtGlobal,2);
      printf("Initial adsorbate's position is X = %g A\tY = %g A\tZ = %g A\n", x,y,z);
    }
      //printf("GetPES2 = %d\n",A->GetPES2);
      
    A->E = (double*)calloc(L->number_of_atoms ,sizeof(double));
    A->R = (double*)calloc(L->number_of_atoms ,sizeof(double));
    A->F = (double*)calloc(L->number_of_atoms ,sizeof(double));
    
    
    // A->AS is an array of pointers to the rows of the adsorbate-substrate distance vectors
    A->AS = (gsl_vector_view*)calloc( L->number_of_atoms ,sizeof(gsl_vector_view) );
    for (int k=0; k<L->number_of_atoms; k++)
      {
	A->AS[k] = gsl_matrix_row(A->Adsorbate_Substrate,k);
      }
      
    // allocate memory for other vectors
    A->R1 = gsl_vector_calloc(3);
    A->V1 = gsl_vector_calloc(3);
    A->F1 = gsl_vector_calloc(3);
    A->F2 = gsl_vector_calloc(3);
  
}



void Reset_Adsorbate(struct Adsorbate * A)

{
  // set force and velocity to zero
  gsl_vector_set_zero(A->Ft);
  gsl_vector_set_zero(A->F1);
  gsl_vector_set_zero(A->F2);
  gsl_vector_set_zero(A->Vt);
  gsl_vector_set_zero(A->V1);
  
}

void Free_Adsorbate(struct Adsorbate * A)
{
    // free dynamically allocate members
   gsl_matrix_free(A->Adsorbate_Substrate);
   gsl_vector_free(A->Rt);
   gsl_vector_free(A->RtGlobal);
   //gsl_vector_free(A->dR);
   //gsl_vector_free(A->Vt);
   gsl_vector_free(A->Ft);
   free(A->E);
   free(A->R);
   free(A->F);
   free(A->AS);
   gsl_vector_free(A->R1);
   gsl_vector_free(A->V1);
   gsl_vector_free(A->F1);
   gsl_vector_free(A->F2);
}




void Adsorbate_Substrate_Distance(struct Adsorbate * A, struct Substrate * S)
{

    // Calculate the distance vector between adsorbate and substrate
    // A->A1, A->A2 are temp matrices

    gsl_matrix_set_all(A->Adsorbate_Substrate,0); // reset distance matrix
    
    double x = gsl_vector_get(A->Rt,0); // get x coordinate of adsorbate
    double y = gsl_vector_get(A->Rt,1); // get y coordinate of adsorbate
    double z = gsl_vector_get(A->Rt,2); // get z coordinate of adsorbate

    // do for x,y,z

    A->A1 = gsl_matrix_column (A->Adsorbate_Substrate, 0); // put x coordinate of distances between adsorbate and substrate atoms
    A->A2 = gsl_matrix_column (S->Rt, 0); // put x coordinates of substrate atoms position
    // Add adsorbate x coordinate to A->A1
    gsl_vector_set_all(&A->A1.vector, x);
    // subtract x coordinates of substrate atoms
    gsl_vector_sub (&A->A1.vector, &A->A2.vector);
    //gsl_matrix_fprintf(stdout,A->Adsorbate_Substrate, "%g");

    A->A1 = gsl_matrix_column (A->Adsorbate_Substrate, 1); // put y coordinate of distances between adsorbate and substrate atoms
    A->A2 = gsl_matrix_column (S->Rt, 1); // put y coordinates of substrate atoms position
    // Add adsorbate y coordinate to A->A1
    gsl_vector_set_all(&A->A1.vector, y);
    // subtract y coordinates of substrate atoms
    gsl_vector_sub (&A->A1.vector, &A->A2.vector);


    A->A1 = gsl_matrix_column (A->Adsorbate_Substrate, 2); // put z coordinate of distances between adsorbate and substrate atoms
    A->A2 = gsl_matrix_column (S->Rt, 2); // put coordinates of substrate atoms position
    // Add adsorbate z coordinate to A->A1
    gsl_vector_set_all(&A->A1.vector, z);
    // subtract z coordinates of substrate atoms
    gsl_vector_sub (&A->A1.vector, &A->A2.vector);

}




void Get_Kinetic_Energy_Adsorbate(struct Adsorbate * A)
{

    // find current kinetic energy of adsorbate
    double d;
    gsl_blas_ddot(A->Vt,A->Vt,&d); // dot profuct of velocity vecotr
    A->KE = A->mass*d/2;
}



void Adsorbate_Dipole(register struct Adsorbate * A1, register struct Adsorbate * A2, struct Parameters* P)
{

    /* Calculate the dipole interaction between current adsorbate and all sufficiently near
       adsorbates.

       Interaction is 2*m^2/R^3 where m is the dipole of the adsorbate and r
       is the distance between the pair of interacting adsorbates.
       Therefore the Force is

               F = (6*m^2/|R|^5)*R

    */

    double R = AdsorbateAdsorbateDist(A1,A2);
    if (R<P->DipoleRange) // adsrobate are close enough
       {
            // scale distance vector by 6*m^2/|R|^5
        gsl_vector_scale(A1->InterAdDistance,6*A1->AdsorbateDipole*A2->AdsorbateDipole/gsl_pow_5(R));
            // add force to first adsorbate and subtract from the second
            gsl_vector_add(A1->Ft,A1->DipoleForce);
            gsl_vector_sub(A2->Ft,A1->DipoleForce);
       }
}



double AdsorbateAdsorbateDist(register struct Adsorbate * A1, register struct Adsorbate * A2)
{
    /* Calculate and return the distance R between adsorbates. */


    // copy position of first adsorbate to InterAdDistance
    gsl_vector_memcpy(A1->InterAdDistance, A1->RtGlobal);
    // subtract second adsorbate's position
    gsl_vector_sub(A1->InterAdDistance, A2->RtGlobal);
    // scale distance by R^-4
    return gsl_blas_dnrm2(A1->InterAdDistance);


}


void Adsorbate_Morse(register struct Adsorbate * A, struct Substrate * S, struct Lattice *L)
{

  /* Calculate the substrate-adsorbate interaction
   *
   *
   * The potential and force term is given by
   * 
   * V[r] = A*( exp(2*b*r0)*exp(-2*b*r[j]) - 2*exp(b*r0)*exp(-b*r[j]) )
   * 
   * Fj[r] = -2*b*A*( exp(2*b*r0)*exp(-2*b*r[j]) + exp(b*r0)*exp(-b*r[j]) )*r[j]/|r|
   *
   * where |r| is the magnitude of the distance between the adsorbtae and a substrate atom, and r[j] is the j'th
   * component of the distance vector.
   *
   * For the sake of speed, 2bA and exp(b*r0) are predefined as constants.
   *
   */

   A->PE = 0; // reset current potential energy
   register size_t l = L->number_of_atoms;
   register size_t N = L->number_of_bulk_atoms;
   //double register  r, F, expr;
   double register Amp = A->Force_Amplitude; // load constant to register
   double register preexp = A->preexp; // load constant to register
   double * restrict E = A->E; // will hold exp(-b*R)
   double * restrict Ra = A->R; // will hold R (magnitude of the vector R = A->AS)
   double * restrict f = A->F; // will hold the force vector
   double b = A->b;
   register double tmp;
   gsl_vector_view * restrict R = A->AS ; // pointer to first row of Adsorbate_Substrate distance matrix
   gsl_vector_view * restrict Ft = S->BulkForce; // pointer to first row of the force acting on the substrate atoms
   register gsl_vector_view *  end = &A->AS[l]; // pointer to last element
   gsl_vector_set_all(A->Ft,0); // set force acting on adsorbate to zero
   int register m = 0;

   // Norm of the distance
   for (; R<end; R++)
     {
      Boundary_shift(L,&R->vector);
     }

   // b*R - Norm of the distance times b
   for (R = A->AS; R<end; R++, Ra++) // loop over all the adsorbate-substrate distance vectors
     {
      *Ra = b*gsl_blas_dnrm2 (&R->vector);
     }
    

   //exp(-b*R)
   register double *END = Ra; // pointer to last element of Ra
   for ( Ra = A->R; Ra<END; E++,Ra++)
      {	
         *E = exp(-1*(*Ra));
      }


   // calculate the force (and potential energy if needed)
   register double *EEND = E;
   for (E = A->E, Ra = A->R; E<EEND ;E++,Ra++,f++)  // loop over all substrate atoms (loop over rows of A->Adsorbate_Substrate)
   {	
      tmp = (*E)*preexp; //exp(-b*(R-r0))
      // recall Amp=2*A*b*b, Ra = b|r|
      *f = -Amp*(tmp*tmp - tmp)/(*Ra); // -2*b*b*A*( exp(-2*b*(R-r0)) - exp(-b*(R-r0)) )/(|r|*b)
      if (A->savepot) // if we wish to track the potential energy along the adsorbate's trajectory
      {
	// notice we use here A->Amp and not Amp
	A->PE += A->Amp*(tmp*tmp - 2*tmp); // V = A( exp(2*b*r0)*exp(-2*b*R) - 2*exp(b*r0)*exp(-b*R) ) 
      }	  
   }

   
   // apply the force
   for (m=0,R = A->AS, f = A->F ; R<end ; R++, f++, Ft++)
   {
      double tmp = *f;
      //printf("Morse force is %g\n", *f);
      gsl_blas_daxpy (-tmp, &R->vector, A->Ft) ; // multiply vector r by F and assign to the force acting on the adsorbate
      if (m < N) {
      gsl_blas_daxpy (tmp, &R->vector, &Ft->vector) ; // assign the minus of the force to the forces acting on the substrate atom
      }
      m++;
      if (A->GetPES2 != 0) // quenching with z free.  keep just z component for the adsorbate
	{
	  gsl_vector_set(A->Ft,0,0);
	  gsl_vector_set(A->Ft,1,0);
	}
	  
   }
   
}


void Update_Velocity_Beeman_Adsorbate(struct Adsorbate* restrict A)
{
        // use force vector to calculate velocity using the Beeman algorithm

        // V(t+dt) = v(t) + 1/6 * (2*F(t+dt)+5*F*a(t)-F(t-dt))*dt/m

        // add to V(t) the other 3 terms on the right hand side

        gsl_blas_daxpy ((A->dt_over_m*2.0), A->Ft, A->Vt);
        gsl_blas_daxpy (A->dt_over_m*5.0, A->F1, A->Vt);
        gsl_blas_daxpy (-A->dt_over_m, A->F2, A->Vt);

}



void Quench_Velocity_Adsorbate(struct Adsorbate* restrict A)
{
  // if the scalar product between force acting on substrate atom and its
  // velocity is negative, set velocity to zero

	double b;
	// Scalar product between force and velocity
	gsl_blas_ddot (A->Ft, A->Vt, &b);
	if (b<0) gsl_vector_set_zero(A->Vt);
}






void Update_Position_Beeman_Adsorbate(struct Adsorbate* restrict A, struct Lattice* restrict L, struct Parameters* P)
{
        // use velocity vector to calculateposition using the Beeman algorithm

        // R(t+dt) = R(t) + v(t)*dt + 1/6 * (4*a(t)-a(t-dt))*dt^2/m
        // dR(t) = v(t)*dt + 1/6 * (4*a(t)-a(t-dt))*dt^2/m
        // integration over dR(t') from t'=0 to t'=t  == Displacement from equilibrium (since all atoms start at equilibrium)

        // set dR to 0
	gsl_vector_set_all(A->dR,0);
        gsl_blas_daxpy (P->Timestep, A->V1, A->dR);
        gsl_blas_daxpy (4.0*A->dt2_over_m, A->F1, A->dR);
        gsl_blas_daxpy (-A->dt2_over_m, A->F2, A->dR);

	// update global position and position over lattice  
	gsl_blas_daxpy (1, A->dR, A->Rt);
	gsl_blas_daxpy (1, A->dR, A->RtGlobal);

	// Apply boundary conditions to Rt

	Boundary_shift(L, A->Rt);
}

void Update_Temporaries_Adsorbate(struct Adsorbate* restrict A)
{
        // update R1,V1,F1,F2

        gsl_vector_memcpy(A->F2,A->F1);
        gsl_vector_memcpy(A->F1,A->Ft);
        gsl_vector_memcpy(A->V1,A->Vt);
        gsl_vector_memcpy(A->R1,A->Rt);
}


void Follow_Adsorbate(struct Adsorbate* restrict A)
{
   // print the attributes of adsorbate  and  track its dynamics

   // get potision ,velocity and force
  
   double x,y,z,vx,vy,vz,fx,fy,fz;
   static double Z=0,V=0;
   static size_t k=1;
   fx = gsl_vector_get(A->Ft,0);
   fy = gsl_vector_get(A->Ft,1);
   fz = gsl_vector_get(A->Ft,2);
   vx = gsl_vector_get(A->Vt,0);
   vy = gsl_vector_get(A->Vt,1);
   vz = gsl_vector_get(A->Vt,2);
   x = gsl_vector_get(A->RtGlobal,0);
   y = gsl_vector_get(A->RtGlobal,1);
   z = gsl_vector_get(A->RtGlobal,2);
   V = V+vz;
   Z = Z+z;
   
 

   printf("Adsorbate's force is  X = %g 	Y = %g		Z = %g \n", fx,fy,fz);
   printf("Adsorbate's velocity is  X = %g m/sec		Y = %g m/sec		Z = %g m/sec\n", 100*vx,100*vy,100*vz);
   printf("Average Z position  is %g\n" , Z/k);
   printf("Average Z velocity  is %g\n" , V/k); k++; 
   printf("Adsorbate's  position is  X = %f		Y = %f		Z = %f\n\n",x,y,z);
}



void Total_CM(struct Substrate * S,struct Adsorbate* restrict A, struct Lattice* restrict L)
{
    // CM motion of substrate+adsorbate
    
    double X,Y,Z,ax,ay,az;
    double cmx,cmy,cmz;
    gsl_vector_view x,y,z;
    ax = gsl_vector_get(A->Rt,0);
    ay = gsl_vector_get(A->Rt,1);
    az = gsl_vector_get(A->Rt,2);    
    
    gsl_vector * temp2 = gsl_vector_calloc(L->number_of_atoms);
    gsl_vector_set_all(temp2,1);

    x = gsl_matrix_subcolumn (S->Rt, 0,0,L->number_of_atoms); // velocities in the x coordinate
    gsl_blas_ddot(&x.vector,temp2,&X); // dot product with vectors of one gives sum of x velocities

    y = gsl_matrix_subcolumn (S->Rt, 1,0,L->number_of_atoms); // velocities in the y coordinate
    gsl_blas_ddot(&y.vector,temp2,&Y); // dot product with vectors of one gives sum of y velocities

    z = gsl_matrix_subcolumn (S->Rt, 2,0,L->number_of_atoms); // velocities in the z coordinate
    gsl_blas_ddot(&z.vector,temp2,&Z); // dot product with vectors of one gives sum of z velocities

    X = X/L->number_of_atoms;
    Y = Y/L->number_of_atoms;
    Z = Z/L->number_of_atoms;    
 
    cmx = (X*L->number_of_atoms*L->mass + ax*A->mass)/(L->number_of_atoms*L->mass + A->mass);
    cmy = (Y*L->number_of_atoms*L->mass + ay*A->mass)/(L->number_of_atoms*L->mass + A->mass);
    cmz = (Z*L->number_of_atoms*L->mass + az*A->mass)/(L->number_of_atoms*L->mass + A->mass);
    
    printf("\n***********\n");
    printf("CM position is X = %g A	Y = %g A	Z = %g A\n", cmx, cmy ,cmz);
    printf("***********\n");
    
    gsl_vector_free(temp2); // release memory
    
}



