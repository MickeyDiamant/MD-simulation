#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "Lattice.h"

#include <omp.h>
#include <glib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_rng.h>
#define PI (3.141592653589793)


void Init_Substrate(struct Substrate* restrict S, struct Lattice * restrict L, struct Parameters* P)
{
	 
        // initialize substrate's data
        S->Rt = gsl_matrix_calloc(L->number_of_atoms,3); // position matrix
        S->R1 = gsl_matrix_calloc(L->number_of_atoms,3); // position matrix
        gsl_matrix_memcpy (S->Rt, L->LatticeRp); // initial position in equilibrium
        S->Displacement = gsl_matrix_calloc(L->number_of_pairs,3); // displacement matrix
	S->Temperature = gsl_vector_calloc((size_t)(P->Simulation_Time*P->Timestep)); // vector to hold temperature every picosec
        // allocate force matrices and set to zero
        S->Ft = gsl_matrix_calloc(L->number_of_atoms,3);
        S->F1 = gsl_matrix_calloc(L->number_of_atoms,3);
        S->F2 = gsl_matrix_calloc(L->number_of_atoms,3);
        // allocate velocity matrices and set to zero
        S->Vt = gsl_matrix_calloc(L->number_of_atoms,3);
        S->V1 = gsl_matrix_calloc(L->number_of_atoms,3);
	// center of mass position
	S->CMPosition = gsl_vector_calloc(3);
	// center of mass velocity
	S->CMVelocity = gsl_vector_calloc(3);
	// matrix view of the submatrix Ft which holds the force acting on the bottom atoms
	S->BottomForce = gsl_matrix_submatrix (S->Ft, L->number_of_bulk_atoms, 0, L->number_of_atoms-L->number_of_bulk_atoms, 3);
	// holds the products of the force and velocity components
	S->Quenchtemp = gsl_vector_calloc(3);
        // initiate velocity
	printf("Setting initial velocity\n");	
        Initial_Velocity(S, L) ;
	printf("Initial velocity set\n");

  
}



// void Reset_Substrate(struct Substrate* restrict S, struct Lattice * restrict L)
// {
//     // reset velocity and force of substrate atoms
//     gsl_matrix_memcpy (S->Rt, L->LatticeRp);
//     //gsl_matrix_set_zero(S->Displacement);
//     //gsl_matrix_set_zero(S->Vt);
//     //gsl_matrix_set_zero(S->V1);
//     //gsl_matrix_set_zero(S->Ft);
//     //gsl_matrix_set_zero(S->F1);
//     //gsl_matrix_set_zero(S->F2);
//     Initial_Velocity(S,  L, gsl_rng * Gauss) ;
//     Initial_Temperature_Rescale(S,L);
// }


void Free_Substrate(struct Substrate* restrict S)
{
  // free dynamically allocate members
  gsl_matrix_free(S->Rt);
  gsl_matrix_free(S->R1);
  gsl_matrix_free(S->Displacement);
  gsl_vector_free(S->Temperature);
  gsl_matrix_free(S->Ft);
  gsl_matrix_free(S->F1);
  gsl_matrix_free(S->F2);
  gsl_matrix_free(S->Vt);
  gsl_matrix_free(S->V1);
  gsl_vector_free(S->Quenchtemp);
  free(S->NN);
  free(S->Forceatoms);
  free(S->Force);
  free(S->ForceN);
  free(S->BulkForce);
  free(S->Atomic_Force_Matrix_Row);
  free(S->Disp);
  free(S->Velview);
  free(S->Forview);
  free(S->neighbour_index);  
}



void Garray_To_Array(struct Lattice *L, struct Substrate* restrict S)
{
        // accessing GArray is slow.
        // Replace some of the GArrays with regular C arrays which point to position of lattice atoms.
        // that is, arrays whose elements are pointers to rows of matrices.
        
        int m;
        size_t n;
	// pointer to position vectors of a substrate pair's neighbour atom
        S->NN = (gsl_vector_view*)calloc( L->number_of_pairs ,sizeof(gsl_vector_view) );
	// pointer to position vectosr of substrate pair's atom
        S->Forceatoms = (gsl_vector_view*)calloc( L->number_of_pairs ,sizeof(gsl_vector_view) );
	// pointer to force vectors of substrate pair's atom
        S->Force = (gsl_vector_view*)calloc( L->number_of_pairs ,sizeof(gsl_vector_view) );
	// pointer to force vectors of substrate pair's neighbour atom
        S->ForceN = (gsl_vector_view*)calloc( L->number_of_pairs ,sizeof(gsl_vector_view) );
	// pointer to rows of force acting on atoms
	S->BulkForce = (gsl_vector_view*)calloc( L->number_of_bulk_atoms ,sizeof(gsl_vector_view) );
	// pointer to a rows from which the Atomic force matrix can be constructed by outer product
        S->Atomic_Force_Matrix_Row = (gsl_vector_view*)calloc( L->number_of_pairs ,sizeof(gsl_vector_view) );
	// pointer to rows of displacemnt matrix
	S->Disp = (gsl_vector_view*)calloc( L->number_of_pairs ,sizeof(gsl_vector_view) );
	// pointer to rows of velocity matrix
	S->Velview = (gsl_vector_view*)calloc( L->number_of_atoms ,sizeof(gsl_vector_view) );
	// pointer to rows of force matrix
	S->Forview = (gsl_vector_view*)calloc( L->number_of_atoms ,sizeof(gsl_vector_view) );	
	// array of k items, where item k is the index of the neighbour atoms in pair k
	S->neighbour_index = (size_t*)calloc( L->number_of_pairs ,sizeof(size_t) );

	for (m=0;m<L->number_of_bulk_atoms; m++)
	{
	    S->BulkForce[m]  = gsl_matrix_row (S->Ft,m);
	      // if quenching is later requested these pointer point to rows of the force and velocity matrix for each individual atom (no pairs involved)
	      S->Velview[m] = gsl_matrix_row (S->Vt, m); // velocity vector of atom
	      S->Forview[m] = gsl_matrix_row (S->Ft, m); // velocity vector of atom	    
	} 

        for (m=0; m<L->number_of_pairs; m++) { // loop over all pairs. In each pair assign pointers to force and position vecotr of both atoms in pair
                n = g_array_index(L->GNN, size_t, m); // index of neighbour atom in pair m
		S->neighbour_index[m] = n;
                S->NN[m] =  gsl_matrix_row (S->Rt, n); // Position vector of neighbour. Neighbour could be a bottom atom!
                S->ForceN[m] = gsl_matrix_row (S->Ft, n); // Force vector of neighbour. Neighbour could be a bottom atom!
                n = g_array_index(L->GForceatoms, size_t, m); // index of atom in pair m
                assert(n<L->number_of_bulk_atoms);
		S->Forceatoms[m] = gsl_matrix_row (S->Rt,  n); // Position vector of atom
                S->Force[m] = gsl_matrix_row (S->Ft, n); // Force vector of atom 
                S->Atomic_Force_Matrix_Row[m] = gsl_matrix_row (L->NN_Distance_Vectors,m);
		S->Disp[m] = gsl_matrix_row (S->Displacement,m);

        }
}



gsl_rng * GaussianDist(uint seed)
{
	//generates and returns a Gaussian random number generator
        const gsl_rng_type * T;
        gsl_rng * r;
	
	/* create a generator chosen by the
          environment variable GSL_RNG_TYPE */
        gsl_rng_env_setup();
	gsl_rng_set(r,seed);
        T = gsl_rng_default;
        r = gsl_rng_alloc(T); 
	
	printf("Generating Gaussian random distribution\n");
	return r;
}


void Initial_Velocity(struct Substrate* restrict S, struct Lattice* restrict L)
{
        // Initiate velocity using a Boltzmann distribution
        // working in units of eV for energy, psec for time, Kelvin for temperature and AMU for mass
        // The velocity is given in units of 98.17614 * Angstrom / Picosec, so we divide by this numerical factor.


        // generate seed
	
        const gsl_rng_type * T;
        gsl_rng * r;
	//printf("Temperature is %d\n", S->temperature);
        size_t i, s,  n = L->number_of_bulk_atoms;
        double sigma = sqrtf(L->mass/(PI*L->Kb*S->temperature));

	/* create a generator chosen by the
          environment variable GSL_RNG_TYPE */

        gsl_rng_env_setup();

        T = gsl_rng_default;
        r = gsl_rng_alloc (T);
        gsl_rng_set(r,S->seed);

	/* print n random variates chosen from
          the poisson distribution with mean
          parameter mu */

	for (i = 0; i < n; i++) {  // run over atoms
                for (s=0; s<3; s++)  {  // run ove directions
                        gsl_matrix_set( S->Vt, i, s, gsl_ran_gaussian(r, sigma)/sigma ); // Get a Gaussian distribution
                }
        }  

      Zero_Velocity_CM(S, L); // move to center of mass system
      gsl_rng_free(r);
}





void Zero_Velocity_CM(struct Substrate* restrict S, struct Lattice* restrict L)
{
  // Moves to cnter of mass
  gsl_vector_view x,y,z;
  gsl_vector * temp = gsl_vector_calloc(L->number_of_bulk_atoms);
  gsl_vector_set_all(temp,1);
  double Vx=0, Vy=0, Vz=0;
  
  x = gsl_matrix_subcolumn (S->Vt, 0,0,L->number_of_bulk_atoms); // velocities in the x coordinate
  gsl_blas_ddot(&x.vector,temp,&Vx); // dot product with vectors of one gives sum of x velocities

  y = gsl_matrix_subcolumn (S->Vt, 1,0,L->number_of_bulk_atoms); // velocities in the y coordinate
  gsl_blas_ddot(&y.vector,temp,&Vy); // dot product with vectors of one gives sum of y velocities

  z = gsl_matrix_subcolumn (S->Vt, 2,0,L->number_of_bulk_atoms); // velocities in the z coordinate
  gsl_blas_ddot(&z.vector,temp,&Vz); // dot product with vectors of one gives sum of z velocities

  Vx = Vx/L->number_of_bulk_atoms;
  Vy = Vy/L->number_of_bulk_atoms;
  Vz = Vz/L->number_of_bulk_atoms;
  
  // change velocity by subtracting the components of the CM velocity

  gsl_vector_add_constant (&x.vector, -Vx);
  gsl_vector_add_constant (&y.vector, -Vy);
  gsl_vector_add_constant (&z.vector, -Vz);

  // in order to verify we check the Vcm is indeed 0

  x = gsl_matrix_subcolumn (S->Vt, 0,0,L->number_of_bulk_atoms); // velocities in the x coordinate
  gsl_blas_ddot(&x.vector,temp,&Vx); // dot product with vectors of one gives sum of x velocities

  y = gsl_matrix_subcolumn (S->Vt, 1,0,L->number_of_bulk_atoms); // velocities in the y coordinate
  gsl_blas_ddot(&y.vector,temp,&Vy); // dot product with vectors of one gives sum of y velocities

  z = gsl_matrix_subcolumn (S->Vt, 2,0,L->number_of_bulk_atoms); // velocities in the z coordinate
  gsl_blas_ddot(&z.vector,temp,&Vz); // dot product with vectors of one gives sum of z velocities

  gsl_vector_free(temp); // release memory
}


void Initial_Temperature_Rescale(struct Substrate* restrict S, struct Lattice* restrict L)
{
	// set initial substrate atoms velocities to a Bolztmann distribution
        Get_Kinetic_Energy_Substrate(S,L);
        Get_Temperature(S,L);
        gsl_matrix_scale(S->Vt,sqrtf(2*S->temperature/S->Current_Temp)); // Scale Velocity. Factor 2 because half of KE goes to PE.
        Get_Kinetic_Energy_Substrate(S,L);
        Get_Temperature(S,L);
}




void Harmonic_Force_Index2(register struct Substrate*  S, register struct Lattice*  L)
{
        // Calculate the harmonic force between nearest neighbours:
        // 1. Get the displacement between pair of nearest neighbours
        // 2. Scalar product between displacement and row of atomic force matrix
        // 3. Multiply row of atomic force matrix with above scalar reuslt to get force vector

        
        // Set all forces to zero
        gsl_matrix_set_all(S->Ft,0);
        // set potential energy to zero
        S->PE = 0;
	register size_t l = L->number_of_bulk_atoms;
	register unsigned int k;
        // set displacement to equilibrium values    
        gsl_matrix_memcpy(S->Displacement,L->LatticeDp); // negative of equilibrium distance vector
 
       
	for (k = 0; k<L->number_of_pairs; k++){Force(S, L,  k, l);}
}	  



void Harmonic_Force_Index1(register struct Substrate*  S, register struct Lattice*  L)
{
        // Calculate the harmonic force between nearest neighbours:
        // 1. Get the displacement between pair of nearest neighbours
        // 2. Scalar product between displacement and row of atomic force matrix
        // 3. Multiply row of atomic force matrix with above scalar reuslt to get force vector

        double b;
        // Set all forces to zero
        gsl_matrix_set_all(S->Ft,0);
        // set potential energy to zero
        S->PE = 0;
	register size_t l = L->number_of_bulk_atoms;
	register unsigned int k;
       // set displacement to equilibrium values    
       gsl_matrix_memcpy(S->Displacement,L->LatticeDp); // negative of equilibrium distance vector
       //gsl_matrix_scale(S->Displacement,-1);
       // Using pointer arithmetic instead of indices to run over arrays
       gsl_vector_view  nn;
       gsl_vector_view  Fatoms;
       gsl_vector_view  afmr;
       gsl_vector_view  f;
       gsl_vector_view  fN;
       gsl_vector_view   disp;
       register size_t * nindex = S->neighbour_index;
       register gsl_vector_view *  end = &S->NN[L->number_of_pairs];
    
       // loop over all pairs

	  for (k = 0; k<L->number_of_pairs; k++)
         
         { 	  
               // Displacement

               Fatoms = S->Forceatoms[k];
	       disp = S->Disp[k];
	       nn = S->NN[k];
               gsl_blas_daxpy (1, &Fatoms.vector, &disp.vector); // add position of atom to displacement vector
               gsl_blas_daxpy (-1, &nn.vector, &disp.vector); // subtract position of neighbour from displacement vector
  	       
	       // Scalar product

	       afmr = S->Atomic_Force_Matrix_Row[k];
               gsl_blas_ddot (&afmr.vector, &disp.vector, &b);

	       // Update force

	       f = S->Force[k];
	       fN = S->ForceN[k];
               gsl_blas_daxpy (-b, &afmr.vector, &f.vector); // update atom's force
	       if (*nindex < l){ // if neighbour atom in the pair is not a bottom atom, apply force on it
               gsl_blas_daxpy (b, &afmr.vector, &fN.vector); // update neighbour's force
	       }
 	}

}



void Force(register struct Substrate*  S, register struct Lattice*  L, size_t k, size_t l)
       {
	 	 
	  gsl_vector_view  nn;
	  gsl_vector_view  Fatoms;
	  gsl_vector_view  afmr;
	  gsl_vector_view  f;
	  gsl_vector_view  fN;
	  gsl_vector_view   disp;
	  register size_t * nindex = S->neighbour_index;
	  double b;
	  Fatoms = S->Forceatoms[k];
	  disp = S->Disp[k];
	  nn = S->NN[k];
	  gsl_blas_daxpy (1, &Fatoms.vector, &disp.vector); // add position of atom to displacement vector
	  gsl_blas_daxpy (-1, &nn.vector, &disp.vector); // subtract position of neighbour from displacement vector
	  
	  // Scalar product

	  afmr = S->Atomic_Force_Matrix_Row[k];
	  gsl_blas_ddot (&afmr.vector, &disp.vector, &b);

	  // Update force

	  f = S->Force[k];
	  fN = S->ForceN[k];
	  gsl_blas_daxpy (-b, &afmr.vector, &f.vector); // update atom's force
	  if (*nindex < l)
	      { // if neighbour atom in the pair is not a bottom atom, apply force on it
	            gsl_blas_daxpy (b, &afmr.vector, &fN.vector); // update neighbour's force
	      }	 
	 
       }



void Harmonic_Force(register struct Substrate*  S, register struct Lattice*  L)
{
	// Calculate the harmonic force between nearest neighbours:
	// 1. Get the displacement between pair of nearest neighbours
	// 2. Scalar product between displacement and row of atomic force matrix
	// 3. Multiply row of atomic force matrix with above scalar reuslt to get force vector

	double b;
	// Set all forces to zero
	gsl_matrix_set_all(S->Ft,0);
	// set potential energy to zero
	S->PE = 0;
	register size_t l = L->number_of_bulk_atoms;
       // set displacement to equilibrium values
       gsl_matrix_memcpy(S->Displacement,L->LatticeDp); // negative of equilibrium distance vector
       //gsl_matrix_scale(S->Displacement,-1);
       // Using pointer arithmetic instead of indices to run over arrays
       register gsl_vector_view * nn = S->NN;
       register gsl_vector_view * Fatoms = S->Forceatoms;
       register gsl_vector_view * afmr = S->Atomic_Force_Matrix_Row;
       register gsl_vector_view * f = S->Force;
       register gsl_vector_view * fN = S->ForceN;
       register gsl_vector_view * disp = S->Disp;
       register size_t * nindex = S->neighbour_index;
       register gsl_vector_view *  end = &S->NN[L->number_of_pairs];
       
       // loop over all pairs
       for (; nn<end; nn++,Fatoms++,afmr++,f++,fN++,disp++,nindex++)
        {
               // Displacement
	       // displacement = (x-y)-(x0-y0)
	       // disp is -(x0-y0), thats why L->LatticeDp was multiplyed by -1 is Lattice.c
               gsl_blas_daxpy (1, &Fatoms->vector, &disp->vector); // add position of atom to displacement vector
               gsl_blas_daxpy (-1, &nn->vector, &disp->vector); // subtract position of neighbour from displacement vector
	       
	       // Scalar product between displacement from eq. distance vector and a vecotr whose outer product with itself gives the atomic force matrix
               gsl_blas_ddot (&afmr->vector, &disp->vector, &b);
	       
	       // add to potential energy
	       S->PE += b*b; 

	       // Update force
               gsl_blas_daxpy (-b, &afmr->vector, &f->vector); // update atom's force
	       if (*nindex < l) // if neighbour atom in the pair is not a bottom atom, apply force on it
	       { 
		  gsl_blas_daxpy (b, &afmr->vector, &fN->vector); // update neighbour's force
	       }
	}
}

void Get_Temperature(struct Substrate* restrict S, struct Lattice* restrict L )
{
        // find the temperature from kinetic energy per atom
        S->Current_Temp = 2*(S->KE_Per_Atom)/( 3*L->Kb); // T = 2E/3Kb [Kb] = eV/K
}

void Get_Kinetic_Energy_Substrate(struct Substrate* restrict S, struct Lattice* restrict L)
{
        //Take the current velocity matrix and find the KE of the substrate
        double register d=0,KE=0;
        gsl_vector_view x,y,z;

        // get columns of matrix
        x = gsl_matrix_subcolumn (S->Vt, 0,0,L->number_of_atoms); // velocities in the x coordinate
        d = gsl_blas_dnrm2 (&x.vector); // d is the L2-norm of x
        KE += d*d; // add square of norm
        y = gsl_matrix_subcolumn (S->Vt, 1,0,L->number_of_atoms);
        d = gsl_blas_dnrm2 (&y.vector);
        KE += d*d;
        z = gsl_matrix_subcolumn (S->Vt, 2,0,L->number_of_atoms);
        d = gsl_blas_dnrm2 (&z.vector);
        KE += d*d;

        // multiply a proportionality factor [A^2][amu]/[psec^2] = 1.0364269E-04 eV
        S->KE = KE*((L->mass)/2);// units of [E] 
        S->KE_Per_Atom = S->KE/L->number_of_bulk_atoms;
	

}

void Scale_Velocity(struct Substrate* restrict S, struct Lattice* restrict L, struct Parameters* P)
{
        // Scale the velocity using a Berendsen thermostat to reach a desired temperature
        double Scale_factor = sqrtf(1 + (P->Timestep / P->Vel_Scale_factor)*(S->temperature/S->Current_Temp - 1)); // Calculate scale factor
        //printf("Scale factor is %g\n", Scale_factor);
        if (S->Current_Temp/S->temperature < 1) {
                assert(Scale_factor>=1);
        } else {
                assert(Scale_factor<=1);
        }
        gsl_matrix_scale(S->Vt, Scale_factor); // Scale velocity
}



void Print_CM(struct Substrate* restrict S, struct Lattice* restrict L)
{
    // Prints current center of mass

	size_t k;
        double Xsum = 0;
        double Ysum = 0;
        double Zsum = 0;

        for (k=0; k<L->number_of_atoms; k++) {
                Xsum += gsl_matrix_get(S->Rt,k,0);
                Ysum += gsl_matrix_get(S->Rt,k,1);
                Zsum += gsl_matrix_get(S->Rt,k,2);

        }

        printf("Xcm = %f\tYcm = %f\tZcm = %f\n", Xsum/(L->number_of_atoms),Ysum/(L->number_of_atoms),Zsum/(L->number_of_atoms));      
}



void Linear_Sum_Matrices(double a, gsl_matrix * restrict  M2,  gsl_matrix * restrict M1, struct Substrate* restrict S)
{
        // due to lack of a routine to sum matrices.
        // computes M1 => M1+a*M2 (row by row)


        size_t k;

        //get number of rows
        size_t row1 = M1->size1;
        size_t row2 = M2->size1;
        //get number of columns
        size_t col1 = M1->size2;
        size_t col2 = M2->size2;
        // assert row1=row2, col1=col2
        assert(col1==col2);
        assert(row1==row2);

        for (k=0; k<col1; k++) {

		assert (k<3);
                // get row k of both matrices
                S->Vtemp2 = gsl_matrix_column (M2, k);
                S->Vtemp1 = gsl_matrix_column (M1, k);
                // sum the two rows
                gsl_blas_daxpy (a, &S->Vtemp2.vector, &S->Vtemp1.vector); // m1 => m1+a*m2
        }

}



void Update_Velocity_Beeman(struct Substrate* restrict S, struct Lattice* restrict L)
{
        // use force matrices to calculate velocity using the Beeman algorithm

        // V(t+dt) = v(t) + 1/6 * (2*F(t+dt)+5*F*a(t)-F(t-dt))*dt/m

        // add to V(t) the other 3 terms on the right hand side

        Linear_Sum_Matrices ((L->dt_over_m*2.0), S->Ft, S->Vt,S);
        Linear_Sum_Matrices (L->dt_over_m*5.0, S->F1, S->Vt,S);
        Linear_Sum_Matrices (-1*L->dt_over_m, S->F2, S->Vt,S);

}


void Quench_Velocity_Substrate(struct Substrate* restrict S, register struct Lattice*  L)
{
  // if the scalar product between force acting on substrate atom and its
  // velocity is negative, set velocity to zero

	//double b;

       // Using pointer arithmetic instead of indices to run over arrays
	
       register gsl_vector_view * F = S->Forview;
       register gsl_vector_view * V = S->Velview;
       register gsl_vector_view * end = &S->Velview[L->number_of_bulk_atoms];
       
       // loop over all pairs
       for (; V<end ; V++, F++)
        {
	       // product between cmponents of velocity an force
               //gsl_blas_ddot (&F->vector, &V->vector, &b);
               gsl_vector_set_all(S->Quenchtemp,1); // set to a vector of 1
	       gsl_vector_mul(S->Quenchtemp,&F->vector);
	       gsl_vector_mul(S->Quenchtemp,&V->vector);
	       for (int j=0; j<3; j++) // for each vecotr component
	       {  // if product with force component is negative set to zero
		  if ( gsl_vector_get(S->Quenchtemp,j) < 0) gsl_vector_set(&V->vector,j,0);
	       }
	}
}  



void Update_Position_Beeman(struct Substrate* restrict S, struct Lattice* restrict L, struct Parameters* P)
{
        // use velocity matrices to calculate position using the Beeman algorithm

        // R(t+dt) = R(t) + v(t)*dt + 1/6 * (4*a(t)-a(t-dt))*dt^2/m


        Linear_Sum_Matrices (P->Timestep, S->V1, S->Rt,S);
        Linear_Sum_Matrices (4.0*L->dt2_over_m, S->F1, S->Rt,S);
        Linear_Sum_Matrices (-L->dt2_over_m, S->F2, S->Rt,S);

}

void Update_Temporaries(struct Substrate* restrict S)
{
        // update R1,V1,V2,F1,F2

        gsl_matrix_memcpy(S->F2,S->F1);
        gsl_matrix_memcpy(S->F1,S->Ft);
        gsl_matrix_memcpy(S->V1,S->Vt);
        gsl_matrix_memcpy(S->R1,S->Rt);
        //gsl_matrix_memcpy(S->dR1,S->dR);

        //printf("Updated temporaries\n");
}

 void Follow_Random_Atom(struct Substrate* restrict S,  struct Lattice* restrict L, int pair)
 {
   // print the attributes of a single pair of atoms and to track its dynamics

   // get potision and velocity
   //gsl_vector_view atomF = S->Forceatoms[pair]; // atom's force
   gsl_vector_view atomP = gsl_matrix_row (S->Rt,pair); // atom's position
   gsl_vector_view atomV = gsl_matrix_row (S->Vt,pair); // atom's position
   //gsl_vector_view atom1 = S->NN[pair]; // atom's first neighbour's position
   //int n=0;
   //double x,y,z,vx,vy,vz,fx,fy,fz;
   double x,y,z,vx,vy,vz;
   static double X=0,V=0;
   static size_t k=1;
 //   x = gsl_matrix_get(S->Rt, atom,0);
 //   y = gsl_matrix_get(S->Rt, atom,1);
 //   z = gsl_matrix_get(S->Rt, atom,2);
    vx = gsl_vector_get(&atomV.vector,0);
    vy = gsl_vector_get(&atomV.vector,1);
    vz = gsl_vector_get(&atomV.vector,2);
   x = gsl_vector_get(&atomP.vector,0);
   y = gsl_vector_get(&atomP.vector,1);
   z = gsl_vector_get(&atomP.vector,2);
   V = V+vx;
   X = X+x;
   
   // Magnitude of veloicty
//    fx = gsl_vector_get(&atomF.vector,0);
//    fy = gsl_vector_get(&atomF.vector,1);
//    fz = gsl_vector_get(&atomF.vector,2);
   //if (fx*fxold<0) {n++; printf("completed %d periods\n", n/2);}
   //fxold = fx;
 //   fy = gsl_matrix_get(S->Ft, atom,1);
 //   fz = gsl_matrix_get(S->Ft, atom,2);

   //printf("\nAtom %d's position is  X = %g		Y = %g		Z = %g\n", atom,x,y,z);
   printf("Atom %d's velocity is  X = %g m/sec		Y = %g m/sec		Z = %g m/sec\n", pair,100*vx,100*vy,100*vz);
   printf("Average X position  is %g\n" , X/k);
   printf("Average X velocity  is %g\n" , V/k); k++; 
    printf("Atom %d's  position is  X = %f		Y = %f		Z = %f\n\n", pair,x,y,z);
    //printf("Atom %d's force is  FX = %g		FY = %g		FZ = %g\n\n", pair,fx,fy,fz);
 //
 //   x = gsl_matrix_get(S->Rt, atom1,0);
 //   y = gsl_matrix_get(S->Rt, atom1,1);
 //   z = gsl_matrix_get(S->Rt, atom1,2);
 //   vx = gsl_matrix_get(S->Vt, atom1,0);
 //   vy = gsl_matrix_get(S->Vt, atom1,1);
 //   vz = gsl_matrix_get(S->Vt, atom1,2);
 //   fx = gsl_matrix_get(S->Ft, atom1,0);
 //   fy = gsl_matrix_get(S->Ft, atom1,1);
 //   fz = gsl_matrix_get(S->Ft, atom1,2);
 //
 //   //printf("\nAtom %d's position is  X = %g		Y = %g		Z = %g\n", atom1,x,y,z);
 //   //printf("Atom %d's velocity is  X = %g		Y = %g		Z = %g\n", atom1,vx,vy,vz);
 //   printf("Atom %d's force is  X = %g		Y = %g		Z = %g\n\n", atom1,fx,fy,fz);
 }


 void CM_Velocity(struct Substrate* restrict S,  struct Lattice* restrict L)
 {

    // Prints out the position and velocity of the center of mass

  gsl_vector_view x,y,z;
  gsl_vector * temp1 = gsl_vector_calloc(L->number_of_bulk_atoms);
  gsl_vector_set_all(temp1,1);
  double Vx=0, Vy=0, Vz=0;

  x = gsl_matrix_subcolumn (S->Vt, 0,0,L->number_of_bulk_atoms); // velocities in the x coordinate
  gsl_blas_ddot(&x.vector,temp1,&Vx); // dot product with vectors of one gives sum of x velocities

  y = gsl_matrix_subcolumn (S->Vt, 1,0,L->number_of_bulk_atoms); // velocities in the y coordinate
  gsl_blas_ddot(&y.vector,temp1,&Vy); // dot product with vectors of one gives sum of y velocities

  z = gsl_matrix_subcolumn (S->Vt, 2,0,L->number_of_bulk_atoms); // velocities in the z coordinate
  gsl_blas_ddot(&z.vector,temp1,&Vz); // dot product with vectors of one gives sum of z velocities

  Vx = Vx/L->number_of_bulk_atoms;
  Vy = Vy/L->number_of_bulk_atoms;
  Vz = Vz/L->number_of_bulk_atoms;

  gsl_vector_set(S->CMVelocity,0,Vx);
  gsl_vector_set(S->CMVelocity,1,Vy);
  gsl_vector_set(S->CMVelocity,2,Vz);
  
  gsl_vector_free(temp1);
 }
 
void CM_Position(struct Substrate* restrict S,  struct Lattice* restrict L)
{
  
  // ***************************

  gsl_vector_view x,y,z;
  gsl_vector * temp2 = gsl_vector_calloc(L->number_of_atoms);
  gsl_vector_set_all(temp2,1);
  double X=0, Y=0, Z=0;

  x = gsl_matrix_subcolumn (S->Rt, 0,0,L->number_of_atoms); // velocities in the x coordinate
  gsl_blas_ddot(&x.vector,temp2,&X); // dot product with vectors of one gives sum of x velocities

  y = gsl_matrix_subcolumn (S->Rt, 1,0,L->number_of_atoms); // velocities in the y coordinate
  gsl_blas_ddot(&y.vector,temp2,&Y); // dot product with vectors of one gives sum of y velocities

  z = gsl_matrix_subcolumn (S->Rt, 2,0,L->number_of_atoms); // velocities in the z coordinate
  gsl_blas_ddot(&z.vector,temp2,&Z); // dot product with vectors of one gives sum of z velocities

  X = X/L->number_of_atoms;
  Y = Y/L->number_of_atoms;
  Z = Z/L->number_of_atoms;
  
  gsl_vector_set(S->CMPosition,0,X);
  gsl_vector_set(S->CMPosition,1,Y);
  gsl_vector_set(S->CMPosition,2,Z);  

//   static double V;
//   static size_t k=1;
//   V = V + Vx;
//   printf("***********\n");
//   printf("\nCM velocity is Vx = %g m/sec	Vy = %g m/sec	Vz = %g m/sec\n", 100*Vx,100*Vy,100*Vz);
//   printf("\nAverage X velocity of CM is %g\n", 100*V/k); k++;
//   printf("\nCM position is X = %g A	Y = %g A	Z = %g A\n", X ,Y ,Z);
//   printf("***********\n");
  gsl_vector_free(temp2); // release memory
  
}


void Sample_Temperature(struct Substrate* restrict S, struct Lattice* restrict L, size_t t)
{
  
  // Writes to an array the current temperature every picosec
  Get_Kinetic_Energy_Substrate(S,L);
  Get_Temperature(S,L);  
  gsl_vector_set(S->Temperature, t, S->Current_Temp) ; 

}



void Temperature_Dist(struct Substrate* restrict S, int p)

{
  // calculates the mean and variance of the temperature as a function of time
  
  // find mean and square of values
  
  size_t *numTsteps = &S->Temperature->size; // number of samples
  double mean = 0, mean2 = 0; // mean and square mean
  for (size_t t=0; t<*numTsteps; t++)
  {
    mean += gsl_vector_get(S->Temperature,t);
  }
  mean = mean/(*numTsteps);
  mean2 = mean*mean;
  
  double var = 0, std = 0;
  
  for (size_t t=0; t<*numTsteps; t++)
  {
    var += gsl_vector_get(S->Temperature,t)*gsl_vector_get(S->Temperature,t)-mean2  ;      
  }
  
  std = sqrtf(var/(*numTsteps-1));
  
  printf("%dK: The average temperature for part %d was %.6fK with standard deviation of %.6fK\n",S->temperature, p, mean, std);
  
  
  
}