#pragma once


#include <stdlib.h>
#include <glib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "Simulation.h"
#include "Lattice.h"
#include "Substrate.h"
#include "Adsorbate.h"
#include "Quench.h"



typedef struct Adsorbate Adsorbate;
typedef struct Lattice Lattice;
typedef struct Parameters Parameters;




typedef struct Substrate  {

        // Temperature
        int temperature;
        // seed for random number generator
        uint seed;
        // pointer to the current location of atoms
        gsl_matrix  * Rt;
        // pointer to the previous location of atoms
        gsl_matrix  * R1;
        // change in position matrix
        gsl_matrix  * dR;
        // pointer to displacement matrix
        //gsl_matrix  * Displacement;
        // pointer to the current force matrix
        gsl_matrix  * Ft;
        // pointers to auxillary force matrices
        gsl_matrix  * F1;
        gsl_matrix  * F2;
        // pointer to the current velocity matrix
        gsl_matrix  * Vt;
        // change in position matrix
        gsl_matrix  * dV;
	// Displacement
	gsl_matrix * Displacement;
        // pointers to auxillary velocity matrices
        gsl_matrix  * V1;
	// vectors of center of mass position and velocity
	gsl_vector * CMPosition;
	gsl_vector * CMVelocity;
	gsl_vector * Quenchtemp;
	// pointer to temperature vector
        gsl_vector  * Temperature;
        // pointer to kinetic energy vector
        gsl_vector  * Kinetic_Energy;
        // pointer to potential energy vector
        gsl_vector  * Potential_Energy;
        // current kinetic energy
        double KE;
        // current kinetic energy per atom
        double KE_Per_Atom;
        // current potential energy
        double PE;
        // curent temperature
        double Current_Temp;
	// commulative temperature
	double Com_Temp;

	FILE * AdTraj; // pointer file which holds adsorbates trajectories
	FILE * SUBSfile; // pointer to file which holds substrate atoms position
	FILE * SUBSVEL;	// pointer to file which holds substrate atoms velocities
	
	
        // pointer to array which holds pointers to position of atoms (as the number of neighbours they have)
        gsl_vector_view * restrict Forceatoms;
        // pointer to array which holds pointers to position of nearest neighbours
        gsl_vector_view * restrict NN;
        // pointer to array which holds pointers to force on atoms (as the number of neighbours they have)
        gsl_vector_view * restrict Force;
        // pointer to array which holds pointers to force on nearest neighbours
        gsl_vector_view * restrict ForceN;
        // pointer to vector of atomic force matrix
        gsl_vector_view * restrict Atomic_Force_Matrix_Row;
	// pointer to rows of displacement matrix
	gsl_vector_view * restrict Disp;
// 	// pointer to rows of force acting on bulk and top atoms (no bottom atoms)
	gsl_vector_view * BulkForce;
	// pointer to rows of individual atoms velocities and forces for use with quenching
	gsl_vector_view * Forview;
	gsl_vector_view * Velview;
	// Pointer to force acting on bottom atoms
	gsl_matrix_view BottomForce;
	// array which holds index number of the neighbour atom in each pair
	size_t * neighbour_index;

	
        // Temporary vectors and matrices

      
        gsl_vector_view View1;
        gsl_vector_view View2;
        gsl_vector_view View3;
        gsl_vector_view View4;
        gsl_vector_view Vtemp;
        gsl_vector_view Vtemp1;
        gsl_vector_view Vtemp2;
        gsl_vector_view Vtemp3;
        gsl_vector_view VF;
        gsl_vector_view VF1;
        gsl_vector_view VF2;
        gsl_matrix * Quench;

} Substrate;


void Init_Substrate(struct Substrate* restrict S, struct Lattice * restrict L, struct Parameters* P);

void Pair_Displacement(struct Substrate* restrict S, struct Lattice * restrict L);

gsl_rng * GaussianDist(uint seed);

void Harmonic_Force(struct Substrate* restrict S, struct Lattice* restrict L);

void Harmonic_Force_Index1(register struct Substrate*  S, register struct Lattice*  L);

void Harmonic_Force_Index2(register struct Substrate*  S, register struct Lattice*  L);

void Calculate_Force(struct Substrate * restrict S, struct Lattice * restrict L, size_t p,int register n);

void Initial_Velocity(struct Substrate* restrict S, struct Lattice* restrict L);

void Get_Kinetic_Energy_Substrate(struct Substrate* restrict S, struct Lattice* restrict L);

void Get_Temperature(struct Substrate* restrict S, struct Lattice* restrict L);

void Scale_Velocity(struct Substrate* restrict S, struct Lattice* restrict L, struct Parameters* P);

void Update_Velocity_Beeman(struct Substrate* restrict S, struct Lattice* restrict L);

void Update_Position_Beeman(struct Substrate* restrict S, struct Lattice* restrict L, struct Parameters* P);

void Update_Temporaries(struct Substrate* restrict S);

void Linear_Sum_Matrices(double a, gsl_matrix * M2,  gsl_matrix * M1, struct Substrate* restrict S);

void Follow_Random_Atom(struct Substrate* restrict S,  struct Lattice* restrict L, int pair);

void Initial_Temperature_Rescale(struct Substrate* restrict S, struct Lattice* restrict L);

void Garray_To_Array(struct Lattice *L,struct Substrate* restrict S);

void Zero_Velocity_CM(struct Substrate* restrict S, struct Lattice* restrict L);

void CM_Velocity(struct Substrate* restrict S,  struct Lattice* restrict L);

void CM_Position(struct Substrate* restrict S,  struct Lattice* restrict L);

void Force(register struct Substrate*  S, register struct Lattice*  L, size_t k, size_t l);

void Print_CM(struct Substrate* restrict S, struct Lattice* restrict L);

void Sample_Temperature(struct Substrate* restrict S, struct Lattice* restrict L, size_t t);

void Temperature_Dist(struct Substrate* restrict S, int p);

void Quench_Velocity_Substrate(struct Substrate* restrict S, register struct Lattice*  L);

void Free_Substrate(struct Substrate* restrict S);

//void Reset_Substrate(struct Substrate* restrict S, struct Lattice * restrict L);