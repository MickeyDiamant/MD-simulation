#pragma once
#ifndef LATTICE_INCLUDED
#define LATTICE_INCLUDED
#include <stdlib.h>
#include <glib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "Lattice.h"
#include "Substrate.h"
#include "Adsorbate.h"
#include "Simulation.h"
#include "Parameters.h"
#include "Quench.h"




typedef struct Adsorbate Adsorbate;
typedef struct Substrate Substrate;
typedef struct Parameters Parameters;
typedef struct Parameters Parameters;


typedef struct Lattice {



        // lattice constant in x direction
        double latconstx;
        // lattice constant in y direction
        double latconsty;
        // lattice constant in z direction
        double latconstz;
        // x coordinate boundary
        double xboundary;
        // y coordinate boundary
        double yboundary;
        // z coordinate boundary
        double zboundary;
        // max x coordinate
        double xmax;
        // min x coordinate
        double xmin;
        // max y coordinate
        double ymax;
        // min y coordinate
        double ymin;
        // max z coordinate
        double zmax;
        // min z coordinate
        double zmin;
        // force constant
        double Force_constant;
        // Boltzmann conts
        double Kb;

        // Lattice atom mass in AMU
        double mass;
	//
        // Velocity scaling related factor
        double Vel_Scale_factor;
        // time step divided by mass by 6
        double dt_over_m;
        // time step square divided by mass by 6
        double dt2_over_m;
        // do we wish to monitor the potential energy
        int Get_PE; // if 1, calculate potential energy each time step
        // number of lattice points in x direction
        int numcellx;
        // number of lattice points in y direction
        int numcelly;
        // number of lattice points in z direction
        int numcellz;
        // number of atoms in a single layer
        size_t num_of_atoms_in_layer;
        // number og layers
        size_t number_of_layers;
        // total number of atoms
        size_t number_of_atoms;
        // number of top layer atoms
        size_t number_of_top_atoms;
        // number of bottom layer atoms
        size_t number_of_bottom_atoms;
        // number of bulk atoms
        size_t number_of_bulk_atoms;
        // number of distinct neighbours pairs
        size_t number_of_pairs;
        // pointer to a matrix of the perfect lattice points
        gsl_matrix * LatticeRp;
        // pointer to distance vectors between nearest neighbours in a perfect lattice (no periodic boundary conditions)
        gsl_matrix * LatticeDp;
        // matrix which will hold the distance vectors between nearest neighbours
        gsl_matrix * NN_Distance_Vectors;
        // pointer to a vector which holds the indices of nearest neighbours
        GArray * GNN;
        // pointer to vector which holds repeated indices of atoms (as the number of neighbours they have)
        GArray * GForceatoms;
        // pointer to vector which holds indices of bottom layer atoms
        GArray * Bottom;
        // pointer to vector which holds indices of top layer atoms
        GArray * Top;
        // pointer to vector which holds indices of bulk atoms
        GArray * Bulk;
        // pointer to array to keep record of what type each atom is (8=top,12=bulk,0=bottom)
        GArray * Type;

	FILE * QuenchTracking; // pointer to file which hold tracking of the quenching process


        //size_t *Pairs;
        // temporary matrice
        gsl_matrix *DMtemp;
        // temporary vectors
        gsl_vector *Dtemp;	
        gsl_vector *Vtemp;
        gsl_vector *Vtemp1;
        gsl_vector_view Vtemp2;


} Lattice ;


void Pes2FCC(struct Lattice *Lattice);

void Distance_Bewteen_Atoms (struct Lattice* L);

void Test_NN_Distance(struct Lattice* L);

double Pair_Distance(const struct Lattice* L, const  size_t a, const  size_t b, size_t boundary_condition) ;

void Boundary_shift(const struct Lattice* L, gsl_vector * temp);

void Find_Top_Bottom_Bulk(struct Lattice* L);

void Find_Lattice_Extent(struct Lattice* L);

void Move_To_Center_Of_Mass(struct Lattice* L);

void GenFCC (struct Lattice *Lattice);

void Construct_Force_Matrix_Array(struct Lattice* L);

void Pair_Equilibrium_Distance(struct Lattice *L);

void Print_Neighbours_List(struct Lattice *L);

void Init_Lattice(struct Lattice* L, struct Parameters* P);

void Sort_Distance(struct Lattice* L, gsl_vector * distances, size_t k);


#endif // LATTICE_INCLUDED
