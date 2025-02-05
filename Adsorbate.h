#pragma once

#include "Lattice.h"
#include "Simulation.h"

#include "Substrate.h"
#include <stdlib.h>
#include <glib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "Parameters.h"
#include "Quench.h"



typedef struct Lattice Lattice;
typedef struct Substrate Substrate;
typedef struct Parameters Parameters;
typedef struct Parameters Parameters;



typedef struct Adsorbate {

     // vector of exponential
    double * E;
    double * F;
    gsl_matrix * restrict Trajectory; // array to hold trajectory
    // vector of R
    double *R;


    // Morse interaction related parameters
    double Force_Amplitude;
    double preexp;
    double b;
    double dt_over_m;
    double dt2_over_m;
    double init; // (x,y) initial coordinate
    double Amp;
    double R0;

    // Dipole interaction
    double AdsorbateDipole;
    // pointer to force and distance vectors
    gsl_vector  * DipoleForce;
    gsl_vector  * InterAdDistance;


    double mass;
    // current kinetic energy
    double KE;
    // current potential energy
    double PE;      

    // number of times the trajectory is sampled
    size_t num_of_samples;        

    // do we save the potential due to substrate
    int savepot;
    // do we want to extract the PES by freezing the adsorbate at place
    int GetPES2;        
    
    // pointer to vector of current position
    gsl_vector  * Rt; // position in over lattice
    gsl_vector  * RtGlobal; // position in space
    // pointer to the previous location of atom
    gsl_vector  * R1;
    // change in position vector
    gsl_vector  * dR;
    // pointer to the current force vector
    gsl_vector  * Ft;
    // pointers to auxillary force matrices
    gsl_vector  * F1;
    gsl_vector  * F2;
    // pointer to the current velocity vector
    gsl_vector  * Vt;
    // pointers to auxillary velocity matrices
    gsl_vector  * V1;
    // Distance matrix from substrate atoms
    gsl_matrix * Adsorbate_Substrate;
    // pointer to rows of Adsorbate_Substrate matrix
    gsl_vector_view * AS;
    
    
    // File to save potential to
    FILE * Potential_File;
    // Auxillary vector views
    gsl_vector_view A1;
    gsl_vector_view A2;
    gsl_vector_view A3;
    gsl_vector_view A4;

} Adsorbate;

// Generate adsorbate
void Gen_Adsorbate(struct Adsorbate * A, struct Lattice *L, struct Parameters *P);

void Adsorbate_Init_position(struct Adsorbate * A, size_t seed, double xmax, double ymax);

// Initiate Adosrbate
void Init_Adsorbate(struct Adsorbate * A, struct Lattice *L, size_t seed);

// get distance matrix between adsorbate and substrate
void Adsorbate_Substrate_Distance(struct Adsorbate * A, struct Substrate * S);

// The Adsorbate-Substrate Force_Amplitude
void Adsorbate_Morse(register struct Adsorbate * A, struct Substrate * S, struct Lattice *L);
// update adsorbate's velocity
void Update_Velocity_Beeman_Adsorbate(struct Adsorbate* restrict A);

// update adsorbat's position
void Update_Position_Beeman_Adsorbate(struct Adsorbate* restrict A, struct Lattice* restrict L, struct Parameters* P);

void Update_Temporaries_Adsorbate(struct Adsorbate* restrict A);

void Follow_Adsorbate(struct Adsorbate* restrict A);

void Total_CM(struct Substrate * S,struct Adsorbate* restrict A, struct Lattice* restrict L);

void Adsorbate_Fixed_Position(struct Adsorbate * A, double x, double y, double z);

void Free_Adsorbate(struct Adsorbate * A);

void Quench_Velocity_Adsorbate(struct Adsorbate* restrict A);

void Get_Kinetic_Energy_Adsorbate(struct Adsorbate * A);

void Reset_Adsorbate(struct Adsorbate * A);

void Adsorbate_Dipole(register struct Adsorbate * A1, register struct Adsorbate * A2, struct Parameters* P);

double AdsorbateAdsorbateDist(register struct Adsorbate * A1, register struct Adsorbate * A2);

