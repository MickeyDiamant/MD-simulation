#pragma once

#include "Lattice.h"
#include "Substrate.h"
#include "Adsorbate.h"
#include <stdlib.h>
#include <glib.h>
#include "Parameters.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

typedef struct Adsorbate Adsorbate;
typedef struct Lattice Lattice;
typedef struct Substrate Substrate;
typedef struct Parameters Parameters;



void GetPES2(int temp, struct Lattice* restrict L, struct Parameters* P);

void Get_PES_XYZ(int k, int temp, gsl_vector *R, gsl_vector *PEZ, struct Lattice* restrict L, struct Parameters* P);

double Get_PES_XY(int k, int temp, gsl_vector *R, struct Lattice* restrict L, struct Parameters* P);

void Assign_Min_Potnetial(size_t j, size_t k, size_t l, gsl_matrix *PE, gsl_vector *PEZ);

void FreeZ(int k, int temp, struct Lattice* restrict L, struct Parameters* P, gsl_matrix *PE, gsl_matrix * Height );

void FixedZ(int temp, struct Lattice* restrict L, struct Parameters* P, gsl_matrix *PE );

int QuenchPE(struct Parameters* P, double TotalPE, double * PreviousTotalPE);

int QuenchKE(struct Parameters* P, double TotalKE, double * PreviousTotalKE);
