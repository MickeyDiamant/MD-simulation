#pragma once

#include "Lattice.h"
#include "Substrate.h"
#include "Adsorbate.h"
#include <stdlib.h>
#include <glib.h>
#include "Parameters.h"
#include "Quench.h"


typedef struct Adsorbate Adsorbate;
typedef struct Lattice Lattice;
typedef struct Substrate Substrate;
typedef struct Parameters Parameters;



void Adsorbate_Timestep(struct Adsorbate * Adatoms[], struct Substrate * S, struct Lattice* restrict L, struct Parameters* P);

void Substrate_Timestep(struct Adsorbate ** Adatoms, struct Substrate * S, struct Lattice* restrict L, struct Parameters* P);

void Substrate_Timestep_No_Adsorbates(struct Substrate * S, struct Lattice* restrict L, struct Parameters* P);

void Init_Simulation(struct Lattice* restrict L, struct Parameters* P, size_t argc, char * paramfile);

char *  Gen_Filename(char * name, int temp, int part, char * suffix);

void MD(int temp, int part, struct Lattice* restrict L, struct Parameters* P);

void Openfiles(int part, struct Substrate * S, struct Parameters* P, struct Adsorbate ** Adatoms, struct Lattice* restrict L);

void Track(struct Substrate * S, struct Parameters * P,   struct Lattice* restrict L, struct Adsorbate ** Adatoms, size_t j);

void Closefiles(struct Substrate * S, struct Parameters* P, struct Adsorbate ** Adatoms, int part);

void Relax(struct Substrate * S, struct Lattice* restrict L, struct Parameters* P,struct Adsorbate ** Adatoms );

char * Get_Date();

void Track_Temp(struct Substrate * S, struct Lattice* restrict L, size_t k);

//void Show_Time( size_t j, size_t display, struct Parameters* P, int temperature, int part, int temp);
void Show_Time( size_t t, struct Parameters* P);

int Check_Datafile_exist(struct Parameters * P, int temp, int part);

int PartsCompleted(struct Parameters * P);



