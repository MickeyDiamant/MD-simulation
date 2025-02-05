#pragma once

#include "Lattice.h"
#include "Substrate.h"
#include "Adsorbate.h"
#include <stdlib.h>
#include <glib.h>
#include "Quench.h"


typedef struct Adsorbate Adsorbate;
typedef struct Lattice Lattice;
typedef struct Substrate Substrate;


typedef struct Parameters {

    double Amass;
    double b;
    double R0;
    double Amp;
    double Dipole; // do we use dipole-dipole interaction
    double AdsorbateDipole; // dipole interaction parameter
    double DipoleRange; // cutoff radius for dipole-dipole interaction
    double Vel_Scale_factor;

    //double QuenchEnergy; // relax the system so each substrate atoms has no more than QuenchEnergy[E] = QuenchEnergy*1.03584E-4 eV

    int cpus; // how many cores to use
    
    int Track_Adsorbate;

    int Track_Adsorbate_Potential;

    int Track_Substrate;

    int QuenchTrack;
     
    int Track; // do we track something

    int Parts; // how many times run for each temperature

    int GetPES2; // do we want to extract the PES by freezing the adsorbate at place

    int GridXY; // how many grid points in the x,y direction

    int GridZ; // how many posints in the z direction
    
    int zinit; // initial z coordinate when quenching with free z

    int PESIterations; // over how many runs to average the PES

    GDateTime * Start; //start time of the simulation

    size_t templen; // how many tmeperatures to run over

    size_t Mins_To_Average; // average over this number of minimal values to get potential at point (x,y,z)

    size_t Temperature_Rescale_Rate; // number of time steps between temperature rescaling

    double Timestep; // Time step in unit of the used time unit

    double Elimit; // energy level that should be reached during quenching
    
    double shiftx; // shift of lattice to center origin over hollow site
    double shifty;

    size_t Simulation_Time; // number of time steps

    size_t Substrate_Relax; // number of time steps to relax substrate

    size_t Number_of_Adatoms; // How many adatoms

    size_t Substrate_Samplerate; // How often we sample the substrate position/velocity

    size_t Adsorbate_Samplerate; // How often we sample the adsorbates position/velocity

    size_t number_of_adatoms; // How many adatoms

    uint Seed;

    int Completedparts;

    int PDOS; // Do we want to calculate the phonon density of states (track substrate atoms velocities)

    char * Simtype; // Type of simulation, MD or Langevin

    char * Savedir; // Folder where data files will be saved

    GKeyFile * Datafile;
    // Holds the temperatures to run the simulation with
    int * Temperatures;
    // pointer to file to write parameters and remarks to
    FILE * Summary;

    // Timer
    GTimer *StartGlobalTime;


} Parameters;



void Get_parameters(char * datafile, struct  Parameters * P, struct Lattice * L);

void Get_Data_File(char * datafile, struct  Parameters * P);

void Copy_Param_File(char * datafile, struct  Parameters * P);

void Get_Time(GKeyFile * Key, struct  Parameters * P);

void Get_Geometry(GKeyFile * Key,struct Lattice * L);

void Get_Physics(GKeyFile * Key, struct Lattice * L, struct  Parameters * P);

void Get_Files(GKeyFile * Key, struct  Parameters * P);

void Get_Simulation(GKeyFile * Key, struct Lattice * L, struct  Parameters * P);
