#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Lattice.h"
#include "Substrate.h"
#include "Adsorbate.h"
#include "Simulation.h"
#include "Quench.h"
#include "Parameters.h"
#include <glib.h>
#include <glib/gstdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <assert.h>
#include <omp.h>
#include <unistd.h>


/*  Units are:
 *	[x] = Angstrom
 *	[m] = amu
 *	[t] = psec
 *	[T] = K
 *
 * 	1[x]/[t] = 100 m/sec	1m/sec = 0.01 [x]/[t]
 *	1[E] = 1.03584E-4 eV	1eV = 9654 [E]
 * 	1[F] = 1.66E-13 N       1N = 602.41E10 [F]
 *
 * 	k = 9036.1446 [F]/[x]
 *	Kb = 0.83192 [E]/[T]
 * 	A = 1303.3[E]
 * 	b = 0.875 1/[x]
 * 	R0 = 3.3 [x]
 *
 */

int main(int argc, char **argv)
{
	
	// Start simulation
	//int numcores = omp_get_num_procs(); // get number of available processors
	//omp_set_num_threads(numcores-1); // leave one core free
	struct Parameters * P; // declare parameters struct
	struct Lattice * restrict  lfcc ; // declaer lattice strut
	lfcc = (struct Lattice*)malloc(sizeof(struct Lattice)); // allocate memory
	P = (struct Parameters*)malloc(sizeof(struct Parameters)); // allocate memory
	if (argc<2) {fprintf(stderr,"\n*****\tPlease supply in the command line a Perameters file\t*****\n"); exit(1);}
	char * pfile = argv[1]; // pointer to name of parameters file
	
	
	// initialization: get parameters from file, generate lattice, find nearest neighbours...
	Init_Simulation(lfcc, P, argc, pfile);
	
	int temp, part; // counter of the list of temperatures and parts
	char cwd[1024];
	getcwd(cwd, sizeof(cwd)); // find the folder we work in
        if (P->Track) Copy_Param_File(pfile, P); // copy the parameters file to the folder where the output is saved, iכ output is saved
	

	if (strcmp(P->Simtype, "MD")==0) printf("\n*****\tRunning MD simulation\t*****\n");// do we do MD simulation>

        //Start timer

        P->StartGlobalTime = g_timer_new ();
        g_timer_start(P->StartGlobalTime);


	// loop over temperatures
	P->Start = g_date_time_new_now_local (); // time when simulation started
	for (temp = 0 ; temp < P->templen; temp++)
	{ 	
    
	  if (P->GetPES2) // extracting PES 
	  {
	    printf("Getting PES\n");
	    GetPES2(temp, lfcc, P);
	  }  
	  
	  
	  else // tracking adsorbates/substrate
	  {
	  
	    if (P->Track) P->Completedparts = PartsCompleted(P); // parts completed so far if we keep track
	    printf("Completed %d parts out of %d for temperature %dK\n",P->Completedparts, P->Parts, P->Temperatures[temp]);
	    #pragma omp parallel for schedule(static,1) private(part) shared(P,lfcc,temp)
	    for (part = 0; part < P->Parts; part++)
	    {		      
		 printf("Entering folder %s\n", P->Savedir);
		  g_chdir(P->Savedir); // enter the folder where data will be saved			
		  
		  // dont overwrite exisiting data file
		  // if datafile exist, skip to next part.
		  if (Check_Datafile_exist(P, temp, part)==1) continue;			
		  else MD(temp, part, lfcc, P); // run simulation

	     }		    
	  }
	} 
	 
	// return to original dir
	g_chdir(cwd); // return to original folder
        printf("\n\n***Finished***\n");
        return 0;
}
