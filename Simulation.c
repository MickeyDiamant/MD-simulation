
#include <glib.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "Simulation.h"
#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <fnmatch.h>
#include <omp.h>
#include <gsl/gsl_sort_vector.h>




void Init_Simulation(struct Lattice* restrict L, struct Parameters* P, size_t argc, char * paramfile)

{
    // Get parameters and set up lattice

      if (argc>1) Get_Data_File(paramfile, P);
      printf("starting\n");

      Get_parameters(paramfile, P, L);
	      
      printf("initializing lattice\n");

      Init_Lattice(L, P);
      printf("Force constant is %g\n", L->Force_constant);
      printf ("Number of atoms is %zu\n",L->number_of_atoms);
  
}




void Adsorbate_Timestep(struct Adsorbate * Adatoms[], struct Substrate * S, struct Lattice* restrict L, struct Parameters* P)
{
	    // Finds the force acting on the adsorbate and calculate its new position and velocity
	    for (int i=0; i<P->number_of_adatoms; i++)
	    {
	      Update_Position_Beeman_Adsorbate(Adatoms[i],L,P);
	      Adsorbate_Substrate_Distance(Adatoms[i], S);
	      Adsorbate_Morse(Adatoms[i], S, L);
              if (P->AdsorbateDipole) // do we include adsorbate dipole-dipole interaction
                    {
                        for (int j=i+1; j<P->number_of_adatoms; j++) // loop over adsorbates pairs
                        {
                            Adsorbate_Dipole(Adatoms[i],Adatoms[j], P);
                        }
                    }
	      Update_Velocity_Beeman_Adsorbate(Adatoms[i]);
	      Update_Temporaries_Adsorbate(Adatoms[i]);
	    }	     
 
}   
   
   
char *  Gen_Filename(char * name, int temp, int part, char * suffix)
{
    // returns pointer to a string which serves as the name of the file
    // into which results of type 'name' at temperatures 'temp' are saved
  

    char * filename = malloc(100);
    if (part <0)
    {
      int len = sprintf(filename, "%s-%d.%s",name,temp,suffix);
      free(filename);
      filename = malloc(len+10);
      len = sprintf(filename, "%s-%d.%s",name,temp,suffix);
    }
    else
    {
      int len = sprintf(filename, "%s-%dK-%d.%s",name,temp,part,suffix);
      free(filename);
      filename = malloc(len+10);
      len = sprintf(filename, "%s-%dK-%d.%s",name,temp,part,suffix);
    }
    return filename;
}  


   
void Substrate_Timestep(struct Adsorbate ** Adatoms, struct Substrate * S, struct Lattice* restrict L, struct Parameters* P)
{
	// Dynamic of substrate and adsorbates
	// if P->GetPES2 == 1 quench system with adsorbate fixed
	// if P->GetPES2 == 2 quench, but adsorbate free to move in z direction. Use PE to determine if we reached minimum.
	// if P->GetPES2 == 3 quench, but adsorbate free to move in z direction. Use KE to determine if we reached minimum.

	Update_Position_Beeman(S,L,P);
	Harmonic_Force(S,L); // calculate the force between substrate atoms	
	// Adsorbates dynamics
	if (P->GetPES2 == 1) // if we want to extract the PES by freezing the adsorbate at place
	    
	    {
		for (int i=0; i<P->number_of_adatoms; i++)
		  {
		      Adsorbate_Substrate_Distance(Adatoms[i], S);
		      Adsorbate_Morse(Adatoms[i], S, L); // also calculates the force the adsorbate applys on the substrate atoms		
		  }
	    }
	    
	else // free dynamics or quenching with free z coordinate of adsorbate
            {
		for (int i=0; i<P->number_of_adatoms; i++)
		  {
		      Update_Position_Beeman_Adsorbate(Adatoms[i],L,P);
		      Adsorbate_Substrate_Distance(Adatoms[i], S);
		      Adsorbate_Morse(Adatoms[i], S, L); // also calculates the force the adsorbate applys on the substrate atoms
		      Update_Velocity_Beeman_Adsorbate(Adatoms[i]);
		      if (P->GetPES2 > 1) Quench_Velocity_Adsorbate(Adatoms[i]); // is asking for quenching
		      Update_Temporaries_Adsorbate(Adatoms[i]);
		  }    
	    }
	    
		
        Update_Velocity_Beeman(S,L); 
	if (P->GetPES2 > 0) Quench_Velocity_Substrate(S,L); // quench (set velocity to zero if its scalar product with the force is negative)
	Update_Temporaries(S);
}


void Substrate_Timestep_No_Adsorbates(struct Substrate * S, struct Lattice* restrict L, struct Parameters* P)
{
	// Dynamic of substrate only (if no adsorbate is found)

	Update_Position_Beeman(S,L,P);
	Harmonic_Force(S,L); // calculate the force (and potential energy is asked)

	Update_Velocity_Beeman(S,L);
	Update_Temporaries(S);
}


void MD(int temp, int part, struct Lattice* restrict L, struct Parameters* P)
{
	// generate adsrobate trajectories using MD
      
	// Generate substrate trajectories using MD if required
	//gsl_rng * Gauss = GaussianDist();
	struct Substrate * subs ;
	subs = (struct Substrate*)malloc(sizeof(struct Substrate));
	subs->temperature = P->Temperatures[temp];
	printf("Simulation temperature is %d K\n", subs->temperature);
	printf("Initializing substrate\n");
	subs->seed = P->Seed + (uint)P->Temperatures[temp] + (uint)part;
	printf("Seed is %u\n", subs->seed);
	
	Init_Substrate(subs, L, P);
	Garray_To_Array(L,subs);
	
	
	// Generate Adsorbates if required
	struct Adsorbate * Adatoms[P->number_of_adatoms];	  
	if (P->number_of_adatoms>0)
	{	
	      
	      for (int i=0; i<P->number_of_adatoms; i++) {
		  // allocate memory for adsorbate 
		  printf("Initializing adsrobate #\%d\n",i);
		  Adatoms[i] = (struct Adsorbate*)malloc(sizeof(struct Adsorbate));
		  // Generata adsorbate
		  Gen_Adsorbate(Adatoms[i],L,P);
		  // Initialize data structure		  
		  Init_Adsorbate(Adatoms[i], L, P->Seed+i);
	      }
	    
	      if (P->number_of_adatoms>1) {
		    assert(Adatoms[0]!=Adatoms[1]);
		    assert(&Adatoms[0]->Trajectory != &Adatoms[1]->Trajectory);
	      }
	}

	
	// Set initial temperature to desired value
	Initial_Temperature_Rescale(subs,L);
	//double olddt=P->Timestep;
	//P->Timestep=0.01;
	// Relaxation stage
	Relax(subs, L, P, Adatoms );
	printf("Finished initial relaxing\n");
	//P->Timestep=olddt;
	
/* Dynamics */	
	  
	// Start taking data

	//Open files to save adsorbate's trajectory and other data
	#pragma omp barrier
	Openfiles(part, subs, P, Adatoms, L);

	// TIMELOOP
	
	size_t display = 1000000/L->number_of_atoms; // how often to update info about time of simulation
	size_t track_temp = (size_t)(1/P->Timestep); // how often to take temperature
	size_t track_time = 0; // how many times tracked temperature 
	size_t total_temp_times = (size_t)(P->Simulation_Time*P->Timestep); // how many times to track temperautre
	printf("\n\r\r");
	#pragma omp barrier
	if (P->number_of_adatoms>0) // include adsorbate in simulation
	{   	    
	    // main loop
	    for (size_t j=0; j < P->Simulation_Time; j++) 
	    {	      
		  
	        if (j%display==1)
		{
		  if (omp_get_thread_num()==0)
		  {
                    Show_Time( j/display, P );
		    // Follow_Adsorbate(Adatoms[0]); // for debuggin purpose
		    //printf("\r  %c ","|/-\\"[(j/display)%4]);
	
		    fflush(stdout);  
		  }

		    
		}
		
 		Substrate_Timestep(Adatoms, subs, L, P); // dynamics
		
		// track temperature
		if ( (j%track_temp == 0) & (track_time < total_temp_times) ) 
		{
		  track_time += 1;
		  Sample_Temperature(subs, L, j/track_temp); // take temperature every picosec
		}  
		
		// Tracking trajectories
		if (P->Track==1)   Track(subs, P, L, Adatoms, j)  ;
		

		// Track adsorbate's position and velocity, CM position and velocity, total kinetic and total potental energy
	    }		
	}
	
	else // no adsorbate in simulation
	  
	{
	    // the following lines are to ensure the time related output of the threads don't overlap
	    int steps = 0;
	    while (steps < omp_get_thread_num())
	    {
	      printf("\n");
	      steps++;
	    }
	    
	    // main loop
	    for (size_t j=0; j < P->Simulation_Time; j++) 
	    {
		//if (j%display==0) Show_Time(j, display, P, subs->temperature, part, temp);
	    
		Substrate_Timestep_No_Adsorbates(subs, L, P);
		// Tracking
		if (P->Track==1)   Track(subs, P, L, NULL, j)  ;
		//Track_Temp(subs, L);
	    }	  
	}

	// Temperature statistics. done serialy cause we write to same file
	printf("\nTemperature Statistics:\n");
	Temperature_Dist(subs, part);

	// Close open files if any
	Closefiles(subs, P, Adatoms, part) ;
	#pragma omp barrier
}  





//void Show_Time( size_t j, size_t display, struct Parameters* P, int temperature, int part , int temp)
//{
//	// finds the elapsed time, and time left (based on how many time steps were completed so far)
//	GDateTime * Current;
//	GTimeSpan Elapsed;
//	size_t hour, minute, sec, fhour, fminute, fsec;
//	//double finish, rate, perpart, totaltime, totalpassed, timeleft; // in sec;
//	//double ratio;
//	Current = g_date_time_new_now_local();
//	Elapsed = g_date_time_difference (Current,P->Start);
//	Elapsed/=1000000; // return elapsed time in seconds
//	sec = Elapsed%60;
//	minute = (Elapsed/60)%60;
//	hour = (Elapsed/3600);
//	//rate = (double)j/Elapsed; // timesteps per second
//	//perpart = (double)P->Simulation_Time/rate; // sec
//	//totaltime = ( (double)(P->templen*P->Parts) )*perpart; // total time is the time for one part times the total number
//										    // of parts times the total number of temperatures

//	//totalpassed = (P->Parts*temp + P->Completedparts)*perpart + Elapsed; // temp is an integer index which starts at 0
//	//timeleft = (totaltime-totalpassed);
//	//ratio = 100*(double)totalpassed/(double)totaltime;

//	//assert(ratio>=0);
	
//// 	fsec = ((size_t)timeleft)%60;
//// 	fminute = (((size_t)timeleft)/60)%60;
//// 	fhour = ((size_t)timeleft)/3600;
	
	
//	printf("\r  %c\tElapsed: %zuh:%zum:%zus\t%c                       ",\
//	"|/-\\"[(j/display)%4],hour,minute,sec, "|/-\\"[(j/display)%4]);
	
//	fflush(stdout);

//}


void Show_Time(size_t t, struct Parameters* P)
{
    /* claculates the elapsed time */

    double T = g_timer_elapsed (P->StartGlobalTime, NULL); // number of elapsed seconds

    // elapsed time
    size_t Telapsed = (size_t)T;
    int seconds = (Telapsed)%60;
    int minutes = (Telapsed/60)%60;
    int hours = (Telapsed/3600)%24;
    int days = Telapsed/(24*3600);

    printf("\r  %c\tElapsed: %dd:%dh:%dm:%ds\t%c                    ",\
    "|/-\\"[t%4],days,hours,minutes,seconds, "|/-\\"[t%4]);

    fflush(stdout);

}

void Relax(struct Substrate * S, struct Lattice* restrict L, struct Parameters* P,struct Adsorbate ** Adatoms )
{
    // Relaxation procedure using Berrendsen thermostat

	double tempratio = 10;
	double relaxtemp = 0; // will hold comulative average temperature
	size_t relax_time = 0;
	size_t j;

	for (j=0; tempratio > 1.001 || tempratio < 0.999 ; j++,relax_time++)
	{
		// keep loop going as long as the ration of average temperature and required temperature is not close enough to 1

		if (P->number_of_adatoms>0)
		      {Substrate_Timestep(Adatoms, S, L, P);}
		else
		      {Substrate_Timestep_No_Adsorbates(S, L, P);}
		      
		Get_Kinetic_Energy_Substrate(S,L);
		Get_Temperature(S,L);
		if (j>P->Substrate_Relax) // start keeping average temperature after initial fluctuation decayed
		{
		  relaxtemp += S->Current_Temp;
		  tempratio = (relaxtemp/(j-P->Substrate_Relax)) / S->temperature; // keep sum of temperatures
		}

		if ( (fabs(S->Current_Temp - S->temperature)>0.0001)  &  (j%P->Temperature_Rescale_Rate == 1)  )
		{
			Scale_Velocity(S,L,P);
		}
		//if (j%1000==1) Follow_Adsorbate(Adatoms[0]); // for debuggin purpose
 	}



	printf("Average temperature after rescaling for %g psec is %f K\n",j*P->Timestep, relaxtemp/(relax_time-1-P->Substrate_Relax) );
	
}  



void Openfiles(int part, struct Substrate * S, struct Parameters* P, struct Adsorbate ** Adatoms, struct Lattice* restrict L)
{
    // open all required files to save data
    
    if (P->Track_Adsorbate > 0)
    {			
	    printf("Opening adsorbate trajectory file part %d for temperature %dK\n",part,S->temperature);
	    char * filename = Gen_Filename("MD-TRAJ", S->temperature, part, "tmp");
	    //printf("Opening adsorbate trajectory file %s.\n",filename);
	    S->AdTraj = fopen (filename,"w");
	    // print some relevant parameters to file
	    fprintf(S->AdTraj,"%d\n",7); // later will be used to know at which line the adsorbate trajectory starts
	    fprintf(S->AdTraj,"%zu\n",P->number_of_adatoms); // number of adsorbate
	    fprintf(S->AdTraj,"%g\n",P->Timestep*P->Adsorbate_Samplerate); // time between sampling of trajectory
	    fprintf(S->AdTraj,"%zu\n",P->Simulation_Time/P->Adsorbate_Samplerate); // number of trajectory steps
	    fprintf(S->AdTraj,"%d\n",S->temperature); // simulation temperature
	    fprintf(S->AdTraj,"%d\n",P->Parts); // total number of simulation parts
	    fprintf(S->AdTraj,"%u\n",S->seed); // initial random seed
    }				

    
    //Open file to save substrate atom's trajectory
    if (P->Track_Substrate > 0)
    {	
	printf("Opening susbtrate trajectory file for temperature %dK\n",S->temperature);
	char * filename = Gen_Filename("MD-SUBS", S->temperature, -1, "tmp");
	S->SUBSfile = fopen (filename,"w");
    }		

    
    //Open file to save substrate atom's trajectory
    if (P->PDOS > 0)
    {	
	printf("Opening susbtrate velocities file for temperature %dK\n",S->temperature);
	char * filename = Gen_Filename("SUBSVEL", S->temperature, -1, "tmp");
	S->SUBSVEL = fopen (filename,"w");
	fprintf(S->SUBSVEL,"%zu\n",L->number_of_bulk_atoms);
	fprintf(S->SUBSVEL,"%zu\n",L->number_of_atoms);
	fprintf(S->SUBSVEL,"%g\n",P->Timestep*P->Substrate_Samplerate);
	fprintf(S->SUBSVEL,"%zu\n",P->Simulation_Time);
	fprintf(S->SUBSVEL,"%d\n",S->temperature);
    }

      
    if (P->Track_Adsorbate_Potential > 0)
    {
      printf("Opening adsorbate potential file part %d for temperature %dK\n",part, S->temperature);
      char * filename = Gen_Filename("MD-AdsorbatePotential", S->temperature, part,"tmp");
      Adatoms[0]->Potential_File = fopen (filename,"w");
      fprintf(Adatoms[0]->Potential_File,"%d\n",6); // later will be used to know at which line the adsorbate trajectory starts
      fprintf(Adatoms[0]->Potential_File,"%g\n",P->Timestep*P->Adsorbate_Samplerate);
      fprintf(Adatoms[0]->Potential_File,"%zu\n",P->Simulation_Time);
      fprintf(Adatoms[0]->Potential_File,"%d\n",S->temperature);
      fprintf(Adatoms[0]->Potential_File,"%d\n",P->Parts);
      fprintf(Adatoms[0]->Potential_File,"%u\n",S->seed);
    }
    

}



void Track(struct Substrate * S, struct Parameters* P,   struct Lattice* restrict L, struct Adsorbate ** Adatoms, size_t j)
{
  
    // Track adsorbates and substrate atoms
  
    if (P->Track_Adsorbate > 0)
    {	// check it is time to save position again	
	  if (j%P->Adsorbate_Samplerate == 0) {
	    // Keep trajectory of adsorbates
	    for (int a = 0; a<P->number_of_adatoms; a++){
	      // Write to file position of each adsorbate, up to P->Track_Adsorbate adsorbates
	      gsl_vector_fprintf (S->AdTraj, Adatoms[a]->RtGlobal, "%.8f");
	      } 
	  }  
    }
  
    // Substrate atoms trajectories
    
    if (P->Track_Substrate > 0)
    {	// check it is time to save position again	
	    if (j%P->Substrate_Samplerate == 0) {
	      // Keep trajectory of substrate atoms
		gsl_matrix_fprintf (S->SUBSfile, S->Rt, "%.8f");
		    
	    }
    }	  
            
    // Substrate atoms velocities
    
    if (P->PDOS > 0)
    {	// check it is time to save position again	
	    if (j%P->Substrate_Samplerate == 0) {
	      // Keep trajectory of substrate atoms
		gsl_matrix_fprintf (S->SUBSVEL, S->Vt, "%.8f");
		    
	    }
    }	  	      
      
    // Adsorbate's potential energy

    if (P->Track_Adsorbate_Potential > 0)
    {	// check it is time to save position again
	  if (j%P->Adsorbate_Samplerate == 0) {

	      // Write to file position of each adsorbate, up to P->Track_Adsorbate adsorbates
	      fprintf (Adatoms[0]->Potential_File, "%.8f\n", Adatoms[0]->PE );
	      

	      

	  }
    }    
}    
	  
	  
	  
	  
	  
	  
void Closefiles(struct Substrate * S, struct Parameters* P, struct Adsorbate ** Adatoms, int part)

{
      // close open data files and change the suffix from "tmp" to "txt"
      if (P->Track_Adsorbate > 0)
	{
	  fclose(S->AdTraj); // close <dsorbate trajectory file
	  char * oldname = Gen_Filename("MD-TRAJ", S->temperature, part, "tmp");
	  char * newname = Gen_Filename("MD-TRAJ", S->temperature, part, "txt");
	  if (rename(oldname,newname)!=0) printf("Error! Could not rename tmp file %s to %s",oldname, newname);
	}				
	  
      if (P->Track_Substrate > 0)
	{	
	  fclose(S->SUBSfile); // close substrate trajectory file
	  char * oldname = Gen_Filename("MD-SUBS", S->temperature, -1, "tmp");
	  char * newname = Gen_Filename("MD-SUBS", S->temperature, -1, "txt");
	  if (rename(oldname,newname)!=0) printf("Error! Could not rename tmp file %s to %s",oldname, newname);	  
	}		

      if (P->PDOS > 0)
	{	
	  fclose(S->SUBSVEL); // close substrate velocities file
	  char * oldname = Gen_Filename("SUBSVEL", S->temperature, -1, "tmp");
	  char * newname = Gen_Filename("SUBSVEL", S->temperature, -1, "txt");
	  if (rename(oldname,newname)!=0) printf("Error! Could not rename tmp file %s to %s",oldname, newname);	  	  
	}
	
      if (P->Track_Adsorbate_Potential > 0)
	{

	  printf("Closing potential file for part %d\n",part);
	  fclose(Adatoms[0]->Potential_File); // close adsorbate potential file
	  char * oldname = Gen_Filename("MD-AdsorbatePotential", S->temperature, part, "tmp");
	  char * newname = Gen_Filename("MD-AdsorbatePotential", S->temperature, part, "txt");
	  if (rename(oldname,newname)!=0) printf("Error! Could not rename tmp file %s to %s",oldname, newname);	  
	}
	

}	  


void Track_Temp(struct Substrate * S, struct Lattice* restrict L, size_t k)
{
  // tracks comulative average temperature
  if (k==0) S->Com_Temp=0; // initialize
  
  Get_Kinetic_Energy_Substrate(S, L);
  Get_Temperature(S, L);
  S->Com_Temp += S->Current_Temp;
  
  if (k%1000==0) printf("%dK:  Commulative average temperature after relaxation is %g K\n",S->temperature,S->Com_Temp/k);
   
}




int Check_Datafile_exist(struct Parameters * P, int temp, int part)
{  
      char * datafile = malloc(100);
      int len = sprintf(datafile,"MD-TRAJ-%dK-%d.txt",P->Temperatures[temp],part);
      free(datafile);
      datafile = malloc(len);
      sprintf(datafile,"MD-TRAJ-%dK-%d.txt",P->Temperatures[temp],part);
      
      DIR *dir;
      struct dirent *ent;
      dir = opendir (P->Savedir);
      if (dir != NULL) 
      {
	/* print all the files and directories within directory */
	while ((ent = readdir (dir)) != NULL) 
	{
	  if ( fnmatch (datafile, ent->d_name, 0) == 0 ) 
	  {
	    closedir (dir);
	    printf("Data file part %d for temperature %dK already exist. Skipping.\n",part,P->Temperatures[temp]);
	    return 1;
	    
	  } 
	}
	
      }
      return 0; 
}	    


int PartsCompleted(struct Parameters * P)
{
    // Finds out how many parts were completed by counting the number txt files

      DIR *dir = opendir (P->Savedir);
      struct dirent *ent;
      int count = 0;
      if (dir != NULL) 
      {
	/* print all the files and directories within directory */
	while ((ent = readdir (dir)) != NULL) 
	{
	  if ( strcmp("Summary.txt", ent->d_name) == 0) continue; // if file is "Summary.txt" ignore
	  else if ( strstr ( ent->d_name, "txt") != NULL ) count++; // if "txt" is found in the file's name increase counter by 1
	  else continue;
	}
	closedir (dir);
      }
      return count;
}
