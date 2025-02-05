
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
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>










void GetPES2(int temp, struct Lattice* restrict L, struct Parameters* P)
{
    // Find the minimum potential at each (x,y) point over the lattice unit

    int k;
    gsl_matrix * Height = gsl_matrix_calloc(P->GridXY,P->GridXY); // z coordinate of relaed point. global copy.
    gsl_matrix * PE = gsl_matrix_calloc(P->GridXY,P->GridXY); // potential energy grid. global copy.
    
    for (k=0; k<P->PESIterations; k++)
	{
	    if (P->GetPES2 == 1) FixedZ(temp,L,P,PE); // quench with adsorbate's coordinates x,y,z fixed
	    else if (P->GetPES2 > 1) FreeZ(k,temp,L,P,PE, Height); // quench with adsorbate's coordinates x,y fixed
	}
    // divide by number of iterations
    gsl_matrix_scale(PE,1.0/P->PESIterations);
    gsl_matrix_scale(Height,1.0/P->PESIterations);

    // save PES and height in the folder where program was ran from as "PotGrid.##.txt"
    FILE * pot;
    
    if (P->GetPES2 == 2)
	{  
	  pot = fopen("PotGrid.pe.txt","w");
	  printf("\nMinimizing potential energy of adsorbate\n");
	}
    else if (P->GetPES2 == 3)
	{  
	  pot = fopen("PotGrid.ke.txt","w");
	  printf("\nMinimizing kinetic energy of adsorbate\n");
	}
      
    if (pot!=NULL)
    { if (gsl_matrix_fprintf(pot,PE,"%.3f") == 0) printf("Saved PES file\n");}
    else printf("Error. Could not open file to save the PES!\n");
    fclose(pot);

    FILE * height;
    height = fopen("Height.txt","w");
    if (height!=NULL)
    { if (gsl_matrix_fprintf(height,Height,"%.3f") == 0) printf("Saved height file\n");}
    else printf("Error. Could not open file to save the height!\n");
    fclose(height);
        
    if (P->QuenchTrack==1)
      {
	 printf("Closing quenching trakcing file");
	 fclose(L->QuenchTracking);
      }	      	
	
    // free memory
    gsl_matrix_free(Height);
    gsl_matrix_free(PE);
}



void FixedZ(int temp, struct Lattice* restrict L, struct Parameters* P, gsl_matrix *PE )
{

  // Loop over all 3D grid points. Fix for each point the adsorbate's position and call FixedZ. On return assign the potential value at the grid point.
  
  gsl_vector *R = gsl_vector_calloc(3); // position vector
  gsl_vector *PEZ = gsl_vector_calloc(P->GridZ); // potential energy at an (x,y) point vector
  int i,j,k;

   for (i=0; i<P->GridXY; i++) // assign values to gris points
      {gsl_vector_set(R,0,i*L->latconstx/(double)P->GridXY - P->shiftx); // x coordinate. shift by -0.9025A to center above hollow site
      for (j=0; j<P->GridXY; j++)
	{gsl_vector_set(R,1,j*L->latconsty/(double)P->GridXY - P->shifty); // y coordinate
	  gsl_vector_set_zero(PEZ); // set elements to zero for a new x,y point

	  //printf("\t Point X:%g\tY:%g\n",0.7071*i*L->latconstx/(double)P->GridXY,0.7071*j*L->latconsty/(double)P->GridXY);

		#pragma omp parallel for private(k) shared(i,j,R,PE,PEZ,L,P)
		for (k=0; k<P->GridZ; k++) // z coordinate
		    {
		      gsl_vector_set(R,2,L->zmax+1+(double)k*2/(double)P->GridZ); // z coordinate range is 1-3 Ang above lattice
		      Get_PES_XYZ(k,temp, R, PEZ, L, P); // assigns the potential along z coordinate for a fixed (x,y)

		    }
		//if (omp_get_thread_num() == 0)
		#pragma omp barrier
		{
		printf("PE is %g\n",gsl_vector_min(PEZ)/9.654);
		gsl_matrix_set(PE,i,j,gsl_vector_min(PEZ)/9.654); // assign minimum of potential to PE (divide by 9.654 to get meV)
		//Assign_Min_Potnetial(i,j,P->Mins_To_Average, PE, PEZ); // assigns the minimum of the potential at the point (x,y) to PES
		}
	      }

      }
    gsl_vector_free(R);
    gsl_vector_free(PEZ);
}




void FreeZ(int k, int temp, struct Lattice* restrict L, struct Parameters* P, gsl_matrix *PE, gsl_matrix * Height )
{

  // Loop over all 2D grid points. Fix for each point the adsorbate's (x,y) coordinates and
  // set z to 2A above lattice and call Get_PES_XYZ. On return assign the potential value at the grid point.
  // the index k is used to set a uniqe seed.

    
  if (P->GetPES2 == 2)
      {  
	printf("\nMinimizing potential energy of system\n");
      }
  else if (P->GetPES2 == 3)
      {  
	printf("\nMinimizing kinetic energy of system\n");
      }  

  double PotE;

  gsl_matrix * LocalHeight = gsl_matrix_calloc(P->GridXY,P->GridXY); // z coordinate of relaed point. local copy.
  gsl_matrix * LocalPE = gsl_matrix_calloc(P->GridXY,P->GridXY); // potential energy grid. local copy.
  int m,n;
  size_t counter = 0;
  
  // open file to write into the celocity of atoms and potential and kinetc energy
  if (P->QuenchTrack==1)
  {
    char * filename = Gen_Filename("QuenchTracking", temp,-1,"txt");
    L->QuenchTracking = fopen(filename,"w");
    fprintf(L->QuenchTracking,"X\tY\tZ\tAPE\tXcm\tYcm\tZcm\tVxcm\tVycm\tVzcm\tPE\tKE\tTE\n\n\n"); 
    fclose(L->QuenchTracking);
    L->QuenchTracking = fopen(filename,"a");// open to data to file
  }    
   
   
   /*
   in order to calculate the PES above a unit cell of the top layer, we use the following relations between the axis of the lattice system
   (X,Y) and the axis of the top layer (x,y).
   If the grid above the lattice unit cell is MxM then the points in the lattice unit cell grid (i,j)
   which overlap with the top layer unit cell grid (m,n) have the following relation:
   
   [i,j] = [M/2-m+n,m+n]  0<m,n<M/2
   
   so the point [m,n] in the top layer unit cell grid is equal to the point [M/2-m+n,m+n] in the lattice cell grid.
   
   we designate by P->GridXY the number of grid points per dimension in the top layer unit cell.
   */

   double x,y;
 // this will shift the grid so its center is over a hollow site
   for (m=0; m<P->GridXY;m++) // assign values to grid points
      {
	
      #pragma omp parallel for private(n,x,y) shared (k,P,LocalHeight,LocalPE,m,L)

      for (n=0; n<P->GridXY; n++)
        {
	  gsl_vector *R = gsl_vector_calloc(3); // position vector
	  gsl_vector_set(R,2,L->zmax+P->zinit); // adsorbate starts 2 Angs over lattice	  
	  x = 0.5*L->latconstx*(1 + (double)(n-m)/(double)P->GridXY);
	  y = 0.5*(double)(n+m)*L->latconstx/(double)P->GridXY;
          gsl_vector_set(R,0,x+P->shiftx); // x coordinate. 
	  gsl_vector_set(R,1,y+P->shifty); // y coordinate
	  PotE = Get_PES_XY(k,temp,R, L, P); // the index k is used to generate a distinct seed
	  gsl_matrix_set(LocalPE, m, n, PotE/9.654); // assign minimum of potential to PE (divide by 9.654 to get meV)
	  gsl_matrix_set(LocalHeight, m, n , gsl_vector_get(R,2)); // get final height of adsorbate
	  gsl_vector_set(R,2,L->zmax+2); // reset initial z coordinate
	  counter ++;
	  if (omp_get_thread_num() == 0)
	  {
	  printf("\rrun %d %c %zu out of %zu x=%g  y=%g     GetPES = %d",k,"|/-\\"[n%4],counter,(size_t)(P->GridXY)*(P->GridXY), x, y, P->GetPES2);
	  fflush(stdout);
	  }
	   gsl_vector_free(R); 
	 }

      }

    // write local copy to file
    
    FILE * localpot;
    char * filename = malloc(100);
    int len = sprintf(filename, "PotGrid-%d.txt",k);
    free(filename);
    filename = malloc(len);
    len = sprintf(filename, "PotGrid-%d.txt",k);
    localpot = fopen(filename,"w");
    
    if (localpot!=NULL)
    { if (gsl_matrix_fprintf(localpot,LocalPE,"%.3f") == 0) printf("\nSaved PES file\n");}
    else printf("\nError. Could not open file to save the PES!\n");
    fclose(localpot);
    free(filename);

    FILE * localheight;
    filename = malloc(100);
    len = sprintf(filename, "Height-%d.txt",k);
    free(filename);
    filename = malloc(len);
    len = sprintf(filename, "Height-%d.txt",k);
    localheight = fopen(filename,"w");

    if (localheight!=NULL)
    { if (gsl_matrix_fprintf(localheight,LocalHeight,"%.3f") == 0) printf("\nSaved height file\n");}
    else printf("\nError. Could not open file to save the height!\n");
    fclose(localheight);    
    free(filename);

    // add to toal       
    gsl_matrix_add(PE,LocalPE); // add to total
    gsl_matrix_add(Height,LocalHeight); // add to total

    // free memory
    
    gsl_matrix_free(LocalPE);
    gsl_matrix_free(LocalHeight);
}


void Get_PES_XYZ( int k, int temp, gsl_vector *R, gsl_vector *PEZ, struct Lattice* restrict L, struct Parameters* P)
{
	// k is the index of the z coordinate
	// R is the fixed position of the adsorbate
	// PE is the potential vector along z

	//gsl_rng * Gauss = GaussianDist();
	//Reset_Substrate(S,L);
	struct Substrate * S;

	S = (struct Substrate*)malloc(sizeof(struct Substrate));
	S->temperature = P->Temperatures[temp];
	S->seed = P->Seed + P->Temperatures[temp];
	Init_Substrate(S, L, P);
	Garray_To_Array(L,S);
	// Generate Adsorbates
	struct Adsorbate * Adatoms[1];
	Adatoms[0] = (struct Adsorbate*)malloc(sizeof(struct Adsorbate));
	Gen_Adsorbate(Adatoms[0],L,P);
	Init_Adsorbate(Adatoms[0], L, P->Seed);
	gsl_vector_memcpy(Adatoms[0]->Rt,R); // assign position to adsorbate
	gsl_vector_memcpy(Adatoms[0]->R1,R);
	// Set initial temperature to desired value
	Initial_Temperature_Rescale(S,L);


	/* Dynamics */

	// run till kinetic energy is low enough

	while (S->KE > P->Elimit)
	  {
	    Substrate_Timestep(Adatoms, S, L, P);
	    Get_Kinetic_Energy_Substrate(S, L); // update knietic energy value
	    //printf("Kinetic energy is %g [E]\n", S->KE);
	  }
	// finished quenching. Get adsorbate's potential energy.
	gsl_vector_set(PEZ,k,Adatoms[0]->PE);
	//printf("PEZ is %g\n",Adatoms[0]->PE);
	// free memory
	Free_Substrate(S);
	free(S);
	Free_Adsorbate(Adatoms[0]); // free members
	free(Adatoms[0]);

}





double Get_PES_XY(int k, int temp, gsl_vector *R, struct Lattice* restrict L, struct Parameters* P)
{
	// k is the index of the current run
	// R is the fixed position of the adsorbate
	// PE is the potential vector along z

	// Generate substrate and adsorbate. Do quneching. Keep adsorbate's potential at relaxation point and its height.
	  
	//gsl_rng * Gauss = GaussianDist();
	struct Substrate * S;
	S = (struct Substrate*)malloc(sizeof(struct Substrate));
	S->temperature = P->Temperatures[temp];
	S->seed = P->Seed + P->Temperatures[temp]+k;
	Init_Substrate(S, L, P);
	Garray_To_Array(L,S);
	double PotE;
	
	
	// Generate Adsorbates
	struct Adsorbate * Adatoms[1];
	Adatoms[0] = (struct Adsorbate*)malloc(sizeof(struct Adsorbate));
	Gen_Adsorbate(Adatoms[0],L,P);
	Init_Adsorbate(Adatoms[0], L, P->Seed);
	gsl_vector_memcpy(Adatoms[0]->Rt,R); // assign position to adsorbate
	gsl_vector_memcpy(Adatoms[0]->R1,R); // assign position to adsorbate
	gsl_vector_set(Adatoms[0]->Vt,2,5); // set initial z velocity
	gsl_vector_set(Adatoms[0]->V1,2,5); // set initial z velocity
	Adatoms[0]->KE = 1;
	
	// Set initial temperature to desired value
	Initial_Temperature_Rescale(S,L);
  
	double TotalPE,  TotalKE,  PreviousTotalPE = 0 ,  PreviousTotalKE = 0; 
	/* Dynamics */

	// Quenching. run till kinetic energy is minimal
	int nomin = 1;
	if (P->QuenchTrack==1) fprintf(L->QuenchTracking,"**********************************\n"); //separate grid points with this line

	double timestep=0;
	while ( nomin ) // run till we reach a minimum of PE or KE is below set value
	  {
	    
	    Substrate_Timestep(Adatoms, S, L, P); // dynamics
	    Get_Kinetic_Energy_Substrate(S, L); // update substrate kinetic energy
	    Get_Kinetic_Energy_Adsorbate(Adatoms[0]); // update adsorbate's kinetic energy
	    TotalPE = Adatoms[0]->PE + S->PE;
	    TotalKE = Adatoms[0]->KE + S->KE;	 
	    
	    
	    // Track quenching	
	    if (P->QuenchTrack==1)
		{
		  //for (int j=0;j<3;j++) fprintf(L->QuenchTracking,"%g\t",gsl_vector_get(Adatoms[0]->Vt,j)); // print adatom's velocty
// 		  CM_Velocity(S, L);
// 		  for (int j=0;j<3;j++) fprintf(L->QuenchTracking,"%g\t",gsl_vector_get(S->CMVelocity,j)); // print substrate CM velocty
// 		  CM_Position(S,L);
// 		  for (int j=0;j<3;j++) fprintf(L->QuenchTracking,"%g\t",gsl_vector_get(S->CMPosition,j)); // print substrate CM position
		  fprintf(L->QuenchTracking,"%.10f\t\t",timestep*P->Timestep); // print time to file
		  fprintf(L->QuenchTracking,"%.10f\t\t",TotalKE/L->Kb/9654); // print temperature to file
		  fprintf(L->QuenchTracking,"%.10f\t\t",TotalPE/9654); // print total potential energy to file
		  //-fprintf(L->QuenchTracking,"%.g\t",S->PE+Adatoms[0]->PE); // print total potential energy to file
		  fprintf(L->QuenchTracking,"%.10f\t\t",TotalKE/9654+TotalPE/9654); // print total energy energy to file
		  //fprintf(L->QuenchTracking,"%.g\n",S->KE+S->PE+Adatoms[0]->PE); // print substrate kinetic energy to file
		  fprintf(L->QuenchTracking,"\n");
		}	      
		
	    // make sure the adsorbate has not moved in x or y
	    assert(gsl_vector_get(Adatoms[0]->Vt,0)==0);
	    assert(gsl_vector_get(Adatoms[0]->Vt,1)==0);
	    assert(gsl_vector_get(Adatoms[0]->V1,0)==0);
	    assert(gsl_vector_get(Adatoms[0]->V1,1)==0);
	    
	    // check if we reached a minimum using KE or PE
	    if (P->GetPES2 == 2) nomin = QuenchPE(P, TotalPE, &PreviousTotalPE);
	    else if (P->GetPES2 == 3) nomin = QuenchKE(P, TotalKE, &PreviousTotalKE);
	    timestep+=1;
	  }
	  
	
	// finished quenching. Get adsorbate's potential energy.
	PotE = Adatoms[0]->PE;
	if (P->QuenchTrack==1) fprintf(L->QuenchTracking,"\n*** Minimal PE is %.5f ***\n", PotE );
	gsl_vector_set(R,2,gsl_vector_get(Adatoms[0]->Rt,2)); // get final z coordinate
	// free memory
	Free_Substrate(S); // free members
	free(S);
	Free_Adsorbate(Adatoms[0]); // free members
	free(Adatoms[0]);

	return PotE;
}




int QuenchPE(struct Parameters* P, double TotalPE, double * PreviousTotalPE)
{
  
      // compare the difference between the current and previous potential energy of the entire system.
      // If the difference fell below a preset threshold then the system is declared relaxed.

      double diff = fabs(TotalPE - *PreviousTotalPE);
      *PreviousTotalPE = TotalPE; // update previous PE to current one
      if (diff < P->Elimit) return 0; // reached relaxation      
      else return 1;	  
}


int QuenchKE(struct Parameters* P, double TotalKE, double * PreviousTotalKE)
{
      // compare the difference between the current and previous kinetic energy of the entire system.
      // If the difference fell below a preset threshold then the system is declared relaxed.

      double diff = fabs(TotalKE - *PreviousTotalKE);
      *PreviousTotalKE = TotalKE; // update previous KE to current one
      if (diff < P->Elimit) return 0; // reached relaxation      
      else return 1;	  
}




void Assign_Min_Potnetial(size_t j, size_t k, size_t l, gsl_matrix *PE, gsl_vector *PEZ)
{
  // sorts PEZ and finds the minimum potential at (x,y)
  // l is the number of minimal values in the sorted vector we average on to get minimum
  // j,k are the indices of the potential matrix PE.

  assert( (l>0) & (l<PEZ->size)); // make sure l is not zero or greater than the number of grid points in the z direction
  double p,P[l]; // array to hold l minimal values
  gsl_sort_vector_smallest (P, l, PEZ); // sort minimal values and place in P
  size_t i;
  for (i=0,p=0; i<l; i++)
  {
    p += P[i]; // sum minimal values
  }
  printf("PE is %g\n",p/l/9.654);
  gsl_matrix_set(PE,j,k,p/l/9.654); // assign average minimal value to PE
}
