
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




void Get_parameters(char * datafile, struct  Parameters * P, struct Lattice * L)
{
    // top level routine which calls other routines to parse the parameters file and assign the values to the simulation's parameters
    printf("\n\n*****\n\n");
    Get_Data_File(datafile, P);
    Get_Time(P->Datafile, P);
    Get_Geometry(P->Datafile, L);
    Get_Physics(P->Datafile, L, P);
    Get_Simulation(P->Datafile, L, P);
    // create a folder to save data in only if needed
    Get_Files(P->Datafile, P);
    char * Date = Get_Date();
    printf("%s\n",Date);
    printf("\n\n*****\n\n");
}


void Copy_Param_File(char * datafile, struct  Parameters * P)
{
    // copy the parameters file 'datafile' to the data folder


    char * filename = malloc(100), ch;

    int len = sprintf(filename, "%s/Summary.txt",P->Savedir);
    free(filename);
    filename = malloc(len+1);    
    len = sprintf(filename, "%s/Summary.txt",P->Savedir);

    FILE * Source;
    Source = fopen(datafile,"r");
    
    printf("Saving parameters file <%s> to folder %s\n",datafile,filename);
    
    
    P->Summary = fopen(filename,"a");
    assert(P->Summary!=NULL);
    assert(Source!=NULL);
    
    char * date = Get_Date();

    // put starting time of simulation in copy of Parameters file
    fprintf(P->Summary,"#  %s\n", date);    
    while( ( ch = fgetc(Source) ) != EOF )
	fputc(ch, P->Summary);

    printf("File copied successfully.\n");
    free(filename);
    fclose(Source);
    fclose(P->Summary);
}

char * Get_Date()
{
    // get the date and time when the program began
    GDateTime * date = g_date_time_new_now_local();
    // get year,month,day
    gint year,month,day;
    g_date_time_get_ymd(date,&year,&month,&day);
    // get hour,minute
    gint hour,minute;
    hour =  g_date_time_get_hour (date);
    minute =  g_date_time_get_minute (date);
    GDate * Date = g_date_new_dmy(day,month,year);
    GDateMonth Month = g_date_get_month (Date);
    GDateYear Year = g_date_get_year (Date);
    char * buff = malloc(512);    
    if (minute<10) sprintf(buff,"Simulation began on %d/%d/%d at %d:0%d\n",day, month, year, hour, minute);
    else sprintf(buff,"Simulation began on %d/%d/%d at %d:%d\n",day, month, year, hour, minute);
    return buff;

}




void Get_Data_File(char * datafile, struct  Parameters * P)

{
      // Loads the parameters file
      P->Datafile = g_key_file_new ();
      g_key_file_load_from_file (P->Datafile, datafile, 0, NULL);

}



void Get_Time(GKeyFile * Key, struct  Parameters * P)
{

    // Get time related parameters (timestep, simulation duration, etc.) from the parameters file
    printf("\n\nGetting time parameters\n");
    // Set time related parameters
    char *Group = "Time";
    // Get keys of group
    char **Keys;
    size_t chainlen;
    Keys = g_key_file_get_keys (Key, Group, &chainlen, NULL);
    char *key;

    key = Keys[0];

    P->Timestep = g_key_file_get_double (Key, Group, key, NULL);
    size_t inverse = (size_t)(1/P->Timestep); // convert timsteps to picosec
    printf("Simulation Timestep is %g Picosec\n", P->Timestep);

    key = Keys[1];

    double temptime = (double) g_key_file_get_uint64 (Key, Group, key, NULL);
    P->Simulation_Time = (size_t)(temptime*inverse); 
    printf("Simulation duration is %zu timesteps or %g Picosec\n", P->Simulation_Time, ((double)P->Simulation_Time)*(P->Timestep));

    key = Keys[2];

    P->Substrate_Relax = (size_t)g_key_file_get_uint64 (Key, Group, key, NULL);

    key = Keys[3];

    P->Substrate_Samplerate = (size_t)g_key_file_get_uint64 (Key, Group, key, NULL);

    key = Keys[4];

    P->Adsorbate_Samplerate = (size_t)g_key_file_get_uint64 (Key, Group, key, NULL);

    key = Keys[5];

    P->Temperature_Rescale_Rate = (size_t)g_key_file_get_uint64 (Key, Group, key, NULL);

    key = Keys[6];

    P->Vel_Scale_factor = g_key_file_get_double (Key, Group, key, NULL);


}


void Get_Geometry(GKeyFile * Key,struct Lattice * L)
{
    // Get geometry related parameters (lattice constant, number of lattice cells)
    printf("\n\nGetting geometry parameters\n");
    // Get the shape of the lattice
    char *Group = "Geometry";
    // Get keys of group
    char **Keys;
    size_t chainlen;
    Keys = g_key_file_get_keys (Key, Group, &chainlen, NULL);
    char *key;


    key = Keys[0];
    int * Cells;
    size_t numcell;
    Cells = g_key_file_get_integer_list (Key, Group, key, &numcell, NULL);
    assert(numcell == 3);
    L->numcellz = Cells[2];
    assert(L->numcellz>0);
    L->numcelly = Cells[1];
    assert(L->numcelly>0);
    L->numcellx = Cells[0];
    assert(L->numcellx>0);
    printf("There are %d, %d and %d lattice unit cells in the X,Y and Z directions respectively\n",L->numcellx,L->numcelly,L->numcellz);

    key = Keys[1];
    double * Consts;
    size_t conlen;
    Consts = g_key_file_get_double_list (Key, Group, key, &conlen, NULL);
    L->latconstz = Consts[2]; // Angstrom
    assert(L->latconstz>0);
    L->latconsty = Consts[1];
    assert(L->latconsty>0);
    L->latconstx = Consts[0];
    assert(L->latconstx>0);
    printf("The lattice constants in the X,Y and Z directions are %gA, %gA and %gA respectively\n",L->latconstx,L->latconsty,L->latconstz);
}



void Get_Physics(GKeyFile * Key, struct Lattice * L, struct  Parameters * P)
{
    // Get physical parameters. some parameters are read into P instead of A which has not been initialized yet.
    printf("\n\nGetting physical parameters\n");
    // Get the shape of the lattice
    char *Group = "Physics";
    // Get keys of group
    char **Keys;
    size_t chainlen;
    Keys = g_key_file_get_keys (Key, Group, &chainlen, NULL);
    char *key;
    GError * error = NULL;

    key = Keys[0];
    L->Force_constant = g_key_file_get_double (Key, Group, key, NULL);
    L->Force_constant *= 602.41; // multipled by factor to convert from N/m to [F]/A

    key = Keys[1];
    size_t conlen;
    double * Morse;
    Morse = g_key_file_get_double_list (Key, Group, key, &conlen, &error);
    P->b = Morse[1];
    P->Amp = Morse[0]*9654; // multipled by a factor to convert from eV to [E]
    P->R0 = Morse[2];


    key = Keys[2];
    size_t templen;
    int * Temp = g_key_file_get_integer_list (Key, Group, key, &templen, NULL);
    //if (error != NULL) {fprintf (stderr, "ERROR: %s\n", error->message);}
    P->Temperatures = (int *)calloc(templen,sizeof(int));
    printf("Simulation will run for the following Temperatures:  ");
    for (size_t k=0; k<templen; k++) {P->Temperatures[k]=Temp[k]; printf("%dK  ",P->Temperatures[k]);}
    printf("\n");
    P->templen = templen;

    key = Keys[3];
    L->Kb = g_key_file_get_double (Key, Group, key, NULL); // Boltzmann const
    L->Kb *= 9654; // multipled by a factor to convert from eV/K to [E]/K
    printf("Boltzmann constant is %g [E]/T\n",L->Kb);
    key = Keys[4];
    L->mass = g_key_file_get_double (Key, Group, key, NULL);
    printf("Mass of substrate atom is %g AMU\n",L->mass);
    key = Keys[5];
    P->Amass = g_key_file_get_double (Key, Group, key, NULL);
    printf("Mass of adatom is %g AMU\n",P->Amass);
    key = Keys[6];
    P->shiftx = g_key_file_get_double (Key, Group, key, NULL);
    key = Keys[7];
    P->shifty = g_key_file_get_double (Key, Group, key, NULL);

    key = Keys[8];
    double * D;
    D = g_key_file_get_double_list (Key, Group, key, &conlen, &error);
    P->Dipole = D[0];
    P->AdsorbateDipole = D[1]*sqrt(602.41); // convert dipole parameter to simulation units,
                                            // so to give units of [F] in the total dipole interaction term
    P->DipoleRange = D[2];
    if (D[0]) printf("Using dipole-dipole interaction with cutoff range %gA\n", D[2]);

}


void Get_Simulation(GKeyFile * Key, struct Lattice * L, struct  Parameters * P)
{
    // Get some simulation related parameters
    printf("\n\nGetting simulation parameters\n");
    // Get the shape of the lattice
    char *Group = "Simulation";
    // Get keys of group
    char **Keys;
    size_t chainlen;
    Keys = g_key_file_get_keys (Key, Group, &chainlen, NULL);
    char *key;
    int k = 0;
    GError * error = NULL;


    key = Keys[k];
    P->cpus = g_key_file_get_integer (Key, Group, key, &error);
    // It is the user's responsibility to choose a suitable number of threads
    printf("Using %d threads in parallel\n", P->cpus);
    omp_set_num_threads(P->cpus);

    
    k++;
    key = Keys[k];
    P->Parts = g_key_file_get_integer (Key, Group, key, &error);
    if (error != NULL) {fprintf (stderr, "ERROR\n");}
    printf("Running %d times per configuration\n", P->Parts);


    k++;
    key = Keys[k];
    P->Seed = g_key_file_get_integer (Key, Group, key, NULL);
    


    k++;
    key = Keys[k];
    P->number_of_adatoms = g_key_file_get_integer (Key, Group, key, NULL);



    k++;
    key = Keys[k];
    P->Simtype = g_key_file_get_string (Key, Group, key, NULL);


    
    k++;
    key = Keys[k];
    size_t peslen;
    int * GetPES2;
    GetPES2 = g_key_file_get_integer_list (Key, Group, key, &peslen, NULL);
    // if we want to extract the PES using method 2, we run just for 300K.
    if ( (peslen>1) & (GetPES2[0]!=0) )
      {
	P->GetPES2 = GetPES2[0]; // type of quenching
	P->GridXY = GetPES2[1]; // grid resolution in x,y
	if (P->GetPES2 == 2 || P->GetPES2 == 3) P->zinit = GetPES2[2]; // initial z coordinate above substrate in angstrom when quenching with the vertical coordinate free
	else if (P->GetPES2 == 1) P->GridZ = GetPES2[2]; // grid resolution in z if quenching with z fixed also.
	P->Elimit = pow(10,GetPES2[3])*L->Kb; // energy in [E] the system has to reach during quenching
	printf("KElimit = %g\n",P->Elimit);
	P->Mins_To_Average = GetPES2[4]; // if P->GetPES2 = 1, average on this number of lowest potential values for a give x,y point
	P->PESIterations = GetPES2[4]; // if P->GetPES2 = 2, average the PES on this number of runs
	if (P->GetPES2 == 1) printf("Running to generate PES by fixing adsorbate over substrate\n");
	else if ( (P->GetPES2 > 1) & (P->GetPES2 < 4) ) printf("Running to generate PES by fixing only the adsorbate's x,y coordinates\n");
	else { printf ("Unknown option (GetPES2>3)!\n"); exit(1);}
	printf("Setting number of parts to 1!\n");
	printf("Grid is %d points in the x,y direction and %d points in the z direction.\n",P->GridXY,P->GridZ);
	P->Parts = 1;
	if (P->templen>1)
	  {
	    printf("You have chosen more than 1 temperature to extract the potential with.\n");
	    printf("Setting temperature to 300K.\n");
	    P->templen = 1;
	    P->Temperatures[0] = 300;
	}
	if (P->number_of_adatoms!=1)
	{
	   printf("Number of adsorbate is incorrect. Setting to 1.\n");
	   P->number_of_adatoms = 1;
	}
      }


    
    k++;
    key = Keys[k];
    if ( (P->Track_Adsorbate = g_key_file_get_integer (Key, Group, key, NULL))!=0) printf("Will keep adsorbates trajectories\n");

    // Check that if we track adsorbates, that there is at least one
    if ( (P->Track_Adsorbate>0) & (P->number_of_adatoms==0) )
    {
      printf("\n!!!\tYou chose to track adsorbate's trajectory, but adsorbates number is set to 0\t!!!\n");
      exit(1);
    }



    k++;
    key = Keys[k];
    if ( (P->Track_Adsorbate_Potential = g_key_file_get_integer (Key, Group, key, NULL))!=0)
    {
      assert(P->number_of_adatoms==1); // for now track potential due to sybstrate alone
      printf("Will keep adsorbates potential\n");
    }
    if ( ((P->Track_Adsorbate_Potential = g_key_file_get_integer (Key, Group, key, NULL))==0) & (P->GetPES2>0) )
    {
      printf("You asked to generate a PES byt forgot to set the switch \"Track_Adsorbate_Potential\". Seting it to 1.\n");
      P->Track_Adsorbate_Potential = 1;
    }



    k++;
    key = Keys[k];
    if ( (P->Track_Substrate = g_key_file_get_integer (Key, Group, key, NULL))!=0) printf("Will keep substrate atoms trajectories\n");



    k++;
    key = Keys[k];
    if ( (P->PDOS = g_key_file_get_integer (Key, Group, key, NULL))!=0) printf("Will keep substrate atoms velocities\n");
    
    
    k++;
    key = Keys[k];
    if ( (P->QuenchTrack = g_key_file_get_integer (Key, Group, key, NULL))!=0) printf("Will track quenching\n");    


    // check if we wish to track something
    for( int j = 4; j<chainlen; j++)
    {
       key = Keys[k-j]; // count from last key
       if (g_key_file_get_integer (Key, Group, key, NULL)==1) {P->Track=1; break;}
    }


}


void Get_Files(GKeyFile * Key, struct  Parameters * P)
{
    // Get parameters related to files and directories
    printf("Getting Files parameters\n");
    // Get the shape of the lattice
    char *Group = "Files";
    // Get keys of group
    char **Keys;
    size_t chainlen;
    Keys = g_key_file_get_keys (Key, Group, &chainlen, NULL);
    char *key;

    struct passwd *pw = getpwuid(getuid());

    char *homedir = pw->pw_dir;

    key = Keys[0];
    char * folder;
    folder =  g_key_file_get_string (Key, Group, key, NULL);
    P->Savedir = malloc(100);
    int len = sprintf(P->Savedir,"%s%s/Output",homedir,folder);
    free(P->Savedir);
    P->Savedir = malloc(len);
    sprintf(P->Savedir,"%s%s/Output",homedir,folder);

    if (P->Track>0 ) // if we track trajectories
    {
	// check dir does not exist already
	struct stat st;
	if(stat(P->Savedir,&st) == 0)
	{
	    printf("Directory %s is already present. Appending new files.\n", P->Savedir);
	}

	else
	{
	  // Create dir
	  if (g_mkdir_with_parents(P->Savedir,0755)==-1) {printf("Failed to generate the folder %s\n", P->Savedir); exit(1);}
	  else printf("Data will be save in the folder %s\n", P->Savedir);
	}
    }
}
