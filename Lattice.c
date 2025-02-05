#include <stdio.h>
#include <stdlib.h>
#include <glib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sort_vector.h>
#include <assert.h>
#include "Lattice.h"



/****************************************************


void Pes2FCC(struct Lattice *Lattice)
{
  
  // version of FCC lattice whose axis are parallel to the the x,y,z axes of the lab so we can extract the potential using method 2
  
  int register l,m,n;
  size_t a = 0; //counter of substrate atoms position vectors
  
  // calculate number of atoms (second term on right hand side is the number of atoms in the bottom layer
  Lattice->number_of_atoms = 4*Lattice->numcellx*Lattice->numcelly*Lattice->numcellz + 2*Lattice->numcellx*Lattice->numcelly;

  printf("N = %zu\n", Lattice->number_of_atoms);

  // alloctae matrix to hold position of lattice atoms in a prefect lattice

  assert(Lattice->number_of_atoms>0);
  Lattice->LatticeRp = gsl_matrix_calloc ( (size_t)Lattice->number_of_atoms, 3);
  
  
  gsl_vector * FCCvec[4]; // pointer to array of the lattice base vectors
  for (int i=0; i<4; i++) FCCvec[i] = gsl_vector_calloc(3); // allocate memory for lattice base vectors
  for (size_t i=0; i<3; i++) // [0.5,0.5,0]; [0,0.5,0.5]; [0.5,0,0.5]
    {
	gsl_vector_set(FCCvec[i+1],i%3,0.5); 
	gsl_vector_set(FCCvec[i+1],(i+1)%3,0.5); 
    }	  
  gsl_vector *ones = gsl_vector_calloc(3);
  gsl_vector_set_all(ones,1); // (1,1,1)


  // scale elements by lattice constants
  gsl_vector * latconst = gsl_vector_calloc(3);
  gsl_vector_set(latconst, 0, Lattice->latconstx);
  gsl_vector_set(latconst, 1, Lattice->latconsty);
  gsl_vector_set(latconst, 2, Lattice->latconstz);
  for (int i=1; i<4; i++) gsl_vector_mul(FCCvec[i],latconst); // scale the element of each vector by the lattice constant
  gsl_vector_mul(ones,latconst);
  
  // 45 degree rotation matrix
  
  gsl_matrix *Rot45 = gsl_matrix_calloc(3,3); 
  gsl_matrix_set(Rot45,0,0,1/sqrt(2));
  gsl_matrix_set(Rot45,0,1,1/sqrt(2));
  gsl_matrix_set(Rot45,1,0,-1/sqrt(2));
  gsl_matrix_set(Rot45,1,1,1/sqrt(2));
  gsl_matrix_set(Rot45,2,2,1);
  

  
  gsl_vector * vec; // template for the position vector
  vec = gsl_vector_calloc(3); // allocate memory for lattice base vectors
  gsl_vector * temp = gsl_vector_calloc(3); //temp vector
  printf("Generating lattice atoms positions\n");
  for (n=Lattice->numcellz; n>0; n--) 
    {
	  for(l=Lattice->numcellx; l>0; l--) 
	    {
		  for(m=Lattice->numcelly; m>0; m--)
		    {  
 
			// add base vectors and mulitiply by rotation matrix
			for (int i=0; i<4; i++) 
			  {	
				gsl_vector_set(vec,0,(double)l);
				gsl_vector_set(vec,1,(double)m);
				gsl_vector_set(vec,2,(double)n);					

				gsl_vector_mul(vec,ones); // this give the vector (l,m,n)
				gsl_vector_add(vec,FCCvec[i]); // this gives one of the 4 base vector e.g. (l,m+0.5,n)
				gsl_vector_set_zero(temp); // set temp vector to zero
				gsl_blas_dgemv(CblasNoTrans, 1, Rot45, vec, 0, temp); // put result of rotation into temp
				gsl_vector_fprintf(stdout,vec,"%g");
				// set potision of substrate atom number a
				gsl_matrix_set_row (Lattice->LatticeRp, a, temp);
				a++;
				//printf("a=%zu, l=%d, m=%d, n=%d\n",a,l,m,n);  
			  }      
		    } 
	    }
  }	
  
  n=0.5;  // generate the bottom layer
  for(l=Lattice->numcellx; l>0; l--) 
    {
	for(m=Lattice->numcelly; m>0; m--) 
	  {

	    
		for (int i=2; i<4; i++) // only use [0,0.5,0.5]; [0.5,0,0.5]
		  {	
			gsl_vector_set(vec,0,(double)l);
			gsl_vector_set(vec,1,(double)m);
			gsl_vector_set(vec,2,(double)n);			    
			    
			gsl_vector_mul(vec,ones); 
			gsl_vector_add(vec,FCCvec[i]); 
			gsl_vector_set_zero(temp); // set temp vector to zero
			gsl_blas_dgemv(CblasNoTrans, 1, Rot45, vec, 0, temp); // put result of rotation into temp
			gsl_vector_fprintf(stdout,vec,"%g");
			// set potision of substrate atom number a
			gsl_matrix_set_row (Lattice->LatticeRp, a, temp);
			//printf("a=%zu, l=%d, m=%d, n=%d\n",a,l,m,n);  
			a++;
		  }      
	    }
     }  
  
   
  printf("a=%zu\n",a);  
  assert (a == Lattice->number_of_atoms);
  printf("Finished generating lattice\n");
  FILE * Lat;
  Lat = fopen("lattice", "w");
  gsl_matrix_fprintf(Lat,Lattice->LatticeRp,"%g");
}

*******************************************************************/



void GenFCC (struct Lattice *Lattice)
{

        // Generates an FCC (001) lattice

        size_t register l,m,n;
        size_t a=0; // a runs over all N atoms
        // calculate number of atoms (second term on right hand side is the number of atoms in the bottom layer
        Lattice->number_of_atoms = 4*Lattice->numcellx*Lattice->numcelly*Lattice->numcellz + 2*Lattice->numcellx*Lattice->numcelly;

        printf("N = %zu\n", Lattice->number_of_atoms);

        // alloctae matrix to hold position of lattice atoms in a prefect lattice

	assert(Lattice->number_of_atoms>0);
        Lattice->LatticeRp = gsl_matrix_calloc ( (size_t)Lattice->number_of_atoms, 3);

        // Assign initial position to each atom and append to list
        printf("Generating lattice atoms positions\n");
	
        for (n=Lattice->numcellz; n>0; n--) {
                for(l=Lattice->numcellx; l>0; l--) {
                        for(m=Lattice->numcelly; m>0; m--) {

                                gsl_matrix_set(Lattice->LatticeRp, a,0, (Lattice->latconstx)*l);
                                gsl_matrix_set(Lattice->LatticeRp, a,1, (Lattice->latconsty)*m);
                                gsl_matrix_set(Lattice->LatticeRp, a,2, (Lattice->latconstz)*n);
                                a++;

                                // second lattice base vector
                                gsl_matrix_set(Lattice->LatticeRp, a,0, Lattice->latconstx*(l+0.5));
                                gsl_matrix_set(Lattice->LatticeRp, a,1, Lattice->latconsty*(m+0.5));
                                gsl_matrix_set(Lattice->LatticeRp, a,2, (Lattice->latconstz)*n);
                                a++;

                                // third lattice base vector
                                gsl_matrix_set(Lattice->LatticeRp, a,0, Lattice->latconstx*(l+0.5));
                                gsl_matrix_set(Lattice->LatticeRp, a,1, (Lattice->latconsty)*m);
                                gsl_matrix_set(Lattice->LatticeRp, a,2, Lattice->latconstz*(n+0.5));
                                a++	;


                                // forth lattice base vector
                                gsl_matrix_set(Lattice->LatticeRp, a,0, (Lattice->latconstx)*l);
                                gsl_matrix_set(Lattice->LatticeRp, a,1, Lattice->latconsty*(m+0.5));
                                gsl_matrix_set(Lattice->LatticeRp, a,2, Lattice->latconstz*(n+0.5));
                                a++;

                        }
                }
        }
        // Add bottom layer which will be frozen
	for(l=Lattice->numcellx; l>0; l--) {
	    for(m=Lattice->numcelly; m>0; m--) 
	    {
		    gsl_matrix_set(Lattice->LatticeRp, a,0, (Lattice->latconstx)*(l+0.5));
		    gsl_matrix_set(Lattice->LatticeRp, a,1, (Lattice->latconsty)*m);
		    gsl_matrix_set(Lattice->LatticeRp, a,2, (Lattice->latconstz)*(0.5));
		    a++;

		    // second lattice base vector
		    gsl_matrix_set(Lattice->LatticeRp, a,0, Lattice->latconstx*l);
		    gsl_matrix_set(Lattice->LatticeRp, a,1, Lattice->latconsty*(m+0.5));
		    gsl_matrix_set(Lattice->LatticeRp, a,2, (Lattice->latconstz)*(0.5));
		    a++;
	    }
        }
        // make sure number of generated atoms is equal to lenght of array
        assert (a == Lattice->number_of_atoms);
        printf("Finished generating lattice\n");

}


void Init_Lattice(struct Lattice* L, struct Parameters* P)
{
	// initiate lattice
	GenFCC(L) ;       
        Move_To_Center_Of_Mass(L);
        Find_Lattice_Extent(L);
        Find_Top_Bottom_Bulk(L);
        Distance_Bewteen_Atoms(L);
        Construct_Force_Matrix_Array(L);
        Pair_Equilibrium_Distance(L);
        Test_NN_Distance(L);
	// set some parameters
	L->Get_PE = 0;
	L->dt_over_m = P->Timestep/(6*L->mass);
        L->dt2_over_m = P->Timestep*P->Timestep/(6*L->mass);
	L->Vel_Scale_factor = P->Temperature_Rescale_Rate;
	// save the lattice configuration to file
 	FILE * Lat;
 	Lat = fopen("lattice", "w");
 	gsl_matrix_fprintf(Lat,L->LatticeRp,"%g");  
	
}



void Move_To_Center_Of_Mass(struct Lattice* L)
{

        // Find center of mass and shift it to x=0,y=0,z=0

        size_t k,l;
        double Xsum = 0;
        double Ysum = 0;
        double Zsum = 0;

        printf("Finding center of mass\n");
        for (k=0; k<L->number_of_atoms; k++) {
                Xsum += gsl_matrix_get(L->LatticeRp,k,0);
                Ysum += gsl_matrix_get(L->LatticeRp,k,1);
                Zsum += gsl_matrix_get(L->LatticeRp,k,2);

        }

        gsl_vector *CM = gsl_vector_alloc(3);
        gsl_vector_set(CM,0,-Xsum/(L->number_of_atoms));
        gsl_vector_set(CM,1,-Ysum/(L->number_of_atoms));
        gsl_vector_set(CM,2,-Zsum/(L->number_of_atoms));

        printf("Xcm = %f\tYcm = %f\tZcm = %f\n", Xsum,Ysum,Zsum);

        // allocate temporary vector
        gsl_vector *temp = gsl_vector_alloc(L->number_of_atoms);

        for (l=0; l<3; l++) {

                gsl_matrix_get_col(temp,L->LatticeRp,(size_t)l);
                gsl_vector_add_constant (temp, gsl_vector_get(CM,l));
                gsl_matrix_set_col(L->LatticeRp,l, temp);

        }

        printf("Moved to center of mass system\n");
        Xsum = Ysum = Zsum = 0;
        for (k=0; k<L->number_of_atoms; k++) {
                Xsum += gsl_matrix_get(L->LatticeRp,k,0);
                Ysum += gsl_matrix_get(L->LatticeRp,k,1);
                Zsum += gsl_matrix_get(L->LatticeRp,k,2);
        }
        printf("Xcm = %f\tYcm = %f\tZcm = %f\n", Xsum,Ysum,Zsum);
        // free temporary vectors
        gsl_vector_free(CM);
        gsl_vector_free(temp);
}


void Find_Lattice_Extent(struct Lattice* L)
{

        // Find minimum and maximum coordinate and lenghts of lattice

        // Get column vectors
        gsl_vector *V = gsl_vector_alloc_col_from_matrix(L->LatticeRp, 0);
        // get the max and min values
        gsl_vector_minmax (V, &(L->xmin), &(L->xmax));
        // set the boundaries of the lattice
        L->xboundary = L->xmax - L->xmin + 0.5*L->latconstx;
        gsl_vector_free(V);

        printf("Lattice x edges are at %f and %f\n", L->xmin, L->xmax);


        gsl_vector *T = gsl_vector_alloc_col_from_matrix(L->LatticeRp, 1);
        gsl_vector_minmax (T, &(L->ymin), &(L->ymax));
        L->yboundary = L->ymax - L->ymin + 0.5*L->latconsty;
        gsl_vector_free(T);

        printf("Lattice y edges are at %f and %f\n", L->ymin, L->ymax);


        gsl_vector *W = gsl_vector_alloc_col_from_matrix(L->LatticeRp, 2);
        gsl_vector_minmax (W, &(L->zmin), &(L->zmax));
        L->zboundary = L->zmax - L->zmin;
        gsl_vector_free(W);

        printf("Lattice z edges are at %f and %f\n", L->zmin, L->zmax);

}




void Find_Top_Bottom_Bulk(struct Lattice* L)
{


        // Find indices of atoms in the bottom and top layers, and bulk atoms(all but bottom atoms).

        size_t k, top,bottom, type;
        // initialize array
        L->Top = g_array_new(TRUE, TRUE, sizeof( size_t));
        L->Bottom = g_array_new(TRUE, TRUE, sizeof( size_t));
        L->Bulk = g_array_new(TRUE, TRUE, sizeof( size_t));
        L->Type = g_array_new(TRUE, TRUE, sizeof( size_t));
        L->number_of_top_atoms = 0;
        L->number_of_bottom_atoms = 0;
        L->number_of_bulk_atoms = 0;
        // loop over all atoms
        for (k=0; k<L->number_of_atoms; k++) {

                //
                bottom = 0;
                top = 0;
                if (gsl_matrix_get(L->LatticeRp,k,2) > 0.99*L->zmax) { // test for top atom
                        top = 1;
                        //printf("Atom number %d is a top atom\n", k);
                        // append to list of top atoms
                        g_array_append_val(L->Top,k);
                        // increase number of top atoms is it's a top atom
                        L->number_of_top_atoms += 1;
                        // write type of atom k
                        type = 8; // number of neighbours
                        g_array_append_val(L->Type,type);
                }

                if (gsl_matrix_get(L->LatticeRp,k,2) < 0.99*L->zmin) { // test for bottom atom
                        bottom = 1;
                        //printf("Atom number %d is a bottom atom\n", k);
                        g_array_append_val(L->Bottom,k);
                        L->number_of_bottom_atoms += 1;
                        type = 0; // bottom atoms are not affetced by neighbours
                        g_array_append_val(L->Type,type);
                }

                if (bottom == 0) {
                        g_array_append_val(L->Bulk,k);
                        L->number_of_bulk_atoms += 1;
                        if (top == 0) {
                                //printf("Atom number %d is a bulk atom\n", k);
                                type = 12;
                                g_array_append_val(L->Type,type);
                        }
                }

        }
        printf("There are %zu top atoms\n", L->number_of_top_atoms);
        printf("There are %zu bulk atoms\n", L->number_of_bulk_atoms);
        printf("There are %zu bottom atoms\n", L->number_of_bottom_atoms);


}



void Boundary_shift(const struct Lattice* L, gsl_vector * temp)
{
	// shift according to periodic boundary conditions
        // shift x distance
        double xdist = gsl_vector_get(temp,0);
        if (xdist > L->xboundary/2) {
                xdist = xdist - L->xboundary;
        } else if (xdist < -(L->xboundary/2)) {
                xdist = xdist + L->xboundary;
        }
        gsl_vector_set(temp, 0, xdist);

        // shift y distance
        double ydist = gsl_vector_get(temp,1);
        if (ydist > L->yboundary/2) {
                ydist = ydist - L->yboundary;
        } else if (ydist < -(L->yboundary/2)) {
                ydist = ydist + L->yboundary;
        }
        gsl_vector_set(temp, 1, ydist);

}


double Pair_Distance(const struct Lattice* L, const  size_t a, const  size_t b, size_t boundary_condition)
{

        // Finds the distance between atoms a and b

        double dR;
        gsl_matrix_get_row (L->Vtemp, L->LatticeRp, a);
        gsl_matrix_get_row (L->Vtemp1, L->LatticeRp, b);
        // subtract position of atom b from position of atom a
        gsl_vector_sub (L->Vtemp, L->Vtemp1 );
        // Shift according to periodic boundary conditions if boundary_condition == 1
        if (boundary_condition == 1) {
                Boundary_shift(L, L->Vtemp);
        }
        dR = gsl_blas_dnrm2 (L->Vtemp);
        //printf("Distance between atom %d and atom %d is %f\n", a, b, dR);
        return dR;


}




void Test_NN_Distance(struct Lattice* L)
{

        // prints the distance to the nearest neighbours

        size_t k,n,m,p;
        printf("Total number of %zu pairs\n", L->number_of_pairs);
        // loop over all pairs
        for (k=0; k<L->number_of_pairs; k++) {
                m = g_array_index(L->GForceatoms, size_t, k); // index of k'th bulk atom
                p = g_array_index(L->GNN,size_t,k); // index of neighbour
                n = g_array_index(L->Type, size_t, m); // number of neighbours for k'th bulk atom
                assert(n!=0); // can not be a bottom layer atom
                //printf("Distance between atom %d and atom %d of pair %d is %f\n", m,p,k, gsl_matrix_get(L->Mtemp,m,p));
        }
        // free memory
        gsl_vector_free(L->Dtemp);
}



void Sort_Distance(struct Lattice* L, gsl_vector * distances, size_t k)
{
    // For each atom (with index number k) sorts "distances", the distances vector 
    // to other atoms (including self)
    size_t n,m;
    n = 1+g_array_index(L->Type, size_t, k) ; // add one to count the atom itself
    if (n!=1) // Not a bottom atom. Participate in dynamics part.
	{
	    // alloctae array to hold indices of nearest atoms
	    size_t *p = (size_t*)calloc( n ,sizeof(size_t) );

	    gsl_sort_vector_smallest_index (p, n, distances); // fill p with indices of nearest atoms

	    // now just collect the first 13 (or 9) indices and drop the first
	    assert(p[0]==k); // make sure the first index is indeed the atom itself

	    // append the neighbours indices to the NN array, IF their index is greater than k
	    for (m=0; m<n; m++) 
	    {
		if (p[m]>k) 
		    {
			    L->GNN = g_array_append_val (L->GNN, p[m]); // append neighbour index to nearest
			    // neighbours array
			    //printf("p[%d] is %d and its distance is %g\n", m,p[m],gsl_matrix_get(L->Mtemp,k,p[m]));
			    L->GForceatoms = g_array_append_val (L->GForceatoms, k);// append atom's index to Forceatoms array
			    L->number_of_pairs +=1 ; // increase number of pairs by one
		    }
	    }
	}    
}




void Distance_Bewteen_Atoms (struct Lattice* L)
{

        // Find distance magnitude between lattice atoms
        size_t a,b; // hold idices of atoms
        // allocate temporary matrix and vectors to hold the values
        L->Dtemp = gsl_vector_calloc(L->number_of_atoms); // holds distances between one atom and all others
        L->Vtemp = gsl_vector_calloc(3); // temporary vector
        L->Vtemp1 = gsl_vector_calloc(3); // temporary vector

        // initialize array to hold indices of nearest neighbours for a specific atom
        L->GNN = g_array_new(FALSE, TRUE, sizeof( size_t));
        // initialize array to hold indices of atoms with neighbours
        L->GForceatoms = g_array_new(FALSE, TRUE, sizeof( size_t));
        // loop over all atom pairs
        for (a=0; a<L->number_of_atoms; a++)
	  {
		// do for atom a
                for (b=0; b<L->number_of_atoms; b++) 
		    {
                        // get distance between atom a and some atom b
                        double dist = Pair_Distance(L, a, b, 1);
                        // set vector value for distance between atoms a and b
                        gsl_vector_set(L->Dtemp,b,dist);
		    }	
		assert(gsl_vector_get(L->Dtemp,a)==0); // distance between atom a and itself is 0
		Sort_Distance(L, L->Dtemp, a); // Find the nearest neighbours
		gsl_vector_set_all(L->Dtemp,0); // reset the vector of distances for the next atom
                
	  }
        // free memory
        gsl_vector_free(L->Vtemp);
	gsl_vector_free(L->Vtemp1);
}



void Construct_Force_Matrix_Array(struct Lattice* L)
{
        // here we construct a matrix whose rows are the distance vectors between nearest negihbours
        // (with periodic boundary conditions). The atomic force matrix is proportional to the outer products
        // of such vectors


        L->NN_Distance_Vectors = gsl_matrix_calloc(L->number_of_pairs,3);
        size_t k,m,p;
        double d,b;
        L->Vtemp = gsl_vector_calloc(3);
        L->Vtemp1 = gsl_vector_calloc(3);
        // loop over all pairs
        for (k=0; k<L->number_of_pairs; k++) {
                m = g_array_index(L->GForceatoms, size_t, k); // index of k'th bulk atom
                p = g_array_index(L->GNN,size_t,k); // index of neighbour
                // get current atom's position
                gsl_matrix_get_row (L->Vtemp, L->LatticeRp, m);
                // get current atom's neighbour's position
                gsl_matrix_get_row (L->Vtemp1, L->LatticeRp, p);
                // subtract position of neighbour atom from position of current atom
                gsl_vector_sub (L->Vtemp, L->Vtemp1 );
                // Adjust distances according to periodic boundary conditions
                Boundary_shift(L, L->Vtemp);
                // scale the distance vector so all elemnets are 1,-1 or 0
                if ((d=gsl_vector_max(L->Vtemp))==0) d=gsl_vector_min(L->Vtemp);
                gsl_vector_scale(L->Vtemp,fabs(1/d));
                // some tests
                assert(gsl_vector_max(L->Vtemp)<1.1);
                assert(gsl_vector_min(L->Vtemp)>-1.1);
                gsl_blas_ddot (L->Vtemp, L->Vtemp, &b);
                assert(b>1);
                // append to matrix
                gsl_matrix_set_row(L->NN_Distance_Vectors,k,L->Vtemp);

        }
        // free memory
        gsl_vector_free(L->Vtemp);
        gsl_vector_free(L->Vtemp1);
	
	/***************************************************************************************
         multiply the vectors by sqrt of the force constant since its their outer product that
         generates the atomic force matrix
        ***************************************************************************************/
        gsl_matrix_scale(L->NN_Distance_Vectors,sqrt(L->Force_constant));
}



void Pair_Equilibrium_Distance(struct Lattice *L)
{

        // Distance vectors between nearest neighbours pair without using periodic boundary conditions
        L->LatticeDp = gsl_matrix_calloc(L->number_of_pairs,3);
        size_t k,m,p;
        //double d;
        L->Vtemp = gsl_vector_calloc(3);
        L->Vtemp1 = gsl_vector_calloc(3);
        // loop over all pairs
        for (k=0; k<L->number_of_pairs; k++) {
                m = g_array_index(L->GForceatoms, size_t, k); // index of first atom in pair k
                p = g_array_index(L->GNN,size_t,k); // index of second atom in pair k (neighbour)

                gsl_matrix_get_row (L->Vtemp, L->LatticeRp, m); // get new atom equilibrium position
                gsl_matrix_get_row (L->Vtemp1, L->LatticeRp, p); //get current atom's neighbour's equilibrium position
                gsl_vector_sub (L->Vtemp, L->Vtemp1 ); // subtract position of neighbour from position of atom 
                //Boundary_shift(L,L->Vtemp); // update according to periodic boundary conditions
                //assert(gsl_blas_dnrm2 (L->Vtemp)<2.6); // assert distance between neighbours is smaller than 2.6 (correct for Copper)
                gsl_matrix_set_row(L->LatticeDp,k,L->Vtemp);  // append to matrix

        }
        // displacement = (x-y)-(x0-y0)
	// disp is -(x0-y0), thats why L->LatticeDp is multiplyed by -1 
        gsl_matrix_scale(L->LatticeDp,-1); // multiply by -1
	gsl_vector_free(L->Vtemp);
	gsl_vector_free(L->Vtemp1);
}




void Print_Neighbours_List(struct Lattice *L)
{
        // print a list of atoms and their neighbours

        int *a = L->GForceatoms,*b=L->GNN;
        int k=0;
        for (; a<&(L->GForceatoms[L->number_of_pairs]); a++,b++) {
                k++;
                printf("Atom: %d	Neighbour: %d\n",*a,*b);
        }
        printf("Counted %d pairs\n",k);
}
