/*
Copyright (C) 2005 Jens Kunstmann

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

*/

#include	"lm_all.h"
#include	"data.h"
#include	"read_RHO.h"




/************************** read_RHO ****************************
Calls read_field() to read data from 'RHO' files.
Converts the data to Angstrom units.
Finally submits the data to the central data class.

It construct a local instance of the cental data class
writes DIRECTLY data INTO it and finally submits this data to the 
central data object. 

filename:	filename to read the data from
gd:		pointer of the central data class
*/

int read_RHO(char *filename, char *progname, struct data *gd)
{
   int		i,
   		err = 0;

   struct data	*ld;	/* local  */
   

   if(ld = data_init())  /* construct the LOCAL data class */
   {
      if( !(err = read_field(filename, &ld->rho)) )
      {
	 /* make identifier string*/
	 sprintf(ld->rho.identifier, "%s_%s", filename, progname);
	 ld->rho.identifier[ID_SIZE-1] = 0; /* if the buffer is too small, terminate by hand */
	 
	 /* conversion of dimensions vom r_B to Angstrom */
	 for(i=0; i < 3; i++) ld->rho.origin[i] = ld->rho.origin[i] * CONV;
	 for(i=0; i < ld->rho.npoints; i++) ld->rho.field[i] = ld->rho.field[i] * 1/(CONV*CONV*CONV);

	 debug( print_somedata(&ld->rho) );

	 /* save the data to the global data object */
	 data_put(gd, ld, PUT_RHO);         
      }
      /* error message is already printed in read_field() */         

      data_free(ld); /* is the the destructor of the 'data' class */
   }
   else err = print_err(err_cam, 0L);
   
   return(err);

}



/************************** read_ELF ****************************
Calls read_field() to read data from 'ELF' files.
NO conversion is needed.
Finally submits the data to the central data class.

It construct a local instance of the cental data class
writes DIRECTLY data INTO it and finally submits this data to the 
central data object.

filename:	filename to read the data from
gd:		pointer of the central data class
*/

int read_ELF(char *filename, char *progname, struct data *gd)
{
   int		i,
   		err = 0;
		
   struct data	*ld;	/* local  */   
      
   
   if(ld = data_init())  /* construct the LOCAL data class */
   {
      if( !(err = read_field(filename, &ld->elf)) )
      {
	 /* make identifier string*/
	 sprintf(ld->elf.identifier, "%s_%s", filename, progname);
	 ld->elf.identifier[ID_SIZE-1] = 0; /* if the buffer is too small, terminate by hand */
	 
	 /* NO conversion of dimensions for ELF*/
	 debug( print_somedata(&ld->elf) );

	 /* save the data to the global data object */
	 data_put(gd, ld, PUT_ELF);         
      }         

      data_free(ld); /* the the destructor of the 'data' class 
      			it will also free all sub-arrays */
   }
   else err = print_err(err_cam, 0L);
   
   return(err);
   
}





/********************* read_field *************************
Read field data from 'RHO' or 'ELF' files.

It sequentially reads the data from the rho_file,
reorders the data from (y,z,x) [lmto] to (x,y,z) [xcrysden] format
and stores the data in the data_field given via the pointer 'rho'.

!!!It allocates memory for the datafield!!!
Don't forgat to free it afterwards.

filename:	name of the file to open
rho:		the central data object
*/

int read_field(char *filename, struct data_field *rho)
{
   FILE		*fh;	
   	
   char		buf[90];   	/* line buffer */
   
   float	*tmp_field,	/* yet another field, needed for data conversion */
   		delta[3][3];	/* the three delta vectors are used to construct the 
				   'datacell', index ordering: [A,B,C][x,y,z]*/

   int		err = 0,	/* error value */
   		i, ix, iy, iz, 	/* loop counters */
		ind,		/* index for the map ind = f(ix,iy,iz) */
		nx, ny, nz,	/* grid values */
		cnt_RHO = 0;	/* count the number of lines that are read in successfully */

		
   if (fh = fopen(filename, "r"))
   {       
      /*** first read the grid, origin and the deltas ***/
      
      fgets(buf, sizeof(buf), fh); 
      if(sscanf(buf, "object 1 class gridpositions counts %d %d %d", 
         &rho->grid[0], &rho->grid[1], &rho->grid[2] ) == 3) cnt_RHO++;

      fgets(buf, sizeof(buf), fh); 
      if(sscanf(buf, "origin %f %f %f", 
         &rho->origin[0], &rho->origin[1], &rho->origin[2] ) == 3) cnt_RHO++;
      
      for(i=0; i < 3; i++)
      {
         fgets(buf, sizeof(buf), fh);
	 if(sscanf(buf, "delta %f %f %f", &delta[i][0], &delta[i][1], &delta[i][2] ) == 3) cnt_RHO++;
      }	 
	 
      
      /* check for correctness */
      if(cnt_RHO == 5)    	 
      {
	 /* this will improve the readablility of the code later on */
	 nx = rho->grid[0];
	 ny = rho->grid[1];
	 nz = rho->grid[2];	 
       	 
	 rho->npoints = nx*ny*nz;
	 
	 /* construct the 'datacell' from the 'delta's */
	 for(i=0; i < 3; i++) rho->datacell.lat_A[i] = CONV * (nx-1) * delta[0][i];
	 for(i=0; i < 3; i++) rho->datacell.lat_B[i] = CONV * (ny-1) * delta[1][i];
	 for(i=0; i < 3; i++) rho->datacell.lat_C[i] = CONV * (nz-1) * delta[2][i];
	 
	 /* skip the next 2 lines */
	 for(i=0; i < 2; i++) fgets(buf, sizeof(buf), fh);
	 
	 /* allocate the sub-array for the field */
	 if(rho->field) free(rho->field);
	 if(rho->field = (float *) calloc(sizeof(float), rho->npoints))
	 { 
	    /* allocate temporary field, it will hold the raw data */ 
	    if(tmp_field = (float *) calloc(sizeof(float), rho->npoints))
	    {
	       /*** read the field data line by line into 'tmp_field' ***/

	       for(i=0; i < rho->npoints; i++)
	       {
		  fgets(buf, sizeof(buf), fh);
        	  if(sscanf(buf, " %f ", &tmp_field[i] ) == 1) cnt_RHO++;	          
	       }
	       	       
	       /* check for correctness again */
	       if(cnt_RHO == rho->npoints + 5)
	       {
		  /*** reorder the data from (y,z,x) to (x,y,z) order 
	               (first is the fastest) and save the result in 'rho->field' 
		       'i' is the general counter that increases continuously
		  ***/
		  for (i=0, iz=0; iz < nz; iz++)
		  {
	             for (iy=0; iy < ny; iy++)
		     {
	        	for (ix=0; ix < nx; ix++, i++)
			{
	        	   /* ind = f(ix,iy,iz)
			      maps the 3D-skalar field to the corresponding buffer index 
			      in the (y,z,x) ordering given by the lmto output
			   */
			   ind = ix*(ny*nz) + iy*nz + iz; 
			   rho->field[i] = tmp_field[ind];
			}
	             }
		  }
     
	       }
	       else err = print_err(err_crd, filename);
	       	       
	       free(tmp_field);
	    }
	    else err = print_err(err_cam,0L); /* 'tmp_field' allocation */   	    	    
	 }
	 else err = print_err(err_cam,0L); /* 'rho->field' memory allocation*/ 	          		 	 
      }
      else err = print_err(err_crd, filename);	  	 

      fclose(fh);
   }
   else err = print_err(err_cof, filename); /* open RHO file*/ 
  
  

   return(err);   

}



/********************** print_somedata **************************
   prints the processed data 
*/

void print_somedata(struct data_field *rho)
{
   int i;
   
   printf("Grid=%d %d %d\n", rho->grid[0], rho->grid[1], rho->grid[2]);
   printf("npoints=%d\n", rho->npoints);
   printf("Origin=%f %f %f\n", rho->origin[0], rho->origin[1], rho->origin[2]);
   
   /* print the first and the last 10 datapoints */
   for(i=0; i < 10; i++) printf("%f\n", rho->field[i]);  
   printf("\n"); 
   for(i = (rho->npoints-10); i < rho->npoints; i++) printf("%f\n", rho->field[i]);
      
} 


