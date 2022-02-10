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

#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	"error.h"
#include	"data.h"



/****************** write_xsf ****************************
writes the data stored in the 'gd' object out to a xsf-file
  
conf:		configuration data
progname:	name of the executing program, 
gd:		global data
*/

int write_xsf(char *out_filename, struct data *gd)
{
   
   FILE		*fh;
   char 	identifier[20];
   int 		i,
   		err = 0;
		
   
   if(fh = fopen(out_filename, "w"))
   {
      /* STRUCT part */
      fprintf(fh, "CRYSTAL\nPRIMVEC\n");
      write_unitcell(fh, &gd->primcell);
      
      fprintf(fh, "CONVVEC\n");
      write_unitcell(fh, &gd->convcell);      
      
      fprintf(fh, "PRIMCOORD\n%d  1\n", gd->natoms);
      for(i=0; i < gd->natoms; i++)
      {	       	         
	 fprintf(fh, "%d    %f    %f    %f\n",
	        gd->list[i].atnum, gd->list[i].coord[0], gd->list[i].coord[1], gd->list[i].coord[2]);
      }
      
      /* RHO part */
      if(gd->rho.field) /* write only if data exists */
	 write_field(fh, &gd->rho);
      
      /* ELF part */
      if(gd->elf.field) /* write only if data exists */
	 write_field(fh, &gd->elf);
      
      fflush(fh); /* otherwise the compression doesn't work */
      fclose(fh);
   }
   else err = print_err(err_cof, out_filename);
  
  
   return(err); 

}



/******************* write_field *****************************
write field data (RHO or ELF) to the output file handle

for the RHO and ELF we use different BLOCK_DATAGRIDs in the xsf-file because
they can have a different number of datapoints, which is not possible 
for sets that are written within the same BLOCK

fh:		output file handle
rho:		field sub-class ('rho' or 'elf')
identifier:	string that serves a identifier for the data field
*/

void write_field(FILE *fh, struct data_field *rho)
{
   int	i;   
   
   /* write BEGIN_BLOCK */
   fprintf(fh,"\n");
   fprintf(fh, "BEGIN_BLOCK_DATAGRID_3D\n");
   fprintf(fh, "%s\n", rho->identifier); /* one-word identifier, rho->identifier*/
   fprintf(fh, "   BEGIN_DATAGRID_3D_%s\n", rho->identifier);

   /* grid size */
   fprintf(fh, "   %d    %d    %d\n", rho->grid[0], rho->grid[1], rho->grid[2]);
   
   /* origin */
   fprintf(fh, "   %f    %f    %f\n", rho->origin[0], rho->origin[1], rho->origin[2]);
   
   /* unit vectors */   
   write_unitcell(fh, &rho->datacell);
   
   /* write field data */
   for(i=0; i < rho->npoints; i++)
   {
      if( !(i % 5) ) fprintf(fh, "\n"); /* newline after 5 data points */
      fprintf(fh, "    %f", rho->field[i]);   
   }
   
   /* write END_BLOCK */
   fprintf(fh, "\n   END_DATAGRID_3D\n");
   fprintf(fh, "END_BLOCK_DATAGRID_3D\n");
      
}



/************************ write_unitcell ********************
write a specific unit cell data out to the fh

fh:	file handle
gc:	unit cell
*/

void write_unitcell(FILE *fh, struct data_unitcell *gc)
{
 
   fprintf(fh, "   %f    %f    %f\n", gc->lat_A[0], gc->lat_A[1], gc->lat_A[2]);   
   fprintf(fh, "   %f    %f    %f\n", gc->lat_B[0], gc->lat_B[1], gc->lat_B[2]);
   fprintf(fh, "   %f    %f    %f\n", gc->lat_C[0], gc->lat_C[1], gc->lat_C[2]);          

}



/********************** data_put **************************
  copies die data defined by the 'token' from the local to
  the global data object, see 'data.h'
  
  !! Allocates data lists !!!
  Don't forget to free the memory afterwards!
    
ld:	local data
gd:	global data 
what:	token that defines what will be read from the local data  
*/

int data_put(struct data *gd, struct data *ld, int token)
{
   int   i,		/* loop counter */
   	 err = 0;	/* error code */
   
   switch(token)
   {      
      case PUT_STRUCT: /* copy the atomic structure data */
	 data_put_unitcell(&gd->primcell, &ld->primcell);	 
	 data_put_unitcell(&gd->convcell, &ld->convcell);	 
	 gd->natoms = ld->natoms;
	 
	 if(gd->list) free(gd->list); /* for the case that it will be called many times */
	 if(gd->list = (struct data_alist *) calloc(sizeof(struct data_alist), gd->natoms))
	 {
	    memcpy(gd->list, ld->list, sizeof(struct data_alist) * gd->natoms);		   	       
	 }
	 else err = print_err(err_cam, 0L);         
	 break;
      
      
      case PUT_RHO: /* copy the charge density */
	 err = data_put_field(&gd->rho, &ld->rho);         
         break;
            
      case PUT_ELF: /* copy the ELF */
         err = data_put_field(&gd->elf, &ld->elf);
         break;
   }   
   
   return(err);   
   
}



/****************** data_put_field ********************************************
copies ELF or RHO data from a local to a global buffer

gf:	globel field
lf:	local field
*/

int data_put_field(struct data_field *gf, struct data_field *lf)
{
   int	err = 0,
   	i;
   
   strcpy(gf->identifier, lf->identifier);
   for(i=0; i < 3; i++) gf->grid[i] = lf->grid[i];   
   for(i=0; i < 3; i++) gf->origin[i] = lf->origin[i];   
   data_put_unitcell(&gf->datacell, &lf->datacell);   
   gf->npoints = lf->npoints;

   if(gf->field) free(gf->field); /* for the case that it will be called many times */	 	 	 	 
   if(gf->field = (float *) calloc(sizeof(float), gf->npoints))
   {
      memcpy(gf->field, lf->field, sizeof(float) * gf->npoints);		   	       
   }
   else err = print_err(err_cam, 0L); 

   return(err);
}



/************************ data_put_unitcell ********************
simply copies the three vectors of a unit cell from 
a local to a global data_unitcell instance

gc:	global unit cell
lc:	local unit cell
*/

void data_put_unitcell(struct data_unitcell *gc, struct data_unitcell *lc)
{
   int i;
   
   for(i=0; i < 3; i++) gc->lat_A[i] = lc->lat_A[i];
   for(i=0; i < 3; i++) gc->lat_B[i] = lc->lat_B[i];
   for(i=0; i < 3; i++) gc->lat_C[i] = lc->lat_C[i];
         
}



/********************** data_init ********************
   initializes the central or global data structure,
   e.g. it allocates memory for it
   
   cannot return an error code because it is returning
   the object pointer ;)
*/

struct data *data_init(void)
{
   struct data	*gd;
   
   if(gd = (struct data *) calloc(sizeof(struct data), 1))
   {
      gd->list = 0L;	/* no sub-arrays allocated */
      gd->rho.field = 0L;
      gd->elf.field = 0L;
      return(gd);		   	       
   }
   else return(0L);
   
}



/********************** data_free ********************
   frees the the global data structure and the sub-arrays
*/

void data_free(struct data *gd)
{
   if(gd->list) free(gd->list);
   if(gd->rho.field) free(gd->rho.field);
   if(gd->elf.field) free(gd->elf.field);
   free(gd);
}
