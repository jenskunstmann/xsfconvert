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

#include	<string.h>
#include	"lm_all.h" 
#include	"data.h"
#include	"read_CTRL.h"



/********************* read_CTRL *************************
Read data from 'CTRL' file.

It reads the class DIM to get the dimension of the data,
checks for success, then it reads STRUC, CLASS and SITE, checks 
for success again, then combines the data of the CLASS and SITE
to a single data set and finally sends the data to the 
central data object.
 
conf:		configuration data 
gd:		the central data object
*/

int read_CTRL(struct config *conf, struct data *gd)
{

   FILE		*fh;	
   	
   char		buf[80];   	/* line buffer */

   int		err = 0,	/* error value */
   		i, j,		/* loop counters */
   		rd_DIM = 0,	/* tags the classes that are read in successfully */
		rd_STRUC = 0,
		rd_CLASS = 0,
		rd_SITE = 0;

   struct data_CTRL	*ld;		/* local data */  
   struct data		*cntr_dat;	/* local instance of the central data */
   
   
   	
   if(ld = CTRL_init()) /* allocate and initialize the LOCAL data object */
   {
      if (fh = fopen(conf->ctrl_file, "r"))
      {
	 /*** read class DIM  ***/       

	 /* read the file line by line and store the line in 'buf' */
	 debug( printf("Read DIM class\n") );
	 while(fgets(buf, sizeof(buf), fh))
	 {
             if(sscanf(buf, "DIM NBAS=%d NCLASS=%d", &ld->nsites, &ld->nclass) == 2)
	     {
		rd_DIM = 1;
		printf("   Found %d atomic site(s) and %d atomic class(es).\n", ld->nsites, ld->nclass);
        	break;                  	  	  
	     }
	 }

	 /*** go back to beginning of the file and read the other classes ***/	

	 fseek(fh, 0L, SEEK_SET);          
	 if(rd_DIM)
	 {
	    /*** now read only the FIRST appearance of the classes: 
	         STRUC, CLASS, SITE ***/

	    /* read the file line by line and store the line in 'buf' */
	    while( fgets(buf, sizeof(buf), fh) && (err == 0) )
	    {    
	       /* strstr(): find a sub-string within another string */
	       if(strstr(buf, "STRUC ") && (rd_STRUC == 0)) 
	       {
		  debug( printf("call read_STRUC\n") );
		  err = read_STRUC(fh, conf->ctrl_file, buf, ld);
		  rd_STRUC = 1; /* in the case of an error the loop will abort
		  		   that's why we don't check the read_STRUC() */		  
		  continue;
	       }   		  
	       if(strstr(buf, "CLASS ") && (rd_CLASS == 0))
	       {
		  debug( printf("call read_CLASS\n") );
		  err = read_CLASS(fh, conf->ctrl_file, buf, ld);
		  rd_CLASS = 1;     		  
		  continue;
	       }		  
	       if(strstr(buf, "SITE ") && (rd_SITE == 0))
	       {
		  debug( printf("call read_SITE\n") );
		  err = read_SITE(fh, conf->ctrl_file, buf, ld);
		  rd_SITE = 1;
		  continue;
	       } 
	    } 
	    debug( print_rawdata(ld) );

	    /*** check if everything was read in correctly and 
	       continue with combining the CLASS and SITE data to 'ld->at_list' ***/

	    if( !err && rd_STRUC && rd_CLASS && rd_SITE) 
	    {	             		     
	       if( !(err = CTRL_SortOut(conf->noempty, ld)) )
	       {		 		  	       
		  /*** Now save all data to the global structure ***/

		  debug( print_finaldata(ld) );
		  if(cntr_dat = data_init())  /* construct the local data class properly */
        	  {
		     /* write data to it */	     
		     for(i=0; i < 3; i++) cntr_dat->primcell.lat_A[i] = ld->lat_A[i];
		     for(i=0; i < 3; i++) cntr_dat->primcell.lat_B[i] = ld->lat_B[i];
		     for(i=0; i < 3; i++) cntr_dat->primcell.lat_C[i] = ld->lat_C[i];
		     
		     /*** generate the conventinal unit cell.  
		          it will be different from the prim. cell, 
			  if a,b,c are defined in 'conf' ***/
		     
		     if((conf->a != 0) && (conf->b != 0) && (conf->c != 0))
		     {
		        cntr_dat->convcell.lat_A[0] = conf->a * ld->latconst;
			cntr_dat->convcell.lat_A[1] = 0;
			cntr_dat->convcell.lat_A[2] = 0;
			
			cntr_dat->convcell.lat_B[0] = 0;
			cntr_dat->convcell.lat_B[1] = conf->b * ld->latconst;
			cntr_dat->convcell.lat_B[2] = 0;
			
			cntr_dat->convcell.lat_C[0] = 0;
			cntr_dat->convcell.lat_C[1] = 0;
			cntr_dat->convcell.lat_C[2] = conf->c * ld->latconst;
		     }
		     else /* conventional = primitive (default) */
		     {
			for(i=0; i < 3; i++) cntr_dat->convcell.lat_A[i] = ld->lat_A[i];
			for(i=0; i < 3; i++) cntr_dat->convcell.lat_B[i] = ld->lat_B[i];
			for(i=0; i < 3; i++) cntr_dat->convcell.lat_C[i] = ld->lat_C[i];		     
		     }
		     
		     cntr_dat->natoms = ld->natoms;

		     cntr_dat->list = ld->at_list;                    
		     err = data_put(gd, cntr_dat, PUT_STRUCT);
		     cntr_dat->list = 0L; /* otherwise data_free() would free that buffer */
		     
		     data_free(cntr_dat); /* is the the destructor of the 'data' class */      
		  }
		  else err = print_err(err_cam, 0L);  	
	       }
	    }
	    else err = print_err(err_crd, conf->ctrl_file);	       	  
	    debug( printf("\n%s\n",buf) );

	 }
	 else err = print_err(err_crd, conf->ctrl_file); /* endif(rd_DIM)*/
	 
	 fclose(fh);
      }
      else err = print_err(err_cof, conf->ctrl_file); /* open CTRL file*/ 
      
      CTRL_free(ld);   
   }
   else err = print_err(err_cam, 0L);  /* construct the ld object */

  

   return(err);
}




/********************** CTRL_init ********************
   initializes the CTRL data structure,
   e.g. it allocates memory for it
   
   cannot return an error code because it is returning
   the object pointer ;)
*/

struct data_CTRL *CTRL_init(void)
{
   struct data_CTRL	*ld;
   
   if(ld = (struct data_CTRL *) calloc(sizeof(struct data_CTRL), 1))
   {
      ld->at_list = 0L;	/* no sub-arrays allocated */
      ld->cl_list = 0L;
      ld->si_list = 0L;
      return(ld);		   	       
   }
   else return(0L);
   
}



/********************** CTRL_free ********************
   frees the the CTRL data structure and the sub-arrays
*/

void CTRL_free(struct data_CTRL *ld)
{      
   if(ld->at_list) free(ld->at_list);
   if(ld->cl_list) free(ld->cl_list);
   if(ld->si_list) free(ld->si_list);
   free(ld);
}



/**************************  CTRL_SortOut ****************************************
combinines the CLASS and SITE data to 'ld->at_list'
converts the positions to cartesian coordinates in Angstrom units
depending on the format (pos_mod) they are given in

noempty:	1 = sort out empty spheres,
		0 = leave empty spheres
ld:		local data class
*/

int CTRL_SortOut(int noempty, struct data_CTRL *ld)
{
   int	err = 0,
   	i,j,k;
  
   /* data block might actually be smaller if empty speres are excluded */
   if(ld->at_list) free(ld->at_list);
   if(ld->at_list = (struct data_alist *) calloc(sizeof(struct data_alist), ld->nsites))
   {    		
      /* sort out empty spheres by simply erasing the 'name' label */

      ld->ntypes = ld->nclass;
      if(noempty)
      {	 
	 for(i=0; i < ld->nclass; i++)
	 {
	    if(ld->cl_list[i].atnum == 0)
	    {
	       strcpy(ld->cl_list[i].name, " ");
	       ld->ntypes--;
	    }   
	 }      
      }
      else printf("   Don't remove empty sphere(s).\n");


      /* match CLASS with SITE to 'ld->at_list' */

      ld->natoms = 0;
      for(i=0; i < ld->nsites; i++)
      {
	 for(j=0; j < ld->nclass; j++)
	 {
            /* strcmp(): returns 0 if the two strings are equal */
	    if(strcmp(ld->si_list[i].name, ld->cl_list[j].name) == 0)
	    {			      
	       ld->at_list[ld->natoms].atnum = ld->cl_list[j].atnum;
	       
	       /* convert to cartesian angstrom depending on whether 
 	          the positions are given as 'X=' or 'POS=' */
	       if(ld->pos_mod) /* 'X=' : units of the lattice vectors */
	       {	          
		  for(k=0; k < 3; k++)
		     ld->at_list[ld->natoms].coord[k]	= ld->si_list[i].coord[0] * ld->lat_A[k]
		  					+ ld->si_list[i].coord[1] * ld->lat_B[k]
							+ ld->si_list[i].coord[2] * ld->lat_C[k];
	       }
	       else /* 'POS=' : cartesian units */
	       {		  
		  for(k=0; k < 3; k++)
		     ld->at_list[ld->natoms].coord[k] = ld->latconst * ld->si_list[i].coord[k];       
	       }

               ld->natoms++;
	       ld->cl_list[j].numat++; /* numat is initialized to 0 during allocation */

	       break;
	    }
	 }
      }		     

      if(noempty)
      {
         if(ld->ntypes != ld->nclass)
	    printf("   Removed empty sphere(s):\n   %d site(s) and %d class(es) remain.\n", ld->natoms, ld->ntypes);	
	 else printf("   There are no empty spheres.\n");   
      } 

      /* print how many atoms of each atmic class there are */		     

      for(i=0; i < ld->nclass; i++)
      {
	 if(strcmp(ld->cl_list[i].name, " ") != 0)
	    printf("   There are %d sites of class '%s'.\n", ld->cl_list[i].numat, ld->cl_list[i].name);
      } 
      
   }
   else err = print_err(err_cam,0L);
   
     	     
   return(err);

}



/**************************** read_STRUC *****************************************
Read the STRUC class in the CTRL file and store in the local data class.
These are 4 lines.
	
fh:		filehandle for CTRL file
filename:	filename of CTRL file
abuf:		"given" buffer that stores the first line of the class
ld:		local data class
*/

int read_STRUC(FILE *fh, char *filename, char *abuf, struct data_CTRL *ld)
{
   
   char		buf[80];   	/* local line buffer */     
   int		err = 0,
   		i, 		/* loop counter */
   		cnt_STRUC = 0;	/* line count for successfully read lines */
   float	vec[3];		/* temporal vector */
		
   
   /* read the first line from the given buffer 'abuf' */
   if(sscanf(abuf, "STRUC ALAT=%f", &ld->latconst) == 1) cnt_STRUC++;
   ld->latconst = CONV * ld->latconst;  /* convert to Angstrom */

   /* continue reading from the local buffer 'buf' */
   fgets(buf, sizeof(buf), fh);  
   if(sscanf(buf, " PLAT=%f %f %f", &vec[0], &vec[1], &vec[2]) == 3)
   {      
      cnt_STRUC++;
      for(i=0; i < 3; i++) ld->lat_A[i] = ld->latconst * vec[i];
   }
    
   fgets(buf, sizeof(buf), fh);  
   if(sscanf(buf, " %f %f %f", &vec[0], &vec[1], &vec[2]) == 3)
   {      
      cnt_STRUC++;
      for(i=0; i < 3; i++) ld->lat_B[i] = ld->latconst * vec[i];
   }   
     
   fgets(buf, sizeof(buf), fh);  
   if(sscanf(buf, " %f %f %f", &vec[0], &vec[1], &vec[2]) == 3)
   {      
      cnt_STRUC++;
      for(i=0; i < 3; i++) ld->lat_C[i] = ld->latconst * vec[i];
   }
   
   if(cnt_STRUC != 4) err = print_err(err_crd, filename);
   
   return(err);
   
}



/**************************** read_CLASS *****************************************
Read the CLASS class in the CTRL file to the local data class. 
These are 'ld->nclass' lines.

fh:		filehandle for CTRL file
filename:	filename of CTRL file
abuf:		"given" buffer that stores the first line of the class
ld:		local data class
*/

int read_CLASS(FILE *fh, char *filename, char *abuf, struct data_CTRL *ld)
{
   
   char		buf[80];   	/* local line buffer */     
   int		cnt_CLASS = 0,	/* line count for successfully read lines */
   		err = 0,
   		i;		/* loop counter */
		
   if(ld->cl_list) free(ld->cl_list);
   if(ld->cl_list = (struct data_CLASS *) calloc(sizeof(struct data_CLASS), ld->nclass))
   {
      /* read the first line from the given buffer 'abuf' */    
      if(sscanf(abuf, "CLASS ATOM=%s Z=%d", ld->cl_list[0].name, &ld->cl_list[0].atnum) == 2) cnt_CLASS++;

      /* continue reading from the local buffer 'buf' */
      for(i=1; i < ld->nclass; i++)
      {
	 fgets(buf, sizeof(buf), fh);
	 if(sscanf(buf, " ATOM=%s Z=%d", ld->cl_list[i].name, &ld->cl_list[i].atnum) == 2) cnt_CLASS++;   
      }
      
      if(cnt_CLASS != ld->nclass) err = print_err(err_crd, filename);
   }
   else err = print_err(err_cam, 0L);   

   
   return(err);
   
}



/**************************** read_SITE *****************************************
Read the SITE class in the CTRL file to the local data class.
These are 'ld->nsites' lines.
Disciminate between positions given in cartesion coordinates (POS=)
and positions given in units of the translations vectors (X=),
depending on what is given in the first line of the class

It reads the first line from the given buffer 'abuf' 
and continues reading from the local buffer 'buf'
	
fh:		filehandle for CTRL file
filename:	filename of CTRL file
abuf:		"given" buffer that stores the first line of the class
ld:		local data class

pos_mod	 	positions given in: 
		0 = cartesian coord. (POS=) 
		1 = lattice vectors  (X=) 
*/

int read_SITE(FILE *fh, char *filename, char *abuf, struct data_CTRL *ld)
{
   
   char		buf[80],   	/* local line buffer */     
   		*frmt[] = {" ATOM=%s POS=%f %f %f", 
			   " ATOM=%s X=%f %f %f"}; 
		
   int		err = 0,
   		cnt_SITE = 0,	/* line count for successfully read lines */
   		i;		/* loop counter */


   
   if(ld->si_list) free(ld->si_list);
   if(ld->si_list = (struct data_SITE *) calloc(sizeof(struct data_SITE), ld->nsites)) 
   {         
      /*** determine if positions are gives as 'POS='  or 'X=' ***/
      
      if( sscanf(abuf, "SITE ATOM=%s POS=%f %f %f", ld->si_list[0].name, 
		 &ld->si_list[0].coord[0], 
		 &ld->si_list[0].coord[1], 
		 &ld->si_list[0].coord[2]) == 4 )
      {
         debug( printf("Positions given in cartesian coordinates.\n") );
	 ld->pos_mod = 0;
	 cnt_SITE++;	 		 
      }
      else
      {           
	 if( sscanf(abuf, "SITE ATOM=%s X=%f %f %f", ld->si_list[0].name, 
		    &ld->si_list[0].coord[0], 
		    &ld->si_list[0].coord[1], 
		    &ld->si_list[0].coord[2]) == 4 )
	 {
            debug( printf("Positions given in units of the transl. vectors.\n") );
	    ld->pos_mod = 1;
	    cnt_SITE++;	 		 
	 }
      }
      
      /*** read the positions ***/
      
      for(i=1; i < ld->nsites; i++)
      {
	 fgets(buf, sizeof(buf), fh);
	 if( sscanf(buf, frmt[ld->pos_mod], ld->si_list[i].name, 
	    	    &ld->si_list[i].coord[0], 
		    &ld->si_list[i].coord[1], 
		    &ld->si_list[i].coord[2]) == 4 ) 
	     cnt_SITE++;  
      }	 

      
      if(cnt_SITE != ld->nsites) err = print_err(err_crd, filename);
   }
   else err = print_err(err_cam, 0L);   
   
   return(err);
   
}



/********************** print_rawdata **************************
   prints the read in and unprocessed data
*/

void print_rawdata(struct data_CTRL *ld)
{
   int i;
   
   /* STRUC */
   printf("latconst=%f\nlat_A=[%f %f %f]\nlat_B=[%f %f %f]\nlat_C=[%f %f %f]\n\n", 
       ld->latconst,
       ld->lat_A[0], ld->lat_A[1], ld->lat_A[2],
       ld->lat_B[0], ld->lat_B[1], ld->lat_B[2],
       ld->lat_C[0], ld->lat_C[1], ld->lat_C[2]);
       
   /* CLASS */
   for(i=0; i < ld->nclass; i++)
   {
      printf("name='%s', atnum=%d\n", ld->cl_list[i].name, ld->cl_list[i].atnum);
   }
   
   /* SITE */
   for(i=0; i < ld->nsites; i++)
   {	       	         
      printf("name=%s, coord=[%f %f %f]\n",
      ld->si_list[i].name, ld->si_list[i].coord[0], ld->si_list[i].coord[1], ld->si_list[i].coord[2]);
   }
} 

 

/********************** print_finaldata **************************
   prints the processed data 
*/

void print_finaldata(struct data_CTRL *ld)
{
   int i;
   

   printf("lat_A=[%f %f %f]\nlat_B=[%f %f %f]\nlat_C=[%f %f %f]\n\n",  
       ld->lat_A[0], ld->lat_A[1], ld->lat_A[2],
       ld->lat_B[0], ld->lat_B[1], ld->lat_B[2],
       ld->lat_C[0], ld->lat_C[1], ld->lat_C[2]);
       
   for(i=0; i < ld->natoms; i++)
   {	       	         
      printf("%d: atnum=%d,  coord=[%f %f %f]\n",
      i, ld->at_list[i].atnum, ld->at_list[i].coord[0], ld->at_list[i].coord[1], ld->at_list[i].coord[2]);
   }
} 
