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


#include	<math.h>
#include	<string.h>
#include	"v_all.h"
#include	"data.h"
#include	"read_CAR.h"


/********************* read_CAR *************************
Read data from 'CAR' file.

First read the crystal structure and save in 'gd' 
then, if present, read the ELF or charge field and save.

conf:		configuration data 
progname:	name of the program 
gd:		the central data object
*/


int read_CAR(struct config *conf, char *progname, struct data *gd)
{

   FILE		*fh;	   
   	   
   int		err = 0,	/* error value */
   		i; 		/* loop counters */	
		
   float	scale;		/* the scaling factor for the field data */
   
   char		*comment = 0L;		
			
   struct data		*ld;	/* local data */     
   struct data_CAR	ld_add;	/* additional data, than is necessary
   				   but cannot be saved in 'ld' */

   	
   if(ld = data_init()) /* allocate and initialize the LOCAL data object */
   {
      if (fh = fopen(conf->in_file, "r"))
      {	 
	 ld_add.species = conf->species; /* connect species pointer with data */
	  
	 printf("  Read crystal structure: ");
	 if( !(err = read_POS(fh, conf->in_file, ld, &ld_add)) )
	 {	    
	    /* generate the conventinal unit cell from the configuration data */
	    if((conf->a != 0) && (conf->b != 0) && (conf->c != 0))
	    {
	       ld->convcell.lat_A[0] = conf->a;
	       ld->convcell.lat_A[1] = 0;
	       ld->convcell.lat_A[2] = 0;

	       ld->convcell.lat_B[0] = 0;
	       ld->convcell.lat_B[1] = conf->b;
	       ld->convcell.lat_B[2] = 0;

	       ld->convcell.lat_C[0] = 0;
	       ld->convcell.lat_C[1] = 0;
	       ld->convcell.lat_C[2] = conf->c;
	    }
	    else /* conventional = primitive (default) */
	    {
	       for(i=0; i < 3; i++) ld->convcell.lat_A[i] = ld->primcell.lat_A[i];
	       for(i=0; i < 3; i++) ld->convcell.lat_B[i] = ld->primcell.lat_B[i];
	       for(i=0; i < 3; i++) ld->convcell.lat_C[i] = ld->primcell.lat_C[i];		     
	    }
	    
	    /* save structure data to global data class */
	    data_put(gd, ld, PUT_STRUCT);
	    
	    if(conf->struct_only)
	    {
	       printf("  Read structure data, only.\n");
	    }
	    else
	    {
	       /* configure the scaling of the field values */
	       if(conf->scale_field)
	       {
	          scale = 1/ld_add.volume;
		  comment = "    Scale field values by 1/Vol_cell.\n";
	       }	  
	       else scale = 1;
	       
	       /* err=10: error occurred, 
		  err=1 : there is no field data
		  err=0 : field data was read in correctly */
	       if( !(err = read_CHG(fh, conf->in_file, &ld->rho, scale)) )
	       {
		  printf("  Read field data:\n");
		  if(comment) printf(comment);

		  /* make identifier string*/
		  sprintf(ld->rho.identifier, "%s_%s", conf->in_file, progname);
		  ld->rho.identifier[ID_SIZE-1] = 0; /* if the buffer is too small, terminate by hand */

		  /* Origin is (0,0,0) for VASP */
		  for(i=0; i < 3; i++) ld->rho.origin[i] = 0.0;

		  /* data cell is the same as the unit cell for VASP */
		  for(i=0; i < 3; i++) ld->rho.datacell.lat_A[i] = ld->primcell.lat_A[i];
		  for(i=0; i < 3; i++) ld->rho.datacell.lat_B[i] = ld->primcell.lat_B[i];
		  for(i=0; i < 3; i++) ld->rho.datacell.lat_C[i] = ld->primcell.lat_C[i];

		  /* send this field to global data class */
		  data_put(gd, ld, PUT_RHO);

	       }
	       if(err == 1)
	       {
	          err = 0; /* having no data field doesn't mean that this is an error */
		  printf("   File does not contain field data.\n");
	       }	  
	    }
	    
	 }
	 	 	 
	 fclose(fh);
      }
      else err = print_err(err_cof, conf->in_file); /* open LMFS file*/ 
      
      data_free(ld);   
   }
   else err = print_err(err_cam, 0L);  /* construct the ld object */

  

   return(err);
}



/************************* read_POS **********************************
Read the crystal structure part from the input file

fh:		filehandle to read from
filename:	filename of filehandle for error msges
ld:		local data to write to
ld_add: 	additional local data 
*/

int read_POS(FILE *fh, char *filename, struct data *ld, struct data_CAR	*ld_add)
{

   char		buf[LINEBUF_SIZE],
   		*pos;		/* points to a position/character in buf  */
		
   int		err = 0,
   		isdirect,	/* =1 direct coordinates, =2 cartesion coordinates, =0 error */
   		i,j,k,s, 	/* loop counters */
   		cnt_POS = 0,	/* line count for successfully read lines */
		skiped,		/* when skipping white space in buf, this holds the number of
				   character skipped */
		vasp52 = 0;	/* =0 (VASP 4.6 or older) ,=1 (VASP 5.2 or higher) */
		
   float	A[3],		/* temporary lattice vectors */ 
   		B[3], 
		C[3],
		X[3],		/* temporary position vector */	
   		tmp_vol;	/* temporary volume */

   atom periodic_table[110] = 
      {{"Ac",89}, {"Al",13},{"Am",95},{"Sb",51},{"Ag",47},{"Ar",18},{"As",33},{"At",85},{"N",7},{"Ba",56},
      {"Bk",97}, {"Be",4},{"Bi",83},{"Bh",107},{"B",5},{"Br",35},{"Cd",48},{"Ca",20},{"Cf",98},{"C",6},
      {"Ce",58},{"Cs",55},{"Cl",17},{"Cr",24},{"Co",27},{"Cu",29},{"Cm",96},{"Ds",110},{"Db",105},{"Dy",66},
      {"Es",99},{"Er",68},{"Sn",50},{"Eu",63},{"Fe",26},{"Fm",100},{"F",9},{"Fr",87},{"Gd",64},{"Ga",31},
      {"Ge",32},{"Hf",72},{"Hs",108},{"He",2},{"Ho",67},{"H",1},{"In",49},{"I",53},{"Ir",77},{"Kr",36},
      {"La",57},{"Lw",103},{"Li",3},{"Lu",71},{"Mg",12},{"Mn",25},{"Mt",109},{"Md",101},{"Hg",80},{"Mo",42},
      {"Nd",60},{"Ne",10},{"Np",93},{"Ni",28},{"Nb",41},{"No",102},{"Au",79},{"Os",76},{"O",8},{"Pd",46},
      {"P",15},{"Pt",78},{"Pb",82},{"Pu",94},{"Po",84},{"K",19},{"Pr",59},{"Pm",61},{"Pa",91},{"Ra",88},
      {"Rn",86},{"Re",75},{"Rh",45},{"Rb",37},{"Ru",44},{"Rf",104},{"Sm",62},{"Sc",21},{"Sg",106},{"Se",34},
      {"Si",14},{"Na",11},{"S",16},{"Sr",38},{"Ta",73},{"Tc",43},{"Te",52},{"Tb",65},{"Tl",81},{"Th",90},
      {"Tu",69},{"Ti",22},{"W",74},{"U",92},{"V",23},{"Xe",54},{"Yb",70},{"Y",39},{"Zn",30},{"Zr",40}};

   
   /*** read in the header sequentially  ***/
   
   /* skip the comment line */
   fgets(buf, sizeof(buf), fh);  
   
   /* read lattice constant/-volume */
   fgets(buf, sizeof(buf), fh);     
   if(sscanf(buf, " %f ", &ld_add->scale) == 1) cnt_POS++;
   
   /* read the three lattice vectors */
   fgets(buf, sizeof(buf), fh);  
   if(sscanf(buf, " %f %f %f", &A[0], &A[1], &A[2]) == 3)
      cnt_POS++;
   fgets(buf, sizeof(buf), fh);  
   if(sscanf(buf, " %f %f %f", &B[0], &B[1], &B[2]) == 3)
      cnt_POS++;
   fgets(buf, sizeof(buf), fh);  
   if(sscanf(buf, " %f %f %f", &C[0], &C[1], &C[2]) == 3)
      cnt_POS++;

   /* parse the atom names, this lines was added in VASP 5.2 */
   fgets(buf, sizeof(buf), fh);
   pos = buf;
   for(s=0; s < SPECIES_SIZE; s++)
   {
      if( !(pos = strpbrk(pos, "HLBCNOFMASPKTVZGRYIXWULEH")) ) break; 	/* find first character of first/next number in buf  */
      sscanf(pos, "%s ", &ld_add->species[s].name);	/* read letter */
      if( !(pos = strpbrk(pos, " \t\n")) ) break;	/* find position in buf where this number is over */
   }
   if(s)
   {
      cnt_POS++;
      vasp52 = 1;
   }   
       
  
   /* read the species information line
     from VASP 5.2 on this line has to be read in,
     for earlier versions the string 'buf' is simply parsed again */
   if(vasp52) fgets(buf, sizeof(buf), fh);          
   pos = buf;
   for(s=0; s < SPECIES_SIZE; s++)
   {         
      /* strpbrk(s, tokens) returns a pointer to the FIRST position in the string 's' 
         where ONE of the characters given in 'tokens' is situated
	 and NULL if none of these characters can be found.  */
      if( !(pos = strpbrk(pos, "123456789")) ) break; 	/* find first character of first/next number in buf  */
      sscanf(pos, "%d ", &ld_add->species[s].num);	/* read number */
      if( !(pos = strpbrk(pos, " \t\n")) ) break;	/* find position in buf where this number is over */
   }
   if(s)
   {
      cnt_POS++;
      ld_add->nspecies = s;	 
   }              

   /* read next line */
   fgets(buf, sizeof(buf), fh);     
   skiped = strspn(buf, " \t"); /* skip white space */

   /* skip line if "selective dynamics" is switched on*/
   if(buf[skiped] == 'S' || buf[skiped] == 's')
   {
      fgets(buf, sizeof(buf), fh);     
      skiped = strspn(buf, " \t"); /* also skip white space here */     
   } 

   /* figure out if position mode is "cartesian" or "direct" */
   switch(buf[skiped])
   {
      case 'c':
      case 'C':
      case 'k':
      case 'K':
	 isdirect = 2;
	 break;

      case 'd':
      case 'D':	   
	 isdirect = 1;
	 break;

      default:
	 isdirect = 0; /* means: read error */
	 break;      
   }
   if(isdirect) cnt_POS++;
      
   
   /*** check if header was read in correctly 
        and do necessary computations;
        VASP 5.2: 7 lines are read 
        prior to VASP 5.2 only 6 lines are read ***/      
	
   if( ((cnt_POS == 7) && vasp52) || ((cnt_POS == 6) && (!vasp52)) ) 
   {                 
      /* figure out how the unit cell has to be scaled */

      /* calculate volume of unit cell = Spatprodukt*/
      tmp_vol = C[0]*(A[1]*B[2]-A[2]*B[1]) + C[1]*(A[2]*B[0]-A[0]*B[2]) + C[2]*(A[0]*B[1]-A[1]*B[0]);
      debug( printf("tmp_vol=%f\n", tmp_vol) );

      /* negative scaling factor means volume of unit cell */
      if(ld_add->scale < 0)
      {
         debug( printf("negative scaling \n") );
	 ld_add->volume = (-1)*ld_add->scale;	 
	 ld_add->scale = powf(ld_add->volume/tmp_vol,1.0/3.0); 
	 /* scaling factor is third root of ratio V_soll/V_ist */	 	 
      }
      else ld_add->volume = ld_add->scale * ld_add->scale * ld_add->scale * tmp_vol;
 
      /* now do the scaling*/
      for(i=0; i < 3; i++) ld->primcell.lat_A[i] = ld_add->scale * A[i];
      for(i=0; i < 3; i++) ld->primcell.lat_B[i] = ld_add->scale * B[i]; 	 
      for(i=0; i < 3; i++) ld->primcell.lat_C[i] = ld_add->scale * C[i];

      
      /* sum up 'species[j].num' to 'natoms' */
      ld->natoms = 0;
      for(j=0; j < ld_add->nspecies; j++) ld->natoms += ld_add->species[j].num;	 
      
      printf("\n    Found %d atom(s) and %d type(s).\n", ld->natoms, ld_add->nspecies);

      /* match the atomic symbols with the atomic numbers (VASP 5.2 and beyond) 
         or assign serial atomic number (prior to VASP 5.2)
         do NOT overrrive the assignments via the command line  */
      for(j=0; j < ld_add->nspecies; j++)
      {
         if (ld_add->species[j].atnum == 0)  /* no assignments via the command line specified */
         {
            ld_add->species[j].atnum = j+1;  /* the serial index = field index + 1  */
            if(vasp52)
            {
               for(i=0; i < 110; i++)
               {
                  if ( strcmp((char *) &ld_add->species[j].name, (char *) periodic_table[i].nom) == 0 )
                  { 
                     ld_add->species[j].atnum = periodic_table[i].Z;
                     break;
                  }
               }
            }
         }
      }

      /* write out informations about the atomic types */
      for(i=0; i < ld_add->nspecies; i++)
      {
         if (vasp52)
         {
	    printf("    Type %d: %s: %d atom(s), assigned atomic number Z=%d\n", i+1,
                   &ld_add->species[i].name, ld_add->species[i].num, ld_add->species[i].atnum);      
         }
         else
         {
	    printf("    Type %d: %d atom(s), assigned atomic number Z=%d\n", i+1,
                   ld_add->species[i].num, ld_add->species[i].atnum);      
         }
      }


      /*** now read the atomic positions ***/
      
      /* allocate memory for atoms data */
      if(ld->list) free(ld->list);
      if(ld->list = (struct data_alist	*) calloc(sizeof(struct data_alist), ld->natoms))
      {
         /* loop to read the atomic positions	    
	    	s: loop over species
	    	j: loop over all atoms of ONE species
	    	i: total loop counter 
	 */
	 for(i=0, s=0; s < ld_add->nspecies; s++)
	 {
	    for(j=0; j < ld_add->species[s].num; j++, i++)
	    {
	       fgets(buf, sizeof(buf), fh);
	       if(sscanf(buf, " %f %f %f", &X[0], &X[1], &X[2]) == 3) cnt_POS++;

              ld->list[i].atnum = ld_add->species[s].atnum;

	       if(isdirect == 1)	/* direct coordinates */
	       {
		  for(k=0; k < 3; k++)
		     ld->list[i].coord[k] = X[0] * ld->primcell.lat_A[k]
		  			  + X[1] * ld->primcell.lat_B[k]
					  + X[2] * ld->primcell.lat_C[k];	       
	       }
	       else 		/* isdirect == 2: cartesian coordinates*/
	       { 
	          for(k=0; k < 3; k++)
		     ld->list[i].coord[k] = ld_add->scale * X[k]; 		    
	       }	  
		  
	    }
	 }
	
	 
	 /* last check if everything was read in correctly */	 
	 if( ((cnt_POS == 7 + ld->natoms) && vasp52) || ((cnt_POS == 6 + ld->natoms) && (!vasp52)) ) 
	 {
	    debug (printf("isdirect=%d\n", isdirect) );
            debug( print_data(ld, ld_add) );      
	 }
	 else err = print_err(err_crd, filename);
	 
      }
      else err = print_err(err_cam, 0L);
         
   }
   else err = print_err(err_crd, filename);
   
   
   return(err);
}



/************************* read_CHG **********************************
Read the crystal structure part from the input file

fh:		filehandle to read from
filename:	filename of filehandle for error msges
rho:		local data_field to write to
scale: 		scaling factor for the field data: either 1 or 1/Vol 
*/

int read_CHG(FILE *fh, char *filename, struct data_field *rho, float scale)
{

   char		buf[LINEBUF_SIZE],
   		*pos;
		
   int		err = 0,
   		count,
   		i,ix,iy,iz, 	/* loop counter and indices*/
		jx,jy,jz,
		nx,ny,nz,
		ind,
		npnt_tmp,	/* size of the temporary field */
   		cnt_POS = 0;	/* line count for successfully read lines */
		
   float	tmp_f,		/* temporary float buffer */
   		*tmp_field;	/* temporary field */			
		
		
   /* skip the empty line */
   fgets(buf, sizeof(buf), fh);  		
   
   /* read grid size */
   fgets(buf, sizeof(buf), fh);
   sscanf(buf, " %d %d %d", &rho->grid[0], &rho->grid[1], &rho->grid[2]);

   /* check carefully whether file continues with charge information */
   if((rho->grid[0] > 0) && (rho->grid[1] > 0) && (rho->grid[2] > 0))
   {      

      /* the XCrysDen grid must be extended by the redundant boundary points */
      rho->grid[0]++;
      rho->grid[1]++;
      rho->grid[2]++;

      /* this will improve the readablility of the code later on */
      nx = rho->grid[0];
      ny = rho->grid[1];
      nz = rho->grid[2];      

      /* number of points on the XCrysDen grid */
      rho->npoints = nx * ny * nz;
      debug( printf("rho->npoints=%d\n", rho->npoints) );

      /* number of points on the VASP grid */
      npnt_tmp = (nx-1) * (ny-1) * (nz-1);
      debug( printf("npnt_tmp=%d\n", npnt_tmp) );

      /* allocate the sub-array for the field */    
      if(rho->field) free(rho->field);
      if(rho->field = (float *) calloc(sizeof(float), rho->npoints))
      {            
         /* allocate temporary field, it will hold the raw data from VASP */ 
	 if(tmp_field = (float *) calloc(sizeof(float), npnt_tmp))
	 {
	    
	    /*** read the raw field data ***/
	    
	    count = 0;
	    while((count < npnt_tmp) && (err == 0))
	    {
	       if(fgets(buf, sizeof(buf), fh))
	       {
		  /* read all floats from that line */
		  pos = buf;
		  while(1)
		  {         
	        	/* find first character of first/next number in buf  */
		     if( !(pos = strpbrk(pos, "1234567890.-")) ) break; 		          
		     if(sscanf(pos, "%f ", &tmp_f) == 1) /* read number */
		     {	          
			tmp_field[count] = scale * tmp_f;
			count++;		  
		     }
	        	/* find position in buf where this number is over */
		     if( !(pos = strpbrk(pos, " \t\n")) ) break;	
		  }	    
	       }
	       else err = print_err(err_crd, filename);	     
	    }
	    debug( printf("number of field data points read in: count=%d\n",count));
	    
   	    /* check if all numbers are read */
	    if(npnt_tmp == count)
	    {
	       
	       /*** add redundant boundary points ***/
	       
	       for (i=0, iz=0; iz < nz; iz++)
	       {
	          for (iy=0; iy < ny; iy++)
		  {
	             for (ix=0; ix < nx; ix++, i++)
		     {
		        /* ix,iy,iz run over the (bigger) CyrsDen grid and 
			   jx,jy,jz run over the (smaller) VASP grid */
			jx = ix;
			jy = iy;
			jz = iz;
			
			/* if we have a boundary point */
			if(ix == nx-1) jx = 0;
			if(iy == ny-1) jy = 0;
			if(iz == nz-1) jz = 0;

			ind = jz*((ny-1)*(nx-1)) + jy*(nx-1) + jx; 
			rho->field[i] = tmp_field[ind];
			
		     }
	          }
	       }
   
	    }
	    else err = print_err(err_crd, filename);	    	    

	    debug( print_somedata(rho) );  
	    free(tmp_field);	 	 
	 }
	 else err = print_err(err_cam, 0L);	/* 'tmp_field' allocation */ 
      }
      else err = print_err(err_cam, 0L);	/* 'rho->field' memory allocation*/        
   }
   else err = 1;
   /* error code 1 means that ther is simply no field written in this file,
      but there is no serious error */
   
   return(err);
		
}	



/********************** print_data **************************
   prints the structural data 
*/

void print_data(struct data *ld, struct data_CAR *ld_add)
{
   int i;
   
   printf("scaling factor=%f\nvolume=%f\n", ld_add->scale, ld_add->volume);
   
   printf("lat_A=[%f %f %f]\nlat_B=[%f %f %f]\nlat_C=[%f %f %f]\n",  
       ld->primcell.lat_A[0], ld->primcell.lat_A[1], ld->primcell.lat_A[2],
       ld->primcell.lat_B[0], ld->primcell.lat_B[1], ld->primcell.lat_B[2],
       ld->primcell.lat_C[0], ld->primcell.lat_C[1], ld->primcell.lat_C[2]);   
       
   for(i=0; i < ld->natoms; i++)
   {	       	         
      printf("%d: atnum=%d,  coord=[%f %f %f]\n",
      i+1, ld->list[i].atnum, ld->list[i].coord[0], ld->list[i].coord[1], ld->list[i].coord[2]);
   }

} 



/********************** print_somedata **************************
   prints the some field data 
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
	
