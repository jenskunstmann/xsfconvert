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


#include <stdio.h>
#include <stdlib.h>
#include "error.h"


#define		debug(a)		/* simple removal for additional debugging output */
#define		FILENAME_SIZE	256	/* maximum length of filename strings */
#define		SPECIES_SIZE	50	/* maximum number of atomic species that can be handled */
#define		LINEBUF_SIZE	130	/* size of the line buffer */


struct data_species
{
   int		num,	/* number of atoms of this species */
   		atnum;	/* atomic number of that species */
   char		*name;   /* atomic symbol of that species */
};

typedef struct
{
   char        *nom;
   int          Z;
} atom;


struct config
{
   char			*in_file,
   			*xsf_file; 
		
   int			struct_only,	/* 1: only read structure data
   				   	0: read structure AND field data */
   			scale_field,	/* 1: scale field data by 1/Vol
   				   	0: don't scale field values */
   			compress;	/* 1: compress the output file using gzip; 
				   	0: don't compress  */
			   
   float		a, b, c;	/* lattice constants of conv. unit cell 
   				   	in Angstroms,
				   	if(a,b,c != 0) {conv. != prim.}
				   	else {conv. = prim.}  */
				   
   struct data_species	species[SPECIES_SIZE];	/* we need to read the atomic number from command line 
   						   that's why it is defined here as part of 'config' */
   			   
};
