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


#define		debug(a)			/* simple removal for additional debugging output */
#define		CONV		0.52917721	/* conversion factor from Bohr radii to Angstrom */
#define		FILENAME_SIZE	256


struct config
{
   char		*ctrl_file,
   		*rho_file, 
		*elf_file, 
   		*xsf_file; 
		
   int		noempty,	/* 1 = sort out empty spheres
    				   0 = leave empty spheres */
		compress;	/* 1: compress the output file using gzip; 
				   0: don't compress  */
			   
   float	a, b, c;	/* lattice constants of conv. unit cell 
   				   in units of the lmto lattice constant,
				   if(a,b,c != 0) {conv. != prim.}
				   else {conv. = prim.}  */
};
