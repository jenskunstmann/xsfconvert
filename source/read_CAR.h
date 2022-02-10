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

/* this holds all the data that has to be stored in 
addition to the 'data' class */
struct data_CAR
{
     
   int			nspecies;	/* number of atomic species present */
	
   float		scale,		/* lattice constant or scaling factor */
   			volume;		/* volume of unit cell */
   
   struct data_species	*species;	/* this will just point to the structure that is 
   					   actually defined as part 'config' in v_data.h */
};



int read_CAR(struct config *conf, char *progname, struct data *gd);
int read_POS(FILE *fh, char *filename, struct data *ld, struct data_CAR	*ld_add);
int read_CHG(FILE *fh, char *filename, struct data_field *rho, float scale);
void print_data(struct data *ld, struct data_CAR *ld_add);
void print_somedata(struct data_field *rho);
