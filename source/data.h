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

/* define the tokens for data_put() */
#define		PUT_STRUCT	(1)
#define		PUT_RHO		(2)
#define		PUT_ELF		(3)

/* the length of the field identifier string */
#define		ID_SIZE		64


/* self explaining */
struct data_unitcell
{
   float	lat_A[3], 
   		lat_B[3], 
		lat_C[3];      
};


/* Will be allocated as array. Each element holds 
the atomic number and the coordinates of each atom */
struct data_alist
{
   int		atnum;
   float	coord[3];
};


/* holds the data for the ELF or RHO */
struct data_field
{
   char			identifier[ID_SIZE];
   int			grid[3];	/* dimension of data grid */
   float		origin[3];	/* position of first data point */
   struct data_unitcell	datacell;   	/* unit cell where the grid is defined */
   int			npoints; 	/* = grid[1]*grid[2]*grid[3] */   
   float		*field;		/* the actual data points */
};


/* holds all the data necessary to write an xsf-file */
struct data
{   
   /* STRUCT part */
   /*float		lat_A[3], lat_B[3], lat_C[3]; */
   struct data_unitcell	primcell;
   struct data_unitcell	convcell;   
   int			natoms;   
   struct data_alist	*list;	
   
   /* RHO part */
   struct data_field	rho;
   
   /* ELF part */
   struct data_field	elf;			
};



int write_xsf(char *out_filename, struct data *gd);
void write_field(FILE *fh, struct data_field *rho);
void write_unitcell(FILE *fh, struct data_unitcell *gc);
int data_put(struct data *gd, struct data *ld, int token);
int data_put_field(struct data_field *gf, struct data_field *lf);
void data_put_unitcell(struct data_unitcell *gc, struct data_unitcell *lc);
struct data *data_init(void);
void data_free(struct data *gd);


