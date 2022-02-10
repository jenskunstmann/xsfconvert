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

#define		NAMESIZE	11 /* size of the array that stores the name string for: 'ATOM=%s' */



struct data_CLASS
{
   char		name[NAMESIZE];	/* name for the class, usually the atomic symbol */
   int		atnum;		/* atomic number Z */
   int		numat;		/* number of atoms of that class */
};

struct data_SITE
{
   char		name[NAMESIZE];
   float	coord[3];
};

struct data_CTRL 
{
   int			natoms, ntypes,
   			nsites,		/* number of sites (from DIM) */   
			nclass,		/* number of classes (from DIM) */
			pos_mod;	/* positions given in: 
					   0 = cartesian coord. (POS=) 
					   1 = lattice vectors  (X=) */
			
   float		latconst,	/* in Angstrom */
   			lat_A[3], lat_B[3], lat_C[3];   
			
   struct data_alist	*at_list;	/* 3 Lists */
   struct data_CLASS	*cl_list;
   struct data_SITE	*si_list;
};



int read_CTRL(struct config *conf, struct data *gd);
struct data_CTRL *CTRL_init(void);
void CTRL_free(struct data_CTRL *ld);
int CTRL_SortOut(int noempty, struct data_CTRL *ld);
int read_STRUC(FILE *fh, char *filename, char *abuf, struct data_CTRL *ld);
int read_CLASS(FILE *fh, char *filename, char *abuf, struct data_CTRL *ld);
int read_SITE(FILE *fh, char *filename, char *abuf, struct data_CTRL *ld);
void print_rawdata(struct data_CTRL *ld);
void print_finaldata(struct data_CTRL *ld);
