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
#include        <ctype.h>
#include	"v_all.h"
#include	"data.h"
#include	"read_CAR.h"



#ifndef PROGNAME
  #define	PROGNAME	"v2xsf"	/* this program */
#endif  
#define		VERSION		"1.05"
#define		RELDATE		"12 May 2017"
#define		AUTHOR		"Jens Kunstmann"


int parse_args(int argc, char **argv, struct config *conf);
int print_help(void);



/********************* main ****************************** 
converts VASP output files for structures and skalar fields
to the xsf format. Possible input files:
POSCAR, CONTCAR, CHG, CHGCAR, LOCPOT
*/

int main(int argc, char **argv)
{		
   int			err = 0,
   			i;
			
   char			outfile[FILENAME_SIZE],
   			cmd[FILENAME_SIZE];			
		
   struct data  	*gd;
   struct config	conf;
   

   
   /* initialize the configuration 
      to the default values */
      
   conf.in_file     = "CONTCAR";
   conf.xsf_file    = 0L;
   conf.struct_only = 0;	/* read structure AND field data */
   conf.scale_field = -1;	/* -1 means default treatment */
   conf.compress    = 1;	/* compress output file, using gzip */
   conf.a 	    = 0;	/* by default: conv.unit cell = prim.unit cell */
   conf.b 	    = 0;
   conf.c 	    = 0;
   /* initialize the atomic numbers */   
   for(i=0; i < SPECIES_SIZE; i++) conf.species[i].atnum = 0;


   /* standard output */
   printf(" +++ %s version %s (%s) written by %s +++\n", PROGNAME, VERSION, RELDATE, AUTHOR);
   
   if( !(err = parse_args(argc, argv, &conf)) )
   {      	   
      if(gd = data_init())  /* construct the data class */
      {         
         /* if no output filename was specified via command line,
	    simply add *.xsf to the input filename */
	 if(conf.xsf_file == 0L)
	 {
	    sprintf(outfile, "%s.xsf", conf.in_file);
	    conf.xsf_file = outfile;
	 }
	 
	 /* if no -c command was specified, automatically switch 'scale_field' 
	    on if the filename contains the "CHG" string, 
	    indicating one of the charge density files CHG or CHGCAR*/
	 if(conf.scale_field == -1)
	 {
	    if(strstr(conf.in_file, "CHG") || strstr(conf.in_file, "chg")) 
	       conf.scale_field = 1;
	    else conf.scale_field = 0;
	 }

	 debug( printf("conf.a=%f conf.b=%f conf.c=%f\n", conf.a, conf.b, conf.c) );
	 debug( printf("outfile='%s'\n", conf.xsf_file) );	 
	 debug( printf("scale_field=%d\nstruct_only=%d\n",conf.scale_field, conf.struct_only) );
	 
	 printf("Read from file '%s':\n", conf.in_file);
	 if( !(err = read_CAR(&conf, PROGNAME, gd)) )
	 {
	    printf("Write xsf-output to file '%s':\n", conf.xsf_file);
	    err = write_xsf(conf.xsf_file, gd); 
	    
	    if( !err && conf.compress )
	    {
	       printf("Compress xsf-output to '%s.gz':\n", conf.xsf_file);
	       sprintf(cmd, "gzip %s", conf.xsf_file);	       
	       if (system(cmd)) err = print_err(err_cf, cmd);
	    }	    	    
	 }	 

	 data_free(gd); /* is the the destructor of the 'data' class */      
      }
      else err = print_err(err_cam, 0L);   
   }
   
   
   return(err);
   
}



/**************************** pars_args *********************************
parse the command line options 

argc:	argument counter
argv:	agrument value
conf:	configuation data (see all.h)
*/

int parse_args(int argc, char **argv, struct config *conf)
{
   int	err = 0,
   	def_infile = 0,	/* remember if 'in_file' was already defined */
   	i,
	num;		/* for the -1,-2,-3,... options it will be 1,2,3,... */
   
   if(argc > 0)
   {     
      /* argv[0] is the name of the program */
      for(i=1; (i < argc) && (err == 0); i++)
      {
	 if(argv[i][0] == '-')
	 {
	    	switch(argv[i][1])
		{
		   case 'C':
		   case 'c':
		      i++;
		      if(i < argc)
		      {
			 switch(argv[i][0])
			 {
		            case 'Y':
			    case 'y':
			       conf->scale_field = 1;
			       break;

			    case 'N':
			    case 'n':
			       conf->scale_field = 0;
			       break;

			    default:
			       err = print_help();
			       break;         
			 }			 
		      }
		      break;
		   
		   case 'S':
		   case 's':
		      conf->struct_only = 1;
		      break;
		   
		   case 'D':
		   case 'd':
		      conf->compress = 0;
		      break;		      	      
		      
		   case 'O':
		   case 'o':
		      i++;
		      conf->xsf_file = argv[i];
		      break;
		       
		      
		   case 'U':
		   case 'u':
		      i++;
		      if(i < argc) /* this is necessary because the sscanf 
		                      will enter the string argv[i]*/
			 if(sscanf(argv[i], "%f", &conf->a) != 1)
			 {
		            err = print_help();
			    break;
			 } 
		      
		      i++;
		      if(i < argc)
			 if(sscanf(argv[i], "%f", &conf->b) != 1)
			 {
		            err = print_help();
			    break;
			 }
		      
		      i++;
		      if(i < argc)
			 if(sscanf(argv[i], "%f", &conf->c) != 1)
			 {
		            err = print_help();
			    break;
			 }		      		      		      
		      break;   
		      
		   default: /* here the -1, -2 -3, ..., -n  are read */
		      if(isdigit(argv[i][1])) 
		      {
			 sscanf(&argv[i][1], "%d", &num);
			 if((num > 0) && (num <= SPECIES_SIZE)) /* don't accept -0 or too big values */
			 {
			    i++;
			    if(i < argc)
			    {
			       if(sscanf(argv[i], "%d", &conf->species[num-1].atnum) != 1)
		        	  err = print_help();
			    }	  
			 }
			 else err = print_help();
		      }
		      else err = print_help();
		      break;   
		               
		} /* end switch() */		
		if(i >= argc) err = print_help();

		
	 }	/* end: if(argv[i][0] == '-')*/
	 else	/* define first item without '-' as 'in_file' */
	 {
	    if(def_infile) err = print_help();
	    else
	    {
	       conf->in_file = argv[i];
	       def_infile = 1;
	    }   
	 } 

	 	   
      } /* end loop over i */
   } /* end if(argc > 0) */
   
   return(err);
}



/******************************* print_help ********************************
print the help text
*/

int print_help(void)
{
   printf("Syntax error!\n");
   printf("\nUsage: %s\t[<infile>] [-1 <#1>] [-2 <#2>] .. [-n <#n>]\n\t\t[-c <y|n>] [-d] [-o <outfile>] [-s] [-u <a> <b> <c>]\n\n", PROGNAME);
   printf("infile\t\tinput file (e.g. POSCAR, CHGCAR, LOCPOT)\n\t\tdefault: CONTCAR\n");
   printf("options:\n");
   printf("  -1    \tspecifies atomic number of species 1\n"); 
   printf("  -2    \tspecifies atomic number of species 2\n");
   printf("   :    \t:                                  :\n");   
   printf("  -n    \tspecifies atomic number of species n\n\n");
   printf("  -c, -C\ty[es]: scale field data by 1/Vol_cell\n\t\tn[o]:  do not scale field data\n");        
   printf("  -d, -D\tdo not compress the xsf output file\n");   
   printf("  -o, -O\txsf output file (without .gz extension)\n\t\tdefault: <infile>.xsf\n");         
   printf("  -s, -S\tdo not read field data, convert crystal structure, only\n");   
   printf("  -u, -U\tdefine lattice constants a,b,c of conventional unit cell\n\t\t(a,b,c are in Angstrom)\n\n");       
   return(10);
}
