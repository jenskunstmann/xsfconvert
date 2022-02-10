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
#include	"read_CTRL.h"
#include	"read_RHO.h"


#ifndef PROGNAME
  #define	PROGNAME	"lm2xsf"	/* this program */
#endif  
#define		VERSION		"1.03"
#define		RELDATE		"12 September 2014"
#define		AUTHOR		"Jens Kunstmann"


int parse_args(int argc, char **argv, struct config *conf);
int print_help(void);



/********************* main ****************************** 
lm2xsf is a conversion tool for 3D structures and
two real space scalar fields.

The two fields RHO and ELF differ in the fact that 
RHO has a dimension, which has to be converted by the lm2xsf,
while ELf is dimensionless and doesn't need to be converted.

*/

int main(int argc, char **argv)
{		
   int			err = 0,
   			i;
			
   char			cmd[FILENAME_SIZE];			
		
   struct data  	*gd;
   struct config	conf;
   
   
   /* initialize the configuration to the default values */
   conf.ctrl_file = "CTRL";
   conf.rho_file  = "RHO";
   conf.elf_file  = "ELF";
   conf.xsf_file  = "this.xsf";
   conf.compress  = 1;	/* compress output file, using gzip */
   conf.noempty   = 1;
   conf.a 	  = 0;	/* by default: conv.unit cell = prim.unit cell */
   conf.b 	  = 0;
   conf.c 	  = 0;


   /* standard output */
   printf(" +++ %s version %s (%s) written by %s +++\n", PROGNAME, VERSION, RELDATE, AUTHOR);
   
   if( !(err = parse_args(argc, argv, &conf)) )
   {      	   
      if(gd = data_init())  /* construct the data class */
      {
         debug( printf("conf.a=%f conf.b=%f conf.c=%f\n", conf.a, conf.b, conf.c) );
	 
	 printf("Read CTRL file named '%s':\n", conf.ctrl_file);
	 if( !(err = read_CTRL(&conf, gd)) ) 
	 { 
	    /* additional to reading, the units are converted;
	       'err' is not used because these files may simply be not present
	       which is not an error in the sense of the program */
	    printf("Read RHO file named '%s':\n", conf.rho_file);	 
	    read_RHO(conf.rho_file, PROGNAME, gd);

	    /* dimensionless units, no conversion necessary */
	    printf("Read ELF file named '%s':\n", conf.elf_file);
	    read_ELF(conf.elf_file, PROGNAME, gd);

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
   	i;
   
   if(argc > 0)
   {     
      /* argv[0] is the name of the program */
      for(i=1; (i < argc) && (err == 0); i++)
      {
	 if(argv[i][0] == '-')
	 {
	    	switch(argv[i][1])
		{
		   case 'D':
		   case 'd':
		      conf->compress = 0;
		      break;
		      
		   case 'C':
		   case 'c':		      
		      i++;
		      conf->ctrl_file = argv[i];
		      break;
		      
		   case 'R':
		   case 'r':
		      i++;
		      conf->rho_file = argv[i];
		      break;
		      
		   case 'E':
		   case 'e':
		      i++;
		      conf->elf_file = argv[i];
		      break;
		      
		   case 'O':
		   case 'o':
		      i++;
		      conf->xsf_file = argv[i];
		      break;
		      
		   case 'S':
		   case 's':
		      conf->noempty = 0;
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
		      
		   default:
		      err = print_help();   
		               
		}
		if(i >= argc) err = print_help();
	 }
	 else err = print_help();
      }
   }
   
   return(err);
}



/******************************* print_help ********************************
print the help text
*/

int print_help(void)
{
   printf("Syntax error!\n");
   printf("\nUsage: %s [-s] [-d] [-u a b c] [<option> file <option> file ...]\n\n", PROGNAME);
   printf("options:\n");
   printf("  -s, -S\tdo not remove empty spheres\n"); 
   printf("  -d, -D\tdo not compress the xsf output file\n");   
   printf("  -u, -U\tdefine lattice constants a,b,c of conventional unit cell\n\t\t(a,b,c are in units of the LMTO lattice constant 'ALAT')\n");  
   printf("  -c, -C\tCTRL file\n");
   printf("  -r, -R\tRHO or RHOS file (density is converted to (1/Ang)^3)\n");   
   printf("  -e, -E\tELF file (dimensionless, no conversion)\n");
   printf("  -o, -O\txsf output file (without .gz extension)\n\n");
   return(10);
}
