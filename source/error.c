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
#include "error.h"

/* declare global error codes */

char *err_cof = "cannot open file '%s'.";
char *err_crd = "cannot read data in file '%s'.";
char *err_cam = "cannot allocate memory.";
char *err_cf  = "command failed: '%s'.";



/**************************** print_err ******************************
   prints a formatted error message to stderr and return a error code
*/

int print_err(char *formatstr, char *inputstr)
{
   fprintf(stderr, "   ERROR: ");
   if(inputstr) fprintf(stderr, formatstr, inputstr);
   else fprintf(stderr, formatstr);
   fprintf(stderr, "\n");
   
   return(10);

}
