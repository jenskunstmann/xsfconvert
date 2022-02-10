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

int read_RHO(char *filename, char *progname, struct data *gd);
int read_ELF(char *filename, char *progname, struct data *gd);
int read_field(char *filename, struct data_field *rho);
void print_somedata(struct data_field *rho);
