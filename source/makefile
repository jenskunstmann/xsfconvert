#### simple definitions ####
#compiler options
CC = gcc
C_OPTS =

#linker options
LN = gcc
LN_OPTS = -lm

#preprocessor
PP = cpp
PP_OPTS =

#program name and location
PR1 = lm2xsf
PR2 = v2xsf
BINPATH = ../bin/

#source files
SOURCE1 = lm_main.o data.o error.o read_CTRL.o read_RHO.o
SOURCE2 =  v_main.o data.o error.o read_CAR.o

#### the rules ####
# the syntax for rules:
# targets : dependencies
#        <tab> command [;command]
#        <tab> command 
#	...

# calling 'make' without arguments will always make 
# the FIRST target 
all: lmto vasp

vasp: $(BINPATH)$(PR2)

lmto: $(BINPATH)$(PR1) 

$(BINPATH)$(PR1): $(SOURCE1)
	$(LN) -o $(BINPATH)$(PR1) $(SOURCE1) $(LN_OPTS)
		
$(BINPATH)$(PR2): $(SOURCE2)
	$(LN) -o $(BINPATH)$(PR2) $(SOURCE2) $(LN_OPTS)	

# this is an IMPLICIT target, telling make how to make
# the .o files from .c files. It's automatically invoked 
# when the source files in target '$(PR1)' are checked
%.o: %.c
	$(CC) -c $(C_OPTS) $*.c

# calls the preprocessor to define the macro PROGNAME	
lm_main.o: lm_main.c
	$(PP) $(PP_OPTS) -DPROGNAME='"$(PR1)"' lm_main.c  tmp.c
	$(CC) -c $(C_OPTS) -o lm_main.o tmp.c
	rm -f tmp.c
	
v_main.o: v_main.c
	$(PP) $(PP_OPTS) -DPROGNAME='"$(PR2)"' v_main.c  tmp.c
	$(CC) -c $(C_OPTS) -o v_main.o tmp.c
	rm -f tmp.c	

# 'make clean' will always work because no 
# dependencies are given
clean:
	rm -f *.o
	touch *.c	
