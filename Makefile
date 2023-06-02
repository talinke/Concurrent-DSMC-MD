CC=gcc
CFLAGS=-I.

DSMC-MD: DSMC-MD-Coupling.c DSMC_func.c MD_func.c DSMC-MD.h
	$(CC) -o DSMC-MD DSMC-MD-Coupling.c DSMC_func.c MD_func.c
