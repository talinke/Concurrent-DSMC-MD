/**************************************************************************************************
 * Author:      Tim Linke (talinke@ucdavis.edu)
 * Date:        2023-06-02
 *************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* Define molecular properties */
#define MASS 2.6566962e-26 // mass of oxygen in kg
#define BOLTZMANN 1.380649e-23
#define CUTOFF 0.5

/* Define simulation parameters */
#define NUM_TIME_STEPS 100
#define NUM_CELLS_X 10
#define NUM_CELLS_Y 10
#define NUM_CELLS (NUM_CELLS_X * NUM_CELLS_Y)
#define SIMPAR 1 //ratio of real to simulated particles
#define PARTICLES_X 20
#define PARTICLES_Y 20
#define PARTICLES (PARTICLES_X * PARTICLES_Y)

typedef struct {double x, y;} Vec_2;
typedef struct {Vec_2  r, v, f; double mass;} Particle;
typedef struct {int counter; Particle* cell_atoms[4*(int)(PARTICLES/NUM_CELLS)];} Cell;

/* Declare molecular dynamics simulation functions */
void initialize_positions(Particle* atom, Vec_2 Dim);
void initialize_velocities(Particle* atom);
void update_particles(Cell cells[][NUM_CELLS_Y], int index_x, int index_y, double dt, Vec_2 Dim);
Vec_2 force_LJ(Vec_2 dr);
double energy_LJ(Vec_2 dr);



/* Declare DSMC simulation functions */
void drift_particles(Particle* atom, double dt);
void apply_BC(Particle* atom, Vec_2 Dim);
void initialize_cells(Cell cells[][NUM_CELLS_Y]);
void assign_particles_to_cells(Particle* atom, Cell cells[][NUM_CELLS_Y], Vec_2 Cell_Size);
void collide_particles_in_cells(Particle* atom, Cell cells[][NUM_CELLS_Y], Vec_2 Cell_Size, double dt);
void uplift_particles(Cell cells[][NUM_CELLS_Y], int index_x, int index_y);
void average_flow(Particle* atom, Cell cells[][NUM_CELLS_Y], Vec_2 velocity_field[][NUM_CELLS_Y]);