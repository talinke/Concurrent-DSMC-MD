/**************************************************************************************************
 * DSMC-MD-Coupling.c
 * Run a concurrent Direct Simulation Monte Carlo (DSMC) and Molecular Dynamics (MD) simulation.
 * This code is written for educational purposes and has not been validated for physical accuracy.
 * 
 * Author:      Tim Linke (talinke@ucdavis.edu)
 * Date:        2023-06-02
 * 
 * Compilation: gcc DSMC-MD-Coupling.c DSMC_func.c MD_func.c -o DSMC-MD-Coupling -lm
 * Execution:   ./DSMC-MD-Coupling
 *************************************************************************************************/


#include "DSMC-MD.h"

int main() {

    double     dt, egy_ave = 0;
    Vec_2      Dim, Num_Particles, Cell_Size;

    char    *rundat = "MD_Data.dat";
    char    datfile[32];
    FILE    *f_dat;
    sprintf(datfile, "%s", rundat);
    if ((f_dat = fopen(datfile, "w")) == NULL) {
        fprintf(stdout, "Can't open Dat file\n");
        fflush(stdout);
        exit(1);
    }
    char    *egydat = "Energy.dat";
    char    egyfile[32];
    FILE    *egy_dat;
    sprintf(egyfile, "%s", egydat);
    if ((egy_dat = fopen(egyfile, "w")) == NULL) {
        fprintf(stdout, "Can't open Dat file\n");
        fflush(stdout);
        exit(1);
    }
    
    /* Simulation Domain */
    Dim.x = 10;
    Dim.y = 10;
    Cell_Size.x = Dim.x / NUM_CELLS_X;
    Cell_Size.y = Dim.y / NUM_CELLS_Y;
    Particle   *atom;
    atom = (Particle *)malloc(PARTICLES* sizeof(Particle));

    /* Initialize results */
    Vec_2 velocity_field[NUM_CELLS_X][NUM_CELLS_Y] = {0};


    /* Initialize particles */
    initialize_positions(atom, Dim);
    initialize_velocities(atom);
    
    /* Initialize DSMC simulation */
    Cell cells[NUM_CELLS_X][NUM_CELLS_Y];
    dt = 1e-4;
    //initialize_particles(particles, NUM_PARTICLES, BOX_SIZE);


    /* Run simulation */
    for (int t = 0; t < NUM_TIME_STEPS; t++) {
        /* DSMC Simulation */
        drift_particles(atom, dt);
        apply_BC(atom, Dim);
        assign_particles_to_cells(atom, cells, Cell_Size);
        collide_particles_in_cells(atom, cells, Cell_Size, dt);
        average_flow(atom, cells, velocity_field);
        fprintf(stdout, "Time step: %d\n", t);
        for (int i = 0; i < (NUM_CELLS_X - 1); i++) {
            for (int j = 0; j < (NUM_CELLS_Y - 1); j++) {
                if ((fabs(velocity_field[i][j].x - velocity_field[i+1][j+1].x) > 5e2) || (fabs(velocity_field[i][j].y - velocity_field[i+1][j+1].y) > 5e2)) {
                    //Gradients are very high -> employ MD for 10 time steps
                    fprintf(stdout, "Gradients are very high -> employ MD for 10 time steps\n");
                    fprintf(stdout, "Velocity field: %f %f\n", velocity_field[i][j].x, velocity_field[i][j].y);
                    //Avoid high-energy bursts
                    uplift_particles(cells, i, j);
                    //Run MD for 10 time steps
                    for (int k = 0; k < 10; k++) {
                        update_particles(cells, i, j, dt/10, Dim);
                        //Print results
                        for (int w = 0; w < PARTICLES; w++) {
                            fprintf(f_dat , "%e %e\n", (atom + w)->r.x, (atom + w)->r.y); fflush(f_dat);
                        }
                    }
                }
                // Obtain energy in flow field
                egy_ave = sqrt(velocity_field[i][j].x * velocity_field[i][j].x + velocity_field[i][j].y * velocity_field[i][j].y);
                //Print results
                fprintf(egy_dat, "%e %e\n", t*dt, egy_ave); fflush(egy_dat); egy_ave = 0;
            }
        }
        for (int i = 0; i < PARTICLES; i++) {
            fprintf(f_dat , "%e %e\n", (atom + i)->r.x, (atom + i)->r.y); fflush(f_dat);
        }
    }

    /* Print results */
    printf("Average velocities: \n");
    for (int j = 0; j < NUM_CELLS_Y; j++) {
        for (int i = 0; i < NUM_CELLS_X; i++) {
            if (cells[i][j].counter == 0){
                printf("[ ---.-- ]\t");
            }
            else {
                printf("[ %.2f ]\t", sqrt(velocity_field[i][j].x * velocity_field[i][j].x + velocity_field[i][j].y * velocity_field[i][j].y));
            }
        }
        printf("\n\n");
    }

    fclose(f_dat);
    free(atom);
    return 0;
}