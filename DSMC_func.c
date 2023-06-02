#include "DSMC-MD.h"

/***********************************************************************************************
 ***************** Direct Simulation Monte Carlo (DSMC) simulation functions *******************
 ***********************************************************************************************
 * Author: Tim Linke (talinke@ucdavis.edu)
 * Date:   2023-06-02
 **********************************************************************************************/


/*------------------------------- Uplift particles --------------------------------------------
 * Check if particles from DSMC simulation result in high-energy collisions. If so, reposition
 * the particles iteratively to prepare for MD simulation.
 * -------------------------------------------------------------------------------------------*/
void uplift_particles(Cell cells[][NUM_CELLS_Y], int index_x, int index_y) {
    Vec_2 dr;
    double distance;
    int particles = cells[index_x][index_y].counter;

    for (int i = 1; i < particles; i++) {
        for (int j = 0; j < i; j++) {
            if (i != j) {
                dr.x = (cells[index_x][index_y].cell_atoms[j])->r.x-(cells[index_x][index_y].cell_atoms[i])->r.x;
                dr.y = (cells[index_x][index_y].cell_atoms[j])->r.y-(cells[index_x][index_y].cell_atoms[i])->r.y;
                distance = sqrt(dr.x*dr.x + dr.y*dr.y);

                if (distance < CUTOFF) {
                    while (distance < CUTOFF) {
                        // Update particle positions randomly
                        // -- Regardless of whether dr.x & dr.y are negative or positive, this update will put distance between them
                        (cells[index_x][index_y].cell_atoms[i])->r.x -= ((double)rand()/RAND_MAX) * dr.x;
                        (cells[index_x][index_y].cell_atoms[i])->r.y -= ((double)rand()/RAND_MAX) * dr.y;
                        
                        dr.x = (cells[index_x][index_y].cell_atoms[j])->r.x-(cells[index_x][index_y].cell_atoms[i])->r.x;
                        dr.y = (cells[index_x][index_y].cell_atoms[j])->r.y-(cells[index_x][index_y].cell_atoms[i])->r.y;
                        distance = sqrt(dr.x*dr.x + dr.y*dr.y);
                    }
                    // Reset checks
                    printf("----> Restarted positioning after atoms (%d, %d)\n", i, j);
                    i = 0;
                    j = 0;
                }
            }
        }
        
    }
}


/*------------------------------- Drift particles ----------------------------------------------
 * Update particle position using Newtonian motion.
 *--------------------------------------------------------------------------------------------*/
void drift_particles(Particle* atom, double dt) {
    // Free motion of particles by a time dt
    for (int i = 0; i < PARTICLES; i++) {
        (atom + i)->r.x += (atom + i)->v.x * dt;
        (atom + i)->r.y += (atom + i)->v.y * dt;
    }
}

/*------------------------------- Apply boundary conditions -------------------------------------
 * Apply periodic boundary conditions to particles.
 *--------------------------------------------------------------------------------------------*/
void apply_BC(Particle* atom, Vec_2 Dim) {
    for (int i = 0; i < PARTICLES; i++) {
        if ((atom + i)->r.x < 0.0) {
            (atom + i)->r.x += Dim.x;
        }
        if ((atom + i)->r.x > Dim.x) {
            (atom + i)->r.x -= Dim.x;
        }
        if ((atom + i)->r.y < 0.0) {
            (atom + i)->r.y += Dim.y;
        }
        if ((atom + i)->r.y > Dim.y) {
            (atom + i)->r.y -= Dim.y;
        }
    }
}

/*------------------------------- Initialize cells --------------------------------------------
 * To run DSMC, the particles must be grouped into cells. The struct CELLS allows the creation
 * of a 2D grid. Each cell contains a counter that keeps track of the number of particles in
 * a cell and allows new cells to be appended. Each cell also contains an array of pointers 
 * to every particle found in the cell. Both the counter and the array are initialized to zero.
 *--------------------------------------------------------------------------------------------*/
void initialize_cells(Cell cells[][NUM_CELLS_Y]) {
    // Initialize all cell values to zero
    for (int i = 0; i < NUM_CELLS_X; i++) {
        for (int j = 0; j < NUM_CELLS_Y; j++) {
            cells[i][j].counter = 0;
            for (int k = 0; k < 4*(int)(PARTICLES/NUM_CELLS); k++) {
                cells[i][j].cell_atoms[k] = NULL;
            }
        }
    }
}

/*------------------------------- Assign particles to cells -----------------------------------
 * This function assigns each particle to a cell. Based on its position, the particle's
 * cell indices are obtained and the particle is added to the corresponding cell. To assure that
 * the particle is appended to the list of existing particles in the cell, the cell counter is
 * used to determine the next available position in the array of pointers to particles. The
 * counter is then incremented to reflect the addition of a new particle to the cell.
 *--------------------------------------------------------------------------------------------*/
void assign_particles_to_cells(Particle* atom, Cell cells[][NUM_CELLS_Y], Vec_2 Cell_Size) {
    int cell_x, cell_y;

    // Clear all cell values
    initialize_cells(cells);
    
    // Assign each particle to a cell
    for (int i = 0; i < PARTICLES; i++) {
        // Determine the cell indices of the particle
        cell_x = (atom + i)->r.x / Cell_Size.x;
        cell_y = (atom + i)->r.y / Cell_Size.y;

        // Check if the cell indices are within bounds
        if (cell_x >= 0 && cell_x < NUM_CELLS_X && cell_y >= 0 && cell_y < NUM_CELLS_Y) {
            // Add particle to corresponding cell and increment the particle count in the cell
            cells[cell_x][cell_y].cell_atoms[cells[cell_x][cell_y].counter] = (atom + i);
            cells[cell_x][cell_y].counter++;
        }
    }
}

/*------------------------------- Collide the particles ----------------------------------------
 * The particles are now collided according to the DSMC method. For each cell, the number of
 * particles in the cell is obtained. Based on the number of particles, the number of collision 
 * candidates is calculated. For the number of candidates, two particles are randomly selected
 * from the cell and the relative velocity is calculated. If the relative velocity passes the
 * acceptance-rejection scheme based on a capped relative velocity, the particles are collided
 * using elastic collisions. The position of both collided particles is then updated with their
 * new velocities.
 *--------------------------------------------------------------------------------------------*/
void collide_particles_in_cells(Particle* atom, Cell cells[][NUM_CELLS_Y], Vec_2 Cell_Size, double dt) {
    double vel_rel;
    int num_collisions = 0, cell_x, cell_y, M_cand, num_particles;
    double v_rel_max = 6.0; // Cap relative velocity
    double r_fac, v_rel, vx_cm, vz_cm, cos_theta, sin_theta, phi, vx_p, vz_p;
    int i_prop, j_prop;
    
    // Iterate over each cell
    for (int i = 0; i < NUM_CELLS_X; i++) {
        for (int j = 0; j < NUM_CELLS_Y; j++) {
            // Get the number of particles in the current cell
            num_particles = cells[i][j].counter;
            
            // If there is only one particle, skip this cell
            if (num_particles <= 1) {
                continue;
            }
            
            // Get the number of candidate collisions
            M_cand = (int)ceil(num_particles * num_particles * M_PI * v_rel_max * SIMPAR * dt / (2 * Cell_Size.x * Cell_Size.y));

            for (int k = 0; k < M_cand; k++) {
                r_fac = (double)rand() / RAND_MAX;
                // Select two particles at random
                i_prop = rand() % num_particles;
                j_prop = rand() % num_particles;

                v_rel = sqrt( pow(cells[i][j].cell_atoms[i_prop]->v.x - cells[i][j].cell_atoms[j_prop]->v.x, 2) + pow(cells[i][j].cell_atoms[i_prop]->v.y - cells[i][j].cell_atoms[j_prop]->v.y, 2));

                // Acceptance-Rejection Scheme
                if (v_rel > r_fac * v_rel_max) {
                    // Collide particles elasticly (TODO: double check this)
                    vx_cm = 0.5 * (cells[i][j].cell_atoms[i_prop]->v.x + cells[i][j].cell_atoms[j_prop]->v.x);
                    vz_cm = 0.5 * (cells[i][j].cell_atoms[i_prop]->v.y + cells[i][j].cell_atoms[j_prop]->v.y);
                    cos_theta = 2 * ((double)rand() / RAND_MAX) - 1;
                    sin_theta = sqrt(1 - pow(cos_theta, 2));
                    phi = 2 * M_PI * ((double)rand() / RAND_MAX);
                    vx_p = v_rel * sin_theta * cos(phi);
                    vz_p = v_rel * cos_theta;
                    cells[i][j].cell_atoms[i_prop]->v.x = vx_cm + 0.5 * vx_p;
                    cells[i][j].cell_atoms[i_prop]->v.y = vz_cm + 0.5 * vz_p;
                    cells[i][j].cell_atoms[j_prop]->v.x = vx_cm - 0.5 * vx_p;
                    cells[i][j].cell_atoms[j_prop]->v.y = vz_cm - 0.5 * vz_p;

                    num_collisions++;
                }
            }
        }
    }
}

/*------------------------------- Average the flow properties ----------------------------------
 * For each cell in the grid, the average flow velocity is computed. The obtained value is then
 * passed to the array velocity_field to reflect the resulting flow field in the entire domain.
 *--------------------------------------------------------------------------------------------*/
void average_flow(Particle* atom, Cell cells[][NUM_CELLS_Y], Vec_2 velocity_field[][NUM_CELLS_Y]) {
    Vec_2 v;
    
    // Iterate over each particle
    for (int i = 0; i < NUM_CELLS_X; i++) {
        for (int j = 0; j < NUM_CELLS_Y; j++) {
            // Reset the flow velocity
            v.x = 0.0;
            v.y = 0.0;

            for (int k = 0; k < cells[i][j].counter; k++) {
                // Add the particle velocity to the flow velocity
                v.x += cells[i][j].cell_atoms[k]->v.x;
                v.y += cells[i][j].cell_atoms[k]->v.y;
            }
            // Compute the average flow velocity
            v.x /= cells[i][j].counter;
            v.y /= cells[i][j].counter;
            velocity_field[i][j] = v;
        }
    }
}
