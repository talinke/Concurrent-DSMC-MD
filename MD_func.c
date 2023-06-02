#include "DSMC-MD.h"

/***********************************************************************************************
 ************************* Molecular dynamics simulation functions *****************************
 ***********************************************************************************************
 * Author: Tim Linke (talinke@ucdavis.edu)
 * Date:   2023-06-02
 **********************************************************************************************/

/*----------------------- Initialize positions of particles ------------------------------------
 * Initializes the positions of each particles using equidistant spacing across the domain. 
 * Force is initialized to zero, and individual masses are defined.
 * -------------------------------------------------------------------------------------------*/
void initialize_positions(Particle* atom, Vec_2 Dim) {
    int k, q = 1;
    Vec_2 dr, offset;
    dr.x = Dim.x/PARTICLES_X;
    dr.y = Dim.y/PARTICLES_Y;
    offset.y = 0.5*dr.y;
    for (int j = 0; j < PARTICLES_Y; j++) {
        offset.x = 0.25*dr.x*(2.0+q);
        q *= -1;
        for (int i = 0; i < PARTICLES_X; i++) {
            k = j*PARTICLES_X+i;
            (atom + k)->r.x   = offset.x + i*dr.x;
            (atom + k)->r.y   = offset.y + j*dr.y;
            (atom + k)->v.x   = 0.0;
            (atom + k)->v.y   = 0.0;
            (atom + k)->f.x   = 0.0;
            (atom + k)->f.y   = 0.0;
            (atom + k)->mass  = 1.0;
        }
    }
}

/*----------------------- Initialize velocities of particles -----------------------------------
 * Initializes the velocities of each particles using a Maxwell-Boltzmann distribution. 
 * -------------------------------------------------------------------------------------------*/
void initialize_velocities(Particle* atom) {
    double x1, x2, v_avg = 0.0;
    double m = MASS;
    double k_B = BOLTZMANN;
    double T = 273.0;
    srand(time(NULL)); // Seed the random number generator

    for (int j = 0; j < PARTICLES; j++) {

        x1 = (double)rand()/(double)RAND_MAX;
        x2 = (double)rand()/(double)RAND_MAX;

        // Box-Mueller Method
        (atom+j)->v.x   = sqrt(-2.0 * log(1.0-x1)) * cos(2.0 * M_PI * x2) * sqrt(k_B*T/m);
        (atom+j)->v.y   = sqrt(-2.0 * log(1.0-x1)) * sin(2.0 * M_PI * x2) * sqrt(k_B*T/m);

        // calculate the average velocity of the domain
        v_avg += sqrt((atom+j)->v.x*(atom+j)->v.x + (atom+j)->v.y*(atom+j)->v.y)/PARTICLES;
    }
    printf("Domain average velocity = %e\n", v_avg);
}


/*----------------------- Update particle positions --------------------------------------------
 * Updates the positions and velocity of each particle using the Verlet algorithm. First, one
 * obtains the acceleration from the current force, then the position is updated using the
 * current velocity and acceleration. The force is then set to zero for the next iteration.
 * The force (which depends on position) is now reevaluated, and the velocity is updated using
 * the average of the current and previous acceleration. The force is obtained using the 
 * Lennard-Jones potential.
 * -------------------------------------------------------------------------------------------*/
void update_particles(Cell cells[][NUM_CELLS_Y], int index_x, int index_y, double dt, Vec_2 Dim){
    Vec_2 dr, frc, a;
    int particles;
    double egy_ave = 0;

    // Get group of particles
    particles = cells[index_x][index_y].counter;


    for (int i=0; i < particles; i++) {
        // Use current force to get acceleration
        a.x = (cells[index_x][index_y].cell_atoms[i])->f.x / (cells[index_x][index_y].cell_atoms[i])->mass;
        a.y = (cells[index_x][index_y].cell_atoms[i])->f.y / (cells[index_x][index_y].cell_atoms[i])->mass;

        // Update particle positions
        (cells[index_x][index_y].cell_atoms[i])->r.x += (cells[index_x][index_y].cell_atoms[i])->v.x * dt + 0.5 * a.x * dt * dt;
        (cells[index_x][index_y].cell_atoms[i])->r.y += (cells[index_x][index_y].cell_atoms[i])->v.y * dt + 0.5 * a.y * dt * dt;

        // Null force for next iteration
        (cells[index_x][index_y].cell_atoms[i])->f.x = 0.0;
        (cells[index_x][index_y].cell_atoms[i])->f.y = 0.0;
    }

    // Calculate the force on each particle
    for (int i=0; i < particles; i++) {
        for (int j=0; j<i; j++) {
            if (i != j) {
                // Calculate the distance between the particles
                dr.x = (cells[index_x][index_y].cell_atoms[j])->r.x-(cells[index_x][index_y].cell_atoms[i])->r.x;
                dr.y = (cells[index_x][index_y].cell_atoms[j])->r.y-(cells[index_x][index_y].cell_atoms[i])->r.y;
                if      (dr.x >= 0.5*Dim.x) dr.x -= Dim.x;
                else if (dr.x < -0.5*Dim.x) dr.x += Dim.x;
                if      (dr.y >= 0.5*Dim.y) dr.y -= Dim.y;
                else if (dr.y < -0.5*Dim.y) dr.y += Dim.y;

                // Calculate the force using LJ potential
                frc = force_LJ(dr);

                (cells[index_x][index_y].cell_atoms[i])->f.x -= frc.x;
                (cells[index_x][index_y].cell_atoms[i])->f.y -= frc.y;
                (cells[index_x][index_y].cell_atoms[j])->f.x += frc.x;
                (cells[index_x][index_y].cell_atoms[j])->f.y += frc.y;

                // Calculate the average energy of the system
                egy_ave += energy_LJ(dr);
            }
        }
    }

    // Update particle velocities using new forces
    for (int i=0; i < particles; i++) {
        (cells[index_x][index_y].cell_atoms[i])->v.x += 0.5 * ((cells[index_x][index_y].cell_atoms[i])->f.x / (cells[index_x][index_y].cell_atoms[i])->mass + a.x) * dt;
        (cells[index_x][index_y].cell_atoms[i])->v.y += 0.5 * ((cells[index_x][index_y].cell_atoms[i])->f.y / (cells[index_x][index_y].cell_atoms[i])->mass + a.y) * dt;
    }
}

/*----------------------- Calculate the force using the Lennard-Jones potential ---------------
 * Calculates the force between two particles using the Lennard-Jones potential.
 * -------------------------------------------------------------------------------------------*/
Vec_2 force_LJ(Vec_2 dr) {
   double r2, r2i, r4i, r6i, fac;
   Vec_2  frc;

   r2 = dr.x*dr.x + dr.y*dr.y;

   r2i = 1.0/r2;
   r4i = r2i*r2i;
   r6i = r2i*r4i;

   fac = 12.0*r6i*r2i*(r6i-1.0);

   frc.x = fac*dr.x;
   frc.y = fac*dr.y;

   return frc;
}

/*----------------------- Calculate the energy using the Lennard-Jones potential ---------------
 * Calculates the energy between two particles using the Lennard-Jones potential.
 * -------------------------------------------------------------------------------------------*/
double energy_LJ(Vec_2 dr) {
   double r2, r2i, r4i, r6i;
   double egy;

   r2 = dr.x*dr.x + dr.y*dr.y;

   r2i = 1.0/r2;
   r4i = r2i*r2i;
   r6i = r2i*r4i;

   egy = r6i*(r6i-2.0);

   return egy;
}