#include <stdio.h>
#include <time.h>
#include <stdint.h>
#include <string.h>
#include "global.h"
#include "cell_simulation.h"
#include "base_types2D.h"
#include "agent_simulation_config.h"
// #include "src/xoshiro_rng.h"


#ifndef LENGTH
#define LENGTH (100)
#endif
int length=LENGTH;
#undef LENGTH
#define LENGTH length 


// Global parameter
int stepLimit;
real stepDuration;
int skipSteps; //only after these steps will tracks be recorded (saves space)
real measurementInterval; //at the end of each interval we measure all observables
int nParticles;
int nGreenParticles;
int nRedParticles;
VecR region;  //Size of the simulation region
real redD, greenD, greenPersistentD; //rotational diffusion constants 
real k;  //potential strength
real tau; //decay time of the persistent state
real Pe; //self-propulsion velocity
real potentialRange; //interaction distance, sigma>1s creates stickyness
int nMeasurements;
VecI initialCellsGreen; //2D vector containing the number of particles in the x and y direction
VecI initialCellsRed;
int LennardJones; //bool that says whether if we use Lennard Jones potential instead of harmonic interactions 
real turnAround; // Number that flags wether the particles move away from each other after colision. If smaller than 1 the cells wil lonly turn by a percentage
real redRedAdhesionMult; // Multiplies the adhesive component of the force between these particle types
real greenGreenAdhesionMutl;
real redGreenAdhesionMult;



// Data structures
uint64_t seed;
struct xoshiro256ss_state rng;
particle * particles;
FILE * tracksFile;
#ifdef TRACK_VELOCITIES
FILE * velocityTracksFile;
#endif
FILE * paramFile;
VecI cells; //Number of cells dividing the space in each coordinate for interaction computation
int * cellList; //"Linked list" of length nParticle+nCells, first entries contain 
//indices of next particle in the cell, the last entries contain indices to first 
//particle in that cell

// Extra observables
real simulationTime; //So that any function can access the current time of the simulation
#ifdef MEASURE_COLLISION_ANGLE
real collisionAngle=0;
int nCollisions=0;
real collisionDuration=0; 
#endif //MEASURE_COLLISION_ANGLE

int main(int argc, char **argv){
    
/* Runs the simulation in a rectangular region [0,region.x]x[0,region.y] with periodic boundary conditions. The agents are refered
to as particles, "cells" only refers to the unit cells for some partition of the region. 
In the array "particles", we first have all the green particles and then all the red particles.
The color scheme for measurements is red=0, green=1, green+ = 2;
*/  

    clock_t t;
    t = clock();
    
    SetParameters(argc, argv);
    SetUpJob();
    //Logs
    printf("Begin Simulation... \n");
    fprintf(paramFile,"Begin Simulation... \n");
    #ifdef DEBUG
    printf("DEBUG MODE\n");
    fprintf(paramFile,"DEBUG MODE\n");
    #endif
    #ifdef turnAroundVariation
    double maxAngle = 360/(2*3.1416)*turnAroundVariation; 
    printf("Randomised turn-around directions (maxAngle = +- %f) \n",maxAngle);
    fprintf(paramFile,"Randomised turn-around directions (maxAngle = +- %f) \n",maxAngle);
    #endif
    
    
    for(int stepIdx = 0; stepIdx <= stepLimit; stepIdx++){
        simulationTime = stepIdx*stepDuration;
        //Iterate simulation, save positions at certain step indices
        SingleStep(stepIdx);

        //Progressbar, one line every 5%
        if(stepIdx % (stepLimit/20)==0){ 
            printf("|");
            fflush(stdout); //Flushes stdout
            fprintf(paramFile,"|");
            fflush(paramFile); //Flushes the buffer
        }
    }


    //Simulation duration
    t = clock() - t;
    double time_taken = ((double)t)/CLOCKS_PER_SEC;
    printf("\nSimulation ended after %f seconds \n",time_taken);
    fprintf(paramFile,"\nSimulation ended after %f seconds \n",time_taken);
    #ifdef STICKY_CONTACTS
    fprintf(paramFile, "Pairs are connected by harmonic springs; ");
    #endif
    fprintf(paramFile, "CIL delay time: %f \n", CIL_DELAY);

    #ifdef MEASURE_COLLISION_ANGLE
    fprintf(paramFile, "Average angle after collision event: %f \n", collisionAngle/nCollisions);
    fprintf(paramFile, "Average collision duration: %f \n", collisionDuration/nCollisions);
    fprintf(paramFile, "Number of collisions: %d \n", nCollisions);
    #endif 

    cleanup(); 
}

