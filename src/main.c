#include <stdio.h>
#include <time.h>
#include <stdint.h>
#include <string.h>
#include "global.h"
#include "cell_simulation.h"
#include "base_types2D.h"
#include "agent_simulation_config.h"
// #include "src/xoshiro_rng.h"


/*
In general
- parallelise with openMP
- write positions at each step, see where it LJ goes wrong
- get HPC to work
- rename/ restructure python files so that they only do one thing at a time
- unit tests
- save non-periodic positions
- do lennart jones with distance2 instead of distance
- have potential named in logfile
*/


#ifndef LENGTH
#define LENGTH (100)
#endif
int length=LENGTH;
#undef LENGTH
#define LENGTH length 


// Global parameter
int stepLimit; //defines the duration of the simulation
int nParticles;
int nGreenParticles;
int nRedParticles;
real stepDuration;
VecR region;  //Size of the simulation region
real redD, greenD, greenPlusD; //rotational diffusion constants 
real k;  //potential strength
real tau; //decay time of the persistent state
real Pe; //self-propulsion velocity
real sigma; //interaction distance, sigma>1 creates stickyness
real measurementInterval; //at the end of each interval we measure all observables
int nMeasurements;
VecI initialCellsGreen; //2D vector containing the number of particles in the x and y direction
VecI initialCellsRed;
int LennardJones; //bool that says whether if we use Lennard Jones potential instead of harmonic interactions 
// Data structures
uint64_t seed;
struct xoshiro256ss_state rng;
particle * particles;
VecR ** positionMeasurements; 
real * measurementTimes;
int ** colorMeasurements;
FILE * tracksFile;
FILE * paramFile;
VecI cells; //Number of cells dividing the space in each coordinate for interaction computation
int * cellList; //"Linked list" of length nParticle+nCells, first entries contain 
//indices of next particle in the cell, the last entries contain indices to first 
//particle in that cell

int main(int argc, char **argv){
    
/* Runs the simulation in a rectangular region [0,region.x]x[0,region.y] with periodic boundary conditions. The agents are refered
to as particles, "cells" only refers to the unit cells for some partition of the region. 
In the array "particles", we first have all the green particles and then all the red particles.
The color scheme for measurements is red=0, green=1, green+ = 2;
*/  

    clock_t t;
    t = clock();
    
    
    SetParameters(argv);
    SetUpJob();
    printf("Begin Simulation... \n");
    fprintf(paramFile,"Begin Simulation... \n");
    for(int stepIdx = 0; stepIdx < stepLimit; stepIdx++){
        SingleStep(stepIdx);
        if(stepIdx % (stepLimit/20)==0){ //Progressbar, one line every 5%
            printf("|");
            fflush(stdout); //Flushes stdout
            fprintf(paramFile,"|");
            fflush(paramFile); //Flushes the buffer
        }
    }
    writeTracks();

    t = clock() - t;
    double time_taken = ((double)t)/CLOCKS_PER_SEC;
    printf("\nSimulation ended after %f seconds \n",time_taken);
    fprintf(paramFile,"\nSimulation ended after %f seconds \n",time_taken);

    cleanup(); 
}

