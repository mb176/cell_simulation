#include <stdio.h>
#include <time.h>
#include <stdint.h>
#include <string.h>
#include "src/global.h"
#include "src/cell_simulation.h"
#include "src/base_types2D.h"
#include "src/xoshiro_rng.h"

/*
In general
- restructure the data analysis files to all be in one folder
- efficient interaction check
- do input files 
- do parallelisation
- do more measurements
- do intermediate write steps
- test integrator for certain analytical results, momementum/ energy conservation?
- unit tests
- optimize log and sin calls in box-mueller
- do updates to stdout
- save non-periodic positions
- random initial conditions -> steepest decent to minimize overlaps
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
// Data structures
uint64_t seed;
struct xoshiro256ss_state rng;
particle * particles;
VecR ** positionMeasurements; 
real * measurementTimes;
int ** colorMeasurements;
FILE * tracksFile;

int main(int argc, char **argv){
    
/* Runs the simulation in a rectangular region [0,region.x]x[0,region.y] with periodic boundary conditions. The agents are refered
to as particles, "cells" only refers to the unit cells for some partition of the region. 
In the array "particles", we first have all the green particles and then all the red particles.
The color scheme for measurements is red=0, green=1, green+ = 2;
*/ 
    clock_t t;
    t = clock();
    printf("Begin Simulation... \n");
    SetParameters(argc, argv);
    SetUpJob();
    for(int stepIdx = 0; stepIdx < stepLimit; stepIdx++){
        SingleStep(stepIdx);
    }
    writeTracks();
    cleanup();
    t = clock() - t;
    double time_taken = ((double)t)/CLOCKS_PER_SEC;
    printf("Simulation ended after %f seconds \n",time_taken);
}

