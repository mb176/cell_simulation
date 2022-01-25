#include <stdio.h>
#include <time.h>
#include <stdint.h>
#include "src/cell_simulation.h"
#include "src/base_types2D.h"
#include "src/global.h"
#include "src/xoshiro_rng.h"


/*
For Tuesday
- get python animation to run it
- fix initial conditions
- fix periodic boundaries
In general
- efficient interaction check
- do input files 
- do parallelisation
- do more measurements
- do intermediate write steps
- test integrator for certain analytical results, momementum/ energy conservation?
- unit tests
- optimize log and sin calls in box-mueller
- do updates to stdout
*/

// Global parameter
int stepLimit; //defines the duration of the simulation
int nParticles;
int nGreenParticles;
int nRedParticles;
real stepDuration;
VecR region;  //Size of the simulation region
real D; //translational diffusion constant
real redRotD, greenRotD, greenRotDPlus; //rotational diffusion constants 
real k;  //potential strength
real tau; //decay time of the persistent state
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
    SetParameters();
    SetUpJob();
    for(int stepIdx = 0; stepIdx < stepLimit; stepIdx++){
        SingleStep(stepIdx);
    }
    writeTracks();
    cleanup();
}

