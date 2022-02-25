#ifndef global_variables
#define global_variables

#define MAX_LINE_LENGTH 200
#define MAX_STRING_LENGTH 50
#include "base_types2D.h"

//Parameters
extern int stepLimit;
extern int nParticles;
extern int nGreenParticles;
extern int nRedParticles;
extern real stepDuration;
extern VecR region;  //Size of the simulation region
extern real redD, greenD, greenPlusD; //rotational diffusion constants 
extern real k;  //potential strength
extern real tau; //decay time of the persistent state
extern real Pe; //self-propulsion velocity
extern real sigma; //interaction distance, sigma>1s creates stickyness
extern real measurementInterval; //at the end of each interval we measure all observables
extern int nMeasurements;
extern VecI initialCellsGreen; //2D vector containing the number of particles in the x and y direction
extern VecI initialCellsRed;

// Data structures
extern uint64_t seed;
extern struct xoshiro256ss_state rng;
extern particle * particles; //First half is green, second half red
extern VecR ** positionMeasurements; 
extern real * measurementTimes;
extern int ** colorMeasurements;
extern FILE * tracksFile;

#endif