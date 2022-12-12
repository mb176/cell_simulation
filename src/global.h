#ifndef global_variables
#define global_variables

#define MAX_LINE_LENGTH 200
#define MAX_STRING_LENGTH 50
#include "base_types2D.h"
#include "agent_simulation_config.h"

//Simulation Parameters
extern int stepLimit;
extern real stepDuration;
extern int skipSteps; //only after these steps will tracks be recorded (saves space)
extern real measurementInterval; //at the end of each interval we measure all observables
extern int nParticles;
extern int nGreenParticles;
extern int nRedParticles;
extern VecR region;  //Size of the simulation region
extern real redD, greenD, greenPersistentD; //rotational diffusion constants 
extern real k;  //potential strength
extern real tau; //decay time of the persistent state
extern real Pe; //self-propulsion velocity
extern real potentialRange; //interaction distance, sigma>1s creates stickyness
extern int nMeasurements;
extern VecI initialCellsGreen; //2D vector containing the number of particles in the x and y direction
extern VecI initialCellsRed;
extern int LennardJones; //bool that says whether if we use Lennard Jones potential instead of harmonic interactions 
extern real turnAround; // Bool that flags wether the particles move away from each other after colision
extern real redRedAdhesionMult; // Multiplies the adhesive component of the force between these particle types
extern real greenGreenAdhesionMutl;
extern real redGreenAdhesionMult;
extern real CIL_delay; // the delay after a contact until CIL kicks in. Also limits the time of the pairing/ harmonic spring

// Data structures
extern uint64_t seed;
extern struct xoshiro256ss_state rng;
extern particle * particles; //First half is green, second half red
extern FILE * tracksFile;
extern FILE * paramFile;
extern VecI cells; //Number of cells dividing the space in each coordinate for interaction computation
extern int * cellList; //"Linked list" of length nParticle+nCells, first entries contain 
//indices of next particle in the cell, the last entries contain indices to first 
//particle in that cell

// Extra observables
extern real simulationTime; //So that any function can access the current time of the simulation
#ifdef MEASURE_COLLISION_ANGLE
extern real collisionAngle;
extern int nCollisions;
extern real collisionDuration;
#endif //MEASURE_COLLISION_ANGLE


#endif