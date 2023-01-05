#define VERSION 1.0

//Technical constants for the rng and the GSL minimizer
#define SEED 1672651852; //time(NULL);
#define STEP_SIZE 0.1; //Initial step size of the GSL minimizer
#define LINE_TOL 1e-4; //Error tolerance for the line minimisation in a given direction
#define TOL 1e-5; //Error tolerance for the overall minimisation
#define MAX_ITER 8000; //Maximum number of iterations for the GSL minimiser


//Enable DEBUG for more output/ extra assert cases
#define DEBUG
#define MAX_STEP_DISPLACEMENT 1.0 //In debug mode we assert that a particle cannot travel further than that in a single step 

// Measure angle between cells after heterotypic collision
#define MEASURE_COLLISION_ANGLE

// Write down the particle orientations at each measurement step
#define TRACK_VELOCITIES

//How long is CIL delayed after contact? (-1) corresponds to no contact delay
#define CIL_DELAY -1 //Inspired from experiments we set the duration of contact to be half the duration of increased persistence

// Should particles be connected via harmonic springs while CIL is delayed?
// #define STICKY_CONTACTS 

// Maximum number of neighbourshood pairs that can be stored; Simulation crashes if the value is exceeded 
#define MAX_NEIGHBOUR_PAIRS 10*nParticles
// Extends the size of computational cells beyond the interaction radius. When particle displacement exceeds this value, the neighbourhood list needs updating
#define CELL_SIZE_EXTENSION 1



// #define turnAroundVariation M_PI/2; //if defined the angle after turnAround will be randomised by a uniformly drawn angle with maximal value turnAroundVariation in radians
