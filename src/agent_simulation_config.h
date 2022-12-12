#define VERSION 1.0

//Technical constants for the rng and the GSL minimizer
#define SEED time(NULL);
#define STEP_SIZE 0.1; //Initial step size of the GSL minimizer
#define LINE_TOL 1e-4; //Error tolerance for the line minimisation in a given direction
#define TOL 1e-5; //Error tolerance for the overall minimisation
#define MAX_ITER 8000; //Maximum number of iterations for the GSL minimiser


//Enable DEBUG for more output/ extra assert cases
#define DEBUG
#define MAX_STEP_DISPLACEMENT 1.0 //In debug mode we assert that a particle cannot travel further than that in a single step 

// Measure angle between cells after heterotypic collision
#define MEASURE_COLLISION_ANGLE

// If defined will cause particles to be attached by a harmonic spring to their last contact
#define STICKY_CONTACTS 




// #define turnAroundVariation M_PI/2; //if defined the angle after turnAround will be randomised by a uniformly drawn angle with maximal value turnAroundVariation in radians
