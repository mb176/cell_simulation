#define VERSION 1.0

//////////////////////////////////////////////////////////// Modifications to the model ////////////////////////////////////////////////////

//If defined CIL only happens in heterotypic contacts, otherwise in all contacts (subject to turnAround being >0)
#define DIFFERENTIAL_CIL

//How long is CIL delayed after contact? (-1) corresponds to no contact delay
#define CIL_DELAY -1.0 //Inspired from experiments we set the duration of contact to be half the duration of increased persistence

// Should particles be connected via harmonic springs while CIL is delayed?
// #define STICKY_CONTACTS 

// #define turnAroundVariation M_PI; //if defined the angle after turnAround will be randomised by a uniformly drawn angle with maximal value turnAroundVariation in radians

// #define CIL_COOLDOWN_DURATION 0.02 // If defined CIL has a (resetting) cooldown (persistence change still possible)

// #define NON_DIFFERENTIAL_PERSISTENCE

//////////////////////////////////////////////////////////// Alternative initialisations ////////////////////////////////////////////////////

// Normally the particles get randomly placed into the box, and then the minimisation algorithm is used to minimise overlaps

// #define INITIAL_BLOB // When set the system will be initialised with all green particles in a blob in the center

// #define ONLY_GREEN_CIL // When set, only green particles will change their direction upon contact 

#define INITIAL_PHASE_SEGREGATION // When defined will initialise in quadratic packing with green on one side and red on the other; will ignore areaFraction

//////////////////////////////////////////////////////////// Numerical configurations ////////////////////////////////////////////////////

//Technical constants for the rng and the GSL minimizer
#define SEED time(NULL); 
#define STEP_SIZE 1e-5; //Initial step size of the GSL minimizer
#define LINE_TOL 1e-4; //Error tolerance for the line minimisation in a given direction
#define TOL 1e-5; //Error tolerance for the overall minimisation
#define MAX_ITER 8000; //Maximum number of iterations for the GSL minimiser
// Maximum number of neighbourshood pairs that can be stored; Simulation crashes if the value is exceeded 
#define MAX_NEIGHBOUR_PAIRS 20*nParticles
// Extends the size of computational cells beyond the interaction radius. When particle displacement exceeds this value, the neighbourhood list needs updating
#define CELL_SIZE_EXTENSION 1

//////////////////////////////////////////////////////////// Optional outputs //////////////////////////////////////////////////////////

//Enable DEBUG for more output/ extra assert cases
// #define DEBUG
// #define MAX_STEP_DISPLACEMENT 1.0 //In debug mode we assert that a particle cannot travel further than that in a single step 

// Measure angle between cells after heterotypic collision
#define MEASURE_COLLISION_ANGLE

// Write down the particle orientations at each measurement step
#define TRACK_VELOCITIES

