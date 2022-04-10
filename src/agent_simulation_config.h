#define VERSION 1.0

//Technical constants for the rng and the GSL minimizer
#define SEED 1;//time(NULL);
#define STEP_SIZE 0.1; //Initial step size of the GSL minimizer
#define LINE_TOL 1e-4; //Error tolerance for the line minimisation in a given direction
#define TOL 1e-5; //Error tolerance for the overall minimisation
#define MAX_ITER 100; //Maximum number of iterations for the GSL minimiser


