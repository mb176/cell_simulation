#include <stdlib.h>  
#include <string.h>  //for strcmp()
#include <stdint.h>  //to use uint64_t
#include <time.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include "global.h"
#include "base_types2D.h"
#include "xoshiro_rng.h"
#include "cell_simulation.h"
#include "agent_simulation_config.h"
#include "stdbool.h"


//Gunnar's malloc macro, gives error if space can't be assigned (e.g. too big -> NULL returned from OS) and keeps taps on the ammount of space used 
long long int total_malloced=0LL;
#define MALLOC(a,n) {if ((a=malloc((n)*sizeof(*a)))==NULL) { fprintf(stderr, "Failed to malloc %i items of size %i (%lli bytes) for variable %s in line %i of %s\n", (int)(n), (int)sizeof(*a), (long long int)((n)*sizeof(*a)), #a, __LINE__, __FILE__); } else { total_malloced+=((long long int)((n)*sizeof(*a)));} }
//The verbose version prints the space used
#define MALLOC_VERBOSE(a,n) {if ((a=malloc((n)*sizeof(*a)))==NULL) { fprintf(stderr, "Failed to malloc %i items of size %i (%lli bytes) for variable %s in line %i of %s\n", (int)(n), (int)sizeof(*a), (long long int)((n)*sizeof(*a)), #a, __LINE__, __FILE__); } else { total_malloced+=((long long int)((n)*sizeof(*a))); printf("# Info: %lli bytes malloced for %s, %lli in total so far.\n", (long long int)((n)*sizeof(*a)), #a, total_malloced); } }


double getNextParameter(FILE * file, char * parameterName){
    char line[MAX_LINE_LENGTH] = {0};
    char name[MAX_STRING_LENGTH];
    char value[MAX_STRING_LENGTH];
    if (fgets(line,MAX_LINE_LENGTH,file)) {
        sscanf(line, "%s : %s", name, value);
    } 
    assert(strcmp(name,parameterName)==0); //Do the names match? If not it's an invalid parameter or parameter out of order or not enough parameters
    char * ptr; //Throughaway, needed to call strtod
    double val = strtod(value,&ptr);
    printf("%s = %f \n",name,val);
    return val;
}

void SetParameters(int argc, char** argv){
    char * PATH;
    char * trackFileSuffix;
    if(argc==2){ //Tracks file name unspecified
        PATH = argv[1];
        trackFileSuffix = "_tracks.csv";
    } else if (argc==3){ //Tracks file name specified
        PATH = argv[1];
        trackFileSuffix = argv[2];
    } else {
        assert(false); // main needs at least one and at most two command line inputs
    }
    
    // Files
    paramFile = fopen(PATH,"r+");
    tracksFile = fopen(strcat(PATH,trackFileSuffix),"w");

    printf("Parameter file: "); printf(PATH); printf("\n");
    printf("Initialise paramters:\n");
    
    //Read out parameters
    stepLimit = getNextParameter(paramFile, "stepLimit");
    stepDuration = getNextParameter(paramFile, "stepDuration");
    skipSteps = getNextParameter(paramFile, "skipSteps"); 
    measurementInterval = getNextParameter(paramFile, "measurementInterval");
    nGreenParticles = getNextParameter(paramFile, "nGreenParticles");
    nRedParticles = getNextParameter(paramFile, "nRedParticles");
    real areaFraction = getNextParameter(paramFile, "areaFraction");
    redD = getNextParameter(paramFile, "redD");
    greenD = getNextParameter(paramFile, "greenD");
    greenPersistentD = getNextParameter(paramFile, "greenPersistentD");
    k = getNextParameter(paramFile, "k");
    tau = getNextParameter(paramFile, "tau");
    Pe = getNextParameter(paramFile, "Pe");
    potentialRange = getNextParameter(paramFile, "potentialRange");
    LennardJones = getNextParameter(paramFile, "LennardJones");
    turnAround = getNextParameter(paramFile, "turnAround");
    redRedAdhesionMult = getNextParameter(paramFile,"redRedAdhesionMult");
    greenGreenAdhesionMutl = getNextParameter(paramFile,"greenGreenAdhesionMutl");
    redGreenAdhesionMult = getNextParameter(paramFile,"redGreenAdhesionMult");

    // Barricades
    //measurementInterval must be a multiple of stepDuration 
    assert((measurementInterval/stepDuration)-round(measurementInterval/stepDuration)< 1e-4);
    //The Simulation duration must be a multiple of measurementInterval
    assert((stepLimit-skipSteps)*(stepDuration/measurementInterval)-round((stepLimit-skipSteps)*(stepDuration/measurementInterval))<1e-4); 
    //Can't skip more steps than we have
    assert(stepLimit>skipSteps);
    //Need it least 20 steps, otherwise the loading bar divides by zero
    assert(stepLimit>=20);

    //Initialise secondary parameters
    nParticles = nGreenParticles + nRedParticles;
    //Setup Square area to fit the area fraction; each particle covers area of size pi*0.5^2 = 0.785398
    real length = sqrt(nParticles*0.785398/areaFraction);
    printf("length=%f \n",length);
    region = (VecR) {.x = length,
                     .y = length};  //Size of the simulation region
    cells.x = region.x/(potentialRange); //Smallest cells that are larger than interaction range
    cells.y = region.y/(potentialRange);

    nMeasurements = round((stepLimit-skipSteps)*(stepDuration/measurementInterval)) + 1; //+1 one because we do the first and the last step 
    seed = SEED;//time(NULL);
    rng = xoshiro256ss_init(seed);
    srand(seed);

    //append the file
    fprintf(paramFile, "nParticles              : %i\n",nParticles);
    fprintf(paramFile, "Length                  : %f\n",length);
    fprintf(paramFile, "Seed                    : %i \n",seed);

    //Set up log file
    fprintf(paramFile, "\n");
    time_t tm;
    time(&tm);
    fprintf(paramFile,"LOGFILE: %s",ctime(&tm));  
}

void AllocArrays(){
    //Particle array
    MALLOC(particles, nParticles); //particles = (particle *) malloc(nParticles*(sizeof(particle)));
    //Linked list for cells
    MALLOC_VERBOSE(cellList, (nParticles + VProd(cells)));//cellList = (int *) malloc((nParticles + VProd(cells))*sizeof(int));
};

//Potential Energy to be minimised
double fHarmonicPotential(const gsl_vector *v, void *params){
    VecR * positions;
    positions = (VecR *) malloc(nParticles*sizeof(VecR));
    for(int particleIdx=0; particleIdx<nParticles; particleIdx++){
        positions[particleIdx].x = gsl_vector_get(v, 2*particleIdx);
        positions[particleIdx].y = gsl_vector_get(v, 2*particleIdx + 1);
    }
    double * para = (double *) params;
    double potentialRange = para[0];

    double V = 0;
    VecR deltaR;
    real distance;
    for(int pIdx1 = 0; pIdx1 < nParticles; pIdx1++){
        for(int pIdx2 = pIdx1+1; pIdx2 < nParticles; pIdx2++){
            VSub(deltaR, positions[pIdx1],positions[pIdx2]);
            VWrapAllTang(deltaR); //Periodic Boundary conditions
            distance = sqrt(VLenSq(deltaR));
            if(distance < potentialRange){
                V += (distance-1)*(distance-1);
            }
        }
    }
    free(positions);
    return V;
}

//Gradient of the potential (For gsl_multimin energy minimisation of initial positons)
void dfHarmonicPotential(const gsl_vector *v, void *params, gsl_vector *df){
    //ToDo: Take care of divergence at 0
    VecR * positions, *gradient;
    positions = (VecR *) malloc(nParticles*sizeof(VecR));
    gradient = (VecR *) malloc(nParticles*sizeof(VecR));
    for(int particleIdx=0; particleIdx<nParticles; particleIdx++){
        positions[particleIdx].x = gsl_vector_get(v, 2*particleIdx);
        positions[particleIdx].y = gsl_vector_get(v, 2*particleIdx + 1);
    }
    double * para = (double *) params;
    double potentialRange = para[0]; //Interaction range
    VecR deltaR;
    real distance;
    //Calculate gradient
    for(int gradientIdx = 0; gradientIdx < nParticles; gradientIdx++){
        VSet(gradient[gradientIdx],0,0);
        for(int particleIdx = 0; particleIdx < nParticles; particleIdx++){
            if(particleIdx != gradientIdx){
                VSub(deltaR, positions[gradientIdx],positions[particleIdx]);
                VWrapAllTang(deltaR); //Periodic Boundary conditions
                distance = sqrt(VLenSq(deltaR));
                if(distance==0){ //Avoid undefined behaviour and just send the particle in some direction
                    VecR direction = {.x = 1, .y = 1};
                    VSAdd(gradient[gradientIdx], gradient[gradientIdx], -2, direction);
                }
                else if(distance < potentialRange){
                    VSAdd(gradient[gradientIdx], gradient[gradientIdx], 2*(1-1/distance), deltaR);
                }
            }
        }
    }
    //Write gradient to gsl vector (For gsl_multimin energy minimisation of initial positons)
    for(int gradientIdx = 0; gradientIdx < nParticles; gradientIdx++){
        gsl_vector_set(df, 2*gradientIdx    , gradient[gradientIdx].x);
        gsl_vector_set(df, 2*gradientIdx + 1, gradient[gradientIdx].y);
    }
    free(positions);
    free(gradient);
}

//Potential and Gradient together (For gsl_multimin energy minimisation of initial positons)
void fdfHarmonicPotential(const gsl_vector *x, void * params, double * f, gsl_vector *df){
    *f = fHarmonicPotential(x, params);
    dfHarmonicPotential(x, params, df);
}

void InitialisePositions(){    
    //Assign random initial positions
    double uniform1, uniform2;
    for(int particleIdx = 0; particleIdx < nParticles; particleIdx++){
        uniform1 = xoshire256ss_uniform(&rng);
        uniform2 = xoshire256ss_uniform(&rng);
        particles[particleIdx].r.x = region.x*uniform1;
        particles[particleIdx].r.y = region.y*uniform2;
    }

    //Minimise Energy of initial positions with GSL conjugate gradient

    //Parameters
    double stepSize = STEP_SIZE; //Initial step size
    double lineTol = LINE_TOL; //Error tolerance for the line minimisation in a given direction
    double tol = TOL; //Error tolerance for the overall minimisation
    int maxIter = MAX_ITER;

    //Initialise function to minimise
    double params[1] = {1.0}; //WARNING: We use a purely repulsive harmonic potential here
    gsl_multimin_function_fdf func;
    func.f      =   &fHarmonicPotential;
    func.df     =  &dfHarmonicPotential;
    func.fdf    = &fdfHarmonicPotential;
    func.n      = 2*nParticles; //Number of variables
    func.params = (void*) params;

    //starting point
    gsl_vector *x;
    x = gsl_vector_alloc(2*nParticles);
    
    for(int particleIdx = 0; particleIdx < nParticles; particleIdx++){
        gsl_vector_set(x,2*particleIdx,     particles[particleIdx].r.x);
        gsl_vector_set(x,2*particleIdx + 1, particles[particleIdx].r.y);
    }

    //Initialise minimizer
    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s;
    T = gsl_multimin_fdfminimizer_conjugate_fr; //Fletcher-Reeves conjugate gradient algorithm
    s = gsl_multimin_fdfminimizer_alloc (T, 2*nParticles); 
    gsl_multimin_fdfminimizer_set (s, &func, x, stepSize, lineTol);

    printf("Initial Potential Energy: %f \n",s->f);
    fprintf(paramFile,"Initial Potential Energy: %f \n",s->f);

    //Minimise potential energy of particle configuration
    size_t iter = 0;
    int status;
    do{
    iter++;
    status = gsl_multimin_fdfminimizer_iterate (s);
    
    if (status)
        break;

    status = gsl_multimin_test_gradient (s->gradient, tol);

    // if (status == GSL_SUCCESS)
    //     printf ("Minimum found at:\n");

    // printf ("%5d (%.3f %.3f) (%.3f %.3f) %10.5f\n", iter,
    //         gsl_vector_get (s->x, 0),
    //         gsl_vector_get (s->x, 1),
    //         gsl_vector_get (s->x, 2),
    //         gsl_vector_get (s->x, 3),
    //         s->f);

    } while (status == GSL_CONTINUE && iter < maxIter);

    printf("Potential energy after %i/%i iterations: %f \n", iter, maxIter, s->f);
    fprintf(paramFile,"Potential energy after %i/%i iterations: %f \n", iter, maxIter, s->f);

    //Copy new starting positions
    VecR pos;
    for(int particleIdx = 0; particleIdx < nParticles; particleIdx++){
        pos.x = gsl_vector_get(s->x,2*particleIdx); 
        pos.y = gsl_vector_get(s->x,2*particleIdx+1); 
        VWrapAll(pos); //Enforce periodic boundary conditions
        particles[particleIdx].r = pos;
    }

    gsl_multimin_fdfminimizer_free (s);
    gsl_vector_free (x);
};

void InitialiseColor(){
    //Green particles
    for(int particleIdx = 0; particleIdx < nGreenParticles; particleIdx++){
        particles[particleIdx].D = greenD;
        particles[particleIdx].decayTimer = 0;
        particles[particleIdx].color = 1;
    }
    //Red particles
    for(int particleIdx = nGreenParticles; particleIdx < nGreenParticles + nRedParticles; particleIdx++){
        particles[particleIdx].D = redD;
        particles[particleIdx].decayTimer = 0;
        particles[particleIdx].color = 0;
    }
};

void InitialiseAngles(){
    //Set theta randomly in [0, 2 pi] and set decay time to 0
    real sum = 0;
    for(int particleIdx=0; particleIdx < nParticles; particleIdx++){
        particles[particleIdx].theta = (float) rand()/(float)(RAND_MAX/(2*M_PI));
        // sum += particles[particleIdx].theta;
    }
    // //set the average of theta to zero
    // for(int particleIdx=0; particleIdx < nParticles; particleIdx++){
    //     particles[particleIdx].theta -= sum/nParticles;
    // }
};


void MeasurePositions(real time){
    //Times
    fprintf(tracksFile, "%f", time);
    for(int particleIdx = 0; particleIdx<nParticles; particleIdx++){
        //Print color
        if(particles[particleIdx].color ==0){
            fprintf(tracksFile,", red");
        } else if (particles[particleIdx].color ==1){
            fprintf(tracksFile,", green");
        } else if (particles[particleIdx].color == 2){
            fprintf(tracksFile,", greenPlus");
        }
        //Print positions
        fprintf(tracksFile,", %f, %f",particles[particleIdx].r.x, \
                                      particles[particleIdx].r.y);
    }
    fprintf(tracksFile,"\n");
};


void SetUpJob(){
    AllocArrays();
    InitialisePositions();
    InitialiseColor();
    InitialiseAngles();
    MeasurePositions(0);
}

//Begin: Support Functions for Single Step
real HarmonicForce(real distance){
    return -k/distance*(distance-1); //WARNING: this needs to be changed if we don't measure in units of potentialRange
}

real LennardJonesForce(real distance){
    //WARNING: needs to be changed if cellRadius!=1
    double invDistance = 1/(distance); 
    double invDistance2 = invDistance*invDistance;
    double invDistance6=invDistance2*invDistance2*invDistance2;
    return 12*k*invDistance6*(invDistance6-1)*invDistance2;
}

real getAngle(VecR r){
    double theta;
    if(r.x>=0){
        theta = atan(r.y/r.x);
    } else {
        theta = atan(r.y/r.x)+M_PI ;
    }
    return theta;
}

void changeDirection(int pIdx1, int pIdx2, VecR deltaR){
    //Particles turn away from each other after contact
    if(turnAround==1){
        real theta = getAngle(deltaR);

        #ifdef turnAroundVariation 
        // //Randomise the turn-around directions
        double randomAngle1 = theta + 2*(xoshire256ss_uniform(&rng)-0.5)*turnAroundVariation;
        double randomAngle2 = theta + M_PI +2*(xoshire256ss_uniform(&rng)-0.5)*turnAroundVariation;
        particles[pIdx1].theta = randomAngle1; //theta + 2*(xoshire256ss_uniform(&rng)-0.5)*turnAroundVariation;
        particles[pIdx2].theta = randomAngle2; //theta + M_PI + 2*(xoshire256ss_uniform(&rng)-0.5)*turnAroundVariation;
        #else
        particles[pIdx1].theta = theta;
        particles[pIdx2].theta = theta + M_PI;
        #endif
    }
    

}

void ComputeInteractions(){
    //Build Linked list
    VecR cellWidth;
    VecI cellIdx;
    int linearCellIdx;
    VDiv(cellWidth, region, cells);

    //Write linked list
    for(int n = nParticles; n < nParticles + VProd(cells); n++) cellList[n] = -1; //Reset list values
    for(int particleIdx = 0; particleIdx < nParticles; particleIdx++){
        VDiv(cellIdx,particles[particleIdx].r,cellWidth);
        linearCellIdx = VLinear(cellIdx, cells) + nParticles;//Linearise matrix index to vector index
        cellList[particleIdx] = cellList[linearCellIdx];
        cellList[linearCellIdx] = particleIdx;
        // Explaination: After this every cellIndex (cellList[nParticle:nParticle+nCells])
        // links to the first particle in that cell, that then links to the next particle,
        // until the particle chain finishes with -1
    }


    //Delete old force values
    for(int particleIdx = 0; particleIdx < nParticles; particleIdx++) VZero(particles[particleIdx].force);

    //Compute interactions:
    VecR deltaR;
    real distance, forceMagnitude;
    VecI cellIdx1, cellIdx2;
    int linearCellIdx1, linearCellIdx2;
    VecI offset[] = {{0,0},{0,1},{-1,0},{-1,1},{1,1}}; // Define cell index offsets that need to be scanned (only half of neighbours to avoid double counting)
    int nOffsets = 5;
    //Go through all cells
    for(int cellXIdx = 0; cellXIdx < cells.x; cellXIdx++){
        for(int cellYIdx = 0; cellYIdx < cells.y; cellYIdx++){
            VSet(cellIdx1,cellXIdx,cellYIdx);
            //Go through this cell and its neighbouring cells
            for(int offsetIdx = 0; offsetIdx < nOffsets; offsetIdx++){
                VSAdd(cellIdx2, cellIdx1, 1, offset[offsetIdx]);
                VCellWrapAll(cellIdx2); //Periodic boundaries for cells
                linearCellIdx1 = VLinear(cellIdx1, cells) + nParticles;
                linearCellIdx2 = VLinear(cellIdx2, cells) + nParticles;
                //Go through all particles in the cells
                for(int pIdx1 = cellList[linearCellIdx1]; pIdx1 >=0; pIdx1 = cellList[pIdx1]){
                    for(int pIdx2 = cellList[linearCellIdx2]; pIdx2 >=0; pIdx2 = cellList[pIdx2]){
                        if(linearCellIdx1!=linearCellIdx2 || pIdx2 < pIdx1){ //Avoid double counting in same cell
                            VSub(deltaR, particles[pIdx1].r,particles[pIdx2].r);
                            VWrapAllTang(deltaR); //Apply periodic boundary condition
                            distance = sqrt(VLenSq(deltaR));
                            if(distance < 1){ //Do particles touch? (Warning: Needs to be changed )
                                //Persistence change
                                if(particles[pIdx1].color==1 && particles[pIdx2].color==0){
                                    particles[pIdx1].color = 2;
                                    particles[pIdx1].D = greenPersistentD;
                                    particles[pIdx1].decayTimer = tau;
                                }
                                else if (particles[pIdx1].color==0 && particles[pIdx2].color==1){
                                    particles[pIdx2].color = 2;
                                    particles[pIdx2].D = greenPersistentD;
                                    particles[pIdx2].decayTimer = tau;
                                }

                                //Direction change (CIL)
                                changeDirection(pIdx1, pIdx2, deltaR);

                                //Repulsive forces
                                if(LennardJones==1){
                                    forceMagnitude = LennardJonesForce(distance);
                                } else {
                                    forceMagnitude = HarmonicForce(distance);
                                }
                                VVSAdd(particles[pIdx1].force,forceMagnitude,deltaR);
                                VVSAdd(particles[pIdx2].force,-forceMagnitude,deltaR);

                            } else if(distance < potentialRange){ //Are particles within attraction range, but not touching?
                                // Attractive forces
                                if(LennardJones==1){
                                    forceMagnitude = LennardJonesForce(distance);
                                } else {
                                    forceMagnitude = HarmonicForce(distance);
                                }
                                
                                //Are these two red particles?
                                if(particles[pIdx1].color==0 && particles[pIdx2].color==0){
                                    VVSAdd(particles[pIdx1].force,redRedAdhesionMult*forceMagnitude,deltaR);
                                    VVSAdd(particles[pIdx2].force,-redRedAdhesionMult*forceMagnitude,deltaR);
                                } //Are these two green particles?
                                else if((particles[pIdx1].color==1 || particles[pIdx1].color==2) && (particles[pIdx2].color==1 || particles[pIdx2].color==2)){
                                    VVSAdd(particles[pIdx1].force,greenGreenAdhesionMutl*forceMagnitude,deltaR);
                                    VVSAdd(particles[pIdx2].force,-greenGreenAdhesionMutl*forceMagnitude,deltaR);
                                } // There are opposite-color particles
                                else {
                                    VVSAdd(particles[pIdx1].force,redGreenAdhesionMult*forceMagnitude,deltaR);
                                    VVSAdd(particles[pIdx2].force,-redGreenAdhesionMult*forceMagnitude,deltaR);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void EulerMaruyamaR(){
    real rootStepDuration = sqrt(stepDuration); //ToDo: Avoid repeat calls
    real root2 = sqrt(2);
    //Updates the positions r of the particle based on the forces calculated in ComputeInteractions
    VecR velocity;
    double noise[2];
    for(int particleIdx=0; particleIdx < nParticles; particleIdx++){
        //Self propulsion
        velocity.x = cos(particles[particleIdx].theta); //ToDo: Can we get less calls to sin/ cos here?
        velocity.y = sin(particles[particleIdx].theta); //WARNING: not 3d compatible
        VVSAdd(particles[particleIdx].r,stepDuration*Pe,velocity);

        //Forces
        VVSAdd(particles[particleIdx].r,stepDuration,particles[particleIdx].force);

        //Noise (ToDo: Vectorise random number generation)
        xoshiro256ss_normal(noise, &rng);
        particles[particleIdx].r.x += rootStepDuration*root2*noise[0];
        particles[particleIdx].r.y += rootStepDuration*root2*noise[1];

        #ifdef DEBUG
        double distanceSq = VLenSq(velocity)*stepDuration*stepDuration*Pe*Pe;
        distanceSq += VLenSq(particles[particleIdx].force)*stepDuration*stepDuration;
        distanceSq += 2*stepDuration*(noise[0]*noise[0]+noise[1]*noise[1]);
        if(distanceSq>MAX_STEP_DISPLACEMENT*MAX_STEP_DISPLACEMENT){
            double max = MAX_STEP_DISPLACEMENT;
            printf("Error: The displacment of %f was larger the the allowed maximum of %f\n",sqrt(distanceSq),max);
            fflush(stdout);
            assert(0);  //Exceeds maximum displacement
        }
        #endif
    }   
}

void EulerMaruyamaTheta(){
    real rootStepDuration = sqrt(stepDuration); //ToDo: Avoid repeat calls
    real root2 = sqrt(2);
    //Updates the angles of the particles
    double noise[2];
    int noiseCount = 0; //Takes care of having to generate random numbers in pairs
    for(int particleIdx=0; particleIdx < nParticles; particleIdx++){
        noiseCount = noiseCount % 2;
        if (noiseCount == 0) xoshiro256ss_normal(noise, &rng);
        particles[particleIdx].theta += rootStepDuration*root2*particles[particleIdx].D*noise[noiseCount];
        noiseCount++;
    }
}

void EnforcePeriodicBoundaries(){
    double xMax = 2*region.x;
    double yMax = 2*region.y;
    for(int particleIdx=0; particleIdx < nParticles; particleIdx++){
        //Check for large jumps
        assert(abs(particles[particleIdx].r.x)<xMax); //Step is larger than the entire box, reduce stepDuration!
        assert(abs(particles[particleIdx].r.y)<yMax);
        //Position
        VWrapAll(particles[particleIdx].r);
    }
}

void UpdatePersistence(){
    for(int particleIdx=0; particleIdx < nParticles; particleIdx++){
        if(particles[particleIdx].color==2){//Is the particle persistent?
            particles[particleIdx].decayTimer -= stepDuration;
            if(particles[particleIdx].decayTimer <= 0){
                particles[particleIdx].color = 1; 
                particles[particleIdx].D = greenD;
                particles[particleIdx].decayTimer = 0;
            }
        }
    }
}

//End: Support Functions for Single Step
void SingleStep (int stepIdx){
    ComputeInteractions();
    EulerMaruyamaR();
    EulerMaruyamaTheta();
    EnforcePeriodicBoundaries();
    UpdatePersistence();
    int measurementSteps = round(measurementInterval/stepDuration); 
    if(((stepIdx > skipSteps) && (stepIdx % measurementSteps == 0)) ){ //Time for measurement? 
        MeasurePositions(stepIdx*stepDuration);
    } else if (stepIdx==(stepLimit-1)){ //Save final step (initial position is measured by SetUpJob()
        MeasurePositions(stepIdx*stepDuration);
    }
}

void cleanup(){
    //Close files
    fclose(tracksFile);
    fclose(paramFile);

    //Free memory
    free(particles);
    free(cellList);
    
}

