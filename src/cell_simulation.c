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

void SetParameters(char** argv){
    char * PATH = argv[1];
    printf("Parameter file: "); printf(PATH); printf("\n");
    printf("Initialise paramters:\n");
    
    // FILE * paramFile;
    paramFile = fopen(PATH,"r+");
    
    //Read out parameters
    stepLimit = getNextParameter(paramFile, "stepLimit");
    nGreenParticles = getNextParameter(paramFile, "nGreenParticles");
    nRedParticles = getNextParameter(paramFile, "nRedParticles");
    stepDuration = getNextParameter(paramFile, "stepDuration");
    measurementInterval = getNextParameter(paramFile, "measurementInterval");
    real areaFraction = getNextParameter(paramFile, "areaFraction");
    redD = getNextParameter(paramFile, "redD");
    greenD = getNextParameter(paramFile, "greenD");
    greenPlusD = getNextParameter(paramFile, "greenPlusD");
    k = getNextParameter(paramFile, "k");
    tau = getNextParameter(paramFile, "tau");
    Pe = getNextParameter(paramFile, "Pe");
    sigma = getNextParameter(paramFile, "sigma");
    LennardJones = getNextParameter(paramFile, "LennardJones");
    skippedSteps = getNextParameter(paramFile, "skippedSteps"); 

    // Barricades
    //measurementInterval must be a multiple of stepDuration
    assert((measurementInterval/stepDuration)==floor(measurementInterval/stepDuration));
    //Can't skip more steps than we have
    assert(stepLimit>skippedSteps);
    //The Simulation duration must be a multiple of measurementInterval
    assert((stepLimit-skippedSteps)*(stepDuration/measurementInterval)==floor((stepLimit-skippedSteps)*(stepDuration/measurementInterval))); 
    //Need it least 20 steps, otherwise the loading bar divides by zero
    assert(stepLimit>=20);

    //Initialise secondary parameters
    nParticles = nGreenParticles + nRedParticles;
    //Setup Square area to fit the area fraction; each particle covers area of size pi*0.5^2 = 0.785398
    real length = sqrt(nParticles*0.785398/areaFraction);
    printf("length=%f \n",length);
    region = (VecR) {.x = length,
                     .y = length};  //Size of the simulation region
    cells.x = region.x/(sigma); //Smallest cells that are larger than interaction range
    cells.y = region.y/(sigma);

    nMeasurements = (stepLimit-skippedSteps)*(stepDuration/measurementInterval) + 1; //+1 one because we do the first and the last step 
    seed = SEED;//time(NULL);
    rng = xoshiro256ss_init(seed);
    srand(seed);

    char  suffix[] = "_tracks.csv";
    // char * name = strcat(PATH,suffix);
    tracksFile = fopen(strcat(PATH,suffix),"w");
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
    particles = (particle *) malloc(nParticles*(sizeof(particle)));
    //Linked list for cells
    cellList = (int *) malloc((nParticles + VProd(cells))*sizeof(int));
    //Measurement arrays
    measurementTimes = (real *) malloc(nMeasurements*sizeof(real));
    positionMeasurements = (VecR **) malloc(nMeasurements*sizeof(VecR *));
    colorMeasurements = (int **) malloc(nMeasurements*sizeof(int*));
    for(int idx = 0; idx < nMeasurements; idx++){
        positionMeasurements[idx] = (VecR *) malloc(nParticles*sizeof(particle));
        colorMeasurements[idx] = (int *) malloc(nParticles*sizeof(particle));
    }   
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
    double sigma = para[0];

    double V = 0;
    VecR deltaR;
    real distance;
    for(int pIdx1 = 0; pIdx1 < nParticles; pIdx1++){
        for(int pIdx2 = pIdx1+1; pIdx2 < nParticles; pIdx2++){
            VSub(deltaR, positions[pIdx1],positions[pIdx2]);
            VWrapAllTang(deltaR); //Periodic Boundary conditions
            distance = sqrt(VLenSq(deltaR));
            if(distance < sigma){
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
    double sigma = para[0]; //Interaction range
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
                else if(distance < sigma){
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
    static int measurementIdx = 0;
    //Times
    measurementTimes[measurementIdx] = time;
    for(int particleIdx=0; particleIdx < nParticles; particleIdx++){
        //Positions
        positionMeasurements[measurementIdx][particleIdx] = particles[particleIdx].r;
        //Colors
        colorMeasurements[measurementIdx][particleIdx] = particles[particleIdx].color;
    }   
    measurementIdx++;
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
    return -k/distance*(distance-1); //WARNING: this needs to be changed if we don't measure in units of sigma
}

real LennardJonesForce(real distance){
    double length = 0.890898718140339; //=0.5^(1/6) This value makes sure the potential decayes to zero for r=1 
    //WARNING: needs to be changed if cellRadius!=1
    double invDistance = length/(distance); 
    double invDistance6=invDistance*invDistance*invDistance*invDistance*invDistance*invDistance;
    return 48*k*invDistance6*(invDistance6-0.5)*invDistance*invDistance;
}

real GetAngle(VecR r){
    double theta;
    if(r.x>=0){
        theta = atan(r.y/r.x);
    } else {
        theta = atan(r.y/r.x)+M_PI ;
    }
    return theta;
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
                            if(distance < sigma){ //Do particles interact?
                                //Forces: 
                                if(LennardJones==1){
                                    forceMagnitude = LennardJonesForce(distance);
                                } else {
                                    forceMagnitude = HarmonicForce(distance);
                                }
                                 //LennardJonesForce(distance); 
                                VVSAdd(particles[pIdx1].force,forceMagnitude,deltaR);
                                VVSAdd(particles[pIdx2].force,-forceMagnitude,deltaR);
                            }
                            if(distance < 1){ //Do particles touch? (Warning: Needs to be changed )
                                //Persistence change (no refreshing)
                                if(particles[pIdx1].color==1 && particles[pIdx2].color==0){
                                    particles[pIdx1].color = 2;
                                    particles[pIdx1].D = greenPlusD;
                                    particles[pIdx1].decayTimer = tau;
                                    //Contact inhibited locomotion: Cells move away after contact
                                    real theta = GetAngle(deltaR);
                                    particles[pIdx1].theta = theta;
                                    particles[pIdx2].theta = theta + M_PI;

                                }
                                else if (particles[pIdx1].color==0 && particles[pIdx2].color==1){
                                    particles[pIdx2].color = 2;
                                    particles[pIdx2].D = greenPlusD;
                                    particles[pIdx2].decayTimer = tau;
                                    //Contact inhibited locomotion: Cells move away after contact
                                    real theta = GetAngle(deltaR);
                                    particles[pIdx1].theta = theta;
                                    particles[pIdx2].theta = theta + M_PI;
                                    
                                } else if (particles[pIdx1].color==2 && particles[pIdx2].color==0){
                                    //Contact inhibited locomotion: Cells move away after contact
                                    real theta = GetAngle(deltaR);
                                    particles[pIdx1].theta = theta;
                                    particles[pIdx2].theta = theta + M_PI;
                                    
                                } else if (particles[pIdx1].color==0 && particles[pIdx2].color==2){
                                    //Contact inhibited locomotion: Cells move away after contact
                                    real theta = GetAngle(deltaR);
                                    particles[pIdx1].theta = theta;
                                    particles[pIdx2].theta = theta + M_PI;
                                    
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
    int measurementSteps = measurementInterval/stepDuration; 
    if(((stepIdx > skippedSteps) && (stepIdx % measurementSteps == 0)) ){ //Time for measurement? 
        MeasurePositions(stepIdx*stepDuration);
    } else if (stepIdx==(stepLimit-1)){ //Save final step (initial position is measured by SetUpJob()
        MeasurePositions(stepIdx*stepDuration);
    }
}


void writeTracks(){
    for(int measurementIdx=0; measurementIdx<nMeasurements; measurementIdx++){
        fprintf(tracksFile, "%f", measurementTimes[measurementIdx]);
        for(int particleIdx = 0; particleIdx<nParticles; particleIdx++){
            //Print color
            if(colorMeasurements[measurementIdx][particleIdx] ==0){
                fprintf(tracksFile,", red");
            } else if (colorMeasurements[measurementIdx][particleIdx] ==1){
                fprintf(tracksFile,", green");
            } else if (colorMeasurements[measurementIdx][particleIdx] == 2){
                fprintf(tracksFile,", greenPlus");
            }
            //Print positions
            fprintf(tracksFile,", %f, %f",positionMeasurements[measurementIdx][particleIdx].x, \
                                          positionMeasurements[measurementIdx][particleIdx].y);
            
        }
        fprintf(tracksFile,"\n");
    }

}

void cleanup(){
    //Close files
    fclose(tracksFile);
    fclose(paramFile);

    //Free memory
    free(particles);
    free(cellList);
    free(measurementTimes);
    for(int idx = 0; idx < nMeasurements; idx++){
        free(positionMeasurements[idx]);
        free(colorMeasurements[idx]);
    }
    free(positionMeasurements);
    free(colorMeasurements);
}

