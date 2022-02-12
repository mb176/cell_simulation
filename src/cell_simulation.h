#ifndef cell_simulation
#define cell_simulation

#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include "base_types2D.h"
#include "xoshiro_rng.h"


double getNextParameter(FILE * file, char * parameterName){
    char line[MAX_LINE_LENGTH] = {0};
    char name[MAX_STRING_LENGTH];
    char value[MAX_STRING_LENGTH];
    if (fgets(line,MAX_LINE_LENGTH,file)) {
        sscanf(line, "%s : %s", name, value);
    } 
    assert(*name==*parameterName); //Invalid parameter or parameter out of order or not enough parameters
    char * ptr; //Throughaway, needed to call strtod
    double val = strtod(value,&ptr);
    printf("%s = %f \n",name,val);
    return val;
}

void SetParameters(int argc, char ** argv){
    char * PATH = argv[1];
    printf(PATH);
    printf("\n");
    //char * NAME = argv[2];

    FILE * file;
    file = fopen(PATH,"r");
    
    //Read out parameters
    stepLimit = getNextParameter(file, "stepLimit");
    nGreenParticles = getNextParameter(file, "nGreenParticles");
    nRedParticles = getNextParameter(file, "nRedParticles");
    stepDuration = getNextParameter(file, "stepDuration");
    measurementInterval = getNextParameter(file, "measurementInterval");
    real areaFraction = getNextParameter(file, "areaFraction");
    redPe = getNextParameter(file, "redPe");
    greenPe = getNextParameter(file, "greenPe");
    greenPlusPe = getNextParameter(file, "greenPlusPe");
    k = getNextParameter(file, "k");
    tau = getNextParameter(file, "tau");
    v0 = getNextParameter(file, "v0");
    sigma = getNextParameter(file, "sigma");

    fclose(file);

    // Barricades
    //measurementInterval must be a multiple of stepDuration
    assert(measurementInterval/stepDuration==floor(measurementInterval/stepDuration));
    //The Simulation duration must be a multiple of measurementInterval
    assert(stepLimit*stepDuration/measurementInterval==floor(stepLimit*stepDuration/measurementInterval)); 
    //The particle numbers must be square of an integer
    assert(sqrt(nGreenParticles)==floor(sqrt(nGreenParticles)));
    assert(sqrt(nRedParticles)==floor(sqrt(nRedParticles)));
    //Particle numbers can't be zero
    assert(nGreenParticles>0);
    assert(nRedParticles >0);

    //Initialise secondary parameters
    initialCellsGreen = (VecI) {.x = sqrt(nGreenParticles), .y = sqrt(nGreenParticles)}; 
    initialCellsRed = (VecI) {.x = sqrt(nRedParticles), .y=sqrt(nRedParticles)};
    nParticles = nGreenParticles + nRedParticles;
    //Setup Square area to fit the area fraction; each particle covers area of size pi*0.5^2 = 0.785398
    real length = sqrt(nParticles*0.785398/areaFraction);
    printf("length=%f \n",length);
    region = (VecR) {.x = length,
                     .y = length};  //Size of the simulation region

    nMeasurements = stepLimit*stepDuration/measurementInterval; 
    seed = time(NULL);
    rng = xoshiro256ss_init(seed);
    srand(seed);

    char  suffix[] = "_tracks";
    tracksFile = fopen(strcat(PATH,suffix),"w");
}

void AllocArrays(){
    //Particle array
    particles = (particle *) malloc(nParticles*(sizeof(particle)));
    //Measurement arrays
    measurementTimes = (real *) malloc(nMeasurements*sizeof(real));
    positionMeasurements = (VecR **) malloc(nMeasurements*sizeof(VecR *));
    colorMeasurements = (int **) malloc(nMeasurements*sizeof(int*));
    for(int idx = 0; idx < nMeasurements; idx++){
        positionMeasurements[idx] = (VecR *) malloc(nParticles*sizeof(particle));
        colorMeasurements[idx] = (int *) malloc(nParticles*sizeof(particle));
    }   
};

void InitialisePositions(){
    // ToDo: Make sure there is no overlap if nGreenParticles!=nRedParticles
    // ToDo: Make sure it works for nRedParticles = 0
    
    int particleIdx = 0;
    //Green particles 
    VecR cellSizeGreen;
    VDiv(cellSizeGreen, region, initialCellsGreen);
    for(int xIdx = 0; xIdx < initialCellsGreen.x; xIdx++){
        for(int yIdx = 0; yIdx < initialCellsGreen.y; yIdx++){
            particles[particleIdx].r.x = cellSizeGreen.x * (xIdx+0.33); //Small offset for green
            particles[particleIdx].r.y = cellSizeGreen.y * (yIdx+0.33);
            particleIdx++;
        }
    }
    //Red particles 
    VecR cellSizeRed;
    VDiv(cellSizeRed, region, initialCellsRed);
    for(int xIdx = 0; xIdx < initialCellsRed.x; xIdx++){
        for(int yIdx = 0; yIdx < initialCellsRed.y; yIdx++){
            particles[particleIdx].r.x = cellSizeRed.x * (xIdx+0.66); //Big offset for red
            particles[particleIdx].r.y = cellSizeRed.y * (yIdx+0.66);
            particleIdx++;
        }
    }

};

void InitialiseColor(){
    //Green particles
    for(int particleIdx = 0; particleIdx < nGreenParticles; particleIdx++){
        particles[particleIdx].Pe = greenPe;
        particles[particleIdx].decayTimer = 0;
        particles[particleIdx].color = 1;
    }
    //Red particles
    for(int particleIdx = nGreenParticles; particleIdx < nGreenParticles + nRedParticles; particleIdx++){
        particles[particleIdx].Pe = redPe;
        particles[particleIdx].decayTimer = 0;
        particles[particleIdx].color = 0;
    }
};

void InitialiseAngles(){
    //Set theta randomly in [0, 2 pi] and set decay time to 0
    real sum = 0;
    for(int particleIdx=0; particleIdx < nParticles; particleIdx++){
        particles[particleIdx].theta = (float) rand()/(float)(RAND_MAX/(2*M_PI));
        sum += particles[particleIdx].theta;
    }
    //set the average of theta to zero
    for(int particleIdx=0; particleIdx < nParticles; particleIdx++){
        particles[particleIdx].theta -= sum/nParticles;
    }
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
void ComputeInteractions(){
    //Delete old force values
    for(int particleIdx = 0; particleIdx < nParticles; particleIdx++) VZero(particles[particleIdx].force);

    //ToDo: Make this more efficient
    VecR deltaR;
    real distance, forceMagnitude;
    for(int pIdx1 = 0; pIdx1 < nParticles; pIdx1++){
        for(int pIdx2 = pIdx1+1; pIdx2 < nParticles; pIdx2++){
                VSub(deltaR, particles[pIdx1].r,particles[pIdx2].r);
                VWrapAllTang(deltaR); //Apply periodic boundary condition
                distance = sqrt(VLenSq(deltaR));
                if(distance < sigma){ //WARNING: this needs to be changed if we don't measure in units of r0
                    //Forces: 
                    forceMagnitude = HarmonicForce(distance);
                    VVSAdd(particles[pIdx1].force,forceMagnitude,deltaR);
                    VVSAdd(particles[pIdx2].force,-forceMagnitude,deltaR);
                    
                    //Persistence change (no refreshing)
                    if(particles[pIdx1].color==1 && particles[pIdx2].color==0){
                        particles[pIdx1].color = 2;
                        particles[pIdx1].Pe = greenPlusPe;
                        particles[pIdx1].decayTimer = tau;
                    }
                    else if (particles[pIdx1].color==0 && particles[pIdx2].color==1){
                        particles[pIdx2].color = 2;
                        particles[pIdx2].Pe = greenPlusPe;
                        particles[pIdx2].decayTimer = tau;
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
        VVSAdd(particles[particleIdx].r,stepDuration*v0,velocity);

        //Forces
        VVSAdd(particles[particleIdx].r,stepDuration,particles[particleIdx].force);

        //Noise (ToDo: Vectorise random number generation)
        xoshiro256ss_normal(noise, &rng);
        particles[particleIdx].r.x += rootStepDuration*root2*noise[0];
        particles[particleIdx].r.y += rootStepDuration*root2*noise[1];
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
        particles[particleIdx].theta += rootStepDuration*root2*particles[particleIdx].Pe*noise[noiseCount];
        noiseCount++;
    }
}

void EnforcePeriodicBoundaries(){
    for(int particleIdx=0; particleIdx < nParticles; particleIdx++){
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
                particles[particleIdx].Pe = greenPe;
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
    if((stepIdx >0) && (stepIdx % measurementSteps == 0)){ //Initial measurement is done be SetUpJob()
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

    //Free memory
    free(particles);
    free(measurementTimes);
    for(int idx = 0; idx < nMeasurements; idx++){
        free(positionMeasurements[idx]);
        free(colorMeasurements[idx]);
    }
    free(positionMeasurements);
    free(colorMeasurements);
}


#endif