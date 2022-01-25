#ifndef cell_simulation
#define cell_simulation

#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "base_types2D.h"
#include "xoshiro_rng.h"
#include "global.h"


void SetParameters(){
    stepLimit = 10000;
    initialCellsGreen = (VecI) {.x = 1, .y = 2}; 
    initialCellsRed = (VecI) {.x = 1, .y=2};
    nGreenParticles = VProd(initialCellsGreen);
    nRedParticles = VProd(initialCellsRed);
    nParticles = nGreenParticles + nRedParticles;
    stepDuration = 0.001;
    region = (VecR) {.x = 10,
               .y = 8 };  //Size of the simulation region
    D =0.1; //translational Diffusion constant
    redRotD = 0.8;
    greenRotD = 0.8;
    greenRotDPlus = 0.1;
    k = 5;  //potential strength
    tau = 2; //decay time of the persistent state
    measurementInterval = 0.1;
    nMeasurements = stepLimit*stepDuration/measurementInterval; //ToDo: Barricade so that this is multiple
    seed = time(NULL);
    rng = xoshiro256ss_init(seed);
    srand(seed);
    tracksFile = fopen("output/tracks.csv","w");
    
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
        particles[particleIdx].rotD = greenRotD;
        particles[particleIdx].decayTimer = 0;
        particles[particleIdx].color = 1;
    }
    //Red particles
    for(int particleIdx = nGreenParticles; particleIdx < nGreenParticles + nRedParticles; particleIdx++){
        particles[particleIdx].rotD = redRotD;
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


void MeasureProperties(real time){
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

void SetUpJob(){
    AllocArrays();
    InitialisePositions();
    InitialiseColor();
    InitialiseAngles();
    MeasureProperties(0);
};
real HarmonicForce(){
    return k;
}
void ComputeInteractions(){
    //Delete old force values
    for(int particleIdx = 0; particleIdx < nParticles; particleIdx++) VZero(particles[particleIdx].force);

    //ToDo: Make this more efficient
    VecR deltaR;
    real distanceSqrt, forceMagnitude;
    for(int pIdx1 = 0; pIdx1 < nParticles; pIdx1++){
        for(int pIdx2 = pIdx1+1; pIdx2 < nParticles; pIdx2++){
                VSub(deltaR, particles[pIdx1].r,particles[pIdx2].r);
                VWrapAllTang(deltaR); //Apply periodic boundary condition
                distanceSqrt = VLenSq(deltaR);
                if(distanceSqrt < 1){ //WARNING: this needs to be changed if we don't measure in units of r0
                    //Forces: 
                    forceMagnitude = HarmonicForce();
                    VVSAdd(particles[pIdx1].force,forceMagnitude,deltaR);
                    VVSAdd(particles[pIdx2].force,-forceMagnitude,deltaR);
                    
                    //Persistence change (no refreshing)
                    if(particles[pIdx1].color==1 && particles[pIdx2].color==0){
                        particles[pIdx1].color = 2;
                        particles[pIdx1].rotD = greenRotDPlus;
                        particles[pIdx1].decayTimer = tau;
                    }
                    else if (particles[pIdx1].color==0 && particles[pIdx2].color==1){
                        particles[pIdx2].color = 2;
                        particles[pIdx2].rotD = greenRotDPlus;
                        particles[pIdx2].decayTimer = tau;
                    }


                }
        }
    }
}
void EulerMaruyamaR(){
    real rootStepDuration = sqrt(stepDuration); //ToDo: Avoid repeat calls
    //Updates the positions r of the particle based on the forces calculated in ComputeInteractions
    VecR velocity;
    double noise[2];
    for(int particleIdx=0; particleIdx < nParticles; particleIdx++){
        //Self propulsion
        velocity.x = cos(particles[particleIdx].theta); //ToDo: Can we get less calls to sin/ cos here?
        velocity.y = sin(particles[particleIdx].theta); //WARNING: not 3d compatible
        VVSAdd(particles[particleIdx].r,stepDuration,velocity);

        //Forces
        VVSAdd(particles[particleIdx].r,stepDuration,particles[particleIdx].force);

        //Noise (ToDo: Vectorise random number generation)
        xoshiro256ss_normal(noise, &rng);
        particles[particleIdx].r.x += rootStepDuration*D*noise[0];
        particles[particleIdx].r.y += rootStepDuration*D*noise[1];
    }   
}

void EulerMaruyamaTheta(){
    real rootStepDuration = sqrt(stepDuration); //ToDo: Avoid repeat calls
    //Updates the angles of the particles
    double noise[2];
    int noiseCount = 0; //Takes care of having to generate random numbers in pairs
    for(int particleIdx=0; particleIdx < nParticles; particleIdx++){
        noiseCount = noiseCount % 2;
        if (noiseCount == 0) xoshiro256ss_normal(noise, &rng);
        particles[particleIdx].theta += rootStepDuration*particles[particleIdx].rotD*noise[noiseCount];
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
                particles[particleIdx].rotD = greenRotD;
                particles[particleIdx].decayTimer = 0;
            }
        }
    }
}

void SingleStep (int stepIdx){
    ComputeInteractions();
    EulerMaruyamaR();
    EulerMaruyamaTheta();
    EnforcePeriodicBoundaries();
    UpdatePersistence();

    int measurementSteps = measurementInterval/stepDuration; //ToDo: Barricades to enforce multiples
    if(stepIdx >0 && stepIdx % measurementSteps == 0){ //Initial measurement is done be SetUpJob()
       MeasureProperties(stepIdx*stepDuration);
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