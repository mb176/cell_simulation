#include <stdio.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include "base_types2D.h"


int nParticles = 2;
real sigma = 1.2;
VecR region = {.x = 2, .y = 2};

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
    // printf("V = %f \n", V);
    return V;
    free(positions);
}

//Gradient of the potential
void dfHarmonicPotential(const gsl_vector *v, void *params, gsl_vector *df){
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
                if(distance < sigma){
                    VSAdd(gradient[gradientIdx], gradient[gradientIdx], 2*(1-1/distance), deltaR);
                }
            }
        }
    }
    //Write gradient to gsl vector
    for(int gradientIdx = 0; gradientIdx < nParticles; gradientIdx++){
        gsl_vector_set(df, 2*gradientIdx    , gradient[gradientIdx].x);
        gsl_vector_set(df, 2*gradientIdx + 1, gradient[gradientIdx].y);
    }
    // printf ("Gradient: (%.3f %.3f) (%.3f %.3f)\n",
        // gsl_vector_get (df, 0),
        // gsl_vector_get (df, 1),
        // gsl_vector_get (df, 2),
        // gsl_vector_get (df, 3));
    free(positions);
    free(gradient);
}

//Potential and Gradient together
void fdfHarmonicPotential(const gsl_vector *x, void * params, double * f, gsl_vector *df){
    *f = fHarmonicPotential(x, params);
    dfHarmonicPotential(x, params, df);
}

int main(){
    //Parameters
    double stepSize = 0.01; //Initial step size
    double lineTol = 1e-4; //Error tolerance for the line minimisation in a given direction
    double tol = 1e-4; //Error tolerance for the overall minimisation
    int maxIter = 100;

    //Initialise function to minimise
    double params[1] = {sigma};
    gsl_multimin_function_fdf func;
    func.f      =   &fHarmonicPotential;
    func.df     =  &dfHarmonicPotential;
    func.fdf    = &fdfHarmonicPotential;
    func.n      = 2*nParticles; //Number of variables
    func.params = (void*) params;

    size_t iter = 0;
    int status;
    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s;

    //starting point
    gsl_vector *x;
    x = gsl_vector_alloc(2*nParticles);

    VecR positions[2] = {{0.2,0},{0.1,0}};
    
    for(int particleIdx = 0; particleIdx < nParticles; particleIdx++){
        gsl_vector_set(x,2*particleIdx,     positions[particleIdx].x);
        gsl_vector_set(x,2*particleIdx + 1, positions[particleIdx].y);
    }

    //Initialise minimizer
    T = gsl_multimin_fdfminimizer_conjugate_fr; //Fletcher-Reeves conjugate gradient algorithm
    s = gsl_multimin_fdfminimizer_alloc (T, 2*nParticles); 
    gsl_multimin_fdfminimizer_set (s, &func, x, stepSize, lineTol);

    printf("Starting position:\n");

    printf ("(%.3f %.3f) (%.3f %.3f) %10.5f\n", 
            gsl_vector_get (s->x, 0),
            gsl_vector_get (s->x, 1),
            gsl_vector_get (s->x, 2),
            gsl_vector_get (s->x, 3),
            s->f);

    do{
    iter++;
    status = gsl_multimin_fdfminimizer_iterate (s);
    
    if (status)
        break;

    status = gsl_multimin_test_gradient (s->gradient, tol);

    if (status == GSL_SUCCESS)
        printf ("Minimum found at:\n");

    printf ("%5d (%.3f %.3f) (%.3f %.3f) %10.5f\n", iter,
            gsl_vector_get (s->x, 0),
            gsl_vector_get (s->x, 1),
            gsl_vector_get (s->x, 2),
            gsl_vector_get (s->x, 3),
            s->f);

    } while (status == GSL_CONTINUE && iter < maxIter);

    gsl_multimin_fdfminimizer_free (s);
    gsl_vector_free (x);
}