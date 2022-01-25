#ifndef base_types
#define base_types

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Data types ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

typedef double real;

//2-dimensional vector
typedef struct { 
    real x,y;
} VecR;

typedef struct { 
    int x,y;
} VecI;

//active brownian motion particle
typedef struct {
    VecR r; // position of particle
    real theta; // unit vector showing velocity direction
    VecR force; //force displacing the particle this step
    real rotD; // rotational diffusion constant
    real decayTimer; //Time it takes to go back to less persistent state
    int color; //red = 0, green =1 , persistent green = 2
} particle ;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Linear algebra ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#define VDot(v1, v2) \
     ((v1).x * (v2).x + (v1).y * (v2).y)
#define VSAdd(v1, v2, s3, v3) \
    ((v1).x = (v2).x + (s3) * (v3).x, \
    (v1).y = (v2).y + (s3) * (v3).y)
#define VSub(v1, v2, v3) \
    (v1).x = (v2).x - (v3).x, \
    (v1).y = (v2).y - (v3).y
#define VSet(v, sx, sy) \
    (v).x = sx, \
    (v).y = sy 
#define VSetAll(v, s)       VSet (v, s, s)
#define VZero(v)            VSetAll (v, 0)
#define VVSAdd(v1, s2, v2)  VSAdd (v1, v1, s2, v2)
#define VLenSq(v)           VDot (v, v)

#define VMul(v1, v2, v3) \
    (v1).x = (v2).x * (v3).x, \
    (v1).y = (v2).y * (v3).y

#define VScalMult(v1,s2,v2) \
    (v1).x = s2*(v2).x, \
    (v1).y = s2*(v2).y

#define VProd(v1) ((v1).x * (v1).y)

#define VDiv(v1, v2, v3) \
    (v1).x = (v2).x / (v3).x, \
    (v1).y = (v2).y /(v3).y


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Periodic boundaries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//Wraps tangent vectors, i.e. those that leave the top right quadrant
//assumes (rectangular) boundary is given as VecR region
#define VWrapTang(v, t) \
    if (v.t >= 0.5 * region.t) v.t -= region.t; \
    else if (v.t < -0.5 * region.t) v.t += region.t 
#define VWrapAllTang(v) \
    {VWrapTang (v, x); \
    VWrapTang (v, y);}
//Wraps position vectors, that only live in the top right quadrant
#define VWrap(v,t)\
    if (v.t >= region.t) v.t -= region.t; \
    else if (v.t < 0) v.t += region.t;
#define VWrapAll(v) \
    {VWrap(v,x); \
     VWrap(v,y);}


#endif