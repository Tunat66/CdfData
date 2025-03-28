#pragma once
#include <complex>
#include <vector>
#include <cuda_runtime.h>
#include <cmath>
#include <stdexcept>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

namespace Cdf
{
    
    //a plain track struct rather than a track class
    typedef struct plain_track {
        double cotTheta;
        double curvature;
        double d0;
        double phi0;
        double z0;
    } plain_track; 

    //overloaded operator to check for two plain tracks
    __host__ __device__ inline bool operator==(const plain_track &trk1, const plain_track &trk2)
    {
        return trk1.cotTheta == trk2.cotTheta &&
               trk1.curvature == trk2.curvature &&
               trk1.d0 == trk2.d0 &&
               trk1.phi0 == trk2.phi0 &&
               trk1.z0 == trk2.z0;
    }

    //standalone methods that mimic the methods within helix
    __host__ __device__ double radius(double curvature);
    __host__ __device__ void center(double h, double d0, double phi0, double* real, double* imag);
    __host__ __device__ double findZFromXY(plain_track trk, double x, double y, int index, int q);
    __host__ __device__ double** intersect(plain_track trk1, plain_track trk2, int* num_intersections);
    
}
