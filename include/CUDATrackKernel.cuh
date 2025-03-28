#pragma once  // Include guard

#include <cuda_runtime.h>
#include <vector>
#include <cmath>
#include <complex>
#include <math.h>
#include <thrust/complex.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdexcept>
#include <stdio.h>
#include <stdlib.h>

#define M_PI 3.14159265358979323846

namespace Cdf
{
    //HELIX DECLERATIONS: these were once part of a helix class


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

    //VERTEX DECLERATIONS: these were once part of a vertex class
    //created some structs as CUDA does not like std::array 
    //type for vertex
    //a plain vertex
    typedef struct plain_vertex {
        double x = -999;
        double y = -999;
        double zbar = -999;
        double deltaz = -999;
        bool isValid = false; // Flag to indicate if the vertex is valid
    } plain_vertex;

    // Standalone method declarations

    // Return the vertex as an array of doubles
    __host__ __device__ plain_vertex return_vertex(plain_track trk1, plain_track trk2);

    // Calculate transverse momentum (pT) for a track
    __host__ __device__ double pt(plain_track trk);

    // Calculate 3-momentum for a track
    __host__ __device__ void p(plain_track trk, plain_vertex vertex, double* result);

    // Calculate total transverse momentum (pT) for the vertex
    __host__ __device__ double pT(plain_track trk1, plain_track trk2);

    // Calculate impact parameter in the transverse plane for a track
    __host__ __device__ double dPV_trk(plain_track trk, double primaryVertex[2]);

    // Calculate invariant mass of the vertex
    __host__ __device__ double mass(double mass1, double mass2, plain_track trk1, plain_track trk2, plain_vertex vertex);

    // Calculate impact parameter of the vertex
    __host__ __device__ double ImpactParameter(plain_track trk1, plain_track trk2, plain_vertex vertex, double primaryVertex[2]);

    // Calculate Lxy (transverse displacement)
    __host__ __device__ double Lxy(plain_track trk1, plain_track trk2, plain_vertex vertex, double primaryVertex[2]);

    // Calculate lifetime in picoseconds
    __host__ __device__ double lifetime(double m, plain_track trk1, plain_track trk2, plain_vertex vertex, double primaryVertex[2]);


    class Vertex {
        public:
        //cutoff values which are used by the CUDA kernel for vertexing
            
            static constexpr double kpc = 0.002116;
            static constexpr double c = 2.99792458; // Speed of light in cgs units
    };
    
}

// Kernel decleration
__global__ void processTracksKernel(double** trackData, int numTracks, double* primaryVertex, double* massArray, double* lifetimeArray, int* massCounter, double m_pion);
__device__ Cdf::plain_vertex track_compare(Cdf::plain_track track_s, Cdf::plain_track track_o, double primaryVertex[2]);
