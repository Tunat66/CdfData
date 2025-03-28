#pragma once

#include <vector>
#include <array>
#include <cmath>
#include <complex>
#include "Helix.cuh"

namespace Cdf
{
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