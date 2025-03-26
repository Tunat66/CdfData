#include <cuda_runtime.h>
#include <vector>
#include <iostream>

// Kernel to process tracks
__global__ void processTracksKernel(const double* trackData, int numTracks, double* massArray, double* lifetimeArray, int* massCounter, double m_pion) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numTracks) return;

    // Process the track (simplified example)
    double track = trackData[idx];

    // Simulate vertex processing and mass/lifetime calculation
    double mass_kaon = track + m_pion; // Replace with actual calculation
    double lifetime_kaon = mass_kaon * 0.1; // Replace with actual calculation

    // Write results to global memory
    massArray[idx] = mass_kaon;
    lifetimeArray[idx] = lifetime_kaon;

    // Increment massCounter (atomic operation to avoid race conditions)
    atomicAdd(massCounter, 1);
}