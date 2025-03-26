#include "K0SAnalysis.hh"
#include <iostream>
#include <cmath>
#include <fstream>
namespace Cdf 
{
// Constructor
K0SAnalysis::K0SAnalysis(double pTcutoff, double TrackImpactParametercutoff, double Lxycutoff, double ImpactParametercutoff)
    : pTcutoff(pTcutoff),
      TrackImpactParametercutoff(TrackImpactParametercutoff),
      Lxycutoff(Lxycutoff),
      ImpactParametercutoff(ImpactParametercutoff) {}

// Start method
void K0SAnalysis::start() {
    // Currently empty, can be filled if needed
}

// Stop method
void K0SAnalysis::stop() {
    // Remove zeroes from the mass and lifetime arrays
    massArray.erase(std::remove(massArray.begin(), massArray.end(), 0), massArray.end());
    lifetimeArray.erase(std::remove(lifetimeArray.begin(), lifetimeArray.end(), 0), lifetimeArray.end());

    //put the data out into Dat files
    // Write massArray to a file
    // Remove zeroes and NaN entries from the mass and lifetime arrays
    massArray.erase(std::remove_if(massArray.begin(), massArray.end(),
                                   [](double value) { return value == 0 || std::isnan(value); }),
                    massArray.end());
    lifetimeArray.erase(std::remove_if(lifetimeArray.begin(), lifetimeArray.end(),
                                       [](double value) { return value == 0 || std::isnan(value); }),
                        lifetimeArray.end());
    std::ofstream massFile("massArray.dat");
    if (massFile.is_open()) {
        for (const auto& mass : massArray) {
            massFile << mass << "\n";
        }
        massFile.close();
        std::cout << "Mass data written to massArray.dat" << std::endl;
    } else {
        std::cerr << "Error: Could not open massArray.dat for writing." << std::endl;
    }

    // Write lifetimeArray to a file
    std::ofstream lifetimeFile("lifetimeArray.dat");
    if (lifetimeFile.is_open()) {
        for (const auto& lifetime : lifetimeArray) {
            lifetimeFile << lifetime << "\n";
        }
        lifetimeFile.close();
        std::cout << "Lifetime data written to lifetimeArray.dat" << std::endl;
    } else {
        std::cerr << "Error: Could not open lifetimeArray.dat for writing." << std::endl;
    }

}

// Event method
void K0SAnalysis::event(Event* ev) {
    // Get the number of tracks
    auto tracks = ev->getTracks();
    int numTracks = tracks.size();
    
    // Allocate memory on the GPU
    double* d_trackData;
    double* d_massArray;
    double* d_lifetimeArray;
    int* d_massCounter;

    cudaMalloc(&d_trackData, numTracks * sizeof(double));
    cudaMalloc(&d_massArray, numTracks * sizeof(double));
    cudaMalloc(&d_lifetimeArray, numTracks * sizeof(double));
    cudaMalloc(&d_massCounter, sizeof(int));
    
    
    for (const auto& track : ev->getTracks()) {
        Helix helicalTrack(track);

        // Get vertices that satisfy cutoff criteria
        auto vertices = vertex_checkOverlaps(helicalTrack, ev);
        if (vertices.empty()) {
            continue;
        }

        // Iterate through vertices and calculate masses and lifetimes
        for (const auto& vertex : vertices) {
            massCounter++;
            double mass_kaon = vertex->mass(m_pion, m_pion);
            massArray.push_back(mass_kaon);

            double lifetime_kaon = vertex->lifetime(mass_kaon);
            lifetimeArray.push_back(lifetime_kaon);
        }
    }
}

// Helper method to check overlaps
std::vector<std::shared_ptr<Vertex>> K0SAnalysis::vertex_checkOverlaps(const Helix& helicalTrack, Event* ev) {
    std::vector<std::shared_ptr<Vertex>> validVertices;

    for (const auto& track : ev->getTracks()) {
        Helix h(track);
        auto vertex = std::make_shared<Vertex>(helicalTrack, h, ev);
        auto dPVarr = vertex->dPV();

        // Skip invalid vertices
        if (vertex->isBroken()) {
            //std::cout << "broken vertex" << std::endl;
            //std::cout << "Broken Vertex Parameters: " << vertex->ImpactParameter() << " " << vertex->pT() << " " << dPVarr[0] << " " << dPVarr[1] << " " << vertex->Lxy() << std::endl;

            continue;
        }
        //std::cout << "Vertex Parameters: " << vertex->ImpactParameter() << " " << vertex->pT() << " " << dPVarr[0] << " " << dPVarr[1] << " " << vertex->Lxy() << std::endl;


        // Apply filtering criteria
        if (std::abs(vertex->ImpactParameter()) > ImpactParametercutoff) {
            continue;
        }
        if (std::abs(vertex->pT()) < pTcutoff) {
            continue;
        }
        if (std::abs(dPVarr[0]) < TrackImpactParametercutoff || std::abs(dPVarr[1]) < TrackImpactParametercutoff) {
            continue;
        }
        if (std::abs(vertex->Lxy()) < Lxycutoff) {
            continue;
        }
        //std::cout << "Accepted Vertex Parameters: " << vertex->ImpactParameter() << " " << vertex->pT() << " " << dPVarr[0] << " " << dPVarr[1] << " " << vertex->Lxy() << std::endl;

        // Add valid vertex
        validVertices.push_back(vertex);
    }

    return validVertices;
}
}
