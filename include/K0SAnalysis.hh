#pragma once

#include <vector>
#include <memory>
#include "Analysis.hh"
#include "Event.hh"
#include "CUDATrackKernel.cuh"
#include <iostream>
#include <cmath>
#include <fstream>

namespace Cdf
{
    class K0SAnalysis : public Analysis {
        private:
            // Adjustable parameters
            double pTcutoff;
            double TrackImpactParametercutoff; // in cm
            double Lxycutoff;
            double ImpactParametercutoff;
        
            // Results
            std::vector<double> massArray;
            std::vector<double> lifetimeArray;
            //std::unique_ptr<Histogram> mass;
            //std::unique_ptr<Histogram> lifetime;
        
            // Constants
            const double GeV = 1e9;
            const double m_pion = 0.13957039; // GeV
            //int massCounter = 0;
            const int nbins = 100;
            const double mmin = 0.4; // GeV
            const double mmax = 0.6; // GeV
            const double tmin = 10;  // ps
            const double tmax = 200; // ps
            static constexpr int maxTracksPerEvent = 500; // Maximum number of tracks per event, adjust as needed
            static constexpr int maxPairs = maxTracksPerEvent * (maxTracksPerEvent - 1) / 2;

            // stuff to allocate on the gpu
            double* d_massArray;
            double* d_lifetimeArray;
            int* d_massCounter;
            double* d_primaryVertex;

            double* d_flattenedTrackData;
            double** d_trackData;

            //stuff to allocate on the CPU
            double* cPrimaryVertex;
        
        public:
            // Constructor
            K0SAnalysis(double pTcutoff, double TrackImpactParametercutoff, double Lxycutoff, double ImpactParametercutoff);
        
            // Overridden methods from Analysis
            void start() override;
            void stop() override;
            void event(Event* ev) override;
        
        private:
            // Helper method to check overlaps
            //std::vector<std::shared_ptr<Vertex>> vertex_checkOverlaps(const Helix& helicalTrack, Event* ev);
        };    
} // namespace Cdf


