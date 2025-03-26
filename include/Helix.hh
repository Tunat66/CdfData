#pragma once
#include <complex>
#include <vector>
#include "Track.hh" // Assuming you have a CdfTrack class
namespace Cdf
{
    class Helix {
        private:
            Track trk; // A CdfTrack object
        
        public:
            // Constructor
            Helix(const Track& t);
        
            // Getter for the track
            Track getTrk() const;
        
            // Calculate the radius
            double radius(double curvature = -1) const;
        
            // Calculate the center
            std::complex<double> center(double h = -1, double d0 = -1, double phi0 = -1) const;
        
            // Parametrize the helix using the points method
            std::vector<double> points(double t, int q = 1) const;
        
            // Find z from xy
            double findZFromXY(const std::complex<double>& xy, int index = 0, int q = 1) const;

            bool operator==(const Helix &other) const;

            // Check for intersection with another helix
            std::vector<std::vector<double>> intersect(const Helix& h) const;
    };
}
