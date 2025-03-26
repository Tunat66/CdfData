#pragma once

#include <vector>
#include <array>
#include <cmath>
#include <complex>
#include "Helix.hh"
#include "Event.hh"
#include "Track.hh"


namespace Cdf
{
    class Vertex {
        private:
            static constexpr double kpc = 0.002116;
            static constexpr double c = 2.99792458; // Speed of light in cgs units
            std::array<double, 4> vertex;           // Vertex properties: x, y, zbar, deltaz
            std::array<double, 2> primaryVertex;    // Primary vertex (proton-antiproton location)
            Track trk1;                          // First track
            Track trk2;                          // Second track
            bool broken_flag = false;               // Indicates if the vertex is invalid
        
        public:
            // Constructor
            Vertex(const Helix& helix1, const Helix& helix2, Event* event);
        
            // Accessors
            bool isBroken() const { return broken_flag; }
            std::array<double, 4> getVertex() const { return vertex; }
            std::array<double, 2> getPrimaryVertex() const { return primaryVertex; }
        
            // Methods
            double pt(const Track& trk) const;
            std::array<double, 3> p(const Track& trk) const;
            double pT() const;
            double dPV_trk(const Track& trk) const;
            std::array<double, 2> dPV() const;
            double mass(double mass1, double mass2) const;
            double ImpactParameter() const;
            double Lxy() const;
            double lifetime(double m) const; // Returns lifetime in picoseconds
    };
}