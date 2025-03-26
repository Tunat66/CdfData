#pragma once
#include <vector>
#include <stdexcept>

namespace Cdf 
{
    class Track
    {
        friend class Helix; // to access the protected variables
        friend class Vertex;
        public:
            Track() = default;
            Track(std::vector<double> TrackParameters);
            bool isValid();
        protected:
            double cotTheta = -999.0;
            double curvature = -999.0;
            double d0 = -999.0;
            double phi0 = -999.0;
            double z0 = -999.0;
    };



}