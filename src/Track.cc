#include "Track.hh"

namespace Cdf
{
    Track::Track(std::vector<double> TrackParameter)
    {
        std::vector<double> &v = TrackParameter; //create a short hand reference
        if(v.size() != 5)
        {
            throw std::invalid_argument("Track: Track constructor requires 5 arguments");    
        }
        cotTheta = v[0];
        curvature = v[1];
        d0 = v[2];
        phi0 = v[3];
        z0 = v[4];
    }

    bool Track::isValid()
    {
        return (cotTheta != -999.0);
    }
}
