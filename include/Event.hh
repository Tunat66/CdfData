#pragma once
#include "Track.hh"
#include <array>

namespace Cdf
{
    class Event 
    {
        public:
            Event();
            bool isValid();
            std::vector<Track> getTracks() const {return tracks;}
            std::array <double, 2> getVertex() const {return vertex;}
        //DataFile is friend class and can access
        friend class DataFile;
        friend class Vertex;
        protected:
            int runNumber = -1;
            int eventNumber = -1;
            std::array <double, 2> vertex = {0,0}; //transverse location xy
            std::vector<Track> tracks; //tracks associated with event
            std::string trim(std::string line);
    };
}

