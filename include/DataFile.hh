#pragma once
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include "Event.hh"

namespace Cdf
{
    class DataFile
    {
        
        public:
            DataFile(std::string filename);
            ~DataFile();
            bool skipLogic(std::string line);
            bool skipLogic();
            const int BundleSize = 256; // 1 event bundle contains 256 events
            std::vector<Event*> next_eventBundle();
            void rewind();
            std::string trim(std::string line);
            Event *next();

        protected:
            int id = -1; //initialize an invalid value
        private:
            std::ifstream fileStream; // Declare fileStream as a member variable    
    };
}