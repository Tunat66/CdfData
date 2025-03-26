#pragma once
#include <string>
#include <vector>
#include <memory>
#include "DataFile.hh" // Assuming you have a CdfDataFile class
#include "Analysis.hh"    // Assuming you have an Analysis base class

namespace Cdf
{
    class Loop {
        private:
            std::string filename; // Name of the data file
            int nevmax;           // Maximum number of events to process
            int report;           // Number of events between each console report
        
        public:
            // Constructor
            Loop(const std::string& filename, int nevmax = 1000, int report = 10);
        
            // Main loop function
            void run(const std::vector<std::shared_ptr<Analysis>>& analyses);
        };
} // namespace Cdf

