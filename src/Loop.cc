#include "Loop.hh"
#include <iostream>
#include <memory>
#include "DataFile.hh"
#include "Analysis.hh"

namespace Cdf
{
// Constructor
Loop::Loop(const std::string& filename, int nevmax, int report)
    : filename(filename), nevmax(nevmax), report(report) {}

// Main loop function
void Loop::run(const std::vector<std::shared_ptr<Analysis>>& analyses) {
    // Check that all the analyses are valid
    for (const auto& analysis : analyses) {
        if (!analysis) {
            std::cerr << "Analysis list contains null or invalid objects." << std::endl;
            return;
        }
    }

    // Open the data file
    DataFile file(filename);

    // Start at the beginning of the data file
    file.rewind();

    // Call start() on each Analysis
    for (const auto& analysis : analyses) {
        analysis->start();
    }
    std::cout << "Begin analyses" << std::endl;

    // Read and analyze events one at a time
    int nev = 0; // Count events
    
    while (true) {
        Event* ev = file.next(); 
        if(ev == nullptr)
        {
            continue;
        }
        if (ev->isValid()) { 
            for (const auto& analysis : analyses) {
                analysis->event(ev);
            }
            nev++;
            if (nev % report == 0) {
                std::cout << nev << " events processed" << std::endl;
            }
            if (nev >= nevmax) {
                delete ev; //delete before exiting to prevent a memory leak
                break;
            }
        } else {
            continue;
        }
        delete ev;
    }

    // Call stop() on each Analysis
    for (const auto& analysis : analyses) {
        analysis->stop();
    }

    std::cout << "End analysis: " << nev << " events processed" << std::endl;
}
}

