#pragma once
#include "Event.hh"
namespace Cdf
{
    class Analysis
    {
    public:
        virtual ~Analysis() = default; // Virtual destructor for proper cleanup of derived classes

        // Pure virtual functions
        virtual void start() = 0; // To be implemented by derived classes
        virtual void stop() = 0;  // To be implemented by derived classes
        virtual void event(Event* ev) = 0; // To be implemented by derived classes
    };
}