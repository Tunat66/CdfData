#include "Event.hh"

namespace Cdf
{
    Event::Event()
    {
        //add stuff to initialize this is not really needed
    }

    bool Event::isValid()
    {
        return (runNumber >= 0);
    }
}