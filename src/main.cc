#include "K0SAnalysis.hh"
#include "Loop.hh"

int main()
{
    //define some cutoff variables
    double pTcutoff = 1.5;
    double TrackImpactParametercutoff = 0.3;
    double Lxycutoff = 2;
    double ImpactParametercutoff = 0.1;
    //create an instance of the K0SAnalysis class
    auto analysis = std::make_shared<Cdf::K0SAnalysis>(pTcutoff, TrackImpactParametercutoff, Lxycutoff, ImpactParametercutoff);

    //create a Loop object
    Cdf::Loop loop("cdf.dat", 1000, 100);
    loop.run({analysis});


    return 0;
}