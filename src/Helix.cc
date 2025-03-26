#include "Helix.hh"
#define M_PI 3.14159265358979323846
#include <cmath>
#include <stdexcept>

//for debugging 
#include <iostream>

namespace Cdf
{
// Constructor
Helix::Helix(const Track& t) : trk(t) {}

// Getter for the track
Track Helix::getTrk() const {
    return trk;
}

// Calculate the radius
double Helix::radius(double curvature) const {
    if (curvature < 0) {
        curvature = trk.curvature;
    }
    return 0.5 / curvature;
}

// Calculate the center
std::complex<double> Helix::center(double h, double d0, double phi0) const {
    if (h < 0) h = trk.curvature;
    if (d0 < 0) d0 = trk.d0;
    if (phi0 < 0) phi0 = trk.phi0;

    double r = radius(h);
    return (r + d0) * std::polar(1.0, phi0 + M_PI / 2);
}

// Parametrize the helix using the points method
std::vector<double> Helix::points(double t, int q) const {
    double r = radius();
    std::complex<double> c = center();
    std::complex<double> radial_point = c - std::complex<double>(0, 1) * r * std::polar(1.0, trk.phi0 + q * t);

    return {
        std::real(radial_point),
        std::imag(radial_point),
        trk.z0 + std::abs(r) * t * trk.cotTheta
    };
}

// Find z from xy
double Helix::findZFromXY(const std::complex<double>& xy, int index, int q) const {
    double r = radius();
    std::complex<double> c = center();
    double t = (-std::imag(std::log(std::complex<double>(0, 1) * (xy - c) / r)) - trk.phi0) / q + 2 * M_PI * index;
    return trk.z0 + std::abs(r) * t * trk.cotTheta;
}

//overloaded the == operator to check if two helix objects are identical
bool Helix::operator==(const Helix& other) const {
    // Compare the track parameters
    return trk.curvature == other.trk.curvature &&
           trk.d0 == other.trk.d0 &&
           trk.phi0 == other.trk.phi0 &&
           trk.z0 == other.trk.z0 &&
           trk.cotTheta == other.trk.cotTheta;
}

// Check for intersection with another helix
std::vector<std::vector<double>> Helix::intersect(const Helix& h) const {
    std::complex<double> c1 = center();
    std::complex<double> c2 = h.center();
    std::complex<double> d = c2 - c1;
    double d2 = std::norm(d); // Square of distance between centers, norm returns square of magnitude
    double r1 = radius();
    double r2 = h.radius();

    if (d2 > std::pow(std::abs(r1) + std::abs(r2), 2)) {
        return {}; // No intersection
    }

    if(*this == h)
    {
        return {}; // Helices are identical, no intersection
    }

    // Calculate transverse positions
    std::complex<double> cos_alpha = (d2 + (r1 * r1) - (r2 * r2)) / (2.0 * std::abs(d) * std::abs(r1));
    std::complex<double> alpha = std::acos(cos_alpha); // Use complex acos to handle edge cases

    if(std::isnan(std::real(alpha)))
    {
        std::cout << "alpha is nan" << std::endl;
        std::cout << "cos(alpha):" << (d2 + r1 * r1 - r2 * r2) / (2 * std::abs(d) * std::abs(r1)) << std::endl;
        std::cout << "Alpha: " << alpha << std::endl;
    }
    double Delta = std::atan2(std::imag(d), std::real(d));
    std::complex<double> cplus = c1 + std::abs(r1) * std::exp(-1 * std::imag(alpha)) * std::polar(1.0, Delta + std::real(alpha));
    std::complex<double> cminus = c1 + std::abs(r1) * std::exp(-1 * std::imag(alpha)) * std::polar(1.0, Delta - std::real(alpha));

    // Calculate z positions
    double z1_plus = findZFromXY(cplus);
    double z1_minus = findZFromXY(cminus);
    double z2_plus = h.findZFromXY(cplus);
    double z2_minus = h.findZFromXY(cminus);

    double zbar_plus = (z1_plus + z2_plus) / 2;
    double zbar_minus = (z1_minus + z2_minus) / 2;
    double deltaz_plus = z1_plus - z2_plus;
    double deltaz_minus = z1_minus - z2_minus;

    // Format the output
    return {
        {std::real(cplus), std::imag(cplus), zbar_plus, deltaz_plus},
        {std::real(cminus), std::imag(cminus), zbar_minus, deltaz_minus}
    };
}
} // namespace cdf


