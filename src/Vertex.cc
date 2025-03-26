#include "Vertex.hh"
#include <stdexcept>
#include <algorithm>
#include <numeric>

//for debugging
#include <iostream>

namespace Cdf
{
// Constructor
Vertex::Vertex(const Helix& helix1, const Helix& helix2, Event* event) {
    // Get intersection points
    auto intersectionPoints = helix1.intersect(helix2);
    if (intersectionPoints.empty()) {
        broken_flag = true;
        return;
    }
    // Choose the intersection point with the smallest deltaz
    double deltaz1 = intersectionPoints[0][3];
    double deltaz2 = intersectionPoints[1][3];
    int choice = (std::abs(deltaz1) < std::abs(deltaz2)) ? 0 : 1;

    // Check for invalid intersection points
    //if (std::any_of(intersectionPoints[choice].begin(), intersectionPoints[choice].end(),
    //                [](double val) { return std::abs(val) > 1e2; })) {
    //    broken_flag = true;
    //    return;
    //}

    if (intersectionPoints[choice][3] > 10) {
        broken_flag = true;
        return;
    }

    // Set the vertex and tracks
    vertex = {intersectionPoints[choice][0], intersectionPoints[choice][1],
              intersectionPoints[choice][2], intersectionPoints[choice][3]};
    trk1 = helix1.getTrk();
    trk2 = helix2.getTrk();

    // Extract the primary vertex from the event
    primaryVertex = event->vertex;
}

// Calculate transverse momentum (pT) for a track
double Vertex::pt(const Track& trk) const {
    return kpc / trk.curvature;
}

// Calculate 3-momentum for a track
std::array<double, 3> Vertex::p(const Track& trk) const {
    double h = trk.curvature;
    double d0 = trk.d0;
    double phi0 = trk.phi0;
    double x = vertex[0];
    double y = vertex[1];

    double ptrans = std::abs(pt(trk));
    double pz = ptrans * trk.cotTheta;
    double py = ptrans * ((1 + 2 * h * d0) * std::sin(phi0) + 2 * h * x);
    double px = ptrans * ((1 + 2 * h * d0) * std::cos(phi0) - 2 * h * y);

    return {px, py, pz};
}

// Calculate total transverse momentum (pT) for the vertex
double Vertex::pT() const {
    auto mom1 = p(trk1);
    auto mom2 = p(trk2);
    double px = mom1[0] + mom2[0];
    double py = mom1[1] + mom2[1];
    return std::sqrt(px * px + py * py);
}

// Calculate impact parameter in the transverse plane for a track
double Vertex::dPV_trk(const Track& trk) const {
    Helix helix(trk);
    std::complex<double> wPV(primaryVertex[0], primaryVertex[1]);
    return std::abs(helix.center() - wPV) - std::abs(helix.radius());
}

// Calculate impact parameters for both tracks
std::array<double, 2> Vertex::dPV() const {
    return {dPV_trk(trk1), dPV_trk(trk2)};
}

// Calculate invariant mass of the vertex
double Vertex::mass(double mass1, double mass2) const {
    auto p1 = p(trk1);
    auto p2 = p(trk2);

    double e1 = std::sqrt(std::inner_product(p1.begin(), p1.end(), p1.begin(), 0.0) + mass1 * mass1);
    double e2 = std::sqrt(std::inner_product(p2.begin(), p2.end(), p2.begin(), 0.0) + mass2 * mass2);

    double px = p1[0] + p2[0];
    double py = p1[1] + p2[1];
    double pz = p1[2] + p2[2];
    double pTot2 = px * px + py * py + pz * pz;
    return std::sqrt((e1 + e2) * (e1 + e2) - pTot2);
}

// Calculate impact parameter of the vertex
double Vertex::ImpactParameter() const {
    std::array<double, 2> displacement = {vertex[0] - primaryVertex[0], vertex[1] - primaryVertex[1]};
    auto mom = p(trk1);
    auto mom2 = p(trk2);
    double px = mom[0] + mom2[0];
    double py = mom[1] + mom2[1];
    double pmag = std::sqrt(px * px + py * py);
    std::array<double, 2> ptrans = {px / pmag, py / pmag};

    return displacement[0] * ptrans[1] - displacement[1] * ptrans[0];
}

// Calculate Lxy (transverse displacement)
double Vertex::Lxy() const {
    std::array<double, 2> displacement = {vertex[0] - primaryVertex[0], vertex[1] - primaryVertex[1]};
    auto mom = p(trk1);
    auto mom2 = p(trk2);
    double px = mom[0] + mom2[0];
    double py = mom[1] + mom2[1];
    double pmag = std::sqrt(px * px + py * py);
    std::array<double, 2> ptrans = {px / pmag, py / pmag};

    return std::abs(displacement[0] * ptrans[0] + displacement[1] * ptrans[1]);
}

// Calculate lifetime in picoseconds
double Vertex::lifetime(double m) const {
    double Lxy_val = Lxy();
    auto mom = p(trk1);
    auto mom2 = p(trk2);
    double px = mom[0] + mom2[0];
    double py = mom[1] + mom2[1];
    double pmag = std::sqrt(px * px + py * py);

    return (Lxy_val * m / (pmag * c)) * 1e2; // Convert to picoseconds
}
}
