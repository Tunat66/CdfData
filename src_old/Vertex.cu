#include "Vertex.cuh"
#include <math.h>
#include <thrust/complex.h>
#include <stdbool.h>
#include <stddef.h>

namespace Cdf
{

__host__ __device__ plain_vertex return_vertex(plain_track trk1, plain_track trk2)
{
    plain_vertex vertex = {0, 0, 0, 0, false};
    double intersectionPoints[2][4];
    int numPoints = 0;
    double** intersections = intersect(trk1, trk2, &numPoints);

    if (numPoints > 0) {
        for (int i = 0; i < numPoints; ++i) {
            for (int j = 0; j < 4; ++j) {
                intersectionPoints[i][j] = intersections[i][j];
            }
            free(intersections[i]);
        }
        free(intersections);
    }

    if (numPoints == 0) {
        return vertex; // vertex by default has isValid = false
    }

    double deltaz1 = intersectionPoints[0][3];
    double deltaz2 = intersectionPoints[1][3];
    int choice = (fabs(deltaz1) < fabs(deltaz2)) ? 0 : 1;

    if (intersectionPoints[choice][3] > 10) {
        return vertex;
    }

    vertex.x = intersectionPoints[choice][0];
    vertex.y = intersectionPoints[choice][1];
    vertex.zbar = intersectionPoints[choice][2];
    vertex.deltaz = intersectionPoints[choice][3];
    vertex.isValid = true;

    return vertex;
}

__host__ __device__ double pt(plain_track trk) {
    return Vertex::kpc / trk.curvature;
}

__host__ __device__ void p(plain_track trk, plain_vertex vertex, double* result) {
    double h = trk.curvature;
    double d0 = trk.d0;
    double phi0 = trk.phi0;
    double x = vertex.x;
    double y = vertex.y;

    double ptrans = fabs(pt(trk));
    double pz = ptrans * trk.cotTheta;
    double py = ptrans * ((1 + 2 * h * d0) * sin(phi0) + 2 * h * x);
    double px = ptrans * ((1 + 2 * h * d0) * cos(phi0) - 2 * h * y);

    result[0] = px;
    result[1] = py;
    result[2] = pz;
}

__host__ __device__ double pT(plain_track trk1, plain_track trk2) {
    plain_vertex vertex = return_vertex(trk1, trk2);
    double mom1[3], mom2[3];
    p(trk1, vertex, mom1);
    p(trk2, vertex, mom2);

    double px = mom1[0] + mom2[0];
    double py = mom1[1] + mom2[1];
    return sqrt(px * px + py * py);
}

__host__ __device__ double dPV_trk(plain_track trk, double primaryVertex[2]) {
    thrust::complex<double> wPV(primaryVertex[0], primaryVertex[1]);
    double centerX, centerY;
    center(trk.curvature, trk.d0, trk.phi0, &centerX, &centerY);
    thrust::complex<double> center(centerX, centerY);
    return thrust::abs(center - wPV) - fabs(radius(trk.curvature));
}

__host__ __device__ double mass(double mass1, double mass2, plain_track trk1, plain_track trk2, plain_vertex vertex) {
    double p1[3], p2[3];
    p(trk1, vertex, p1);
    p(trk2, vertex, p2);

    double e1 = sqrt(p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2] + mass1 * mass1);
    double e2 = sqrt(p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2] + mass2 * mass2);

    double px = p1[0] + p2[0];
    double py = p1[1] + p2[1];
    double pz = p1[2] + p2[2];
    double pTot2 = px * px + py * py + pz * pz;

    return sqrt((e1 + e2) * (e1 + e2) - pTot2);
}

__host__ __device__ double ImpactParameter(plain_track trk1, plain_track trk2, plain_vertex vertex, double primaryVertex[2]) {
    double displacement[2] = {vertex.x - primaryVertex[0], vertex.y - primaryVertex[1]};
    double mom1[3], mom2[3];
    p(trk1, vertex, mom1);
    p(trk2, vertex, mom2);

    double px = mom1[0] + mom2[0];
    double py = mom1[1] + mom2[1];
    double pmag = sqrt(px * px + py * py);
    double ptrans[2] = {px / pmag, py / pmag};

    return displacement[0] * ptrans[1] - displacement[1] * ptrans[0];
}

__host__ __device__ double Lxy(plain_track trk1, plain_track trk2, plain_vertex vertex, double primaryVertex[2]) {
    double displacement[2] = {vertex.x - primaryVertex[0], vertex.y - primaryVertex[1]};
    double mom1[3], mom2[3];
    p(trk1, vertex, mom1);
    p(trk2, vertex, mom2);

    double px = mom1[0] + mom2[0];
    double py = mom1[1] + mom2[1];
    double pmag = sqrt(px * px + py * py);
    double ptrans[2] = {px / pmag, py / pmag};

    return fabs(displacement[0] * ptrans[0] + displacement[1] * ptrans[1]);
}

__host__ __device__ double lifetime(double m, plain_track trk1, plain_track trk2, plain_vertex vertex, double primaryVertex[2]) {
    double Lxy_val = Lxy(trk1, trk2, vertex, primaryVertex);
    double mom1[3], mom2[3];
    p(trk1, vertex, mom1);
    p(trk2, vertex, mom2);

    double px = mom1[0] + mom2[0];
    double py = mom1[1] + mom2[1];
    double pmag = sqrt(px * px + py * py);

    return (Lxy_val * m / (pmag * Vertex::c)) * 1e2; // Convert to picoseconds
}

} // namespace Cdf