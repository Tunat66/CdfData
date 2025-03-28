#include "CUDATrackKernel.cuh"

//cutoff constants copied over from K0SAnalysis
__constant__ double ImpactParametercutoff = 0.1; // cm
__constant__ double pTcutoff = 1.5; // GeV/c
__constant__ double TrackImpactParametercutoff = 0.3; // cm
__constant__ double Lxycutoff = 2; // cm

namespace Cdf
{
//HELIX CLASS METHODS: these were once part of a helix class
// Calculate the radius
__host__ __device__ double radius(double curvature) {
    return 0.5 / curvature;
}

// Calculate the center
__host__ __device__ void center(double h, double d0, double phi0, double* real, double* imag) {
    double r = radius(h);
    *real = (r + d0) * cos(phi0 + M_PI / 2);
    *imag = (r + d0) * sin(phi0 + M_PI / 2);
}

__host__ __device__ double findZFromXY(plain_track trk, double x, double y, int index, int q) {
    double r = radius(trk.curvature);
    double cx, cy;
    center(trk.curvature, trk.d0, trk.phi0, &cx, &cy);
    double dx = x - cx;
    double dy = y - cy;
    double t = (-atan2(dy, dx) - trk.phi0) / q + 2 * M_PI * index;
    return trk.z0 + fabs(r) * t * trk.cotTheta;
}

// Check for intersection with another helix
__host__ __device__ double** intersect(plain_track trk1, plain_track trk2, int* num_intersections) {
    double c1x, c1y, c2x, c2y;
    center(trk1.curvature, trk1.d0, trk1.phi0, &c1x, &c1y);
    center(trk2.curvature, trk2.d0, trk2.phi0, &c2x, &c2y);
    double dx = c2x - c1x;
    double dy = c2y - c1y;
    double d2 = dx * dx + dy * dy;
    double r1 = radius(trk1.curvature);
    double r2 = radius(trk2.curvature);

    if (d2 > (fabs(r1) + fabs(r2)) * (fabs(r1) + fabs(r2))) {
        *num_intersections = 0;
        return NULL; // No intersection
    }

    if (trk1.curvature == trk2.curvature && trk1.d0 == trk2.d0 && trk1.phi0 == trk2.phi0) {
        *num_intersections = 0;
        return NULL; // Helices are identical, no intersection
    }

    double d = sqrt(d2);
    double cos_alpha = (d2 + r1 * r1 - r2 * r2) / (2.0 * d * fabs(r1));
    if (cos_alpha < -1.0 || cos_alpha > 1.0) {
        *num_intersections = 0;
        return NULL; // No valid intersection
    }
    double alpha = acos(cos_alpha);
    double Delta = atan2(dy, dx);

    double cplus_x = c1x + fabs(r1) * cos(Delta + alpha);
    double cplus_y = c1y + fabs(r1) * sin(Delta + alpha);
    double cminus_x = c1x + fabs(r1) * cos(Delta - alpha);
    double cminus_y = c1y + fabs(r1) * sin(Delta - alpha);

    double z1_plus = findZFromXY(trk1, cplus_x, cplus_y, 0, 1);
    double z1_minus = findZFromXY(trk1, cminus_x, cminus_y, 0, 1);
    double z2_plus = findZFromXY(trk2, cplus_x, cplus_y, 0, 1);
    double z2_minus = findZFromXY(trk2, cminus_x, cminus_y, 0, 1);

    double zbar_plus = (z1_plus + z2_plus) / 2;
    double zbar_minus = (z1_minus + z2_minus) / 2;
    double deltaz_plus = z1_plus - z2_plus;
    double deltaz_minus = z1_minus - z2_minus;

    double** result = (double**)malloc(2 * sizeof(double*));
    result[0] = (double*)malloc(4 * sizeof(double));
    result[1] = (double*)malloc(4 * sizeof(double));

    result[0][0] = cplus_x;
    result[0][1] = cplus_y;
    result[0][2] = zbar_plus;
    result[0][3] = deltaz_plus;

    result[1][0] = cminus_x;
    result[1][1] = cminus_y;
    result[1][2] = zbar_minus;
    result[1][3] = deltaz_minus;

    *num_intersections = 2;
    return result;
}


//VERTEX CLASS METHODS: these were once part of a vertex class


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

// Helper function to compare tracks, runs on the device
__device__ Cdf::plain_vertex track_compare(Cdf::plain_track track_s, Cdf::plain_track track_o, double primaryVertex[2]) { 
    Cdf::plain_vertex vertex = Cdf::return_vertex(track_s, track_o);
    Cdf::plain_vertex empty_vertex;
    // Skip invalid vertices
    if (!vertex.isValid) {
        return empty_vertex;
    }
    // Apply filtering criteria
    double impactParameter = Cdf::ImpactParameter(track_s, track_o, vertex, primaryVertex);
    if (std::abs(impactParameter) > ImpactParametercutoff) {
        return empty_vertex;
    }
    double pT = Cdf::pT(track_s, track_o);
    if (std::abs(pT) < pTcutoff) {
        return empty_vertex;
    }
    double dPV_s = Cdf::dPV_trk(track_s, primaryVertex);
    double dPV_o = Cdf::dPV_trk(track_o, primaryVertex);
    if (std::abs(dPV_s) < TrackImpactParametercutoff || std::abs(dPV_o) < TrackImpactParametercutoff) {
        return empty_vertex;
    }
    double Lxy = Cdf::Lxy(track_s, track_o, vertex, primaryVertex);
    if (std::abs(Lxy) < Lxycutoff) {
        return empty_vertex;
    }

    // If we made it to this stage, return the valid vertex
    return vertex;
}

__global__ void processTracksKernel(double** trackData, int numTracks, double* primaryVertexArr, 
    double* massArray, double* lifetimeArray, int* massCounter, double m_pion) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int idy = blockIdx.y * blockDim.y + threadIdx.y;
    //number of pairs
    //int numPairs = numTracks * (numTracks - 1) / 2; //not used in this kernel, but could be useful for debugging
    //we are going to launch as many threads as the number of pairs
    if (idx >= numTracks){
        //printf("exiting thread %d\n", idx); //for debugging purposes
        return;
    } 
    if (idy >= idx){ // equals case to prevent identical tracks from being compared
        //printf("exiting thread %d\n", idy); //for debugging purposes
        return;
    }

    // Process the track
    const double* track1 = trackData[idx];
    const double* track2 = trackData[idy];
    //convert them to a plain_track structs
    Cdf::plain_track track1_p = {track1[0], track1[1], track1[2], track1[3], track1[4]}; //subject track
    Cdf::plain_track track2_p = {track2[0], track2[1], track2[2], track2[3], track2[4]}; //other track
    //convert the primaryVertex to a std::array
    double primaryVertex[2] = {primaryVertexArr[0], primaryVertexArr[1]};
    
    //compare the two tracks
    Cdf::plain_vertex vertex = track_compare(track1_p, track2_p, primaryVertex);
    if (!vertex.isValid) { return; }

    //if we have a valid track, we can now calculate the mass and lifetime
    double mass_kaon = Cdf::mass(m_pion, m_pion, track1_p, track2_p, vertex); 
    double lifetime_kaon = Cdf::lifetime(m_pion, track1_p, track2_p, vertex, primaryVertex);
    int uniqueIndex = idx * (idx - 1) / 2 + idy;

    //printf("Thread (%d, %d): uniqueIndex = %d, mass = %f, lifetime = %f, massCounter = %d\n", idx, idy, uniqueIndex, mass_kaon, lifetime_kaon, *massCounter);
    // Write results to global memory
    massArray[uniqueIndex] = mass_kaon;
    lifetimeArray[uniqueIndex] = lifetime_kaon;
    // Increment massCounter (atomic operation to avoid race conditions)
    atomicAdd(massCounter, 1);

}





