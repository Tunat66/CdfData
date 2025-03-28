#include "Helix.cuh"
#define M_PI 3.14159265358979323846


namespace Cdf
{

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


} // namespace cdf


