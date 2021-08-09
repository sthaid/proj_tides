#ifndef __TIDES_H__
#define __TIDES_H__

//
// defines
//

#define EARTH_MASS         5.972e24         // mass of earth, kg
#define EARTH_RADIUS       6.371e6          // meters
#define MOON_MASS          7.34767309e22    // mass of moon, kg
#define DIST_EARTH_MOON    3.84400e8        // distance between earth and moon, meters
#define SUN_MASS           1.989e30         // mass of sun, kg
#define DIST_EARTH_SUN     149597870000.    // distance between earth and sun, meters
#define G                  6.67408e-11      // gravitational constant
#define EARTH_GRAVITY      9.81             // gravity of earth, m/s/s

#define METERS_TO_MILES(m)          ((m) * 0.000621371)
#define METERS_TO_NAUTICAL_MILES(m) ((m) * 0.000539957)
#define METERS_TO_FEET(m)           ((m) * 3.28084)
#define TWO_PI                      (2 * M_PI)
#define DEG_TO_RAD(d)               ((d) * (M_PI/180))
#define RAD_TO_DEG(r)               ((r) * (180/M_PI))

#define GEO_DISK   0
#define GEO_SPHERE 1

//
// variables
//

// notes:
// - unless otherwise noted x,y,z are relative to the earth-moon barycenter
// - orbits are circular, and in the x,y plane
// - orbit_radius values are the distance to the barycenter of the orbit

struct {
    double x;
    double y;
    double z;
    double es_orbit_radius;
    double es_orbit_w;
} sun;

struct {
    double x;
    double y;
    double z;
    double em_orbit_radius;
    double em_orbit_w;
} moon;

#define MAX_EARTH_SURFACE  50000  // actual max needed is 41258
struct {
    double x;
    double y;
    double z;
    double em_orbit_radius;
    double em_orbit_w;
    double es_orbit_radius;
    double es_orbit_w;
    struct {
        double x;    // x,y,z are relative to the center of the earth
        double y;
        double z;
        double gt;   // the total accel at this location on the earth's surface
        double R;    // the distance to earth center at this location, including the tide level
        vector_t v;  // sum of all accels at this location, excluding earth's gravity
    } surface[MAX_EARTH_SURFACE];
    int max_surface;
} earth;

struct ctrls_s {
    volatile double theta;         // angle position of the moon, relative to x,y origin
    volatile bool   moon_enabled;  // when true, the effect of the moon is included in the simulation
    volatile bool   sun_enabled;   // when true, the effect of the sun is included in the simulation

    volatile bool   motion;        // moon motion is simulatied
    volatile bool   vectors;       // display acceleration vectors (excluding earth gravity) at equator
    volatile int    geometry;      // select simulating the earth as a disk or sphere
} ctrls;

//
// prototypes
//

void tides_init(void);
void tides_get_min_max(double *min, double *max, int *min_idx, int *max_idx);

#endif
