
//
// defines
//

#define EARTH_MASS         5.972e24         // mass of earth, kg
#define EARTH_RADIUS       6.371e6          // meters
#define MOON_MASS          7.34767309e22    // mass of moon, kg
#define DIST_EARTH_MOON    3.84400e8        // distance between earth and moon, meters
#define SUN_MASS           1.989e30         // mass of sun, kg
//#define DIST_EARTH_SUN     151.75e9         // distance between earth and sun, meters
#define DIST_EARTH_SUN     149597870000.
#define G                  6.67408e-11      // gravitational constant
#define EARTH_GRAVITY      9.81             // gravity of earth, m/s/s

#define METERS_TO_MILES(m)          ((m) * 0.000621371)
#define METERS_TO_NAUTICAL_MILES(m) ((m) * 0.000539957)
#define METERS_TO_FEET(m)           ((m) * 3.28084)
#define TWO_PI                      (2 * M_PI)  // xxx is this used enough
#define DEG_TO_RAD(d)               ((d) * (M_PI/180))
#define RAD_TO_DEG(r)               ((r) * (180/M_PI))

#define GEO_DISK   0
#define GEO_SPHERE 1

//
// variables
//

struct {
    double x;
    double y;
    double z;
    double es_orbit_radius;  // xxx add code
    double es_orbit_w;  // xxx add code
} sun;

struct {
    double x;
    double y;
    double z;
    double em_orbit_radius;
    double em_orbit_w;
} moon;

#define MAX_EARTH_SURFACE  50000  // xxx what is needed
struct {
    double x;
    double y;
    double z;
    double em_orbit_radius;
    double em_orbit_w;
    double es_orbit_radius;
    double es_orbit_w;
    struct {
        double x;  // xxx comments
        double y;
        double z;
        double g;
        double r;
        vector_t v;
    } surface[MAX_EARTH_SURFACE];
    int max_surface;
} earth;

struct ctrls_s {
    volatile double theta;
    volatile bool   moon_enabled;
    volatile bool   sun_enabled;

    volatile bool   motion;
    volatile bool   vectors;
    volatile int    geometry;
} ctrls;

//
// prototypes
//

void tides_init(void);
void tides_get_min_max(double *min, double *max, int *min_idx, int *max_idx);
