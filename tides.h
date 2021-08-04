
//
// defines
//

#define EARTH_MASS         5.972e24         // mass of earth, kg
#define EARTH_RADIUS       6.371e6          // meters
#define MOON_MASS          7.34767309e22    // mass of moon, kg
#define DIST_EARTH_MOON    3.84400e8        // distance between earth and moon, meters
#define G                  6.67408e-11      // gravitational constant
#define EARTH_GRAVITY      9.81             // gravity of earth, m/s/s

#define METERS_TO_MILES(m)          ((m) * 0.000621371)
#define METERS_TO_NAUTICAL_MILES(m) ((m) * 0.000539957)
#define TWO_PI                      (2 * M_PI)
#define DEG_TO_RAD(d)               ((d) * (M_PI/180))
#define RAD_TO_DEG(r)               ((r) * (180/M_PI))

//
// variables
//

struct {
    double x;
    double y;
    double z;
    double orbit_radius;
    double mass;
    double w;
} moon;

#define MAX_EARTH_SURFACE  50000  // xxx what is needed
struct {
    double x;
    double y;
    double z;
    double mass;
    double radius;
    double orbit_radius;
    double w;
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

volatile bool reset_request;
volatile bool motion;
volatile bool vectors;

//
// prototypes
//

void tides_init(void);
