// XXX
// - stop using long doulbe

#include <stdio.h>
#include <stdbool.h>
#include <math.h>

// Purpose: To calculate the location and height of the Earth's tidal bulges,
//          caused by the Moon.

// Assumptions:  XXX review
// - the Earth to Moon distance is constant, the orbits are circular
// - calculations are done using 2 dimensions.
// - earth is circular.
// - 24 hour rotation of the Earth is not considered.
// - contributuion to the tidal bulges by the Sun is ignored.

// Diagram
// - the X,Y origin is at the center of mass, 
//   both the Earth and the Moon have circular orbits about this point
//                           Y             
//                           |
//                  xxxxxxxxx|x
//                 x         | x
//                x          |  x
//               x           |   x
//    X  ------  x       *   O   x  ---------------------------------------- Moon
//               x           |   x
//                x          |  x
//                 x         | x
//                  xxxxxxxxx|x
//                           | 
//                       ^                                                      ^
//                    x.earth                                                 moon.x

//
// defines
//

#define EARTH_MASS         5.972e24         // mass of earth, kg
#define EARTH_RADIUS       6.371e6          // meters
#define MOON_MASS          7.34767309e22    // mass of moon, kg
#define DIST_EARTH_MOON    3.84400e8        // distance between earth and moon, meters
#define G                  6.67408e-11      // gravitational constant
#define EARTH_GRAVITY      9.81             // gravity of earth, m/s/s
#define OCEAN_DEPTH        3682.            // average depth of oceans, meters

#define METERS_TO_MILES(m) ((m) * 0.000621371)
#define TWO_PI             (2 * M_PI)
#define DEG_TO_RAD(d)      ((d) * (M_PI/180))
#define RAD_TO_DEG(r)      ((r) * (180/M_PI))

#define DELTA_H            1e-5
#define DELTA_M            1.0

//
// typedefs
//

struct {
    long double x;
    long double y;
    long double mass;
    long double w;
} moon;

struct {
    long double x;
    long double y;
    long double mass;
    long double radius;
    long double w;
    long double g[360];
    long double h[360];
    long double tpe;
} earth;

typedef struct {
    long double a;
    long double b;
} vector_t;

//
// prototypes
//

static long double square(long double x);
static long double magnitude(vector_t *v);
static int set_vector_magnitude(vector_t *v, long double new_magnitude);

// ---------------------------------------------------

int main(int argc, char **argv)
{
    // Init constant fields of the earth and moon structs.
    earth.mass   = EARTH_MASS;
    moon.mass    = MOON_MASS;
    earth.radius = EARTH_RADIUS;
    for (int i = 0; i < 360; i++) {
        earth.h[i] = OCEAN_DEPTH;;
    }

    // Determine the position of the earth and moon relative to the center of mass.
    // Refer to the diagram above, the Earth and Moon have no velocity in the X direction, 
    // so the centrifugal force of the Moon must equal that of the Earth.
    //   moon.mass * moon.w^2 * moon.x = earth.mass * earth.w^2 * earth.x
    earth.x = -(DIST_EARTH_MOON / (1 + earth.mass / moon.mass));
    earth.y = 0;
    moon.x  = DIST_EARTH_MOON + earth.x;
    moon.y  = 0;
    printf("DIST_EARTH_MOON = %8.0f miles\n", METERS_TO_MILES(DIST_EARTH_MOON));
    printf("earth.x         = %8.0Lf miles\n", METERS_TO_MILES(earth.x));
    printf("moon.x          = %8.0Lf miles\n", METERS_TO_MILES(moon.x));

    // Determine the angular velocity of both the earth and the moon; they should be the same.
    // The centrifugal force must equal the gravitational force for a circular orbit
    //   gravity_force = moon.mass * moon.w^2 * moon.x
    long double gravity_force = G * earth.mass * moon.mass / (DIST_EARTH_MOON * DIST_EARTH_MOON);
    moon.w  = sqrt(gravity_force / (moon.mass * fabs(moon.x)));
    earth.w = sqrt(gravity_force / (earth.mass * fabs(earth.x)));
    // xxx sanity check they are within .0000001
    printf("moon.w         = %Le\n", moon.w);
    printf("earth.w        = %Le\n", earth.w);
    printf("orbital period = %Lf days\n", TWO_PI / moon.w / 86400);
    printf("\n");

    // xxx comment
    for (int deg = 0; deg < 360; deg++) {
        long double x, y, d;
        vector_t g, m, c, t;

        x = earth.x + earth.radius * cos(DEG_TO_RAD(deg));
        y = earth.y + earth.radius * sin(DEG_TO_RAD(deg));

        g.a = earth.x - x;
        g.b = earth.y - y;
        set_vector_magnitude(&g, EARTH_GRAVITY);

        m.a = moon.x - x;
        m.b = moon.y - y;
        d = magnitude(&m);
        set_vector_magnitude(&m, G * moon.mass / square(d));

        c.a = -1;
        c.b = 0;
        set_vector_magnitude(&c,  square(earth.w) * -earth.x);

        t.a = g.a + m.a + c.a;
        t.b = g.b + m.b + c.b;
        earth.g[deg] = magnitude(&t);

        if (deg == 0 || deg == 90 || deg == 180 || deg == 270) {
            printf("%3d: earth gravity    = %0.9Lf %0.9Lf\n", deg, g.a, g.b);
            printf("     moon gravity     = %0.9Lf %0.9Lf\n", m.a, m.b);
            printf("     centrigual accel = %0.9Lf %0.9Lf\n", c.a, c.b);
            printf("     total            = %0.9Lf %0.9Lf\n", t.a, t.b);
            printf("     MAGNITUDE        = %0.9Lf\n", earth.g[deg]);
            printf("\n");
        }
    }

    // xxx
    printf("Running ...\n");
    int loops = 0;
    while (true) {
        int num_exchanges = 0;
        for (int i = 0; i < 360; i++) {
            int best_add_idx=-1, best_sub_idx=-1;
            long double min_dpe=1e99, best_dpe=0;
            for (int j = 0; j < 360; j++) {
                long double dpe;
                dpe = (DELTA_M * earth.g[i] * (earth.h[i] + DELTA_H/2)) -
                      (DELTA_M * earth.g[j] * (earth.h[j] - DELTA_H/2));
                if (dpe < 0 && dpe < min_dpe) {
                    best_add_idx = i;
                    best_sub_idx = j;
                    best_dpe = dpe;
                    min_dpe = dpe;
                }
            }
            if (best_add_idx != -1) {
                earth.h[best_add_idx] += DELTA_H;
                earth.h[best_sub_idx] -= DELTA_H;
                earth.tpe += best_dpe;  // xxx  not used ?
                num_exchanges++;
            }
        }
        if ((++loops % 100) == 0) {
            printf("   loops=%d  numex=%d\n", loops, num_exchanges);
        }
        if (num_exchanges <= 2) {
            break;
        }
    }
    printf("   final loops = %d\n", loops);
    printf("\n");
  
    // print results
    long double minh=1e99, maxh=-1e99;
    for (int i = 0; i < 360; i++) {
        long double h = earth.h[i];
        if (i == 0 || i == 90 || i == 180 || i == 270) {
            printf("deg = %3d   height= %0.10Lf\n", i, h);
        }
        if (h < minh) minh = h;
        if (h > maxh) maxh = h;
    }
    printf("\n");
    printf("min = %LF max = %LF  Range = %LF\n", minh, maxh, maxh-minh);

    // done
    return 0;
}

// -----------------  UTILS  ----------------------------------------

static long double square(long double x)
{
    return x * x;
}

static long double magnitude(vector_t *v)
{
    return sqrt( square(v->a) + square(v->b) );
}

static int set_vector_magnitude(vector_t *v, long double new_magnitude)
{
    long double current_magnitude, factor;

    current_magnitude = magnitude(v);
    if (current_magnitude == 0) {
        return -1;  // xxx or abort
    }

    factor = new_magnitude / current_magnitude;

    v->a *= factor;
    v->b *= factor;

    return 0;
}
