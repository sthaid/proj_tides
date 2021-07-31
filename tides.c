// XXX

// - doubles
//   - cast floating point constants to long double
//   - use define for long double

// - spherical earth
//   - does this change the result

// - graphics
//   - motion on/off
//   - switch to remove the centriigul, or change its value
//   - show the vectors at the 4 locations
// - function for the pe

// - write this up in README.md
//   - include assumptions
//   - what were my miconceptions

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Purpose: To calculate the location and height of the Earth's tidal bulges,
//          caused by the Moon.

// Coordinates Diagram
// - The X,Y origin is at the center of mass.   XXX barycenter
//   Both the Earth and the Moon have circular orbits about this point
//
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
//                    earth.x                                                 moon.x

//
// defines
//

#define EARTH_MASS         5.972e24         // mass of earth, kg
#define EARTH_RADIUS       6.371e6          // meters
#define MOON_MASS          7.34767309e22    // mass of moon, kg
#define DIST_EARTH_MOON    3.84400e8        // distance between earth and moon, meters
#define G                  6.67408e-11      // gravitational constant
#define EARTH_GRAVITY      9.81             // gravity of earth, m/s/s

#define METERS_TO_MILES(m) ((m) * 0.000621371)
#define TWO_PI             (2 * M_PI)
#define DEG_TO_RAD(d)      ((d) * (M_PI/180))
#define RAD_TO_DEG(r)      ((r) * (180/M_PI))

#define DELTA_H            1e-3


//
// typedefs
//

typedef struct {
    long double a;
    long double b;
    long double c;
} vector_t;

//
// variables
//

struct {
    long double x;
    long double y;
    long double z;
    long double mass;
    long double w;
} moon;

struct {
    long double x;
    long double y;
    long double z;
    long double mass;
    long double radius;
    long double w;
    long double g[100000];  // xxx rename g_surface
    long double h[100000];  // xxx rename r
    long double tpe;
} earth;

struct {
    long double x;
    long double y;
    long double z;
} xyz[100000];  // xxx malloc ?
int max_xyz;

//
// prototypes
//

static void init_xyz(void);
static long double potential_energy(long double m, long double g_surface, long double r);
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
    for (int i = 0; i < 100000; i++) {
        earth.h[i] = EARTH_RADIUS;
    }

    // Determine the position of the earth and moon relative to the center of mass.
    // Refer to the diagram above, the Earth and Moon have no velocity in the X direction, 
    // so the centrifugal force of the Moon must equal that of the Earth.
    //   moon.mass * moon.w^2 * moon.x = earth.mass * earth.w^2 * earth.x
    earth.x = -(DIST_EARTH_MOON / (1 + earth.mass / moon.mass));
    earth.y = 0;
    earth.z = 0;
    moon.x  = DIST_EARTH_MOON + earth.x;
    moon.y  = 0;
    moon.z  = 0;
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

    //xxx
    init_xyz();
    //max_xyz = 360; //xxx

    // xxx comment
    for (int i = 0; i < max_xyz; i++) {
        long double x, y, z, d;
        vector_t g, m, c, t;

        x = xyz[i].x;
        y = xyz[i].y;
        z = xyz[i].z;

        g.a = earth.x - x;
        g.b = earth.y - y;
        g.c = earth.z - z;
        set_vector_magnitude(&g, EARTH_GRAVITY);

        m.a = moon.x - x;
        m.b = moon.y - y;
        m.c = moon.z - z;
        d = magnitude(&m);
        set_vector_magnitude(&m, G * moon.mass / square(d));

#if 1
        c.a = -1;
        c.b = 0;
        c.c = 0;
        set_vector_magnitude(&c,  square(earth.w) * -earth.x);
#endif
#if 0
        c.a = x;
        c.b = y;
        c.c = z;
        d = magnitude(&c);
        set_vector_magnitude(&c, square(earth.w) * d);
#endif
#if 0
        memset(&c,0,sizeof(c));
#endif

        t.a = g.a + m.a + c.a;
        t.b = g.b + m.b + c.b;
        t.c = g.c + m.c + c.c;
        earth.g[i] = magnitude(&t);

        if (i == 0 || i == 90 || i == 180 || i == 270) {
            printf("%3d: earth gravity     = %0.9Lf %0.9Lf %0.9Lf\n", i, g.a, g.b, g.c);
            printf("     moon gravity      = %0.9Lf %0.9Lf %0.9Lf\n", m.a, m.b, m.c);
            printf("     centrifugal accel = %0.9Lf %0.9Lf %0.9Lf\n", c.a, c.b, c.c);
            printf("     total             = %0.9Lf %0.9Lf %0.9Lf\n", t.a, t.b, t.c);
            printf("     MAGNITUDE         = %0.9Lf\n", earth.g[i]);
            printf("\n");
        }
    }

    printf("Running ...\n");
    int loops = 0;
    long double earth_tpe_last = 0;
    long double m;
    int tpe_unchanged_count = 0;
    while (true) {
        int i = random() % max_xyz;
        int j = random() % max_xyz;
        long double delta_pe;
        if (i == j) continue;

#if 0
        // xxx use this one and clean up the mass calc
        delta_pe = 0;
        m = square((earth.h[i]+(DELTA_H/2)) / earth.radius) * DELTA_H,
        delta_pe += potential_energy(m,
                                     earth.g[i], 
                                     earth.h[i]+(DELTA_H/2));
        m = square((earth.h[j]-(DELTA_H/2)) / earth.radius) * DELTA_H,
        delta_pe -= potential_energy(m,
                                     earth.g[j], 
                                     earth.h[j]-(DELTA_H/2));
#else
        delta_pe = 0;
        delta_pe += potential_energy(square(earth.h[i]+(DELTA_H/2)),
                                     earth.g[i], 
                                     earth.h[i]+(DELTA_H/2));
        delta_pe -= potential_energy(square(earth.h[j]-(DELTA_H/2)),
                                     earth.g[j], 
                                     earth.h[j]-(DELTA_H/2));
#endif

        if (delta_pe < 0) {
            earth.h[i] += DELTA_H;
            earth.h[j] -= DELTA_H;
            earth.tpe += delta_pe;
        }

        loops++;
        if ((loops % 10000) == 0) {
            //printf("   loops = %9d   tpe = %0.10Le   tpe_unchanged_count = %d\n", loops, earth.tpe, tpe_unchanged_count);
            if (earth.tpe == earth_tpe_last) {
                tpe_unchanged_count++;
                // AAAA XXX
                if (tpe_unchanged_count == 10) break;
            } else {
                tpe_unchanged_count = 0;
            }
            earth_tpe_last = earth.tpe;
        }
        if (loops == 50000000) break;
    }
  
    // print results
    printf("Results ...\n");
    long double minh=1e99, maxh=-1e99;
    for (int i = 0; i < max_xyz; i++) {
        long double h = earth.h[i];
        if (i == 0 || i == 90 || i == 180 || i == 270) {
            printf("   deg = %3d   height = %0.10Lf  height-earth_radius = %+0.6Lf\n", 
                   i, h, h-earth.radius);
        }
        if (h < minh) minh = h;
        if (h > maxh) maxh = h;
    }
    printf("\n");
    printf("   min = %LF max = %LF  Range = %LF\n", minh, maxh, maxh-minh);

    // done
    return 0;
}

// -----------------  XXXXXXXXXXXXXXXX  -----------------------------

static void init_xyz(void)
{
//x = earth.x + earth.radius * cos(DEG_TO_RAD(i));
//y = earth.y + earth.radius * sin(DEG_TO_RAD(i));
// XXX deg vs rad
// xxx add earth.x etc   XXX

    double size, latitude, longitude, circ, x, y, z;
    int n, i;

    size = TWO_PI * EARTH_RADIUS / 360;
    printf("size = %.0f\n", size);

    for (latitude = 0; latitude < 90; latitude += 1) {
        circ = TWO_PI * EARTH_RADIUS * cos(DEG_TO_RAD(latitude));
        if (latitude < 90) {
            n = nearbyint(circ / size);
        } else {
            n = nearbyint((M_PI * square(EARTH_RADIUS * sin(DEG_TO_RAD(0.5)))) / square(size));
            printf("XXXXXXXXXXXXXXXXX OOPS\n");
        }
        printf("lat = %.0f  circ = %.0f  n = %d\n", latitude, circ, n);
        for (i = 0; i < n; i++) {
            longitude = i * (360./n);
            x = EARTH_RADIUS * cos(DEG_TO_RAD(latitude)) * cos(DEG_TO_RAD(longitude));
            y = EARTH_RADIUS * cos(DEG_TO_RAD(latitude)) * sin(DEG_TO_RAD(longitude));
            z = EARTH_RADIUS * sin(DEG_TO_RAD(latitude));

            if (latitude == 0) {
                printf("  longitude = %.3f  xyz = %10.0f %10.0f %10.0f\n", longitude, x, y, z);
            }

            // xyz =
            if (latitude == 0) {
                xyz[max_xyz].x = earth.x + x;
                xyz[max_xyz].y = earth.y + y;
                xyz[max_xyz].z = earth.z + z;
                max_xyz++;
            } else {
                xyz[max_xyz].x = earth.x + x;
                xyz[max_xyz].y = earth.y + y;
                xyz[max_xyz].z = earth.z + z;
                max_xyz++;

                xyz[max_xyz].x = earth.x + x;
                xyz[max_xyz].y = earth.y + y;
                xyz[max_xyz].z = earth.z + -z;
                max_xyz++;
            }
        }
    }
    printf("max_xyz = %d\n", max_xyz);
    printf("  %Le  %Le\n",
      max_xyz * square(size),
      4 * M_PI * square(EARTH_RADIUS));
    printf("%20f  %20f\n", M_PI, TWO_PI);
}

// -----------------  POTENTIAL ENERGY  -----------------------------

// xxx cleanup and save old comments
#if 0
The acceleration at the earth surface (due to earth gravity) is
               G * Me
   g_surface = ------
                Re^2

Imagine a tunnel through the center of the earth, the accel of gravity
at distance R from earth center is
       G * (Me * R^3 / Re^3)
   g = --------------------- = (G * Me / Re^3) * R
             R^2

       g_surface
   g = --------- * R
          Re

The amount of Energy to lift an object of mass m from the earth center to R is:
        R
   E = Integral  F * dR
        0

        R
   E = Integral  (m * g) * dR
        0

        R             g_surface
   E = Integral  (m * ---------- * R) * dR
        0                Re

                  g_surface
   E = 1/2 * m  * --------- * R^2
                     Re
#endif
static long double potential_energy(long double m, long double g_surface, long double r)
{
    // xxx change r to R  or change comments above to use 'r'
    // xxx use 0.5L ?
    return 0.5 * m * (g_surface / EARTH_RADIUS) * (r * r);
}

// -----------------  UTILS  ----------------------------------------

static long double square(long double x)
{
    return x * x;
}

static long double magnitude(vector_t *v)
{
    return sqrt( square(v->a) + square(v->b) + square(v->c) );
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
    v->c *= factor;

    return 0;
}
