// XXX

// - make an equation for the total pe of a colum from ctr of earth
//   - verify this equation is correct

// - doubles
//   - cast floating point constants to long double
//   - use define for long double

// - performance
//   - how long to converge if random i,j used instead of full loops

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
#include <string.h>
#include <math.h>

// Purpose: To calculate the location and height of the Earth's tidal bulges,
//          caused by the Moon.

// Coordinates Diagram
// - The X,Y origin is at the center of mass. 
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
//#define EARTH_GRAVITY      9.81             // gravity of earth, m/s/s
#define EARTH_GRAVITY      10.0             // gravity of earth, m/s/s  xxx

#define METERS_TO_MILES(m) ((m) * 0.000621371)
#define TWO_PI             (2 * M_PI)
#define DEG_TO_RAD(d)      ((d) * (M_PI/180))
#define RAD_TO_DEG(r)      ((r) * (180/M_PI))

#define DELTA_H            1e-3

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
    long double tpe_start;
} earth;

typedef struct {
    long double a;
    long double b;
} vector_t;

//
// prototypes
//

static long double potential_energy(long double g_surface, long double R);
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

#if 1
        c.a = -1;
        c.b = 0;
        set_vector_magnitude(&c,  square(earth.w) * -earth.x);
#endif
#if 0
        c.a = x;
        c.b = y;
        d = magnitude(&c);
        set_vector_magnitude(&c, square(earth.w) * d);
#endif
#if 0
        memset(&c,0,sizeof(c));
#endif

        t.a = g.a + m.a + c.a;
        t.b = g.b + m.b + c.b;
        earth.g[deg] = magnitude(&t);

        if (deg == 0 || deg == 90 || deg == 180 || deg == 270) {
            printf("%3d: earth gravity     = %0.9Lf %0.9Lf\n", deg, g.a, g.b);
            printf("     moon gravity      = %0.9Lf %0.9Lf\n", m.a, m.b);
            printf("     centrifugal accel = %0.9Lf %0.9Lf\n", c.a, c.b);
            printf("     total             = %0.9Lf %0.9Lf\n", t.a, t.b);
            printf("     MAGNITUDE         = %0.9Lf\n", earth.g[deg]);
            printf("\n");
        }
    }

    // xxx comment
    for (int i = 0; i < 360; i++) {
        earth.h[i] = earth.radius;
        //earth.tpe += (earth.g[i] * pow(earth.h[i], 5));   // xxx function
        earth.tpe += potential_energy(earth.g[i], earth.h[i]);
    }
    earth.tpe_start = earth.tpe;

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
                //dpe = (earth.g[i] * pow(earth.h[i]+DELTA_H/2, 5)) -
                //      (earth.g[j] * pow(earth.h[j]-DELTA_H/2, 5));
                dpe = potential_energy(earth.g[i], earth.h[i] + DELTA_H/2) -
                      potential_energy(earth.g[j], earth.h[j] - DELTA_H/2);
                //printf("i,j  %d %d   dpe = %0.20Lf\n", i, j, dpe);  //xxx
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
                earth.tpe += best_dpe;
                num_exchanges++;
            }
        }
        if ((++loops % 10) == 0) {
            printf("   loops = %3d  numex = %3d  tpe = %0.10Le  tpe-tpe_start = %0.10Le\n", 
               loops, num_exchanges, earth.tpe,
               earth.tpe - earth.tpe_start);
        }
        if (num_exchanges <= 2) {
            break;
        }
    }
    printf("   final loops = %d\n", loops);
    printf("\n");
  
    // print results
    printf("Results ...\n");
    long double minh=1e99, maxh=-1e99;
    for (int i = 0; i < 360; i++) {
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

// -----------------  POTENTIAL ENERGY  -----------------------------

static long double potential_energy(long double g_surface, long double R)
{
#if 0
XXX recheck this

Define ... XXX
  G
  Me
  Re
  g
  etc

The acceleration at the earth surface (due to earth gravity) is
               G * Me
   g_surface = ------
                Re^2

Imagine a tunnel through the center of the earth, the accel of gravity
at distance R from earth center is
       G * (Me * R^3 / Re^3)
   g = --------------------- = (G * Me / Re^3) * R
             R^2

   K = (G * Me / Re^3)

   g = K * R

The amount of Energy to lift an object of mass m from the earth center to R is:
        R
   E = Integral  F * dR
        0

        R
   E = Integral  (m * g) * dR
        0

        R
   E = Integral  (m * K * R) * dR
        0

   E = 1/2 * K * m * R^2

The total potentail energy of a cone of water where the area of the cone at
the earth surface is A, and the density of water is p, and the height of the
cone is R:

             R^2
   m = (A * ----- * dR) * p
            Re^2

        Re
   E = Integral  1/2 * K * m * R^2
        0

        Re                       R^2
   E = Integral  1/2 * K * (A * ----- * dR * p) * R^2
        0                       Re^2

   K1 = 1/2 * K * A / Re^2 * p

        Re
   E = Integral  K1 * R^4 * dR
        0 

   E = 1/5 * K1 * R^5


   K1 = 1/2 * (G * Me / Re^3) * A / Re^2 * p

              G * Me
   K1 = 1/2 * ------ * A * p
               Re^5

                G * Me          
   E = ( 1/10 * ------ * A * p ) * R^5
                 Re^5           

  The units of the above equation for E are kg * m^2 / s^2.
  This is the correct unit for Energy.

Next, work g_surface into the above equation for total potential energy.
        g_surface * Re^2
   Me = ----------------
             G

                g_surface       
   E = ( 1/10 * --------- * A * p ) * R^5
                 Re^3           

         A * p
   E = --------- * g_surface * R^5
       10 * Re^3

#endif
    //XXX include the constant in the result too?
    //return g_surface * powl(R, 5);
    return g_surface * pow(R, 5);  // xxx this is much faster and with the same result
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
