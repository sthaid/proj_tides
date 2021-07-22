#include <stdio.h>
#include <math.h>

// Purpose: To calculate the location and height of the Earth's tidal bulges,
//          caused by the Moon.

// Assumptions:
// - the Earth to Moon distance is constant, the orbits are circular
// - calculations are done using 2 dimensions.
// - earth is circular.
// - 24 hour rotation of the Earth is not considered.
// - contributuion to the tidal bulges by the Sun is ignored.

// This program has 3 main parts:
//
// Part1 calculates the location of the Earth and Moon relative to the center of mass.
// Also the orbital angular velocity is determined
//
// Part2 program calculates the total acceleration at points around the
// Earth's circumference. The following accelerations are added to arrive at the total:
// - earth's gravity
// - moon's gravity
// - centrifugal, caused by the orbital rotation of the Earth around 
//   the center of mass of the Earth and Moon
//
// Part3 starts with a constant water level around the Earth, and iterates adjusting
// the levels to achieve a total minimum potential energy.

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

#define EARTH_MASS         5.972e24         // mass of earth, kg
#define EARTH_RADIUS       6.371e6          // meters
#define MOON_MASS          7.34767309e22    // mass of moon, kg
#define DIST_EARTH_MOON    3.84400e8        // distance between earth and moon, meters
#define G                  6.67408e-11      // gravitational constant
//#define EARTH_GRAVITY      9.81             // gravity of earth, m/s/s
#define EARTH_GRAVITY      10.0             // gravity of earth, m/s/s

#define METERS_TO_MILES(m) ((m) * 0.000621371)
#define TWO_PI             (2 * M_PI)

struct {
    double x;
    double y;
    double mass;
    double w;
} moon;

struct {
    double x;
    double y;
    double mass;
    double radius;
    double w;
} earth;

typedef struct {
    double a;
    double b;
} vector_t;

static inline double square(double x)
{
    return x * x;
}

double magnitude(vector_t *v)
{
    return sqrt( square(v->a) + square(v->b) );
}

int set_vector_magnitude(vector_t *v, double new_magnitude)
{
    double current_magnitude, factor;

    current_magnitude = magnitude(v);
    if (current_magnitude == 0) {
        return -1;  // xxx or abort
    }

    factor = new_magnitude / current_magnitude;

    v->a *= factor;
    v->b *= factor;

    return 0;
}


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
    printf("earth.x         = %8.0f miles\n", METERS_TO_MILES(earth.x));
    printf("moon.x          = %8.0f miles\n", METERS_TO_MILES(moon.x));

    // Determine the angular velocity of both the earth and the moon; they should be the same.
    // The centrifugal force must equal the gravitational force for a circular orbit
    //   gravity_force = moon.mass * moon.w^2 * moon.x
    double gravity_force = G * earth.mass * moon.mass / (DIST_EARTH_MOON * DIST_EARTH_MOON);
    moon.w  = sqrt(gravity_force / (moon.mass * fabs(moon.x)));
    earth.w = sqrt(gravity_force / (earth.mass * fabs(earth.x)));
    // xxx sanity check they are within .0000001
    printf("moon.w = %e  earth.w = %e\n", moon.w, earth.w);
    printf("orbital period = %f days\n", TWO_PI / moon.w / 86400);
    
    // xxx
    #define DEG_TO_RAD(d)  ((d) * (M_PI/180))
    #define RAD_TO_DEG(r)  ((r) * (180/M_PI))



    printf("\n");
    for (int deg = 0; deg < 360; deg++) {
        double angle, x, y;
        angle = DEG_TO_RAD(deg);
        x = earth.x + earth.radius * cos(angle);
        y = earth.y + earth.radius * sin(angle);

        //if (deg != 0 && deg != 90 && deg != 180 && deg != 270) continue;

        //printf("%.0f  x=%.0f y=%.0f\n", RAD_TO_DEG(angle), METERS_TO_MILES(x), METERS_TO_MILES(y));


        vector_t g, m, c, t;
        double d;

        g.a = earth.x - x;
        g.b = earth.y - y;
        set_vector_magnitude(&g, EARTH_GRAVITY);
        //printf("earth gravity = %0.9f %0.9f\n", g.a, g.b);

        m.a = moon.x - x;
        m.b = moon.y - y;
        d = magnitude(&m);
        set_vector_magnitude(&m, G * moon.mass / square(d));
        //printf("moon gravity = %0.9f %0.9f\n", m.a, m.b);

        c.a = x;
        c.b = y;
        d = magnitude(&c);
        //printf("d = %.0f\n", d);
        set_vector_magnitude(&c, square(earth.w) * d);
        //printf("centrigual accel= %0.9f %0.9f\n", c.a, c.b);

        t.a = g.a + m.a + c.a;
        t.b = g.b + m.b + c.b;
        //printf("total accel= %0.9f %0.9f\n", t.a, t.b);
        printf("%3d  MAGNITUDE = %0.9f\n", deg, magnitude(&t)-10);

        //printf("\n");



    }
        


#if 0


    Wearth = sqrt(F / (Rearth * earth.mass));
    Wmoon  = sqrt(F / (Rmoon  * moon.mass ));

    printf("%e %e\n", Wearth, Wmoon);
    double days = ((2*M_PI) / Wearth) / 86400.;
    printf("%f\n", days);

    double Eradius = 12.742e6 / 2;

    double W = Wearth;
    double a1 = W * W * (Eradius + Rearth);
    printf("a1 = %e\n", a1);

    double a2 = G * moon.mass / (D * D);
    printf("a2 = %e\n", a2);
#endif
    


    return 0;
}
