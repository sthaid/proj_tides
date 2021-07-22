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
#define MOON_MASS          7.34767309e22    // mass of moon, kg
#define DIST_EARTH_MOON    3.84400e8        // distance between earth and moon, meters
#define G                  6.67408e-11      // gravitational constant

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
    double w;
} earth;

int main(int argc, char **argv)
{
    // In the diagram above, the Earth and Moon have no velocity in the X direction, 
    // so the centrifugal force of the Moon must equal that of the Earth.
    //   moon.mass * moon.w^2 * moon.x = earth.mass * earth.w^2 * earth.x

    earth.mass = EARTH_MASS;
    moon.mass  = MOON_MASS;
    earth.x    = -(DIST_EARTH_MOON / (1 + earth.mass / moon.mass));
    moon.x     = DIST_EARTH_MOON + earth.x;
    printf("DIST_EARTH_MOON = %8.0f miles\n", METERS_TO_MILES(DIST_EARTH_MOON));
    printf("earth.x         = %8.0f miles\n", METERS_TO_MILES(earth.x));
    printf("moon.x          = %8.0f miles\n", METERS_TO_MILES(moon.x));

    // Determine the angular velocity of both the earth and the moon;
    // they should be the same.
    // The centrifugal force must equal the gravitational force for a circular orbit
    //   gravity_force = moon.mass * moon.w^2 * moon.x

    double gravity_force = G * earth.mass * moon.mass / (DIST_EARTH_MOON * DIST_EARTH_MOON);
    moon.w  = sqrt(gravity_force / (moon.mass * fabs(moon.x)));
    earth.w = sqrt(gravity_force / (earth.mass * fabs(earth.x)));
    // xxx sanity check they are within .0000001
    printf("moon.w = %e  earth.w = %e\n", moon.w, earth.w);
    printf("orbital period = %f days\n", TWO_PI / moon.w / 86400);
    


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
