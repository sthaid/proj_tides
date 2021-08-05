// xxx
// - ctrls
//   - motion
//   - sun
//   - vectors
//   - maybe disable centrifigul force
// - fill in the blue areafill in the blue areafill in the blue area
// - in right side, print out the min and max tides

// ---------------

// xxx lower prio
// - more comments needed
// - verify the min and max are at latitude 0

// xxx maybe ...
// - option to add the sun ?

// xxx misc,  
// - backup proj and zz dirs

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <assert.h>

#include <vectors.h>
#include <tides.h>

//
// defines
//

//
// typedefs
//

//
// variables
//

static bool moon_enabled;
static bool sun_enabled;

//
// prototypes
//

static void init_earth_surface_xyz(void);
static void *tides_thread(void *cx);
static void set_earth_and_moon_position(double theta);

//
// inline functions
//

static inline double square(double x) { return x * x; }

static inline uint64_t microsec_timer(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC,&ts);
    return  ((uint64_t)ts.tv_sec * 1000000) + ((uint64_t)ts.tv_nsec / 1000);
}

// -----------------  INIT  --------------------------

void tides_init(void)
{
    pthread_t tid;
    double gravity_force;

    // xxx move or delete
    // xxx wait for thread to start
    for (int i = 0; i < MAX_EARTH_SURFACE; i++) {
        earth.surface[i].r = EARTH_RADIUS;
    }

    // determine the orbital radius of the earth and moon, relative to barycenter;
    // the centrifugal force of the Moon must equal that of the Earth.
    //   MOON_MASS * moon.em_orbit_w^2 * moon.x = EARTH_MASS * earth.em_orbit_w^2 * earth.x
    earth.em_orbit_radius = (DIST_EARTH_MOON / (1 + EARTH_MASS / MOON_MASS));
    moon.em_orbit_radius  = DIST_EARTH_MOON - earth.em_orbit_radius;

    // determine the angular velocity of both the earth and the moon; they should be the same;
    // the centrifugal force must equal the gravitational force for a circular orbit
    //   gravity_force = MOON_MASS * moon.em_orbit_w^2 * moon.x  xxx check all comments
    gravity_force = G * EARTH_MASS * MOON_MASS / (DIST_EARTH_MOON * DIST_EARTH_MOON);
    moon.em_orbit_w  = sqrt(gravity_force / (MOON_MASS * moon.em_orbit_radius));
    earth.em_orbit_w = sqrt(gravity_force / (EARTH_MASS * earth.em_orbit_radius));
    assert(fabs(moon.em_orbit_w / earth.em_orbit_w - 1) < 1e-10);



    // xxx check
    earth.es_orbit_radius = (DIST_EARTH_SUN / (1 + EARTH_MASS / SUN_MASS));
    gravity_force = G * EARTH_MASS * SUN_MASS / (DIST_EARTH_SUN * DIST_EARTH_SUN);
    earth.es_orbit_w = sqrt(gravity_force / (EARTH_MASS * earth.es_orbit_radius));
    printf("DIST_EARTH_SUN        = %0.6e m  %8.6f million miles\n", DIST_EARTH_SUN, METERS_TO_MILES(DIST_EARTH_SUN)/1000000);
    printf("earth.es_orbit_radius = %0.6e m  %8.6f million miles\n", earth.es_orbit_radius, METERS_TO_MILES(earth.es_orbit_radius)/1000000);
    printf("orbital period        = %f days\n", TWO_PI / earth.es_orbit_w / 86400);


    // print values
    printf("DIST_EARTH_MOON       = %0.6e m  %8.0f miles\n", DIST_EARTH_MOON, METERS_TO_MILES(DIST_EARTH_MOON));
    printf("EARTH_RADIUS          = %0.6e m  %8.0f miles\n", EARTH_RADIUS, METERS_TO_MILES(EARTH_RADIUS));
    printf("earth.em_orbit_radius = %0.6e m  %8.0f miles\n", earth.em_orbit_radius, METERS_TO_MILES(earth.em_orbit_radius));
    printf("moon.em_orbit_radius  = %0.6e m  %8.0f miles\n", moon.em_orbit_radius, METERS_TO_MILES(moon.em_orbit_radius));
    printf("earth.em_orbit_w      = %0.6e\n", earth.em_orbit_w);
    printf("moon.em_orbit_w       = %0.6e\n", moon.em_orbit_w);
    printf("orbital period        = %f days\n", TWO_PI / moon.em_orbit_w / 86400);
    printf("\n");

    // initialize array of approximately evenly distributed points on the earth surface;
    // these points are each in a square area that is 60 nautical miles per side
    init_earth_surface_xyz();

    // create runtime thread
    pthread_create(&tid, NULL, tides_thread, NULL);

    // xxx wait for thread to start before returning
}

static void init_earth_surface_xyz(void)
{
    #define ADD(_x,_y,_z) \
        do { \
            int idx = earth.max_surface; \
            assert(idx < MAX_EARTH_SURFACE); \
            earth.surface[idx].x = (_x); \
            earth.surface[idx].y = (_y); \
            earth.surface[idx].z = (_z); \
            earth.max_surface++; \
        } while (0)

    double size, latitude, longitude, circ, x, y, z;
    int n, i;

    size = TWO_PI * EARTH_RADIUS / 360;

    for (latitude = 0; latitude <= 90; latitude += 1) {
        circ = TWO_PI * EARTH_RADIUS * cos(DEG_TO_RAD(latitude));
        if (latitude < 90) {
            n = nearbyint(circ / size);
        } else {
            n = nearbyint((M_PI * square(EARTH_RADIUS * sin(DEG_TO_RAD(0.5)))) / square(size));
        }
        //printf("  latitude = %.0f  circ = %.0f  n = %d\n", latitude, circ, n);

        for (i = 0; i < n; i++) {
            longitude = i * (360./n);
            x = EARTH_RADIUS * cos(DEG_TO_RAD(latitude)) * cos(DEG_TO_RAD(longitude));
            y = EARTH_RADIUS * cos(DEG_TO_RAD(latitude)) * sin(DEG_TO_RAD(longitude));
            z = EARTH_RADIUS * sin(DEG_TO_RAD(latitude));

            //if (latitude == 0) {
            //    printf("    longitude = %8.3f  xyz = %10.0f %10.0f %10.0f\n", longitude, x, y, z);
            //}

            if (latitude == 0) {
                ADD(x,y,z);
            } else {
                ADD(x,y,z);
                ADD(x,y,-z);
            }
        }
    }

    //printf("init_earth_surface:\n");
    //printf("  size = %.0f m   %0.6f NM\n", size, METERS_TO_NAUTICAL_MILES(size));
    //printf("  max_surface = %d   area = %e  %e\n",
    //       earth.max_surface, earth.max_surface * square(size), 4 * M_PI * square(EARTH_RADIUS));
    //printf("\n");
}

// -----------------  RUNTIME  -----------------------

// xxx recheck 
static void *tides_thread(void *cx)
{
    int loops = 0;
    double theta;
    uint64_t start_us;

    // set the initial earth and moon positions, at theta=0
    theta = 0;
    for (int i = 0; i < MAX_EARTH_SURFACE; i++) {  // xxx rese tin 2places
        earth.surface[i].r = EARTH_RADIUS;
    }
    set_earth_and_moon_position(theta);

    while (true) {
        // xxx
        start_us = microsec_timer();
        while (true) {
            int i = random() % earth.max_surface;
            int j = random() % earth.max_surface;
            if (i == j) continue;

            double delta_pe = 0;
            double m, g_surface, r;

            double dr_i = (.001 * EARTH_RADIUS * EARTH_RADIUS) / square(earth.surface[i].r);
            double dr_j = (.001 * EARTH_RADIUS * EARTH_RADIUS) / square(earth.surface[j].r);

            m         = 1;
            g_surface = earth.surface[i].g;
            r         = earth.surface[i].r + dr_i/2;
            delta_pe += m * g_surface * square(r);

            m         = 1;
            g_surface = earth.surface[j].g;
            r         = earth.surface[j].r - dr_j/2;
            delta_pe -= m * g_surface * square(r);

            if (delta_pe < 0) {
                earth.surface[i].r += dr_i;
                earth.surface[j].r -= dr_j;
            }

            loops++;

            // every 20000 loops (xxx how long) check:
            // - is it time to advance theta
            // - have ctrls changed, requiring a call to set_earth_and_moon_position
            if ((loops % 20000) == 0) {
                bool call = false;
                if (ctrls.moon_motion_enabled && microsec_timer() > start_us + 10000) {
                    theta += .1;
                    call = true;
                }
                if (ctrls.reset_request) {  // xxx name
                    ctrls.reset_request = false;
                    theta = 0;
                    call = true;
                }
                if (ctrls.moon_enabled != moon_enabled) {
                    call = true;
                }
                if (ctrls.sun_enabled != sun_enabled) {
                    call = true;
                }
                if (call) {
                    set_earth_and_moon_position(theta);
                }
            }

            //if ((loops % 1000000) == 0) {
            //    printf("   loops = %9d\n", loops);
            //}
        }
    }

    return NULL;
}

// xxx better name
static void set_earth_and_moon_position(double theta)
{
    //uint64_t start = microsec_timer();

    moon_enabled = ctrls.moon_enabled;
    sun_enabled = ctrls.sun_enabled;

    // set earth and moon position based on input arg 'theta', 
    // which is the orbit position in degrees
    earth.x = -cos(DEG_TO_RAD(theta)) * earth.em_orbit_radius;
    earth.y = -sin(DEG_TO_RAD(theta)) * earth.em_orbit_radius;
    earth.z = 0;
    moon.x  = cos(DEG_TO_RAD(theta)) * moon.em_orbit_radius;
    moon.y  = sin(DEG_TO_RAD(theta)) * moon.em_orbit_radius;
    moon.z  = 0;

    sun.x = 0;
    sun.y = earth.y - sqrt(square(DIST_EARTH_SUN) - square(earth.x));
    sun.z = 0;
//xxx print

    // compute accel vectors for all earth surface locations just computed
    // - g = earth gravity
    // - m = moon gravity
    // - c = centrifugal accel
    // - t = g + m + c
    for (int i = 0; i < earth.max_surface; i++) {
        vector_t g, m, c, t;
        vector_t s, cs;
        double d;

        vector_init(&g,
                    earth.x - (earth.surface[i].x + earth.x),
                    earth.y - (earth.surface[i].y + earth.y),
                    earth.z - (earth.surface[i].z + earth.z));
        vector_set_magnitude(&g, EARTH_GRAVITY);

        if (moon_enabled) {
            vector_init(&m,
                        moon.x - (earth.surface[i].x + earth.x),
                        moon.y - (earth.surface[i].y + earth.y),
                        moon.z - (earth.surface[i].z + earth.z));
            d = vector_get_magnitude(&m);
            vector_set_magnitude(&m, G * MOON_MASS / square(d));

            vector_init(&c, earth.x, earth.y, earth.z);
            d = vector_get_magnitude(&c);
            vector_set_magnitude(&c,  square(earth.em_orbit_w) * d);
        } else {
            vector_init(&m, 0, 0, 0);
            vector_init(&c, 0, 0, 0);
        }

        if (sun_enabled) {
            vector_init(&s,
                        sun.x - (earth.surface[i].x + earth.x),
                        sun.y - (earth.surface[i].y + earth.y),
                        sun.z - (earth.surface[i].z + earth.z));
            d = vector_get_magnitude(&s);
            vector_set_magnitude(&s, G * SUN_MASS / square(d));

            vector_init(&cs, earth.x-sun.x, earth.y-sun.y, earth.z-sun.z);
            d = vector_get_magnitude(&cs);
            vector_set_magnitude(&cs,  square(earth.es_orbit_w) * d);
        } else {
            vector_init(&s, 0, 0, 0);
            vector_init(&cs, 0, 0, 0);
        }

        vector_init(&t, 0, 0, 0);
        vector_add(&t, &g);
        vector_add(&t, &m);
        vector_add(&t, &c);
        vector_add(&t, &s);
        vector_add(&t, &cs);
        earth.surface[i].g = vector_get_magnitude(&t);

        if (i < 360) {
            vector_t *v = &earth.surface[i].v;
            vector_init(v, 0, 0, 0);
            vector_add(v, &m);
            vector_add(v, &c);
            vector_add(v, &s);
            vector_add(v, &cs);

            //if (i == 0 || i == 90 || i == 180 || i== 270) {
                //printf("vect[%d] = %e %e %e\n", i, v->a, v->b, v->c);
            //}
        }            

#if 1  // xxx check
        if (i == theta+0 || i == theta+90 || i == theta+180 || i == theta+270) {
            if (i == 0) printf("accel vectors:\n");
            printf("%3d: earth gravity     = %+0.9f %+0.9f %+0.9f\n", i, g.a, g.b, g.c);
            printf("     moon gravity      = %+0.9f %+0.9f %+0.9f\n", m.a, m.b, m.c);
            printf("     centrifugal accel = %+0.9f %+0.9f %+0.9f\n", c.a, c.b, c.c);
            printf("     total             = %+0.9f %+0.9f %+0.9f\n", t.a, t.b, t.c);
            printf("     moon+centrifugal  = %e %e %e\n", earth.surface[i].v.a, earth.surface[i].v.b, earth.surface[i].v.c);
            printf("     MAGNITUDE         = %0.9f\n", earth.surface[i].g);
            if (i == 270) printf("\n");
        }
#endif
    }

    //printf("set_earth_and_moon_position DURATION = %ld us\n", (microsec_timer()-start));
}

// -----------------  UNIT TEST  ---------------------

#ifdef UNIT_TEST
int main(int argc, char **argv)
{
    //ctrls.moon_enabled = true;
    ctrls.sun_enabled = true;
    ctrls.moon_enabled = true;

    tides_init();

    while (true) {
        sleep(1);

        printf("Results ...\n");

        double min_r=1e99, max_r=-1e99;
        int min_r_idx=-1, max_r_idx=-1;

        //for (int i = 0; i < earth.max_surface; i++) 
        for (int i = 0; i < 360; i++) 
        {
            double r = earth.surface[i].r;

            if (i == 0 || i == 90 || i == 180 || i == 270) {
                printf("   deg = %3d   height = %0.10f  height-earth_radius = %+0.6f\n", 
                       i, r, r-EARTH_RADIUS);
            }

            if (r < min_r) {
                min_r = r;
                min_r_idx = i;
            }
            if (r > max_r) {
                max_r = r;
                max_r_idx = i;
            }
        }
        printf("   min = %f (idx=%d)   max = %f (idx=%d)   range = %f m  %f ft\n", 
               min_r, min_r_idx, max_r, max_r_idx, max_r-min_r,
               METERS_TO_FEET(max_r-min_r));
        printf("\n");
    }

    return 0;
}
#endif
