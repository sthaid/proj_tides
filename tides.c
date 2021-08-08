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

// the earth.surface array is ordered with the 360 equator values first;
// to simulate the earth as a disk just the first 360 values of earth.surface are used
#define MAX_SURFACE (geometry == GEO_DISK ? 360 : earth.max_surface)

//
// typedefs
//

//
// variables
//

static double        theta;
static bool          moon_enabled;
static bool          sun_enabled;
static int           geometry;
static volatile bool tides_thread_init_complete;

//
// prototypes
//

static void init_earth_surface_xyz(void);
static void *tides_thread(void *cx);
static void set_earth_moon_position_and_surface_values(void);

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

    // determine the following for both the orbit of the earth/moon and the
    // orbit of the earth/sun:
    // - the distance to the barycenter (orbit_radius)
    // - the angular velocity (orbit_w)
    // this calculation assumes circular orbits

    earth.em_orbit_radius = (DIST_EARTH_MOON / (1 + EARTH_MASS / MOON_MASS));
    moon.em_orbit_radius  = DIST_EARTH_MOON - earth.em_orbit_radius;
    gravity_force = G * EARTH_MASS * MOON_MASS / (DIST_EARTH_MOON * DIST_EARTH_MOON);
    moon.em_orbit_w  = sqrt(gravity_force / (MOON_MASS * moon.em_orbit_radius));
    earth.em_orbit_w = sqrt(gravity_force / (EARTH_MASS * earth.em_orbit_radius));
    assert(fabs(moon.em_orbit_w / earth.em_orbit_w - 1) < 1e-10);

    earth.es_orbit_radius = (DIST_EARTH_SUN / (1 + EARTH_MASS / SUN_MASS));
    sun.es_orbit_radius = DIST_EARTH_SUN - earth.es_orbit_radius;
    gravity_force = G * EARTH_MASS * SUN_MASS / (DIST_EARTH_SUN * DIST_EARTH_SUN);
    earth.es_orbit_w = sqrt(gravity_force / (EARTH_MASS * earth.es_orbit_radius));
    sun.es_orbit_w = sqrt(gravity_force / (SUN_MASS * sun.es_orbit_radius));
    assert(fabs(sun.es_orbit_w / earth.es_orbit_w - 1) < 1e-10);

    printf("EARTH_RADIUS          = %0.6e m  %8.0f miles\n", EARTH_RADIUS, METERS_TO_MILES(EARTH_RADIUS));
    printf("\n");

    printf("DIST_EARTH_MOON       = %0.6e m  %8.0f miles\n", DIST_EARTH_MOON, METERS_TO_MILES(DIST_EARTH_MOON));
    printf("earth.em_orbit_radius = %0.6e m  %8.0f miles\n", earth.em_orbit_radius, METERS_TO_MILES(earth.em_orbit_radius));
    printf("moon.em_orbit_radius  = %0.6e m  %8.0f miles\n", moon.em_orbit_radius, METERS_TO_MILES(moon.em_orbit_radius));
    printf("orbital period        = %f days\n", TWO_PI / moon.em_orbit_w / 86400);
    printf("\n");

    printf("DIST_EARTH_SUN        = %0.6e m  %8.6f million miles\n", DIST_EARTH_SUN, METERS_TO_MILES(DIST_EARTH_SUN)/1000000);
    printf("earth.es_orbit_radius = %0.6e m  %8.6f million miles\n", earth.es_orbit_radius, METERS_TO_MILES(earth.es_orbit_radius)/1000000);
    printf("sun.es_orbit_radius   = %0.6e m  %8.6f million miles\n", sun.es_orbit_radius, METERS_TO_MILES(sun.es_orbit_radius)/1000000);
    printf("orbital period        = %f days\n", TWO_PI / earth.es_orbit_w / 86400);
    printf("\n");

    // initialize array of approximately evenly distributed points on the earth surface;
    // these points are each in a square area that is 60 nautical miles per side;
    // along the equator there are 360 of these square areas, at other latitudes there are fewer
    init_earth_surface_xyz();

    // create runtime thread, and
    // wait for thread to have completed its initialization
    pthread_create(&tid, NULL, tides_thread, NULL);
    while (!tides_thread_init_complete) {
        usleep(10000);
    }
}

// This routine initializes the earth.surface[].x,y,z values.
// The 360 positions around the equator (latitude=0, z=0) are added first, 
//  followed by values for latitude=1 and -1, 2 and -2, etc.
// Each location in the earth.surface array is an, approximately, 
//  60 NM x 60 NM square.
// The sum of the areas of these squares divided by the area of the 
//  surface of the earth equals 1.00012 (a close approximation).
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
                // add values for this location both at 'latitude' and '-latitude'
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

// xxx recheck , comments, and descritpiton in notes.txt
static void *tides_thread(void *cx)
{
    int loops = 0;

    // initailze
    for (int i = 0; i < MAX_EARTH_SURFACE; i++) {
        earth.surface[i].r = EARTH_RADIUS;
    }
    set_earth_moon_position_and_surface_values();
    __sync_synchronize();
    tides_thread_init_complete = true;

    // loop forever
    while (true) {
        int i = random() % MAX_SURFACE;
        int j = random() % MAX_SURFACE;
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

        // increment loops
        loops++;
        //if ((loops % 1000000) == 0) printf("   loops = %9d\n", loops);

        // every 10000 loops:
        // - when mtion is enabled, increment theat by 0.1 degrees every 10 ms
        // - if ctrls have changed then call set_earth_moon_position
        if ((loops % 10000) == 0) {
            if (ctrls.motion) {
                uint64_t now_us = microsec_timer();
                static uint64_t last_us;
                if (now_us > last_us + 10000) {
                    ctrls.theta += .1;
                    last_us = now_us;
                }
            }

            if (ctrls.theta != theta ||
                ctrls.moon_enabled != moon_enabled ||
                ctrls.sun_enabled != sun_enabled ||
                ctrls.geometry != geometry)
            {
                set_earth_moon_position_and_surface_values();
            }
        }
    }

    return NULL;
}

// This routine re-calculaes the position of the earth and moon,
// and the earth surface 'g' and 'v' values.
// Where: 
// - g is the magnitude of the sum of the accelerations at a location on earth surface, 
//   including: earth, moon, sun gravity, and centrifigul accels due to 
//   the motion of the earth about the earth-moon and earth-sun barycenters.
// - v is the vector sum of the above, without including earth gravity. These are the
//   vectors that are displayed. Earth gravity is not included in this vector, because
//   the earth gravity is so much larger than the other accels, that all you would see 
//   on the display is the earth gravity vector.

static void set_earth_moon_position_and_surface_values(void)
{
    //uint64_t start = microsec_timer();

    // update the global variables that are the current state of the simulation
    theta        = ctrls.theta;
    moon_enabled = ctrls.moon_enabled;
    sun_enabled  = ctrls.sun_enabled;
    geometry     = ctrls.geometry;

    // set earth, moon and sun's position, 
    // note that :
    // - the origin is the earth-moon barycenter
    // - the locations are in the x,y plane (z is always 0)
    earth.x = -cos(DEG_TO_RAD(theta)) * earth.em_orbit_radius;
    earth.y = -sin(DEG_TO_RAD(theta)) * earth.em_orbit_radius;
    earth.z = 0;
    moon.x  = cos(DEG_TO_RAD(theta)) * moon.em_orbit_radius;
    moon.y  = sin(DEG_TO_RAD(theta)) * moon.em_orbit_radius;
    moon.z  = 0;
    sun.x = 0;
    sun.y = earth.y - sqrt(square(DIST_EARTH_SUN) - square(earth.x));
    sun.z = 0;

    // loop over the surface being simulated, where MAX_SURFACE will be either
    // - 360 when the earth is being simulated as a disk, or
    // - a large number, which is total the number of 60NM x 60NM squares on the 
    //   spherical earth.
    for (int i = 0; i < MAX_SURFACE; i++) {
        vector_t g, m, c, s, cs, t;
        double d;

        // compute the following acceleration vectors, at location earth.surface[i]:
        // - g:  earth gravity
        // - m:  moon gravity
        // - c:  centrifigual accel due to earth circular motion around earth-moon barycenter
        // - s:  sun gravity
        // - cs: centrifigual accel due to earth circular motion around earth-sun barycenter
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

        // add the 5 vectors that were computed above, and store result in 't'; because
        // m,c,s, and cs are tiny the direction and magnitide of 't' are very close to the
        // earth gravity vector 'g'
        vector_init(&t, 0, 0, 0);
        vector_add(&t, &g);
        vector_add(&t, &m);
        vector_add(&t, &c);
        vector_add(&t, &s);
        vector_add(&t, &cs);

        // set earth.surface[i[.g to the magnitude of vector 't';
        // it may be confusing that the vector 'g' and the earth.surface[i].g are 
        //  not the same thing, the vector g is the earth's gravity at this location, and 
        //  the searth_surface[i].g is the magnitude of the sum of all accelerations 
        //  at this location
        earth.surface[i].g = vector_get_magnitude(&t);

        // set earth.surface[i[.v to the sum of all accels, along the equator, 
        //  excluding the accel of earth gravity;
        // this 'v' is used for display purpose only
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

#if 1
        // debug prints
        if (i == theta+0 || i == theta+90 || i == theta+180 || i == theta+270) {
            if (i == 0) printf("accel vectors:\n");
            printf("%3d: earth gravity          = %+0.9f %+0.9f %+0.9f   magnitude = %0.9f\n", 
                   i, g.a, g.b, g.c, vector_get_magnitude(&g));
            printf("     moon gravity           = %+0.9f %+0.9f %+0.9f   magnitude = %0.9f\n", 
                   m.a, m.b, m.c, vector_get_magnitude(&m));
            printf("     moon centrifugal accel = %+0.9f %+0.9f %+0.9f   magnitude = %0.9f\n", 
                   c.a, c.b, c.c, vector_get_magnitude(&c));
            printf("     sun gravity            = %+0.9f %+0.9f %+0.9f   magnitude = %0.9f\n", 
                   s.a, s.b, s.c, vector_get_magnitude(&s));
            printf("     sun centrifugal accel  = %+0.9f %+0.9f %+0.9f   magnitude = %0.9f\n", 
                   cs.a, cs.b, cs.c, vector_get_magnitude(&cs));
            printf("     total                  = %+0.9f %+0.9f %+0.9f   magnitude = %0.9f\n",
                   t.a, t.b, t.c, vector_get_magnitude(&t));
            printf("\n");
        }
#endif
    }

    //printf("set_earth_moon_position_and_surface_values DURATION = %ld us\n", (microsec_timer()-start));
}

// -----------------  xxxxxxxxxxx---------------------

// this routine returns the min and max earth surface radius along the equator
void tides_get_min_max(double *min, double *max, int *min_idx, int *max_idx)
{
    *min =  1e99;
    *max = -1e99;

    for (int i = 0; i < 360; i++) {
        double r = earth.surface[i].r;
        if (r < *min) {
            *min = r;
            if (min_idx) *min_idx = i;
        }
        if (r > *max) {
            *max = r;
            if (max_idx) *max_idx = i;
        }
    }
}

// -----------------  UNIT TEST  ---------------------

#ifdef UNIT_TEST
int main(int argc, char **argv)
{
    double min, max;
    int min_idx, max_idx;

    ctrls.moon_enabled = true;
    tides_init();

    while (true) {
        sleep(1);

        printf("Results ...\n");
        tides_get_min_max(&min, &max, &min_idx, &max_idx);
        printf("   min = %f (idx=%d)   max = %f (idx=%d)   range = %f m  %f ft\n", 
               min, min_idx, max, max_idx, max-min, METERS_TO_FEET(max-min));
        printf("   %3d = %+0.3f\n",   0, METERS_TO_FEET(earth.surface[0].r-EARTH_RADIUS));
        printf("   %3d = %+0.3f\n",  90, METERS_TO_FEET(earth.surface[90].r-EARTH_RADIUS));
        printf("   %3d = %+0.3f\n", 180, METERS_TO_FEET(earth.surface[180].r-EARTH_RADIUS));
        printf("   %3d = %+0.3f\n", 270, METERS_TO_FEET(earth.surface[270].r-EARTH_RADIUS));
        printf("\n");
    }

    return 0;
}
#endif
