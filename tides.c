// xxx more comments needed

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <assert.h>

#include <tides.h>

//
// defines
//

#if 0
#define EARTH_MASS         5.972e24         // mass of earth, kg
#define EARTH_RADIUS       6.371e6          // meters
#define MOON_MASS          7.34767309e22    // mass of moon, kg
#define DIST_EARTH_MOON    3.84400e8        // distance between earth and moon, meters
#define G                  6.67408e-11      // gravitational constant
#define EARTH_GRAVITY      9.81             // gravity of earth, m/s/s
#endif


#if 0
#define METERS_TO_MILES(m)          ((m) * 0.000621371)
#define METERS_TO_NAUTICAL_MILES(m) ((m) * 0.000539957)
#define TWO_PI                      (2 * M_PI)
#define DEG_TO_RAD(d)               ((d) * (M_PI/180))
#define RAD_TO_DEG(r)               ((r) * (180/M_PI))
#endif

//
// typedefs
//

typedef struct {
    double a;
    double b;
    double c;
} vector_t;

//
// variables
//

//
// prototypes
//

void tides_init(void);
static void init_earth_surface_xyz(void);
static void *tides_thread(void *cx);
static double square(double x);
static double magnitude(vector_t *v);
static int set_vector_magnitude(vector_t *v, double new_magnitude);
static uint64_t microsec_timer(void);

// -----------------  INIT  --------------------------

void tides_init(void)
{
    pthread_t tid;

    // Init constant fields of the earth and moon structs.
    earth.mass   = EARTH_MASS;
    moon.mass    = MOON_MASS;
    earth.radius = EARTH_RADIUS;
    for (int i = 0; i < 100000; i++) {
        earth.surface[i].r = EARTH_RADIUS;
    }

    //xxx(0);

    // create runtime thread
    pthread_create(&tid, NULL, tides_thread, NULL);
}

static void xxx(double theta)
{
    uint64_t start = microsec_timer();
    printf("XXXXXXXXXXX  theta   %f\n", theta);

    // Determine the position of the earth and moon relative to the center of mass.
    // Refer to the diagram above, the Earth and Moon have no velocity in the X direction, 
    // so the centrifugal force of the Moon must equal that of the Earth.
    //   moon.mass * moon.w^2 * moon.x = earth.mass * earth.w^2 * earth.x

    double xxx_earth = -(DIST_EARTH_MOON / (1 + earth.mass / moon.mass));
    double xxx_moon = DIST_EARTH_MOON + xxx_earth;

    earth.x = cos(DEG_TO_RAD(theta)) * xxx_earth;
    earth.y = sin(DEG_TO_RAD(theta)) * xxx_earth;
    earth.z = 0;
    moon.x  = cos(DEG_TO_RAD(theta)) * xxx_moon;
    moon.y  = sin(DEG_TO_RAD(theta)) * xxx_moon;
    moon.z  = 0;
    printf("DIST_EARTH_MOON = %0.6e m  %8.0f miles\n", DIST_EARTH_MOON, METERS_TO_MILES(DIST_EARTH_MOON));
    printf("EARTH_RADIUS    = %0.6e m  %8.0f miles\n", EARTH_RADIUS, METERS_TO_MILES(EARTH_RADIUS));
// xxx vvv
    printf("earth.x         = %0.6e m  %8.0f miles\n", earth.x, METERS_TO_MILES(earth.x));
    printf("moon.x          = %0.6e m  %8.0f miles\n", moon.x, METERS_TO_MILES(moon.x));

    // Determine the angular velocity of both the earth and the moon; they should be the same.
    // The centrifugal force must equal the gravitational force for a circular orbit
    //   gravity_force = moon.mass * moon.w^2 * moon.x

double moon_r = sqrt(moon.x*moon.x + moon.y*moon.y);
double earth_r = sqrt(earth.x*earth.x + earth.y*earth.y);
    double gravity_force = G * earth.mass * moon.mass / (DIST_EARTH_MOON * DIST_EARTH_MOON);
    moon.w  = sqrt(gravity_force / (moon.mass * fabs(moon_r)));  // xxx fabs not needed
    earth.w = sqrt(gravity_force / (earth.mass * fabs(earth_r)));
    printf("moon.w          = %0.6e\n", moon.w);
    printf("earth.w         = %0.6e\n", earth.w);
    printf("orbital period  = %f days\n", TWO_PI / moon.w / 86400);
    printf("\n");
    assert(fabs(moon.w / earth.w - 1) < 1e-10);

    // init array of roughly evenly distributed locations on the earth surface
    init_earth_surface_xyz();
    //earth.max_surface = 360; //xxx

    // compute accel vectors for all earth surface locations just computed
    // - g = earth gravity
    // - m = moon gravity
    // - c = centrifugal accel
    // - t = g + m + c
    for (int i = 0; i < earth.max_surface; i++) {
        vector_t g, m, c, t;
        double d;

        if (i == 0) {
            printf("%f %f %f\n", earth.x, earth.y, earth.z);
            printf("%f %f %f\n", earth.surface[i].x, earth.surface[i].y, earth.surface[i].z);
        }

        g.a = earth.x - earth.surface[i].x;
        g.b = earth.y - earth.surface[i].y;
        g.c = earth.z - earth.surface[i].z;
        set_vector_magnitude(&g, EARTH_GRAVITY);

        m.a = moon.x - earth.surface[i].x;
        m.b = moon.y - earth.surface[i].y;
        m.c = moon.z - earth.surface[i].z;
        d = magnitude(&m);
        set_vector_magnitude(&m, G * moon.mass / square(d));

#if 1
        c.a = earth.x;
        c.b = earth.y;
        c.c = earth.z;
        d = magnitude(&c);
        set_vector_magnitude(&c,  square(earth.w) * d);
#else
        memset(&c,0,sizeof(c));
#endif

        t.a = g.a + m.a + c.a;
        t.b = g.b + m.b + c.b;
        t.c = g.c + m.c + c.c;
        earth.surface[i].g = magnitude(&t);

        if (i == theta+0 || i == theta+90 || i == theta+180 || i == theta+270) {
            if (i == 0) printf("accel vectors:\n");
            printf("%3d: earth gravity     = %+0.9f %+0.9f %+0.9f\n", i, g.a, g.b, g.c);
            printf("     moon gravity      = %+0.9f %+0.9f %+0.9f\n", m.a, m.b, m.c);
            printf("     centrifugal accel = %+0.9f %+0.9f %+0.9f\n", c.a, c.b, c.c);
            printf("     total             = %+0.9f %+0.9f %+0.9f\n", t.a, t.b, t.c);
            printf("     MAGNITUDE         = %0.9f\n", earth.surface[i].g);
            if (i == 270) printf("\n");
        }
    }

    printf("DURATION = %ld ms\n", (microsec_timer()-start)/1000);
}

static void init_earth_surface_xyz(void)
{
    #define ADD(_x,_y,_z) \
        do { \
            int idx = earth.max_surface; \
            assert(idx < MAX_EARTH_SURFACE); \
            earth.surface[idx].x = earth.x + (_x); \
            earth.surface[idx].y = earth.y + (_y); \
            earth.surface[idx].z = earth.z + (_z); \
            earth.max_surface++; \
        } while (0)

    double size, latitude, longitude, circ, x, y, z;
    int n, i;

    earth.max_surface = 0;

    size = TWO_PI * EARTH_RADIUS / 360;
    printf("init_earth_surface:\n");
    printf("  size = %.0f m   %0.6f NauticalMiles\n", size, METERS_TO_NAUTICAL_MILES(size));

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

    printf("  earth.max_surface = %d\n", earth.max_surface);
    printf("    %e  %e\n", earth.max_surface * square(size), 4 * M_PI * square(EARTH_RADIUS));
    printf("\n");
}

// -----------------  RUNTIME  -----------------------

static void *tides_thread(void *cx)
{
    int loops __attribute__((unused)) = 0;

    double theta;
    uint64_t time_now, time_last=0;


    while (true) {

        // at 80 ms intervals, advance the position of the earth and moon
        time_now = microsec_timer();
        if (time_now - time_last > 20000 || time_last == 0) {
#ifndef UNIT_TEST
            theta = (time_last == 0 ? 0 : theta+.1);
            if (time_last == 0) {
                xxx(theta);  //xxx
            }
            time_last = time_now;
#else
            theta = 0;
            if (time_last == 0) {
                xxx(theta);
            }
            time_last = time_now;
#endif
        }

        int i = random() % earth.max_surface;
        int j = random() % earth.max_surface;

        if (i == j) continue;

#if 0
        #define DELTA_R            1e-3

        double delta_pe=0, m, g_surface, r;

        m         = square(earth.surface[i].r + (DELTA_R/2));
        g_surface = earth.surface[i].g;
        r         = earth.surface[i].r + (DELTA_R/2);
        delta_pe += m * g_surface * square(r);

        m         = square(earth.surface[j].r - (DELTA_R/2));
        g_surface = earth.surface[j].g;
        r         = earth.surface[j].r - (DELTA_R/2);
        delta_pe -= m * g_surface * square(r);

        if (delta_pe < 0) {
            earth.surface[i].r += DELTA_R;
            earth.surface[j].r -= DELTA_R;
            earth.tpe += delta_pe;  // xxx del ?
        }
#else
        double delta_pe = 0;
        double m, g_surface, r;

        m         = 1;
        g_surface = earth.surface[i].g;
        r         = earth.surface[i].r + .0005;
        delta_pe += m * g_surface * square(r);

        m         = 1;
        g_surface = earth.surface[j].g;
        r         = earth.surface[j].r - .0005;
        delta_pe -= m * g_surface * square(r);

        if (delta_pe < 0) {
            earth.surface[i].r += (.001 * EARTH_RADIUS * EARTH_RADIUS) / square(earth.surface[i].r);
            earth.surface[j].r -= (.001 * EARTH_RADIUS * EARTH_RADIUS) / square(earth.surface[j].r);
            earth.tpe += delta_pe;
        }
#endif


        if ((++loops % 1000000) == 0) {
            printf("   loops = %9d   tpe = %0.9e\n", loops, earth.tpe);
        }
    }

    return NULL;
}

// -----------------  UTILS  -------------------------

static double square(double x)
{
    return x * x;
}

static double magnitude(vector_t *v)
{
    return sqrt( square(v->a) + square(v->b) + square(v->c) );
}

static int set_vector_magnitude(vector_t *v, double new_magnitude)
{
    double current_magnitude, factor;

    current_magnitude = magnitude(v);
    assert(current_magnitude > 0);

    factor = new_magnitude / current_magnitude;

    v->a *= factor;
    v->b *= factor;
    v->c *= factor;

    return 0;
}

static uint64_t microsec_timer(void)
{
    struct timespec ts;

    clock_gettime(CLOCK_MONOTONIC,&ts);
    return  ((uint64_t)ts.tv_sec * 1000000) + ((uint64_t)ts.tv_nsec / 1000);
}

// -----------------  UNIT TEST  ---------------------

#ifdef UNIT_TEST

// gcc -Wall -g -O2 -DUNIT_TEST -I. -o ut tides.c -lm -lpthread

int main(int argc, char **argv)
{
    tides_init();

    while (true) {
        sleep(1);

        printf("Results ...\n");

        double min_r=1e99, max_r=-1e99;
        int min_r_idx=-1, max_r_idx=-1;

        for (int i = 0; i < earth.max_surface; i++) {
            double r = earth.surface[i].r;

            if (i == 0 || i == 90 || i == 180 || i == 270) {
                printf("   deg = %3d   height = %0.10f  height-earth_radius = %+0.6f\n", 
                       i, r, r-earth.radius);
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
        printf("   min = %f (idx=%d)   max = %f (idx=%d)   range = %f\n", 
               min_r, min_r_idx, max_r, max_r_idx, max_r-min_r);
        printf("\n");
    }

    return 0;
}
#endif
