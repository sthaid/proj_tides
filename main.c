#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>

#include <util_sdl.h>
#include <vectors.h>
#include <tides.h>

//
// defines
//

#define REQUESTE_WIN_WIDTH  1500
#define REQUESTE_WIN_HEIGHT 900

//
// variables
//

static const struct ctrls_s default_ctrls = {
    .theta         = 0,
    .moon_enabled  = true,
    .sun_enabled   = false,
    .motion        = false,
    .vectors       = true,
    .geometry      = GEO_DISK,
};

//
// prototypes
//

static int pane_hndlr(pane_cx_t * pane_cx, int request, void * init_params, sdl_event_t * event);
static double sanitize_angle(double x);

// -----------------  MAIN  -------------------------------------------------------

int main(int argc, char **argv)
{
    int win_width  = REQUESTE_WIN_WIDTH;
    int win_height = REQUESTE_WIN_HEIGHT;

    // set control params to their default values
    ctrls = default_ctrls;

    // init tides
    tides_init();

    // init sdl
    if (sdl_init(&win_width, &win_height, false, false) < 0) {
        printf("ERROR: sdl_init %dx%d failed\n", win_width, win_height);
        return 1;
    }
    printf("requested win_width=%d win_height=%d\n", REQUESTE_WIN_WIDTH, REQUESTE_WIN_HEIGHT);
    printf("actual    win_width=%d win_height=%d\n", win_width, win_height);

    // run the pane manger;
    // the sdl_pane_manager is the runtime loop, and it will repeatedly call the pane_hndlr,
    //  first to initialize the pane_hndlr and subsequently to render the display and
    //  process events
    sdl_pane_manager(
        NULL,           // context
        NULL,           // called prior to pane handlers
        NULL,           // called after pane handlers
        1000,           // 0=continuous, -1=never, else us
        1,              // number of pane handler varargs that follow
        pane_hndlr, NULL, 0, 0, win_width, win_height, PANE_BORDER_STYLE_NONE);

    // done
    return 0;
}

// -----------------  PANE HNDLR  -------------------------------------------------

// this routine has 4 secions, which are invoked:
// - INITIALIZE: once
// - RENDER:     periodically to draw the display
// - EVENT:      when an event is available
// - TERMINATE:  when done

static int pane_hndlr(pane_cx_t * pane_cx, int request, void * init_params, sdl_event_t * event)
{
    static double    esf;
    static double    msf;
    static double    ssf;
    static texture_t earth_texture;
    static texture_t moon_texture;
    static texture_t sun_texture;

    #define FONTSZ  40
    #define SDL_EVENT_RESET     (SDL_EVENT_USER_DEFINED + 1)
    #define SDL_EVENT_MOTION    (SDL_EVENT_USER_DEFINED + 2)
    #define SDL_EVENT_MOON      (SDL_EVENT_USER_DEFINED + 3)
    #define SDL_EVENT_SUN       (SDL_EVENT_USER_DEFINED + 4)
    #define SDL_EVENT_VECTORS   (SDL_EVENT_USER_DEFINED + 5)
    #define SDL_EVENT_GEOMETRY  (SDL_EVENT_USER_DEFINED + 6)

    // ----------------------------
    // -------- INITIALIZE --------
    // ----------------------------

    if (request == PANE_HANDLER_REQ_INITIALIZE) {
        printf("PANE x,y,w,h  %d %d %d %d\n",
               pane_cx->pane.x, pane_cx->pane.y, pane_cx->pane.w, pane_cx->pane.h);

        // these scale factors are used to convert from meters to pixels for
        // the earth, moon, and sun 
        esf = 100 / EARTH_RADIUS;
        msf = 350 / moon.em_orbit_radius;
        ssf = 610 / DIST_EARTH_SUN;

        // init textures, containing circles representing the earth, moon, and sun;
        // these textures are used in the RENDER section below
        earth_texture = sdl_create_filled_circle_texture(EARTH_RADIUS * esf, WHITE);
        moon_texture = sdl_create_filled_circle_texture(10, WHITE);
        sun_texture = sdl_create_filled_circle_texture(200, YELLOW);

        return PANE_HANDLER_RET_NO_ACTION;
    }

    // ------------------------
    // -------- RENDER --------
    // ------------------------

    if (request == PANE_HANDLER_REQ_RENDER) {
        rect_t *pane  = &pane_cx->pane;
        int     x_ctr = pane->h/2;
        int     y_ctr = pane->h/2;
        int     w     = pane->h;
        int     h     = pane->h;
        int     moon_texture_width, moon_texture_height;
        int     earth_texture_width, earth_texture_height;
        int     sun_texture_width, sun_texture_height;
        double  min, max;

        // get the min and max radius along the equator, 
        // where this radius is the EARTH_RADIUS + tide height
        tides_get_min_max(&min, &max, NULL, NULL);

        // draw X and Y axis; the origin is the earth-moon barycenter
        sdl_render_line(pane, 0, y_ctr, w-1, y_ctr, WHITE);
        sdl_render_line(pane, x_ctr, 0, x_ctr, h-1, WHITE);

        // draw the earth
        sdl_query_texture(earth_texture, &earth_texture_width, &earth_texture_height);
        sdl_render_texture(pane, 
                           (x_ctr - earth_texture_width/2)  + esf * earth.x,
                           (y_ctr - earth_texture_height/2) - esf * earth.y,
                           earth_texture);

        // draw a point at the earth-moon barycenter (which is same as X,Y origin)
        sdl_render_point(pane, x_ctr, y_ctr, BLACK, 3);

        // if moon simulation is enabled then draw the moon
        if (ctrls.moon_enabled) {
            sdl_query_texture(moon_texture, &moon_texture_width, &moon_texture_height);
            sdl_render_texture(pane, 
                            (x_ctr - moon_texture_width/2)  + msf * moon.x,
                            (y_ctr - moon_texture_height/2) - msf * moon.y,
                            moon_texture);
        }

        // if sun simulation is enabled then draw the sun
        if (ctrls.sun_enabled) {
            sdl_query_texture(sun_texture, &sun_texture_width, &sun_texture_height);
            sdl_render_texture(pane, 
                            (x_ctr - sun_texture_width/2)  + ssf * sun.x,
                            (y_ctr - sun_texture_height/2) - ssf * sun.y,
                            sun_texture);
        }

        // draw a blue area representing the height of the tides at the equator
        int count=0;
        point_t points[361];
        for (int i = 0; i < 360; i++) {
            double r = EARTH_RADIUS + 400000 + (earth.surface[i].R - min) * 4000000;
            points[count].x = x_ctr + nearbyint(esf * (r * cos(DEG_TO_RAD(i)) + earth.x));
            points[count].y = y_ctr - nearbyint(esf * (r * sin(DEG_TO_RAD(i)) + earth.y));
            count++;
        }
        points[360] = points[0];
        count++;
        sdl_render_lines(pane, points, count, LIGHT_BLUE);

        // if vectors are enabled then draw the vectors at 10 degree increments;
        // these vectors represent the sum of the accelerations:
        //  - moon gravity
        //  - earth's centrifigul force around the earth-moon barycenter
        //  - sun gravity
        //  - earth's centrifigul force around the earth-sun barycenter
        // the acceleration of earth gravity is not included in the earth.surface[i].v
        // because it's value is comparitively enourmous
        if (ctrls.vectors) {
            double x1, x2,y1, y2;
            for (int i = 0; i < 360; i += 10) {
                double r = EARTH_RADIUS + 400000 + (earth.surface[i].R - min) * 4000000;
                x1 = x_ctr + nearbyint(esf * (r * cos(DEG_TO_RAD(i)) + earth.x));
                y1 = y_ctr - nearbyint(esf * (r * sin(DEG_TO_RAD(i)) + earth.y));

                vector_t *v = &earth.surface[i].v;
                x2 = x1 + nearbyint(v->a * 100000000.);
                y2 = y1 - nearbyint(v->b * 100000000.);

                sdl_render_line(pane, x1,y1, x2,y2, RED);
                sdl_render_point(pane, x2,y2, RED, 2);
            }
        }

        // draw vertical line marking border between the display and controls areas
        sdl_render_line(pane, w-1, 0, w-1, h-1, WHITE);

        // register for events
        int x_ctrls = w + 10;
        sdl_render_text_and_register_event(
                pane, x_ctrls+COL2X(0,FONTSZ), ROW2Y(0,FONTSZ), FONTSZ, 
                "RESET", 
                LIGHT_BLUE, BLACK,
                SDL_EVENT_RESET, SDL_EVENT_TYPE_MOUSE_CLICK, pane_cx);
        sdl_render_text_and_register_event(
                pane, x_ctrls+COL2X(10,FONTSZ), ROW2Y(0,FONTSZ), FONTSZ, 
                "MOTION",
                LIGHT_BLUE, BLACK,
                SDL_EVENT_MOTION, SDL_EVENT_TYPE_MOUSE_CLICK, pane_cx);
        sdl_render_text_and_register_event(
                pane, x_ctrls+COL2X(0,FONTSZ), ROW2Y(2,FONTSZ), FONTSZ, 
                "MOON",
                LIGHT_BLUE, BLACK,
                SDL_EVENT_MOON, SDL_EVENT_TYPE_MOUSE_CLICK, pane_cx);
        sdl_render_text_and_register_event(
                pane, x_ctrls+COL2X(10,FONTSZ), ROW2Y(2,FONTSZ), FONTSZ, 
                "SUN",
                LIGHT_BLUE, BLACK,
                SDL_EVENT_SUN, SDL_EVENT_TYPE_MOUSE_CLICK, pane_cx);
        sdl_render_text_and_register_event(
                pane, x_ctrls+COL2X(0,FONTSZ), ROW2Y(4,FONTSZ), FONTSZ, 
                "VECTORS",
                LIGHT_BLUE, BLACK,
                SDL_EVENT_VECTORS, SDL_EVENT_TYPE_MOUSE_CLICK, pane_cx);
        sdl_render_text_and_register_event(
                pane, x_ctrls+COL2X(10,FONTSZ), ROW2Y(4,FONTSZ), FONTSZ, 
                "GEOMETRY",
                LIGHT_BLUE, BLACK,
                SDL_EVENT_GEOMETRY, SDL_EVENT_TYPE_MOUSE_CLICK, pane_cx);

        // display status
        sdl_render_printf(pane, x_ctrls, ROW2Y(8,FONTSZ), FONTSZ, WHITE, BLACK,
                          "GEO   : %s", 
                          ctrls.geometry == GEO_DISK ? "DISK" : "SPHERE");
        sdl_render_printf(pane, x_ctrls, ROW2Y(9,FONTSZ), FONTSZ, WHITE, BLACK,
                          "MOON  : %s", 
                          ctrls.moon_enabled ? "ON" : "OFF");
        sdl_render_printf(pane, x_ctrls, ROW2Y(10,FONTSZ), FONTSZ, WHITE, BLACK,
                          "SUN   : %s", 
                          ctrls.sun_enabled ? "ON" : "OFF");
        sdl_render_printf(pane, x_ctrls, ROW2Y(11,FONTSZ), FONTSZ, WHITE, BLACK,
                          "THETA : %0.1f", 
                          ctrls.theta);
        sdl_render_printf(pane, x_ctrls, ROW2Y(12,FONTSZ), FONTSZ, WHITE, BLACK,
                          "RANGE : %0.3f ft", 
                          METERS_TO_FEET(max-min));

        return PANE_HANDLER_RET_NO_ACTION;
    }

    // -----------------------
    // -------- EVENT --------
    // -----------------------

    if (request == PANE_HANDLER_REQ_EVENT) {
        int rc = PANE_HANDLER_RET_NO_ACTION;

        switch (event->event_id) {
        case 'q':
            rc = PANE_HANDLER_RET_PANE_TERMINATE;
            break;
        case SDL_EVENT_RESET: case 'r':
            ctrls = default_ctrls;
            break;
        case SDL_EVENT_MOTION: case 'M':
            ctrls.motion = !ctrls.motion;
            break;
        case SDL_EVENT_VECTORS: case 'v':
            ctrls.vectors = !ctrls.vectors;
            break;
        case SDL_EVENT_MOON: case 'm':
            ctrls.moon_enabled = !ctrls.moon_enabled;
            if (ctrls.moon_enabled == false) {
                ctrls.motion = false;
            }
            break;
        case SDL_EVENT_SUN: case 's':
            ctrls.sun_enabled = !ctrls.sun_enabled;
            break;
        case SDL_EVENT_GEOMETRY: case 'g':
            ctrls.geometry = (ctrls.geometry == GEO_DISK ? GEO_SPHERE : GEO_DISK);
            break;
        case '+': case '=':
            if (ctrls.moon_enabled) {
                ctrls.theta = sanitize_angle(nearbyint(ctrls.theta+1));
            }
            break;
        case '-':
            if (ctrls.moon_enabled) {
                ctrls.theta = sanitize_angle(nearbyint(ctrls.theta-1));
            }
            break;
        }

        return rc;
    }

    // ---------------------------
    // -------- TERMINATE --------
    // ---------------------------

    if (request == PANE_HANDLER_REQ_TERMINATE) {
        return PANE_HANDLER_RET_NO_ACTION;
    }

    // not reached
    assert(0);
    return PANE_HANDLER_RET_NO_ACTION;
}

static double sanitize_angle(double x)
{
    if (x < 0 || x >= 360) {
        while (x < 0) x += 360;
        while (x >= 360) x -= 360;
    }
    return x;
}
