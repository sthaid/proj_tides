#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>

#include <util_logging.h>  // xxx maybe just printf here
#include <util_sdl.h>
#include <vectors.h>
#include <tides.h>

#define REQUESTE_WIN_WIDTH  1500
#define REQUESTE_WIN_HEIGHT 900

static struct ctrls_s default_ctrls = {
    .reset_request           = false,
    .moon_motion_enabled     = false,
    .vectors_display_enabled = true,
    .moon_enabled            = true,
    .sun_enabled             = false,
};

static int pane_hndlr(pane_cx_t * pane_cx, int request, void * init_params, sdl_event_t * event);

// --------------------------------------------------------------------------------

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
        FATAL("sdl_init %dx%d failed\n", win_width, win_height);
    }
    INFO("requested win_width=%d win_height=%d\n", REQUESTE_WIN_WIDTH, REQUESTE_WIN_HEIGHT);
    INFO("actual    win_width=%d win_height=%d\n", win_width, win_height);

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

// --------------------------------------------------------------------------------

static int pane_hndlr(pane_cx_t * pane_cx, int request, void * init_params, sdl_event_t * event)
{
    static rect_t  * pane;
    static double    esf;
    static double    msf;
    static texture_t earth_texture;
    static texture_t moon_texture;

    #define FONTSZ  40
    #define SDL_EVENT_RESET             (SDL_EVENT_USER_DEFINED + 1)
    #define SDL_EVENT_MOON_MOTION       (SDL_EVENT_USER_DEFINED + 2)
    #define SDL_EVENT_DISPLAY_VECTORS   (SDL_EVENT_USER_DEFINED + 3)
    #define SDL_EVENT_MOON              (SDL_EVENT_USER_DEFINED + 4)
    #define SDL_EVENT_SUN               (SDL_EVENT_USER_DEFINED + 5)

    // ----------------------------
    // -------- INITIALIZE --------
    // ----------------------------

    if (request == PANE_HANDLER_REQ_INITIALIZE) {
        DEBUG("PANE x,y,w,h  %d %d %d %d\n", pane->x, pane->y, pane->w, pane->h);

        pane = &pane_cx->pane;
        esf = 120 / 6.4e6;
        earth_texture = sdl_create_filled_circle_texture(EARTH_RADIUS * esf, WHITE);

        //msf = 500 / 4.0e8;
        msf = 430 / 3.797e8;
        moon_texture = sdl_create_filled_circle_texture(20, WHITE);

        return PANE_HANDLER_RET_NO_ACTION;
    }

    // ------------------------
    // -------- RENDER --------
    // ------------------------

    if (request == PANE_HANDLER_REQ_RENDER) {
        //rect_t loc = {pane->w/2-100, pane->h/2-100, 200, 200};
        //sdl_render_rect(pane, &loc, 1, WHITE);

        int x_ctr = pane->h/2;
        int y_ctr = pane->h/2;
        int w     = pane->h;
        int h     = pane->h;

int earth_texture_width, earth_texture_height;
int moon_texture_width, moon_texture_height;
        sdl_render_line(pane, 0, y_ctr, w-1, y_ctr, WHITE);
        sdl_render_line(pane, x_ctr, 0, x_ctr, h-1, WHITE);

        sdl_render_line(pane, w-1, 0, w-1, h-1, WHITE);

        sdl_query_texture(earth_texture, &earth_texture_width, &earth_texture_height);
        sdl_render_texture(pane, 
                           x_ctr - earth_texture_width/2 + esf * earth.x,
                           y_ctr - earth_texture_height/2 - esf * earth.y,
                           earth_texture);

        sdl_render_point(pane, x_ctr, y_ctr, BLACK, 3);

#if 0
        sdl_render_point(pane, 
                         x_ctr + esf * earth.x,
                         y_ctr - esf * earth.y,
                         GREEN, 3);
#endif

        if (ctrls.moon_enabled) {
            sdl_query_texture(moon_texture, &moon_texture_width, &moon_texture_height);
            sdl_render_texture(pane, 
                            x_ctr - moon_texture_width/2 + msf * moon.x,
                            y_ctr - moon_texture_height/2 - msf * moon.y,
                            moon_texture);
        }

        // xxx fill in
        int count=0;
        point_t points[361];
        double k = 2e6;
        for (int i = 0; i < 360; i++) {
            points[count].x = x_ctr + nearbyint(((EARTH_RADIUS + k*(earth.surface[i].r-EARTH_RADIUS+0.5)) * cos(DEG_TO_RAD(i)) + earth.x) * esf);
            points[count].y = y_ctr - nearbyint(((EARTH_RADIUS + k*(earth.surface[i].r-EARTH_RADIUS+0.5)) * sin(DEG_TO_RAD(i)) + earth.y) * esf);
            count++;
        }
        points[360] = points[0];
        count++;
        sdl_render_lines(pane, points, count, LIGHT_BLUE);

        if (ctrls.vectors_display_enabled) {
            double x1, x2,y1, y2;
            for (int i = 0; i < 360; i += 10) {
                x1 = x_ctr + (((EARTH_RADIUS + k*(earth.surface[i].r-EARTH_RADIUS+0.5)) * cos(DEG_TO_RAD(i)) + earth.x) * esf);
                y1 = y_ctr - (((EARTH_RADIUS + k*(earth.surface[i].r-EARTH_RADIUS+0.5)) * sin(DEG_TO_RAD(i)) + earth.y) * esf);

                vector_t *v = &earth.surface[i].v;
                x2 = x1 + (v->a * 100000000.);
                y2 = y1 - (v->b * 100000000.);

                // xxx nearbyint
                sdl_render_line(pane, x1,y1, x2,y2, RED);
                sdl_render_point(pane, x2,y2, RED, 2);
            }
        }
// xxx display theta and tides

        sdl_render_text_and_register_event(
                pane, w+10, ROW2Y(0,FONTSZ), FONTSZ, 
                "RESET", 
                LIGHT_BLUE, BLACK,
                SDL_EVENT_RESET, SDL_EVENT_TYPE_MOUSE_CLICK, pane_cx);

        sdl_render_text_and_register_event(
                pane, w+10, ROW2Y(2,FONTSZ), FONTSZ, 
                ctrls.moon_motion_enabled ? "MOTION IS ON" : "MOTION IS OFF",
                LIGHT_BLUE, BLACK,
                SDL_EVENT_MOON_MOTION, SDL_EVENT_TYPE_MOUSE_CLICK, pane_cx);

        sdl_render_text_and_register_event(
                pane, w+10, ROW2Y(4,FONTSZ), FONTSZ, 
                ctrls.vectors_display_enabled ? "VECTORS ARE ON" : "VECTORS ARE OFF",
                LIGHT_BLUE, BLACK,
                SDL_EVENT_DISPLAY_VECTORS, SDL_EVENT_TYPE_MOUSE_CLICK, pane_cx);

        sdl_render_text_and_register_event(
                pane, w+10, ROW2Y(6,FONTSZ), FONTSZ, 
                ctrls.moon_enabled ? "MOON IS ENABLED" : "MOON IS DISABLED",
                LIGHT_BLUE, BLACK,
                SDL_EVENT_MOON, SDL_EVENT_TYPE_MOUSE_CLICK, pane_cx);

        sdl_render_text_and_register_event(
                pane, w+10, ROW2Y(8,FONTSZ), FONTSZ, 
                ctrls.sun_enabled ? "SUN IS ENABLED" : "SUN IS DISABLED",
                LIGHT_BLUE, BLACK,
                SDL_EVENT_SUN, SDL_EVENT_TYPE_MOUSE_CLICK, pane_cx);

        return PANE_HANDLER_RET_NO_ACTION;
    }

    // -----------------------
    // -------- EVENT --------
    // -----------------------

    if (request == PANE_HANDLER_REQ_EVENT) {
        int rc = PANE_HANDLER_RET_NO_ACTION;

        // first handle events common to all displays
        switch (event->event_id) {
        case 'q':  // quit
            rc = PANE_HANDLER_RET_PANE_TERMINATE;
            break;
        case SDL_EVENT_RESET:
            ctrls = default_ctrls;
            ctrls.reset_request = true;
            break;
        case SDL_EVENT_MOON_MOTION:
            ctrls.moon_motion_enabled = !ctrls.moon_motion_enabled;
            break;
        case SDL_EVENT_DISPLAY_VECTORS:
            ctrls.vectors_display_enabled = !ctrls.vectors_display_enabled;
            break;
        case SDL_EVENT_MOON:
            ctrls.moon_enabled = !ctrls.moon_enabled;
            break;
        case SDL_EVENT_SUN:
            ctrls.sun_enabled = !ctrls.sun_enabled;
            break;
// xxx key ctrls + - m s v r
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

