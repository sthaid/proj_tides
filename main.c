#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>

#include <util_logging.h>
#include <util_sdl.h>
#include <tides.h>

#define REQUESTE_WIN_WIDTH  1200
#define REQUESTE_WIN_HEIGHT 1000

static int pane_hndlr(pane_cx_t * pane_cx, int request, void * init_params, sdl_event_t * event);

// --------------------------------------------------------------------------------

int main(int argc, char **argv)
{
    int win_width  = REQUESTE_WIN_WIDTH;
    int win_height = REQUESTE_WIN_HEIGHT;

    INFO("hello\n");

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
        50000,          // 0=continuous, -1=never, else us
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

    // ----------------------------
    // -------- INITIALIZE --------
    // ----------------------------

    if (request == PANE_HANDLER_REQ_INITIALIZE) {
        DEBUG("PANE x,y,w,h  %d %d %d %d\n", pane->x, pane->y, pane->w, pane->h);

        pane = &pane_cx->pane;
        esf = 120 / 6.4e6;
        earth_texture = sdl_create_filled_circle_texture(EARTH_RADIUS * esf, WHITE);

        //msf = 500 / 4.0e8;
        msf = 500 / 3.797e8;
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

        sdl_render_point(pane, x_ctr, y_ctr, BLACK, 5);
        sdl_render_point(pane, 
                         x_ctr + esf * earth.x,
                         y_ctr - esf * earth.y,
                         GREEN, 5);

        sdl_query_texture(moon_texture, &moon_texture_width, &moon_texture_height);
        sdl_render_texture(pane, 
                           x_ctr - moon_texture_width/2 + msf * moon.x,
                           y_ctr - moon_texture_height/2 - msf * moon.y,
                           moon_texture);


        // xxx use lines, and fill in
        for (int i = 0; i < 360; i++) {
// xxx nearbyint
            int x = x_ctr + ((EARTH_RADIUS + 1e7*(earth.surface[i%360].r-EARTH_RADIUS)) * cos(DEG_TO_RAD(i)) + earth.x) * esf;
            int y = y_ctr - ((EARTH_RADIUS + 1e7*(earth.surface[i%360].r-EARTH_RADIUS)) * sin(DEG_TO_RAD(i)) + earth.y) * esf;
            sdl_render_point(pane, x, y, LIGHT_BLUE, 2);
        }

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

