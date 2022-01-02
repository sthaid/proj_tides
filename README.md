# Tidal Bulges

The Earth's tides are caused by the Moon's gravity.

There is a tidal bulge on both the side of the Earth facing the Moon, and
also on the opposite side of the Earth facing away from the Moon.

It seems that the Moon's gravity should be pulling the water towards the Moon.
This would account for the tidal bulge on the side of the Earth facing the Moon.
But, why is there also a tidal bulge on the opposite side of the Earth?

The answer to that question has to do with the centrifugal force caused by the
rotation of the Earth around the Earth-Moon barycenter. Refer to notes.txt for more.

# Tides Program

The tides program simulates the tidal bulges by minimizing gravitational potential energy.
This is done by repeatedly choosing a pair of random locations on the Earth's
surface and checking if moving a small quantity of water from one location to the other
will reduce potential energy. If so, the height of the water at the two locations
is adjusted to simulate the movement of the water.

# Building Tides

I build and run this on 64 bit Ubuntu.
Required packages include:libsdl2-dev and libsdl2-ttf-dev.

I have not built or tested on 32 bit Linux.

# Program Controls

* MOTION or 'M': enable/disable Moon motion
* VECTORS or 'v': show/hide the acceleration vectors
* MOON or 'm': enable/disable the Moon
* SUN or 's': enable/disable the Sun
* GEOMETRY or 'g': choose between Disk or Spherical Earth simulation
* RESET or 'r': reset to initial state
* '+' or '-':  move Moon by +/- 1 degree
* 'q':  exit program

The vectors displayed include the acceleration from the Moon and Sun's gravity, and
from the centrifugal acceleration of the Earth around the Earth-Moon and Earth-Sun
barycenters. The vectors do not include the Earth's gravity.

# Screenshot

![screenshot.png](/assets/screenshot.png)

