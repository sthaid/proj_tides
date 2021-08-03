
tides: main.c util_sdl.c tides.c
	gcc -g -Wall -O2 -I. -I/usr/include/SDL2 $^ -lm -lpthread -lSDL2 -lSDL2_ttf -o $@

ut: tides.c
	gcc -Wall -g -O2 -DUNIT_TEST -I. $^ -lm -lpthread -o $@

clean:
	rm -f tides ut
