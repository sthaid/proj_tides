
tides: main.c util_sdl.c
	gcc -g -Wall -O2 -I. -I/usr/include/SDL2 $^ -lm -lSDL2 -lSDL2_ttf -o $@

clean:
	rm -f tides
