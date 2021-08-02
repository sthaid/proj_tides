
tides: main.c util_sdl.c tides.c
	gcc -g -Wall -O2 -I. -I/usr/include/SDL2 $^ -lm -lpthread -lSDL2 -lSDL2_ttf -o $@

clean:
	rm -f tides
