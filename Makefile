tides: tides.c
	gcc -g -Wall -O2 -I. $^ -lm -o $@
