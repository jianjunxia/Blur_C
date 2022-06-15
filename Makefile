# make
fastblur_c: main.c fast_gaussian_blur_template.h
	gcc main.c -g -o fastblur_c -lm

all: fastblur_c

clean:
	rm fastblur_c
