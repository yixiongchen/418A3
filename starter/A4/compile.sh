#!/bin/sh
gcc -O2 -g -fopenmp svdDynamic.c RayTracer.c utils.c -lm -o RayTracer
