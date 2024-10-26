#!/bin/bash

gcc -o main bird.c sampa.c spa.c main.c -lm
./main
