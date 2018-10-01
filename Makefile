# A simple Makefile for the project

all: srpt

srpt: srpt_main.c
	gcc -fopenmp -g $< -o $@ -lnlopt -lm

.PHONY: run clean

run: srpt
	./srpt

clean:
	rm -f ./srpt ./inp_semp/* mopac_parameter rms_values e_ab e_srp
