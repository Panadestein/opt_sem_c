# A very simple Makefile for the project

build:
	g++ srpt_main.c -o srpt -lnlopt -lm
run:
	./srpt
clean:
	rm -f ./srpt ./inp_semp/*
