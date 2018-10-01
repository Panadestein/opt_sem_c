# srpt_c

Main program *srpt_main.c* containing the bulk of the optimization procedure. Definition of the MOPAC geometries and automatic generation of them in the header *systems.h*, as this is the most system-dependent part of the process, it can be easily extend with some basic *C programming language* knowledge, so you can add new systems. It should be also noticed that the variable *dim* (which represents the dimension of your problem) has to be manually set up in the same header file. This behavior will be eventually improved to make it more automatic.
