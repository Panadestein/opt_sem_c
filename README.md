# srpt_c

The main program (*srpt_main.c*) contains the bulk of the optimization procedure. Definition of the MOPAC geometries for the system of interest and automatic generation of them can be set up in the header file (*systems.h*), as this is the most system-dependent part of the process, it can be easily extended with some basic *C programming language* knowledge, so you can add new systems. 

It should be also noticed that the variables *dim* (which represents the dimension of your problem) and *method* (a string representing the semiempirical method of choice) have to be manually set up in the same header file. Please make sure that the file *./parameter_ref* contains the default parameters of the semiempirical method. This behavior will be eventually improved to make it more automatic.

This package depends on the non-linear optimization library [NLopt](https://github.com/stevengj/nlopt), so please refer to its reference for information about the algorithms, termination conditions and other fit related questions.
