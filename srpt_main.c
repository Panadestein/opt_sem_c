#include <nlopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*
double opt_me(unsigned n, const double *x, double *grad)
{	
	double e_srp = [1000];

	FILE * f1 = fopen("mopac_parameter", "w");
			
	fclose(f1);

	return target;
}
*/

int main(void)
{
    int i = 0, ch = 0, ndat = 0, idxmin;
    double pdev = 0.7;
    double mineab = HUGE_VAL;
    FILE * fp;
 
    fp = fopen("./inp_ab.txt", "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);
    
    while ( (ch = fgetc(fp)) != EOF) {
    if (ch == '\n') {
    ndat++;
    }
    }

    rewind(fp);

    double data[ndat][3];
    double e_ab[ndat];

    while (i < ndat) {
    fscanf(fp, "%lf %lf %lf %lf", &data[i][0], &data[i][1], &data[i][2], &e_ab[i]);
    if (e_ab[i] < mineab) {
    mineab = e_ab[i];
    idxmin = i;
    }
    ++i;
    }

    fclose(fp);

    return 0;
}
