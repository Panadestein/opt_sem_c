#include <nlopt.h>
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
    double pdev = 0.7;
    FILE * fp;
    int n = 0, ch = 0, ndat = 0;
 
    fp = fopen("./inp_ab.txt", "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);
    
    while ( (ch = fgetc(fp)) != EOF) {
    if (ch == '\n') {
    ndat++;
    }
    }

    rewind(fp);

    double data[ndat][4];

    while (n < ndat) {
    fscanf(fp, "%lf %lf %lf %lf", &data[n][0], &data[n][1], &data[n][2], &data[n][3]);
    printf("%le %le %le %le\n", data[n][0], data[n][1], data[n][2], data[n][3]);
    ++n;
    }

    fclose(fp);

    return 0;
}
