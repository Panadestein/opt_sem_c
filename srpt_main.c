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

int main() {
        double pdev = 0.7;
        int ndat = 0;
        int ch = 0;
        int n = 0;
	
	FILE * f2 = fopen("inp_ab.txt", "r");

	if (!f2)		
	{
        printf("Error: file could not be opened\n");
        return EXIT_FAILURE;
        }

        while (!feof(f2))
        {
        ch = fgetc(f2);
        if (ch == '\n')
        {
        ndat++;
        }
        }
 
	double data[ndat][4];

        while ( (n < ndat) && (!feof(f2)))
        {
        fscanf(f2, "%le %le %le %le", &data[n][0], &data[n][1], &data[n][2], &data[n][3]);
        ++n;
        }

        if ( (!feof(f2)) && (n == ndat))
        {
        printf("Error: file too large for array\n");
        }

	fclose(f2);

        printf("%le\n", data[0][0]);
    
        for (int i = 0; i < n; ++i)
	{
	printf("%d: %d, %d, %d, %d\n", i, data[n][0], data[n][1], data[n][2], data[n][3]);
	}

	return 0;
}

