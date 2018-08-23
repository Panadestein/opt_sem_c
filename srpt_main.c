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
	int i = 0, ch = 0;
	int dim = 3, ndat = 0, idxmin;
	long length;
	double pdev = 0.7;
	double ** data;
	double mineab = HUGE_VAL;
	char * buffer = 0;

	FILE * fn;
	fn = fopen("./inp_ab.txt", "r");
	if (fn == NULL)	exit(EXIT_FAILURE);

	while ( (ch = fgetc(fn)) != EOF) {
		if (ch == '\n') {
			ndat++;
		}
	}

	rewind(fn);

	double e_ab[ndat];
	data = (double **)  malloc(ndat * sizeof(double));
	for (int i = 0; i < ndat; i++) {
		data[i] = (double *) malloc(dim * sizeof(double));
	}

	if (data == NULL) {
		printf("Error: memory not allocated\n");
		exit(0);
	}

	while (i < ndat) {
		fscanf(fn, "%lf %lf %lf %lf", &data[i][0], &data[i][1], &data[i][2],
		       &e_ab[i]);
		if (e_ab[i] < mineab) {
			mineab = e_ab[i];
			idxmin = i;
		}
		++i;
	}

	fclose(fn);

	for (i = 0; i < ndat; i++) {
		e_ab[i] -= mineab;
	}

	FILE * fp = fopen("./naf_geo.xyz", "r");
	if (fp == NULL) exit(EXIT_FAILURE);
	fseek(fp, 0L, SEEK_END);
	length = ftell(fp);
	rewind(fp);
	buffer = (char *) malloc((length+1) * sizeof(char));
	if (buffer) {
		fread(buffer, sizeof(char), length, fp);
	}
	fclose(fp);

	for (i = 0; i < ndat; ++i) {
		char buf[0x100];
		snprintf(buf, sizeof(buf), "./inp_semp/geo_%d.mop", i);
		FILE * fq = fopen(buf, "w");
		fprintf(fq, "pm7 charge=0 1scf EXTERNAL=mopac_parameter\n");
		fprintf(fq, "Dumb title rule\n");
		fprintf(fq, " \n");
		fprintf(fq, "Ar %f %f %f\n", data[i][0], data[i][1], data[i][2]);
		fputs(buffer, fq);
		fprintf(fq, " ");
		fclose(fq);
	}

	for (int i = 0; i < ndat; i++) {
		free(data[i]);
	}

	free(buffer);
	free(data);
	return 0;
}

