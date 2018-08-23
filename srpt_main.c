#include <nlopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double opt_me(unsigned n, const double *x, double *grad, void *func_data)
  {
  double e_srp = [1000];

  FILE * fs;
  fs = fopen("mopac_parameter", "w");
  if (fs == NULL) exit(EXIT_FAILURE);
  fclose(fs);

  return target;
  }

int main(void)
{
	// Input files processing and variable initialization

	int i = 0, ch = 0, pardim = 0;
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

	FILE * fr;
	fr = fopen("./parameter_pm7", "r");
	if (fr == NULL)	exit(EXIT_FAILURE);

	while ( (ch = fgetc(fn)) != EOF) {
		if (ch == '\n') {
			pardim++;
		}
	}

	rewind(fr);

	char param_names[pardim][10];
	char param_atoms[pardim][10];
	double param_values[pardim];
	double value_upper[pardim];
	double value_lower[pardim];

	i = 0;
	while (i < pardim) {
		fscanf(fr, "%s %s %lf", param_names[i], param_atoms[i],
		       &param_values[i]);
		++i;
		if (param_values[i] >= 0) {
			value_upper[i] = param_values[i] * (1.0 + pdev);
			value_lower[i] = param_values[i] * (1.0 - pdev);
		} else {
			value_upper[i] = param_values[i] * (1.0 - pdev);
			value_lower[i] = param_values[i] * (1.0 + pdev);
		}
	}

	fclose(fr);

	// Optimization process

	int maxeval = 2000;
	double minrms = 0.01;
	double tol = 0.001;
	double minf = 0.0;

	nlopt_opt opt = nlopt_create(NLOPT_G_MLSL_LDS, pardim);
	nlopt_set_local_optimizer(opt, nlopt_create(NLOPT_LN_BOBYQA, pardim));

	nlopt_set_lower_bounds(opt, value_lower);
	nlopt_set_upper_bounds(opt, value_upper);

	nlopt_set_min_objective(opt, opt_me, NULL);
	nlopt_set_maxeval(opt, maxeval);
	nlopt_set_stopval(opt, minrms);
	nlopt_set_ftol_abs(opt, tol);

	int dbg = nlopt_optimize(opt, param_values, &minf);

	if (dbg < 0) {
		fprintf(stderr, "%s:%d %s -> Nlopt C function failed: %d expected: %d\n",
		        __FILE__, __LINE__, __FUNCTION__, dbg, NLOPT_SUCCESS);
	} else {
		printf("minimum: f(%lf, %lf) = %lf\n",
		       param_values[0], param_values[1], minf);
	}

	// Cleaning up stuff

	nlopt_destroy(opt);

	for (i = 0; i < ndat; i++) {
		free(data[i]);
	}

	free(buffer);
	free(data);

	return 0;
}

