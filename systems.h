#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <nlopt.h>

// System definition and dimension

int mysys = 0;
int dim = 3;

// Semi-empirical method and parameter variation range

char * method = "pm7";
double pdev = 0.7;

// Algorithm parameters

nlopt_algorithm global_alg = NLOPT_GD_MLSL_LDS;
nlopt_algorithm local_alg = NLOPT_LN_BOBYQA;
int maxeval = 3;
double minrms = 0.01;
double tol = 0.001;

//Define here  your system's input for MOPAC

void gen_srpgeo(int ndat, double ** data) {
	long length;
	char * buffer = 0;
	char procedure[0x100];

	switch (mysys) {
	case 0: {
		// Ar/C10H8
		FILE * fp = fopen("./naf_geo.xyz", "r");
		if (!fp)
			exit(EXIT_FAILURE);
		fseek(fp, 0L, SEEK_END);
		length = ftell(fp);
		rewind(fp);
		buffer = (char *) malloc((length+1) * sizeof(char));
		if (buffer) {
			fread(buffer, sizeof(char), length, fp);
		}
		fclose(fp);

		for (int i = 0; i < ndat; ++i) {
			char buf[0x100];
			snprintf(buf, sizeof(buf), "./inp_semp/geo_%d.mop", i);
			FILE * fq = fopen(buf, "w");
			strcpy(procedure, method);
			strcpy(procedure, " charge=0 1scf EXTERNAL=mopac_parameter\n");
			fprintf(fq, "%s", procedure);
			fprintf(fq, "Dumb title rule\n");
			fprintf(fq, " \n");
			fprintf(fq, "Ar %f %f %f\n", data[i][0], data[i][1], data[i][2]);
			fputs(buffer, fq);
			fprintf(fq, " ");
			fclose(fq);
		}
		free(buffer);
	}

	default:
		break;
	}
}
