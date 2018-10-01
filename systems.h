#include <stdio.h>
#include <stdlib.h>

int dim = 3;
int mysys = 0;

void gen_srpgeo(int ndat, double ** data) {
	long length;
	char * buffer = 0;

	switch (mysys) {
		//Define here  your system's input for MOPAC

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
			fprintf(fq, "pm7 charge=0 1scf EXTERNAL=mopac_parameter\n");
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
