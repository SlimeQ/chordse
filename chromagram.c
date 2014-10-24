#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <omp.h>

double max(double *input, int len) {
	double max = 0;
	int i;
	for (i=0; i<len; i++) {
		if (input[i] > max) max = input[i];
	}
	return max;
}
double min(double *input, int len) {
	double min = input[0];
	int i;
	for (i=1; i<len; i++) {
		if (input[i]  < min) min = input[i];
	}
	return min;
}

int K_thread_count = 2;

void chromagram(short int data[], int len, int f_min, int f_max, int M, int B, int K, float Q, int fs, double * out) {
	printf("chromagram\n");
	printf("----------\n");

	// perform constant Q transform
	int k, n;
	double n_over_Nk, Wkn;
	double alpha = 25/46;
	double *Xcq = malloc(sizeof(double) * K);

	double cq_t1 = omp_get_wtime();

	#pragma omp parallel for num_threads(K_thread_count)
	for(k=0; k<K; k++) {
		float fk = exp2(k/B) * f_min;
		int Nk = ceil(Q * fs/fk);
		if (Nk > len) Nk = len;

		double sum = 0;
		for(n=0; n<Nk; n++) {
			n_over_Nk = n/(double)Nk;
			Wkn = alpha - (1 - alpha) * cos(2 * M_PI * n_over_Nk);
			sum += Wkn * data[n] * exp(-2 * I * M_PI * Q * n_over_Nk);
		}
		Xcq[k] = sum;
		// printf("%d = %f\n", k, sum);
	}

	double cq_t2 = omp_get_wtime();
	double cq_elapsed = cq_t2 - cq_t1;
	printf("constant Q time = %f\n", cq_elapsed);

	// compute chromagram
	int b, m;
	double ch_t1 = omp_get_wtime();
	for(b=0; b<B; b++) {
		double sum = 0;
		for(m=0; m<M; m++) {
			sum += fabs(Xcq[b+m*B]);
		}
		out[b] = sum;
		// printf("%d = %f\n", b, sum);
	}

	// normalize
	double minValue = min(out, B);
	for (b=0; b<B; b++) {
		out[b] -= minValue;
	}
	double maxValue = max(out, B);
	for (b=0; b<B; b++) {
		out[b] /= maxValue;
		if (isnan(out[b])) out[b] = 0;
		printf("%d = %f\n", b, out[b]);
	}

	double ch_t2 = omp_get_wtime();
	double ch_elapsed = ch_t2 - ch_t1;
	printf("chromagram loop = %f\n", ch_elapsed);
}

void multichroma(int * data, int totalLen, int numChunks, int N, int f_min, int f_max, int M, int B, int K, float Q, int fs, double * out) {
	
	printf("in multichroma\n");
	printf("data[0] = %d\n", data[0]);
	short int **chunks = malloc(sizeof(short int *) * numChunks);
	printf("chunks is initialized\n");
	double chroma[numChunks][B];
	printf("stuff is initialized\n");
	int i, n;
	for (i=0; i<numChunks; i++) {
		chunks[i] = malloc(sizeof(short int) * N);
		for (n=0; n<N; n++) {
			if (i * N + n < totalLen) {
				printf("setting values...\n");
				chunks[i][n] = data[i * N + n];
				printf("minor success!\n");
			} else {
				chunks[i][n] = 0;
			}
		}
		chromagram((short int *)chunks[i], N, f_min, f_max, M, B, K, Q, fs, chroma[i]);
	}
}