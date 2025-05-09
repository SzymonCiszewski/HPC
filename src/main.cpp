#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "benchmark/benchmark.h"


void macierz_ukladu(double delta_t, double mi, double h, double* r, int N, double** A);

void wektor_prawych_stron(double a, double T, double delta_t, double t, double* u0, int N, double* b);

void rozwiazanie_calkowite(double* u, int N, double* u_c);

void srednia_predkosc(double* u, double* r, int n, double* u_sr, double h);

void przekroj_rury(double r, double R, double* A);

void wydatek(double u_sr, double A, double* Q);

void trapz(double* x, double h, int n, double *y);

void Thomas(double** A, double* x, double* b, int N);



int main() {
	FILE* fA, * fb, * fu, * fr, * fu_srQ;
	fA = fopen("macierz_ukladu.txt", "w");
	fb = fopen("wektor_prawych_stron.txt", "w");
	fu = fopen("wektor_niewiadomych.txt", "w");
	fr = fopen("wektor_r_wezlow.txt", "w");
	fu_srQ = fopen("czas_srednia_predkosc_wydatek.txt", "w");

	// srednica rury [m]
	double D;
	// srednica rdzenia [m]
	double d;
	//krok miedzy wezlami [m]
	double h = 0.01;
	// pole przekroj rury [m^2]
	double A_r;
	// amplituda sily [N]
	double a;
	// rodzaj cieczy
	unsigned int ciecz;
	// okres ruchu [s]
	double T;
	//czas [s]
	double t = 0;
	//krok czasowy [s]
	double delta_t; 
	// czas trwania symulacji [s]
	double tk;
	// lepkosć dynamiczna cieczy [Pa *s]
	double mi;
	printf("Podaj srednice rury\n");
	scanf("%lf", &D);
	printf("Podaj srednice rdzenia\n");
	scanf("%lf", &d);
	printf("Podaj amplitude sily\n");
	scanf("%lf", &a);
	printf("Podaj okres ruchu\n");
	scanf("%lf", &T);
	printf("Podaj rodzaj cieczy\n");
	printf("\t 1 - woda\n");
	printf("\t 2 - rtec\n");
	printf("\t 3 - powietrze\n");
	printf("\t 4 - wodor\n");
	scanf("%d", &ciecz);
	if (ciecz == 1) {
		mi = 0.89 * 0.001;
		delta_t = pow(10, -5);
	
	}
	if (ciecz == 2) {
		mi = 1.554 * 0.001;
		delta_t = pow(10, -4);


	}
	if (ciecz == 3) {
		mi = 17.08 * 0.000001;
		delta_t = pow(10, -2);
	

	}
	if (ciecz == 4) {
		mi = 21.37 * 0.000001;
		delta_t = pow(10, -2);


	}
	if (ciecz != 1 && ciecz != 2 && ciecz != 3 && ciecz != 4) {
		printf("Blad nieznana ciecz\n");
		exit(-1);
	}
	printf("Podaj czas symulacji\n");
	scanf("%lf", &tk);

	//promien rdzenia wewnetrznego [m]
	double r = d / 2; 
	//promien rury [m]
	double R = D / 2; 
	przekroj_rury(r, R, &A_r);
	//liczba krokow czasowych [-]
	int nt;
	nt = (tk - t) / delta_t;
 	//liczba wezlow na y(r) - liczba wierszy macierzy wezlow [-]
	int n;
	//liczba wezlow z niewiadoma u [-]
	int n1;	
	n = (R - r) / h + 1;
	n1 = n - 2;
	printf("liczba wezlow, krokow czasowych : %d\t%d\n", n, nt);

	double* u, * u1, * u0, * u01, * r_j, * r_j1;
	double u_sr, Q;
	u = (double*)malloc(n * sizeof(double));			//wektor u dla wszytskich wezlow
	u1 = (double*)malloc(n1 * sizeof(double));			//wektor u dla wezlow bez scianek
	u0 = (double*)malloc(n * sizeof(double));			//wektor warunku poczatkowego
	u01 = (double*)malloc(n1 * sizeof(double));			//wektor warunku poczatkowego bez scianek
	r_j = (double*)malloc(n * sizeof(double));			//wektor wsp r dla wszystkich wezlow
	r_j1 = (double*)malloc(n1 * sizeof(double));		//wektor wsp r dla wezlow bez scianek

	double** A, * b;
	A = (double**)malloc(n1 * sizeof(double*));
	b = (double*)malloc(n1 * sizeof(double));
	for (int i = 0; i < n1; ++i) {
		A[i] = (double*)malloc(n1 * sizeof(double));
	}
	for (int i = 0; i < n; ++i) {
		u0[i] = 0;
		r_j[i] = r + i * h;
	}
	for (int i = 0; i < n1; ++i) {
		u01[i] = 0;
		r_j1[i] = r + (i + 1) * h;
	}
	//warunki brzegowe
	u[0] = 0;
	u[n - 1] = 0;

	//Solver
	macierz_ukladu(delta_t, mi, h, r_j1, n1, A);

	for (int j = 0; j < nt; ++j) {

		wektor_prawych_stron(a, T, delta_t, t, u01, n1, b);
		Thomas(A, u1, b, n1);
		rozwiazanie_calkowite(u1, n, u);
		srednia_predkosc(u, r_j, n, &u_sr, h);
		wydatek(u_sr, A_r, &Q);
		//zapis do pliku
		if (j % 100 == 0 && t>0.649) {
			for (int i = 0; i < n1; ++i) {
				fprintf(fb, "%lf\n", b[i]);
			}
			for (int i = 0; i < n; ++i) {
				fprintf(fu, "%lf\n", u[i]);
			}
			for (int i = 0; i < n; ++i) {
				fprintf(fr, "%lf\t%lf\n", r_j[i], t);
			}
			fprintf(fu_srQ, "%lf\t%lf\t%lf\n", t, u_sr, Q);
			printf("%d\n", j);
		}

		//nadpisanie wczesniej chwili
		for (int i = 0; i < n1; ++i) {
			u01[i] = u1[i];
		}
		t += delta_t;
	}
	printf("Done\n");
	for (int i = 0; i < n1; ++i) {
		for (int j = 0; j < n1; ++j) {
			fprintf(fA, "%lf\t", A[i][j]);
		}
		fprintf(fA, "\n");
	}

	for (int i = 0; i < n1; ++i) {
		free(A[i]);
	}
	free(u);
	free(u0);
	free(r_j);
	free(b);
	free(A);
	fclose(fA);
	fclose(fb);
	fclose(fu);
	fclose(fr);
	fclose(fu_srQ);
}


void macierz_ukladu(double delta_t, double mi, double h, double* r, int N, double** A) {
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j)
			A[i][j] = 0;
	}
	A[0][0] = 1. / delta_t + 2. * mi / (h * h);
	A[0][1] = -mi * (r[0] + h / 2.) / (r[0] * h * h);

	for (int i = 1; i < N - 1; ++i) {
		A[i][i - 1] = -mi * (r[i] - h / 2.) / (r[i] * h * h);
		A[i][i] = 1. / delta_t + 2. * mi / (h * h);
		A[i][i + 1] = -mi * (r[i] + h / 2.) / (r[i] * h * h);
	}
	A[N - 1][N - 1 - 1] = -mi * (r[N - 1] - h / 2.) / (r[N - 1] * h * h);
	A[N - 1][N - 1] = 1. / delta_t + 2. * mi / (h * h);

}

void wektor_prawych_stron(double a, double T, double delta_t, double t, double* u0, int N, double* b) {
	double w;
	w = 2. * 4. * atan(1.0) / T;

	for (int i = 0; i < N; ++i) {
		b[i] = -a * sin(w * t) + 1. / delta_t * u0[i];
	}

}

void rozwiazanie_calkowite(double* u1, int N, double* u) {

	for (int i = 1; i < N - 1; ++i)
		u[i] = u1[i - 1];

}

void srednia_predkosc(double* u, double* r, int n, double* u_sr, double h) {
	double u_calka_t;
	trapz(u, h, n, &u_calka_t);
	*u_sr = u_calka_t / (r[n - 1] - r[0]);
}

void przekroj_rury(double r, double R, double* A) {
	*A = 4 * atan(1.) * (R * R - r * r);
}

void wydatek(double u_sr, double A, double* Q) {
	*Q = A * u_sr;
}

void trapz(double* x, double h, int n, double *y) {
    *y = 0;
    for (int i = 0; i < n-1; ++i) {
        *y += (x[i] + x[i + 1]) * h / 2;
    }
}


void Thomas(double** A, double* x, double* b, int N) {
    double *alpha, *beta;
    alpha = (double*)malloc(N * sizeof(double));
    beta = (double*)malloc(N * sizeof(double));
    alpha[0] = -A[0][1] / A[0][0];
    beta[0] = b[0] /  A[0][0];
    for (int i = 0; i < N-2; ++i) {
        alpha[i+1] = -A[i + 1][i +1 + 1] / (A[i+1][i+1] + alpha[i]*A[i+1][i+1-1]);
        beta[i+1] = (b[i+1] -A[i+1][i+1-1] * beta[i]) / (A[i + 1][i + 1] + alpha[i] * A[i + 1][i + 1 - 1]);
    }
    x[N - 1] = (b[N - 1] - A[N - 1][N - 1 -1] * beta[ N - 2]) / (A[N - 1][N - 1]  + alpha[N-2]* A[N - 1][N - 1 - 1]);
    for (int i = N - 2; i >= 0; --i) {
        x[i] = alpha[i] * x[i + 1] + beta[i];
    }
    free(alpha);
    free(beta);
}
