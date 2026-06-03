/*
 * qif_in_micro_ext.c
 *
 * C + OpenMP simulation of QIF neurons with inhibitory interaction,
 * a shared adaptation variable st, and periodic external forcing.
 *
 *
 * Microscopic equations:
 *   dy(1:N)  = ( 1 - cs + (1+cs).*(etam - J*taum*st + Iext) ) / taum
 *   dy(N+1)  = (-st + R) / taua
 *
 * External forcing:
 *   Iext = K_Amp * (1 + sin(K_om*t)/2)^3
 *
 * eta is loaded from file (same as original C code):
 *   ../time_series_rand/eta_N<N>_<run>.txt
 *
 * Compile:
 *   gcc -O3 -march=native -fopenmp -o qif_in_ext qif_in_micro_ext.c -lm
 *
 * Usage:
 *   ./qif_in_ext <N> <run_index>
 *   e.g.: ./qif_in_ext 100000 1
 *
 * Output:
 *   time_series/QIF_IN_N<N>_<run_index>_ext.txt
 *   columns: t  r  v  st
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <omp.h>

/* ------------------------------------------------------------------ */
/*  Physical / simulation parameters                                    */
/* ------------------------------------------------------------------ */
#define J      21.0
#define TAUA   5.0
#define TAUM   10.0

/* External forcing parameters (matching MATLAB K struct) */
#define K_AMP  -0.45
#define K_TPER  28.0
#define K_OM   (2.0 * M_PI / K_TPER)

/* Integration parameters (matching MATLAB) */
#define DT     (1e-3 * TAUM)           /* = 0.01 */
#define TEND   (100.0 * K_TPER)        /* = 2800.0 */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ------------------------------------------------------------------ */
/*  Load eta from file                                                  */
/* ------------------------------------------------------------------ */
static int load_eta(const char *path, double *etam, int N)
{
    FILE *fp = fopen(path, "r");
    if (!fp) { perror(path); return -1; }

    char line[128];
    int count = 0;
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#') continue;           /* skip comment lines */
        if (count >= N) break;
        etam[count++] = atof(line);
    }
    fclose(fp);

    if (count != N) {
        fprintf(stderr, "load_eta: expected %d values, got %d\n", N, count);
        return -1;
    }
    return 0;
}

/* ------------------------------------------------------------------ */
/*  W = (1 - conj(Z)) / (1 + conj(Z))                                  */
/* ------------------------------------------------------------------ */
static inline void Z_to_W(double Zre, double Zim, double *Wre, double *Wim)
{
    double num_re =  1.0 - Zre,  num_im = Zim;
    double den_re =  1.0 + Zre,  den_im = -Zim;
    double den2   = den_re*den_re + den_im*den_im;
    *Wre = (num_re*den_re + num_im*den_im) / den2;
    *Wim = (num_im*den_re - num_re*den_im) / den2;
}

/* ------------------------------------------------------------------ */
/*  Kuramoto order parameter  Z = mean(exp(i*theta))                   */
/* ------------------------------------------------------------------ */
static void kuramoto(const double * restrict theta, int N,
                     double *Zre, double *Zim)
{
    double sre = 0.0, sim = 0.0;
    #pragma omp parallel for reduction(+:sre,sim) schedule(static)
    for (int i = 0; i < N; ++i) {
        sre += cos(theta[i]);
        sim += sin(theta[i]);
    }
    *Zre = sre / N;
    *Zim = sim / N;
}

/* ------------------------------------------------------------------ */
/*  External forcing                                                     */
/*   Iext = K_Amp * (1 + sin(K_om*t)/2)^3                             */
/* ------------------------------------------------------------------ */
static inline double ext_force(double t)
{
    double s = 1.0 + 0.5 * sin(K_OM * t);
    return K_AMP * s * s * s;
}

/* ------------------------------------------------------------------ */
/*  RHS -- inhibitory network + external force, shared adaptation st   */
/*                                                                      */
/*  dy(1:N)  = ( 1-cs + (1+cs).*(etam - J*taum*st + Iext) ) / taum   */
/*  dy(N+1)  = (-st + R) / taua                                        */
/* ------------------------------------------------------------------ */
static void rhs(const double * restrict theta,
                double st,
                const double * restrict etam,
                int N, double R, double t,
                double * restrict dtheta,
                double *dst)
{
    double JtaumSt = J * TAUM * st;
    double Iext    = ext_force(t);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; ++i) {
        double cs = cos(theta[i]);
        double Fn = etam[i] - JtaumSt + Iext;
        dtheta[i] = (1.0 - cs + (1.0 + cs) * Fn) / TAUM;
    }

    *dst = (-st + R) / TAUA;
}

/* ------------------------------------------------------------------ */
/*  RK4 step -- theta[N] + scalar st, time-aware for Iext             */
/* ------------------------------------------------------------------ */
static void rk4_step(double * restrict theta,
                     double *st,
                     const double * restrict etam,
                     int N, double h, double t,
                     double *k1t, double *k2t,
                     double *k3t, double *k4t,
                     double *tmp_t)
{
    double Zre, Zim, Wre, Wim, R;
    double k1s, k2s, k3s, k4s, tmp_s;
    const double h2     = 0.5 * h;
    const double inv6h  = h / 6.0;
    const double two_pi = 2.0 * M_PI;

    /* k1 */
    kuramoto(theta, N, &Zre, &Zim);
    Z_to_W(Zre, Zim, &Wre, &Wim);
    R = Wre / (M_PI * TAUM);
    rhs(theta, *st, etam, N, R, t, k1t, &k1s);

    /* k2 */
    for (int i = 0; i < N; ++i) tmp_t[i] = theta[i] + h2 * k1t[i];
    tmp_s = *st + h2 * k1s;
    kuramoto(tmp_t, N, &Zre, &Zim);
    Z_to_W(Zre, Zim, &Wre, &Wim);
    R = Wre / (M_PI * TAUM);
    rhs(tmp_t, tmp_s, etam, N, R, t + h2, k2t, &k2s);

    /* k3 */
    for (int i = 0; i < N; ++i) tmp_t[i] = theta[i] + h2 * k2t[i];
    tmp_s = *st + h2 * k2s;
    kuramoto(tmp_t, N, &Zre, &Zim);
    Z_to_W(Zre, Zim, &Wre, &Wim);
    R = Wre / (M_PI * TAUM);
    rhs(tmp_t, tmp_s, etam, N, R, t + h2, k3t, &k3s);

    /* k4 */
    for (int i = 0; i < N; ++i) tmp_t[i] = theta[i] + h * k3t[i];
    tmp_s = *st + h * k3s;
    kuramoto(tmp_t, N, &Zre, &Zim);
    Z_to_W(Zre, Zim, &Wre, &Wim);
    R = Wre / (M_PI * TAUM);
    rhs(tmp_t, tmp_s, etam, N, R, t + h, k4t, &k4s);

    /* combine + wrap */
    for (int i = 0; i < N; ++i) {
        theta[i] += inv6h * (k1t[i] + 2.0*k2t[i] + 2.0*k3t[i] + k4t[i]);
        theta[i] -= two_pi * floor(theta[i] / two_pi);
    }
    *st += inv6h * (k1s + 2.0*k2s + 2.0*k3s + k4s);
}

/* ------------------------------------------------------------------ */
/*  Output helper                                                       */
/* ------------------------------------------------------------------ */
static int save_timeseries(const char *path,
                           const double *tm, const double *r,
                           const double *v,  const double *st_arr,
                           int Nt, double dt, int N)
{
    FILE *fp = fopen(path, "w");
    if (!fp) { perror(path); return -1; }
    fprintf(fp, "# QIF inhibitory microscopic simulation with external forcing\n");
    fprintf(fp, "# N=%d  dt=%.4e  Nt=%d\n", N, dt, Nt);
    fprintf(fp, "# J=%.6f  taua=%.6f  taum=%.6f\n", J, TAUA, TAUM);
    fprintf(fp, "# K_Amp=%.6f  K_Tper=%.6f  K_om=%.10f\n", K_AMP, K_TPER, K_OM);
    fprintf(fp, "# columns: t  r  v  st\n");
    for (int i = 0; i < Nt; ++i)
        fprintf(fp, "%.6f %.10f %.10f %.10f\n", tm[i], r[i], v[i], st_arr[i]);
    fclose(fp);
    return 0;
}

/* ------------------------------------------------------------------ */
/*  main                                                                */
/* ------------------------------------------------------------------ */
int main(int argc, char **argv)
{
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <N> <run_index>\n", argv[0]);
        return 1;
    }
    int N   = atoi(argv[1]);
    int fnm = atoi(argv[2]);
    if (N <= 0) { fprintf(stderr, "N must be positive.\n"); return 1; }

    omp_set_num_threads(12);

    mkdir("time_series", 0755);
    printf("QIF-IN-EXT  N=%d  run=%d  threads=%d\n",
           N, fnm, omp_get_max_threads());
    printf("Tend=%.1f  dt=%.4e  K_Amp=%.4f  K_Tper=%.1f\n",
           TEND, DT, K_AMP, K_TPER);

    /* --- load neuron heterogeneity from file (same as original C) --- */
    double *etam = (double*)malloc(N * sizeof(double));
    if (!etam) { fprintf(stderr, "malloc failed\n"); return 1; }

    {
        char eta_path[256];
        snprintf(eta_path, sizeof(eta_path),
                 "../QIF_IN_Noninvasive/time_series_rand/eta_N%d_%d.txt", N, fnm);
        if (load_eta(eta_path, etam, N) != 0) {
            fprintf(stderr, "Failed to load eta from %s\n", eta_path);
            free(etam);
            return 1;
        }
        printf("Loaded eta from %s\n", eta_path);
    }

    /* --- allocate state --- */
    double *theta = (double*)malloc(N * sizeof(double));
    if (!theta) { fprintf(stderr, "malloc failed\n"); return 1; }

    /* --- initial conditions (matching MATLAB: V0=0, r0=0.1) --- */
    const double V0 = 0.0, r0 = 0.1;
    for (int j = 0; j < N; ++j) {
        double v0 = V0 + r0 * M_PI * tan(M_PI / 2.0 * (2.0*(j+1) - N - 1.0) / (N + 1.0));
        theta[j]  = 2.0 * atan(v0);
    }
    double st = 0.005;    /* shared adaptation scalar */

    /* --- RK4 scratch buffers --- */
    double *k1t  = (double*)malloc(N * sizeof(double));
    double *k2t  = (double*)malloc(N * sizeof(double));
    double *k3t  = (double*)malloc(N * sizeof(double));
    double *k4t  = (double*)malloc(N * sizeof(double));
    double *tmpt = (double*)malloc(N * sizeof(double));
    if (!k1t||!k2t||!k3t||!k4t||!tmpt) {
        fprintf(stderr, "malloc failed (scratch)\n"); return 1;
    }

    int Nt = (int)round(TEND / DT) + 1;

    /* --- transient (same duration as recording, matching MATLAB) --- */
    printf("Transient (%d steps)...\n", Nt);
    double wall0 = omp_get_wtime();
    double t = 0.0;
    for (int i = 0; i < Nt; ++i) {
        if (i % (Nt/10) == 0) printf("  %d%%\n", (int)(100.0*i/(Nt-1)));
        rk4_step(theta, &st, etam, N, DT, t, k1t, k2t, k3t, k4t, tmpt);
        t += DT;
    }
    printf("  done in %.2f s\n", omp_get_wtime() - wall0);

    /* --- recording arrays --- */
    double *tm_arr  = (double*)malloc(Nt * sizeof(double));
    double *Zre_arr = (double*)malloc(Nt * sizeof(double));
    double *Zim_arr = (double*)malloc(Nt * sizeof(double));
    double *st_arr  = (double*)malloc(Nt * sizeof(double));
    if (!tm_arr||!Zre_arr||!Zim_arr||!st_arr) {
        fprintf(stderr, "malloc failed (output)\n"); return 1;
    }

    /* --- signal generation --- */
    printf("Signal generation (%d steps)...\n", Nt);
    double wall1 = omp_get_wtime();
    t = 0.0;
    for (int i = 0; i < Nt; ++i) {
        if (i % (Nt/10) == 0) printf("  %d%%\n", (int)(100.0*i/(Nt-1)));

        tm_arr[i] = t;
        kuramoto(theta, N, &Zre_arr[i], &Zim_arr[i]);
        st_arr[i] = st;

        rk4_step(theta, &st, etam, N, DT, t, k1t, k2t, k3t, k4t, tmpt);
        t += DT;
    }
    printf("  done in %.2f s\n", omp_get_wtime() - wall1);

    /* --- derive r, v from Z --- */
    double *r_arr = (double*)malloc(Nt * sizeof(double));
    double *v_arr = (double*)malloc(Nt * sizeof(double));
    if (!r_arr || !v_arr) { fprintf(stderr, "malloc failed\n"); return 1; }

    for (int i = 0; i < Nt; ++i) {
        double Wre, Wim;
        Z_to_W(Zre_arr[i], Zim_arr[i], &Wre, &Wim);
        r_arr[i] = Wre / (M_PI * TAUM);
        v_arr[i] = Wim;
    }

    /* --- save --- */
    {
        char pav[256];
        snprintf(pav, sizeof(pav),
                 "time_series_rand/QIF_IN_N%d_%d_ext.txt", N, fnm);
        if (save_timeseries(pav, tm_arr, r_arr, v_arr, st_arr, Nt, DT, N) == 0)
            printf("Saved: %s\n", pav);
    }

    /* --- cleanup --- */
    free(etam); free(theta);
    free(k1t); free(k2t); free(k3t); free(k4t); free(tmpt);
    free(tm_arr); free(Zre_arr); free(Zim_arr); free(st_arr);
    free(r_arr); free(v_arr);

    printf("Done.\n");
    return 0;
}
