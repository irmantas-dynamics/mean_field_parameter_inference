/*
 * losslib.c
 *
 * Computes an RMSE-based loss for fitting a 3-variable ODE system
 * (a mean-field / Montbrió-type neural-population model) to
 * experimental voltage traces.
 *
 * The ODE state is (r, v, A):
 *   r  – mean firing rate
 *   v  – mean membrane potential
 *   S  – synaptic actyvation 
 *
 * The experimental voltage trace V(t) enters as a driving force via
 * the coupling term K*(Vt − v), interpolated with an Akima spline.
 *
 * Public API
 * ----------
 * init_force_interp()  – build an Akima spline over the (time, voltage) data.
 * force_eval()         – evaluate the spline with clamped extrapolation.
 * compute_rmse()       – integrate the ODE for a single parameter set and
 *                        return 0.5 * MSE(v_sim, V_exp) over a user-specified
 *                        window after an initial transient period.
 * lossPar()            – builds the spline once, then sweeps over n_params
 *                        parameter tuples in parallel (OpenMP), filling
 *                        rmse_out[].
 *
 * Dependencies: GSL (spline), LSODA (ODE integrator), OpenMP.
 */
#include <stdio.h>
#include <stdlib.h>
#include "common.h"
#include "lsoda.h"
#include "math.h"
#include <omp.h>

#include <gsl/gsl_spline.h>

#ifndef M_PI
    #define M_PI 3.1415926535897932384626433
#endif

#define VERB 0


static gsl_spline *force_spline = NULL;
//static gsl_interp_accel *force_acc = NULL;

typedef struct {
    double taum;
    double etab;
    double J;
    double delta;
    double taua;
    double K;
    const double *Texp;
    const double *Vexp;
    int n;
} ode_data_t;


void init_force_interp(const double *t,
                       const double *f,
                       size_t n)
{
   // fprintf(stderr,"force init\n");

    //force_acc = gsl_interp_accel_alloc();
    force_spline = gsl_spline_alloc(gsl_interp_akima, n);

    gsl_spline_init(force_spline, t, f, n);
}

double force_eval(double t)
{
    //fprintf(stderr,"force eval\n");
    double t_min = force_spline->x[0];
    double t_max = force_spline->x[force_spline->size - 1];

    if (t <= t_min) return gsl_spline_eval(force_spline, t_min, NULL);
    if (t >= t_max) return gsl_spline_eval(force_spline, t_max, NULL);

   // return gsl_spline_eval(force_spline, t, force_acc);
    return gsl_spline_eval(force_spline, t, NULL);
}


int fex(double t, double *y, double *ydot, void *data)
{
    ode_data_t *p = (ode_data_t *)data;

    double r = y[0];
    double v = y[1];
    double S = y[2];
    //double Iext =  (1.0 + sin(p->Kom * t) / 2.0);
    //Iext = p->KAmp*Iext*Iext*Iext;

   //double Vt = interp_linear(t, p->Texp, p->Vexp, p->n);
    double Vt = force_eval(t);

    ydot[0] = (p->delta / ( M_PI * p->taum) + 2.0 * r * v) / p->taum;
    ydot[1] = (v * v - (p->taum * M_PI * r) * (p->taum * M_PI * r) + p->etab - p->J  * p->taum * S + p->K*(Vt-v) ) / p->taum;
    ydot[2] = (-S + r ) / p->taua;

    return 0;
}


double compute_rmse(double taum,
                    double etab,
                    double J,
                    double delta,
                    double taua,
                    double K,
                    const double *Texp,
                    const double *Vexp,
                    int n,
                    int nTtr,
                    int nDt)
{
    double y[3];
    double t, tout;
    double rtol[3], atol[3];
    double Rtol = 1e-8, Atol = 1e-10;

    for (int i1 = 0;i1<3;i1++){
        rtol[i1] = Rtol;
        atol[i1] = Atol;
    }

    /* initial condition: match first observation */
    y[0] = 0.1;
    y[1] = Vexp[0];
    y[2] = 0.1;

    t = Texp[0];

    //init_force_interp(Texp, Vexp, n);


    ode_data_t data = {
        .taum = taum,
        .etab = etab,
        .J  = J,
        .delta = delta,
        .taua = taua,
        .K = K,
        .Texp = Texp,
        .Vexp = Vexp,
        .n = n
    };

    struct lsoda_opt_t opt = {0};
    opt.rtol = rtol;
    opt.atol = atol;
    opt.itask = 1;

    struct lsoda_context_t ctx = {
        .function = fex,
        .neq = 3,
        .data = &data,
        .state = 1,
    };

    lsoda_prepare(&ctx, &opt);

    #if VERB==1
        /* allocate arrays to store results */
        double *Tsim = malloc(n * sizeof(double));
        double *Rsim = malloc(n * sizeof(double));
        double *Vsim = malloc(n * sizeof(double));
        double *Asim = malloc(n * sizeof(double));

        if (!Tsim || !Rsim || !Vsim || !Asim) {
            perror("malloc");
            lsoda_free(&ctx);
            return 1e30;
        }

        Tsim[0] = t;
        Rsim[0] = y[0];
        Vsim[0] = y[1];
        Asim[0] = y[2];
    #endif




    // transitional period
    //int nTtr = (int)2000./(Texp[1]-Texp[0]) ;
    // period for calculating loss function
    //int nDt = (int)500./(Texp[1]-Texp[0]) ;

    //fprintf(stderr,"nTtr = %d\n",nTtr);
    //fprintf(stderr,"nDt = %d\n" ,nDt);

    // initial integration
    for (int i = 1; i < nTtr; i++) {
        tout = Texp[i];
        lsoda(&ctx, y, &t, tout);

        if (ctx.state <= 0) {
            lsoda_free(&ctx);
            return 1e30;  /* hard failure penalty */
        }
        #if VERB==1
            Tsim[i] = t;
            Rsim[i] = y[0];
            Vsim[i] = y[1];
            Asim[i] = y[2];
         #endif

    }

    /* RMSE accumulation */
    double err2 = 0.0;
    int count = 0;
    // error estimation
    for (int i = nTtr; i < nTtr+nDt; i++) {
        tout = Texp[i];
        lsoda(&ctx, y, &t, tout);

        if (ctx.state <= 0) {
            lsoda_free(&ctx);
            return 1e30;  /* hard failure penalty */
        }

        #if VERB==1
            Tsim[i] = t;
            Rsim[i] = y[0];
            Vsim[i] = y[1];
            Asim[i] = y[2];
        #endif

        double diff = y[1] - Vexp[i];
        err2 += diff * diff;
        count++;
    }

    lsoda_free(&ctx);
    
    #if VERB==1
          /* write results to file after integration */
        FILE *fout = fopen("results.txt", "w");
        if (!fout) {
            perror("fopen results.txt");
            free(Tsim); free(Rsim); free(Vsim); free(Asim);
            return 1e30;
        }

        fprintf(fout, "# t  r  v  A\n");
        for (int i = 0; i < n; i++) {
            fprintf(fout, "%.10e %.10e %.10e %.10e\n",
                    Tsim[i], Rsim[i], Vsim[i], Asim[i]);
        }

        fclose(fout);


        free(Tsim);
        free(Rsim);
        free(Vsim);
        free(Asim);
        printf("%e\n",err2/count);
    #endif
    //gsl_spline_free(force_spline);
    
    return 0.5*(err2 / count);
}

void lossPar(             const double *taum_vals,
                          const double *etab_vals,
                          const double *J_vals,
                          const double *delta_vals,
                          const double *taua_vals,
                          int n_params,
                          double K,
                          const double *Texp,
                          const double *Vexp,
                          int n,    // length of Texp, Vexp
                          int nTtr, // steps of transitional integration (No of integration steps)
                          int nDt,  // steps for Loss function calculation
                          double *rmse_out)
{
    init_force_interp(Texp, Vexp, n);
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n_params; i++) {
            rmse_out[i] = compute_rmse(
                taum_vals[i],
                etab_vals[i],
                J_vals[i],
                delta_vals[i],
                taua_vals[i],
                K,
                Texp,
                Vexp,
                n,
                nTtr,
                nDt
            );

    }
    gsl_spline_free(force_spline);
    //gsl_interp_accel_free(force_acc);
}


