#include <stdio.h>
#include <math.h>

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_sort_vector.h"
#include "gsl/gsl_sf.h"

#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

#include "R.h"
#include "Rmath.h"

#include "BAFT_LNsurv.h"


/* Updating y */

void BAFT_LNsurv_update_y(gsl_vector *yL,
                          gsl_vector *yU,
                          gsl_vector *yU_posinf,
                          gsl_vector *c0,
                          gsl_matrix *X,
                          gsl_vector *y,
                          gsl_vector *beta,
                          double beta0,
                          double sigSq)
{
    double eta, sample;
    int i;
    int n = y -> size;
    
    gsl_vector *xbeta = gsl_vector_calloc(n);
    gsl_blas_dgemv(CblasNoTrans, 1, X, beta, 0, xbeta);
    
    for(i=0;i<n;i++)
    {
        if(gsl_vector_get(yU, i) != gsl_vector_get(yL, i))
        {
            eta = beta0+gsl_vector_get(xbeta, i);
            c_rtnorm(eta, sqrt(sigSq), gsl_vector_get(yL, i), gsl_vector_get(yU, i), 0, gsl_vector_get(yU_posinf, i), &sample);
            gsl_vector_set(y, i, sample);
        }
        else if(gsl_vector_get(yU, i) == gsl_vector_get(yL, i))
        {
            gsl_vector_set(y, i, gsl_vector_get(yU, i));
        }
    }
    gsl_vector_free(xbeta);
    return;
}



/* Updating beta */

void BAFT_LNsurv_update_beta(gsl_vector *yL,
                             gsl_vector *yU,
                             gsl_vector *yU_posinf,
                             gsl_vector *c0,
                             gsl_vector *c0_neginf,
                             gsl_matrix *X,
                             gsl_vector *y,
                             gsl_vector *beta,
                             double beta0,
                             double sigSq,
                             double beta_prop_var,
                             gsl_vector *accept_beta)
{
    int i, j, u;
    double eta, eta_prop, loglh, loglh_prop, logR;
    
    int n = X -> size1;
    int p = X -> size2;
    
    gsl_vector *beta_prop = gsl_vector_calloc(p);
    gsl_vector *xbeta = gsl_vector_calloc(n);
    gsl_vector *xbeta_prop = gsl_vector_calloc(n);
    
    j = (int) runif(0, p);
    
    loglh = 0;
    loglh_prop = 0;
    
    gsl_vector_memcpy(beta_prop, beta);
    gsl_vector_set(beta_prop, j, rnorm(gsl_vector_get(beta, j), sqrt(beta_prop_var)));
    gsl_blas_dgemv(CblasNoTrans, 1, X, beta, 0, xbeta);
    gsl_blas_dgemv(CblasNoTrans, 1, X, beta_prop, 0, xbeta_prop);
    
    for(i=0;i<n;i++)
    {
        eta = beta0 + gsl_vector_get(xbeta, i);
        eta_prop = beta0 + gsl_vector_get(xbeta_prop, i);
        if(gsl_vector_get(c0_neginf, i) == 0)
        {
            loglh += dnorm(gsl_vector_get(y, i), eta, sqrt(sigSq), 1) - pnorm(gsl_vector_get(c0, i), eta, sqrt(sigSq), 0, 1);
            loglh_prop += dnorm(gsl_vector_get(y, i), eta_prop, sqrt(sigSq), 1) - pnorm(gsl_vector_get(c0, i), eta_prop, sqrt(sigSq), 0, 1);
        }else
        {
            loglh += dnorm(gsl_vector_get(y, i), eta, sqrt(sigSq), 1);
            loglh_prop += dnorm(gsl_vector_get(y, i), eta_prop, sqrt(sigSq), 1);
        }
    }
    
    logR = loglh_prop - loglh;
    u = log(runif(0, 1)) < logR;
    if(u == 1)
    {
        gsl_vector_memcpy(beta, beta_prop);
        gsl_vector_set(accept_beta, j, gsl_vector_get(accept_beta, j) + 1);
    }
    
    gsl_vector_free(beta_prop);
    gsl_vector_free(xbeta);
    gsl_vector_free(xbeta_prop);
    return;
}




/* Updating beta0 */

void BAFT_LNsurv_update_beta0(gsl_vector *yL,
                              gsl_vector *yU,
                              gsl_vector *yU_posinf,
                              gsl_vector *c0,
                              gsl_vector *c0_neginf,
                              gsl_matrix *X,
                              gsl_vector *y,
                              gsl_vector *beta,
                              double *beta0,
                              double sigSq,
                              double beta0_prop_var,
                              int *accept_beta0)
{
    int i, u;
    double eta, eta_prop, loglh, loglh_prop, logR, beta0_prop, logprior, logprior_prop;
    
    int n = X -> size1;
    
    gsl_vector *xbeta = gsl_vector_calloc(n);
    
    loglh = 0;
    loglh_prop = 0;
    beta0_prop = rnorm(*beta0, sqrt(beta0_prop_var));
    gsl_blas_dgemv(CblasNoTrans, 1, X, beta, 0, xbeta);
    
    for(i=0;i<n;i++)
    {
        eta = *beta0 + gsl_vector_get(xbeta, i);
        eta_prop = beta0_prop + gsl_vector_get(xbeta, i);
        if(gsl_vector_get(c0_neginf, i) == 0)
        {
            loglh += dnorm(gsl_vector_get(y, i), eta, sqrt(sigSq), 1) - pnorm(gsl_vector_get(c0, i), eta, sqrt(sigSq), 0, 1);
            loglh_prop += dnorm(gsl_vector_get(y, i), eta_prop, sqrt(sigSq), 1) - pnorm(gsl_vector_get(c0, i), eta_prop, sqrt(sigSq), 0, 1);
        }else
        {
            loglh += dnorm(gsl_vector_get(y, i), eta, sqrt(sigSq), 1);
            loglh_prop += dnorm(gsl_vector_get(y, i), eta_prop, sqrt(sigSq), 1);
        }        
    }
    
    logprior = dnorm(*beta0, 0, pow(10,6)*sqrt(sigSq), 1);
    logprior_prop = dnorm(beta0_prop, 0, pow(10,6)*sqrt(sigSq), 1);
    
    logR = loglh_prop - loglh;
    u = log(runif(0, 1)) < logR;
    if(u == 1)
    {
        *beta0 = beta0_prop;
        *accept_beta0 += 1;
    }
    
    gsl_vector_free(xbeta);
    return;
}












/* Updating sigmaSq */

void BAFT_LNsurv_update_sigSq(gsl_vector *yL,
                              gsl_vector *yU,
                              gsl_vector *yU_posinf,
                              gsl_vector *c0,
                              gsl_vector *c0_neginf,
                              gsl_matrix *X,
                              gsl_vector *y,
                              gsl_vector *beta,
                              double beta0,
                              double *sigSq,
                              double a_sigSq,
                              double b_sigSq,
                              double sigSq_prop_var,
                              int *accept_sigSq)
{
    int i, u;
    double eta, loglh, loglh_prop, logR, gamma_prop, sigSq_prop;
    double logprior, logprior_prop;
    
    int n = X -> size1;
    gsl_vector *xbeta = gsl_vector_calloc(n);
    
    loglh = 0;
    loglh_prop = 0;
    gamma_prop = rnorm(log(*sigSq), sqrt(sigSq_prop_var));
    sigSq_prop = exp(gamma_prop);
    gsl_blas_dgemv(CblasNoTrans, 1, X, beta, 0, xbeta);
    
    for(i=0;i<n;i++)
    {
        eta = beta0 + gsl_vector_get(xbeta, i);
        if(gsl_vector_get(c0_neginf, i) == 0)
        {
            loglh += dnorm(gsl_vector_get(y, i), eta, sqrt(*sigSq), 1) - pnorm(gsl_vector_get(c0, i), eta, sqrt(*sigSq), 0, 1);
            loglh_prop += dnorm(gsl_vector_get(y, i), eta, sqrt(sigSq_prop), 1) - pnorm(gsl_vector_get(c0, i), eta, sqrt(sigSq_prop), 0, 1);
        }else
        {
            loglh += dnorm(gsl_vector_get(y, i), eta, sqrt(*sigSq), 1);
            loglh_prop += dnorm(gsl_vector_get(y, i), eta, sqrt(sigSq_prop), 1);
        }        
    }
    
    logprior = (-a_sigSq-1)*log(*sigSq)-b_sigSq /(*sigSq);
    logprior_prop = (-a_sigSq-1)*log(sigSq_prop)-b_sigSq/sigSq_prop;
    
    logR = loglh_prop - loglh + logprior_prop - logprior + gamma_prop - log(*sigSq);
    
    u = log(runif(0, 1)) < logR;
    
    if(u == 1)
    {
        *sigSq = sigSq_prop;
        *accept_sigSq += 1;
    }
    
    gsl_vector_free(xbeta);
    return;
}






