
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_sort_vector.h"
#include "gsl/gsl_sf.h"

#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

#include "R.h"
#include "Rmath.h"
#include "BAFT_DPscr.h"






/* Updating beta1 */

void BAFT_DPscr_update_beta1(gsl_vector *y1_NA,
                             gsl_vector *c0,
                             gsl_vector *c0_neginf,
                             gsl_matrix *X1,
                             gsl_vector *y1,
                             gsl_vector *y2,
                             gsl_vector *beta1,
                             gsl_vector *gamma,
                             gsl_vector *r1,
                             gsl_vector *mu1_all,
                             gsl_vector *zeta1_all,
                             gsl_vector *r1Uniq,
                             gsl_vector *r1Uniq_count,
                             int *nClass_DP1,
                             double beta1_prop_var,
                             gsl_vector *accept_beta1)
{
    int i, j, uu, k, ii;
    double eta, eta_prop, loglh, loglh_prop, logR;
    double mean, mean_prop, sd;
    
    int n = X1 -> size1;
    int p = X1 -> size2;
    
    int u = *nClass_DP1;
    
    gsl_vector *beta_prop = gsl_vector_calloc(p);
    gsl_vector *xbeta1 = gsl_vector_calloc(n);
    gsl_vector *xbeta1_prop = gsl_vector_calloc(n);
    
    j = (int) runif(0, p);
    
    loglh = 0;
    loglh_prop = 0;
    uu = 0;
    logR = 0;
    
    gsl_vector_memcpy(beta_prop, beta1);
    gsl_vector_set(beta_prop, j, rnorm(gsl_vector_get(beta1, j), sqrt(beta1_prop_var)));
    gsl_blas_dgemv(CblasNoTrans, 1, X1, beta1, 0, xbeta1);
    gsl_blas_dgemv(CblasNoTrans, 1, X1, beta_prop, 0, xbeta1_prop);
    
    for(i=0;i<n;i++)
    {
        eta = gsl_vector_get(xbeta1, i) + gsl_vector_get(gamma, i);
        eta_prop = gsl_vector_get(xbeta1_prop, i) + gsl_vector_get(gamma, i);
        
        for(k = 0; k < u; k++)
        {
            if(gsl_vector_get(r1, i) == gsl_vector_get(r1Uniq, k)) ii = k;
        }
        
        mean = eta + gsl_vector_get(mu1_all, ii);
        mean_prop = eta_prop + gsl_vector_get(mu1_all, ii);
        sd = (double) pow(gsl_vector_get(zeta1_all, ii), -0.5);
        
        if(gsl_vector_get(y1_NA, i) == 0)
        {
            loglh += dnorm(gsl_vector_get(y1, i), mean, sd, 1);
            loglh_prop += dnorm(gsl_vector_get(y1, i), mean_prop, sd, 1);
        }else
        {
            loglh += pnorm(gsl_vector_get(y2, i), mean, sd, 0, 1);
            loglh_prop += pnorm(gsl_vector_get(y2, i), mean_prop, sd, 0, 1);
        }
        
        if(gsl_vector_get(c0_neginf, i) == 0)
        {
            loglh -= pnorm(gsl_vector_get(c0, i), mean, sd, 0, 1);
            loglh_prop -= pnorm(gsl_vector_get(c0, i), mean_prop, sd, 0, 1);
        }
    }
    
    if(loglh_prop != -INFINITY && loglh_prop != INFINITY)
    {
        logR = loglh_prop - loglh;
        uu = log(runif(0, 1)) < logR;
        
        if(uu == 1)
        {
            gsl_vector_memcpy(beta1, beta_prop);
            gsl_vector_set(accept_beta1, j, gsl_vector_get(accept_beta1, j) + 1);
        }
    }
    
    gsl_vector_free(beta_prop);
    gsl_vector_free(xbeta1);
    gsl_vector_free(xbeta1_prop);
    return;
    
}







/* Updating beta2 */

void BAFT_DPscr_update_beta2(gsl_vector *y1_NA,
                             gsl_vector *c0,
                             gsl_vector *c0_neginf,
                             gsl_matrix *X2,
                             gsl_vector *y1,
                             gsl_vector *y2,
                             gsl_vector *beta2,
                             gsl_vector *gamma,
                             gsl_vector *r2,
                             gsl_vector *mu2_all,
                             gsl_vector *zeta2_all,
                             gsl_vector *r2Uniq,
                             gsl_vector *r2Uniq_count,
                             int *nClass_DP2,
                             double beta2_prop_var,
                             gsl_vector *accept_beta2)
{
    int i, j, uu, k, ii;
    double eta, eta_prop, loglh, loglh_prop, logR;
    double mean, mean_prop, sd;
    
    int n = X2 -> size1;
    int p = X2 -> size2;
    
    int u = *nClass_DP2;
    
    gsl_vector *beta_prop = gsl_vector_calloc(p);
    gsl_vector *xbeta2 = gsl_vector_calloc(n);
    gsl_vector *xbeta2_prop = gsl_vector_calloc(n);
    
    j = (int) runif(0, p);
    
    loglh = 0;
    loglh_prop = 0;
    uu = 0;
    logR = 0;
    
    gsl_vector_memcpy(beta_prop, beta2);
    gsl_vector_set(beta_prop, j, rnorm(gsl_vector_get(beta2, j), sqrt(beta2_prop_var)));
    gsl_blas_dgemv(CblasNoTrans, 1, X2, beta2, 0, xbeta2);
    gsl_blas_dgemv(CblasNoTrans, 1, X2, beta_prop, 0, xbeta2_prop);
    
    for(i=0;i<n;i++)
    {
        eta = gsl_vector_get(xbeta2, i) + gsl_vector_get(gamma, i);
        eta_prop = gsl_vector_get(xbeta2_prop, i) + gsl_vector_get(gamma, i);
        
        for(k = 0; k < u; k++)
        {
            if(gsl_vector_get(r2, i) == gsl_vector_get(r2Uniq, k)) ii = k;
        }
        
        mean = eta + gsl_vector_get(mu2_all, ii);
        mean_prop = eta_prop + gsl_vector_get(mu2_all, ii);
        sd = (double) pow(gsl_vector_get(zeta2_all, ii), -0.5);
        
        if(gsl_vector_get(y1_NA, i) == 0)
        {
            loglh += pnorm(gsl_vector_get(y1, i), mean, sd, 0, 1);
            loglh_prop += pnorm(gsl_vector_get(y1, i), mean_prop, sd, 0, 1);
        }else
        {
            loglh += dnorm(gsl_vector_get(y2, i), mean, sd, 1);
            loglh_prop += dnorm(gsl_vector_get(y2, i), mean_prop, sd, 1);
        }
        
        if(gsl_vector_get(c0_neginf, i) == 0)
        {
            loglh -= pnorm(gsl_vector_get(c0, i), mean, sd, 0, 1);
            loglh_prop -= pnorm(gsl_vector_get(c0, i), mean_prop, sd, 0, 1);
        }
    }
    
    if(loglh_prop != -INFINITY && loglh_prop != INFINITY)
    {
        logR = loglh_prop - loglh;
        uu = log(runif(0, 1)) < logR;
        if(uu == 1)
        {
            gsl_vector_memcpy(beta2, beta_prop);
            gsl_vector_set(accept_beta2, j, gsl_vector_get(accept_beta2, j) + 1);
        }
    }
    
    
    
    gsl_vector_free(beta_prop);
    gsl_vector_free(xbeta2);
    gsl_vector_free(xbeta2_prop);
    
    return;
    
}






/* Updating beta3 */

void BAFT_DPscr_update_beta3(gsl_vector *y1_NA,
                             gsl_vector *c0,
                             gsl_vector *c0_neginf,
                             gsl_matrix *X3,
                             gsl_vector *y1,
                             gsl_vector *y2,
                             gsl_vector *beta3,
                             gsl_vector *gamma,
                             gsl_vector *r3,
                             gsl_vector *mu3_all,
                             gsl_vector *zeta3_all,
                             gsl_vector *r3Uniq,
                             gsl_vector *r3Uniq_count,
                             int *nClass_DP3,
                             double beta3_prop_var,
                             gsl_vector *accept_beta3)
{
    int i, j, uu, k, ii;
    double eta, eta_prop, loglh, loglh_prop, logR;
    double mean, mean_prop, sd, yStar;
    
    int n = X3 -> size1;
    int p = X3 -> size2;
    
    int u = *nClass_DP3;
    
    gsl_vector *beta_prop = gsl_vector_calloc(p);
    gsl_vector *xbeta3 = gsl_vector_calloc(n);
    gsl_vector *xbeta3_prop = gsl_vector_calloc(n);
    
    j = (int) runif(0, p);
    
    loglh = 0;
    loglh_prop = 0;
    uu = 0;
    logR = 0;
    
    gsl_vector_memcpy(beta_prop, beta3);
    gsl_vector_set(beta_prop, j, rnorm(gsl_vector_get(beta3, j), sqrt(beta3_prop_var)));
    gsl_blas_dgemv(CblasNoTrans, 1, X3, beta3, 0, xbeta3);
    gsl_blas_dgemv(CblasNoTrans, 1, X3, beta_prop, 0, xbeta3_prop);
    
    for(i=0;i<n;i++)
    {
        if(gsl_vector_get(y1_NA, i) == 0)
        {
            if(gsl_vector_get(y1, i) < gsl_vector_get(y2, i))
            {
                eta = gsl_vector_get(xbeta3, i) + gsl_vector_get(gamma, i);
                eta_prop = gsl_vector_get(xbeta3_prop, i) + gsl_vector_get(gamma, i);
                
                for(k = 0; k < u; k++)
                {
                    if(gsl_vector_get(r3, i) == gsl_vector_get(r3Uniq, k))
                    {
                        ii = (int) k;
                    }
                }
                
                mean = eta + gsl_vector_get(mu3_all, ii);
                mean_prop = eta_prop + gsl_vector_get(mu3_all, ii);
                sd = (double) pow(gsl_vector_get(zeta3_all, ii), -0.5);
                
                yStar = gsl_vector_get(y2, i) + log(1-exp(gsl_vector_get(y1, i)-gsl_vector_get(y2, i)));
                
                loglh += dnorm(yStar, mean, sd, 1);
                loglh_prop += dnorm(yStar, mean_prop, sd, 1);
                
            }
        }
    }
    
    if(loglh_prop != -INFINITY && loglh_prop != INFINITY)
    {
        logR = loglh_prop - loglh;
        uu = log(runif(0, 1)) < logR;
        if(uu == 1)
        {
            gsl_vector_memcpy(beta3, beta_prop);
            gsl_vector_set(accept_beta3, j, gsl_vector_get(accept_beta3, j) + 1);
        }
    }
    
    gsl_vector_free(beta_prop);
    gsl_vector_free(xbeta3);
    gsl_vector_free(xbeta3_prop);
    
    return;
    
}








/* Updating y1 and y2 */

void BAFT_DPscr_update_y(gsl_vector *y1L,
                         gsl_vector *y1U,
                         gsl_vector *y2L,
                         gsl_vector *y2U,
                         gsl_vector *y1L_neginf,
                         gsl_vector *y2L_neginf,
                         gsl_vector *y1U_posinf,
                         gsl_vector *y2U_posinf,
                         gsl_vector *c0_neginf,
                         gsl_vector *y1_NA,
                         gsl_matrix *X1,
                         gsl_matrix *X2,
                         gsl_matrix *X3,
                         gsl_vector *y1,
                         gsl_vector *y2,
                         gsl_vector *beta1,
                         gsl_vector *beta2,
                         gsl_vector *beta3,
                         gsl_vector *r1,
                         gsl_vector *r2,
                         gsl_vector *r3,
                         gsl_vector *mu1_all,
                         gsl_vector *mu2_all,
                         gsl_vector *mu3_all,
                         gsl_vector *mu3_vec,
                         gsl_vector *zeta1_all,
                         gsl_vector *zeta2_all,
                         gsl_vector *zeta3_all,
                         gsl_vector *zeta3_vec,
                         gsl_vector *r1Uniq,
                         gsl_vector *r2Uniq,
                         gsl_vector *r3Uniq,
                         gsl_vector *r1Uniq_count,
                         gsl_vector *r2Uniq_count,
                         gsl_vector *r3Uniq_count,
                         gsl_vector *gamma,
                         int *nClass_DP1,
                         int *nClass_DP2,
                         int *nClass_DP3,
                         double mu03,
                         double sigSq03,
                         double a03,
                         double b03,
                         double tau3,
                         gsl_rng *rr)
{
    double eta1, eta2, eta3, y1_temp, y2_temp, gap, CenNew, CenNew1, CenNew2;
    int i, selInd;
    int n = y1 -> size;
    
    int u1 = *nClass_DP1;
    int u2 = *nClass_DP2;
    int u3 = *nClass_DP3;
    
    gsl_vector *xbeta1 = gsl_vector_calloc(n);
    gsl_vector *xbeta2 = gsl_vector_calloc(n);
    gsl_vector *xbeta3 = gsl_vector_calloc(n);
    gsl_blas_dgemv(CblasNoTrans, 1, X1, beta1, 0, xbeta1);
    gsl_blas_dgemv(CblasNoTrans, 1, X2, beta2, 0, xbeta2);
    gsl_blas_dgemv(CblasNoTrans, 1, X3, beta3, 0, xbeta3);
    
    for(i=0;i<n;i++)
    {
        gap = 0;
        
        eta1 = gsl_vector_get(xbeta1, i)+gsl_vector_get(gamma, i);
        eta2 = gsl_vector_get(xbeta2, i)+gsl_vector_get(gamma, i);
        eta3 = gsl_vector_get(xbeta3, i)+gsl_vector_get(gamma, i);
        
        /*  both end points are right-censored */
        
        if(gsl_vector_get(y1U_posinf, i) == 1 && gsl_vector_get(y2U_posinf, i) == 1)
        {
            selInd = MFunc_BAFT_DP(gsl_vector_get(y1L, i), 100000, gsl_vector_get(y1L_neginf, i), 1, mu1_all, zeta1_all, r1Uniq_count, u1, eta1, rr);
            
            c_rtnorm(eta1 + gsl_vector_get(mu1_all, selInd-1), sqrt(1/gsl_vector_get(zeta1_all, selInd-1)), gsl_vector_get(y1L, i), 100000, gsl_vector_get(y1L_neginf, i), 1, &y1_temp);
            
            selInd = MFunc_BAFT_DP(gsl_vector_get(y2L, i), 100000, gsl_vector_get(y2L_neginf, i), 1, mu2_all, zeta2_all, r2Uniq_count, u2, eta2, rr);
            
            c_rtnorm(eta2 + gsl_vector_get(mu2_all, selInd-1), sqrt(1/gsl_vector_get(zeta2_all, selInd-1)), gsl_vector_get(y2L, i), 100000, gsl_vector_get(y2L_neginf, i), 1, &y2_temp);
            
            if(y1_temp < y2_temp)
            {
                selInd = MFunc_BAFT_DP(-100000, 100000, 1, 1, mu3_all, zeta3_all, r3Uniq_count, u3, eta3, rr);
                
                gap = rnorm(eta3 + gsl_vector_get(mu3_all, selInd-1), sqrt(1/gsl_vector_get(zeta3_all, selInd-1)));
                
                if(y1_temp < gap)
                {
                    y2_temp = gap + log(1+exp(y1_temp - gap));
                }else
                {
                    y2_temp = y1_temp + log(1+exp(gap-y1_temp));
                }
                
                if(y1_temp == y2_temp)
                {
                    y2_temp = log(exp(y1_temp) + exp(gap));
                }
                if(y1_temp == y2_temp)
                {
                    y2_temp = y1_temp + 0.000000000001;
                }
                
                /*                 */
                if(gsl_vector_get(y1_NA, i) == 1)
                {
                    set_r3_mu3_zeta3(r3, mu3_vec, zeta3_vec, mu3_all, zeta3_all, y1_temp, y2_temp, c0_neginf, xbeta3, gamma, r3Uniq, r3Uniq_count, i, u3, mu03, sigSq03, a03, b03, tau3,rr);
                    gsl_vector_set(y1_NA, i, 0);
                    c_uniq_h3(r3, r3Uniq, r3Uniq_count, mu3_all, zeta3_all, mu3_vec, zeta3_vec, y1_NA, &u3);
                }else
                {
                    gsl_vector_set(y1_NA, i, 0);
                }
                
            }else
            {
                y1_temp = y2_temp;
                if(gsl_vector_get(y1_NA, i) == 0)
                {
                    gsl_vector_set(y1_NA, i, 1);
                    gsl_vector_set(r3, i, 0);
                    gsl_vector_set(mu3_vec, i, -100000);
                    gsl_vector_set(zeta3_vec, i, -100000);
                    c_uniq_h3(r3, r3Uniq, r3Uniq_count, mu3_all, zeta3_all, mu3_vec, zeta3_vec, y1_NA, &u3);
                }else
                {
                    gsl_vector_set(y1_NA, i, 1);
                    gsl_vector_set(r3, i, 0);
                    gsl_vector_set(mu3_vec, i, -100000);
                    gsl_vector_set(zeta3_vec, i, -100000);
                }
            }
        }
        
        /*  non-terminal event is right-censored and terminal event is observed (or interval-censored) */
        
        if(gsl_vector_get(y1U_posinf, i) == 1 && gsl_vector_get(y2U_posinf, i) == 0)
        {
            selInd = MFunc_BAFT_DP(gsl_vector_get(y1L, i), gsl_vector_get(y1U, i), gsl_vector_get(y1L_neginf, i), 1, mu1_all, zeta1_all, r1Uniq_count, u1, eta1, rr);
            
            c_rtnorm(eta1 + gsl_vector_get(mu1_all, selInd-1), sqrt(1/gsl_vector_get(zeta1_all, selInd-1)), gsl_vector_get(y1L, i), gsl_vector_get(y1U, i), gsl_vector_get(y1L_neginf, i), 1, &y1_temp);
            
            /*  i) terminal event is interval-censored */
            if(gsl_vector_get(y2L, i) < gsl_vector_get(y2U, i))
            {
                selInd = MFunc_BAFT_DP(gsl_vector_get(y2L, i), gsl_vector_get(y2U, i), gsl_vector_get(y2L_neginf, i), 0, mu2_all, zeta2_all, r2Uniq_count, u2, eta2, rr);
                
                c_rtnorm(eta2 + gsl_vector_get(mu2_all, selInd-1), sqrt(1/gsl_vector_get(zeta2_all, selInd-1)), gsl_vector_get(y2L, i), gsl_vector_get(y2U, i), gsl_vector_get(y2L_neginf, i), 0, &y2_temp);
                
                if(y1_temp < y2_temp)
                {
                    CenNew = gsl_vector_get(y2U, i) + log(1-exp(y1_temp-gsl_vector_get(y2U, i)));
                    
                    selInd = MFunc_BAFT_DP(-100000, CenNew, 1, 0, mu3_all, zeta3_all, r3Uniq_count, u3, eta3, rr);
                    
                    c_rtnorm(eta3 + gsl_vector_get(mu3_all, selInd-1), sqrt(1/gsl_vector_get(zeta3_all, selInd-1)), -100000, CenNew, 1, 0, &gap);
                    
                    if(y1_temp < gap)
                    {
                        y2_temp = gap + log(1+exp(y1_temp - gap));
                    }else
                    {
                        y2_temp = y1_temp + log(1+exp(gap-y1_temp));
                    }
                    if(y1_temp == y2_temp)
                    {
                        y2_temp = log(exp(y1_temp) + exp(gap));
                    }
                    if(y1_temp == y2_temp)
                    {
                        y2_temp = y1_temp + 0.000000000001;
                    }
                    
                    /*                 */
                    if(gsl_vector_get(y1_NA, i) == 1)
                    {
                        set_r3_mu3_zeta3(r3, mu3_vec, zeta3_vec, mu3_all, zeta3_all, y1_temp, y2_temp, c0_neginf, xbeta3, gamma, r3Uniq, r3Uniq_count, i, u3, mu03, sigSq03, a03, b03, tau3,rr);
                        gsl_vector_set(y1_NA, i, 0);
                        c_uniq_h3(r3, r3Uniq, r3Uniq_count, mu3_all, zeta3_all, mu3_vec, zeta3_vec, y1_NA, &u3);
                    }else
                    {
                        gsl_vector_set(y1_NA, i, 0);
                    }
                    
                }else
                {
                    y1_temp = y2_temp;
                    if(gsl_vector_get(y1_NA, i) == 0)
                    {
                        gsl_vector_set(y1_NA, i, 1);
                        gsl_vector_set(r3, i, 0);
                        gsl_vector_set(mu3_vec, i, -100000);
                        gsl_vector_set(zeta3_vec, i, -100000);
                        c_uniq_h3(r3, r3Uniq, r3Uniq_count, mu3_all, zeta3_all, mu3_vec, zeta3_vec, y1_NA, &u3);
                    }else
                    {
                        gsl_vector_set(y1_NA, i, 1);
                        gsl_vector_set(r3, i, 0);
                        gsl_vector_set(mu3_vec, i, -100000);
                        gsl_vector_set(zeta3_vec, i, -100000);
                    }
                } /*  ii) terminal event is observed */
            }else if(gsl_vector_get(y2L, i) == gsl_vector_get(y2U, i))
            {
                y2_temp = gsl_vector_get(y2L, i);
                
                if(y1_temp < y2_temp)
                {
                    /*                 */
                    if(gsl_vector_get(y1_NA, i) == 1)
                    {
                        set_r3_mu3_zeta3(r3, mu3_vec, zeta3_vec, mu3_all, zeta3_all, y1_temp, y2_temp, c0_neginf, xbeta3, gamma, r3Uniq, r3Uniq_count, i, u3, mu03, sigSq03, a03, b03, tau3,rr);
                        gsl_vector_set(y1_NA, i, 0);
                        c_uniq_h3(r3, r3Uniq, r3Uniq_count, mu3_all, zeta3_all, mu3_vec, zeta3_vec, y1_NA, &u3);
                    }else
                    {
                        gsl_vector_set(y1_NA, i, 0);
                    }
                    
                }else
                {
                    y1_temp = y2_temp;
                    if(gsl_vector_get(y1_NA, i) == 0)
                    {
                        gsl_vector_set(y1_NA, i, 1);
                        gsl_vector_set(r3, i, 0);
                        gsl_vector_set(mu3_vec, i, -100000);
                        gsl_vector_set(zeta3_vec, i, -100000);
                        c_uniq_h3(r3, r3Uniq, r3Uniq_count, mu3_all, zeta3_all, mu3_vec, zeta3_vec, y1_NA, &u3);
                    }else
                    {
                        gsl_vector_set(y1_NA, i, 1);
                        gsl_vector_set(r3, i, 0);
                        gsl_vector_set(mu3_vec, i, -100000);
                        gsl_vector_set(zeta3_vec, i, -100000);
                    }
                }
            }
        }
        
        /*  non-terminal event is observed (or interval-censored) and terminal event is right-censored */
        
        if(gsl_vector_get(y1U_posinf, i) == 0 && gsl_vector_get(y2U_posinf, i) == 1)
        {
            /*  i) non-terminal event is interval-censored */
            if(gsl_vector_get(y1L, i) < gsl_vector_get(y1U, i))
            {
                selInd = MFunc_BAFT_DP(gsl_vector_get(y1L, i), gsl_vector_get(y1U, i), gsl_vector_get(y1L_neginf, i), 0, mu1_all, zeta1_all, r1Uniq_count, u1, eta1, rr);
                
                c_rtnorm(eta1 + gsl_vector_get(mu1_all, selInd-1), sqrt(1/gsl_vector_get(zeta1_all, selInd-1)), gsl_vector_get(y1L, i), gsl_vector_get(y1U, i), gsl_vector_get(y1L_neginf, i), 0, &y1_temp);
                
                /*  ii) non-terminal event is observed */
            }else if(gsl_vector_get(y1L, i) == gsl_vector_get(y1U, i))
            {
                y1_temp = gsl_vector_get(y1L, i);
            }
            
            if(y1_temp < gsl_vector_get(y2L, i))
            {
                CenNew = gsl_vector_get(y2L, i) + log(1-exp(y1_temp-gsl_vector_get(y2L, i)));
                selInd = MFunc_BAFT_DP(CenNew, 100000, 0, 1, mu3_all, zeta3_all, r3Uniq_count, u3, eta3, rr);
                c_rtnorm(eta3 + gsl_vector_get(mu3_all, selInd-1), sqrt(1/gsl_vector_get(zeta3_all, selInd-1)), CenNew, 100000, 0, 1, &gap);
            }else
            {
                selInd = MFunc_BAFT_DP(-100000, 100000, 1, 1, mu3_all, zeta3_all, r3Uniq_count, u3, eta3, rr);
                gap = rnorm(eta3 + gsl_vector_get(mu3_all, selInd-1), sqrt(1/gsl_vector_get(zeta3_all, selInd-1)));
            }
            
            if(y1_temp < gap)
            {
                y2_temp = gap + log(1+exp(y1_temp - gap));
            }else
            {
                y2_temp = y1_temp + log(1+exp(gap-y1_temp));
            }
            if(y1_temp == y2_temp)
            {
                y2_temp = log(exp(y1_temp) + exp(gap));
            }
            if(y1_temp == y2_temp)
            {
                y2_temp = y1_temp + 0.000000000001;
            }
            
            /*                 */
            if(gsl_vector_get(y1_NA, i) == 1)
            {
                set_r3_mu3_zeta3(r3, mu3_vec, zeta3_vec, mu3_all, zeta3_all, y1_temp, y2_temp, c0_neginf, xbeta3, gamma, r3Uniq, r3Uniq_count, i, u3, mu03, sigSq03, a03, b03, tau3,rr);
                gsl_vector_set(y1_NA, i, 0);
                c_uniq_h3(r3, r3Uniq, r3Uniq_count, mu3_all, zeta3_all, mu3_vec, zeta3_vec, y1_NA, &u3);
            }else
            {
                gsl_vector_set(y1_NA, i, 0);
            }
        }
        
        
        /*  both events are observed (or interval-censored)*/
        
        if(gsl_vector_get(y1U_posinf, i) == 0 && gsl_vector_get(y2U_posinf, i) == 0)
        {
            /*  non-terminal event is interval-censored */
            if(gsl_vector_get(y1L, i) < gsl_vector_get(y1U, i))
            {
                selInd = MFunc_BAFT_DP(gsl_vector_get(y1L, i), gsl_vector_get(y1U, i), gsl_vector_get(y1L_neginf, i), 0, mu1_all, zeta1_all, r1Uniq_count, u1, eta1, rr);
                
                c_rtnorm(eta1 + gsl_vector_get(mu1_all, selInd-1), sqrt(1/gsl_vector_get(zeta1_all, selInd-1)), gsl_vector_get(y1L, i), gsl_vector_get(y1U, i), gsl_vector_get(y1L_neginf, i), 0, &y1_temp);
                /*  non-terminal event is observed */
            }else if(gsl_vector_get(y1L, i) == gsl_vector_get(y1U, i))
            {
                y1_temp = gsl_vector_get(y1L, i);
            }
            
            /*  terminal event is interval-censored */
            if(gsl_vector_get(y2L, i) < gsl_vector_get(y2U, i))
            {
                if(y1_temp < gsl_vector_get(y2L, i))
                {
                    CenNew1 = gsl_vector_get(y2L, i) + log(1-exp(y1_temp-gsl_vector_get(y2L, i)));
                    CenNew2 = gsl_vector_get(y2U, i) + log(1-exp(y1_temp-gsl_vector_get(y2U, i)));
                    selInd = MFunc_BAFT_DP(CenNew1, CenNew2, 0, 0, mu3_all, zeta3_all, r3Uniq_count, u3, eta3, rr);
                    
                    c_rtnorm(eta3 + gsl_vector_get(mu3_all, selInd-1), sqrt(1/gsl_vector_get(zeta3_all, selInd-1)), CenNew1, CenNew2, 0, 0, &gap);
                    
                }else if(y1_temp < gsl_vector_get(y2U, i) && y1_temp > gsl_vector_get(y2L, i))
                {
                    CenNew2 = gsl_vector_get(y2U, i) + log(1-exp(y1_temp-gsl_vector_get(y2U, i)));
                    selInd = MFunc_BAFT_DP(-100000, CenNew2, 1, 0, mu3_all, zeta3_all, r3Uniq_count, u3, eta3, rr);
                    
                    c_rtnorm(eta3 + gsl_vector_get(mu3_all, selInd-1), sqrt(1/gsl_vector_get(zeta3_all, selInd-1)), -100000, CenNew2, 1, 0, &gap);
                }
                if(y1_temp < gap)
                {
                    y2_temp = gap + log(1+exp(y1_temp - gap));
                }else
                {
                    y2_temp = y1_temp + log(1+exp(gap-y1_temp));
                }
                if(y1_temp == y2_temp)
                {
                    y2_temp = log(exp(y1_temp) + exp(gap));
                }
                if(y1_temp == y2_temp)
                {
                    y2_temp = y1_temp + 0.000000000001;
                }
                
                /*  terminal event is observed */
            }else if(gsl_vector_get(y2L, i) == gsl_vector_get(y2U, i))
            {
                y2_temp = gsl_vector_get(y2L, i);
                gap = log(exp(y2_temp) - exp(y1_temp));
            }
            
            /*                 */
            if(gsl_vector_get(y1_NA, i) == 1)
            {
                set_r3_mu3_zeta3(r3, mu3_vec, zeta3_vec, mu3_all, zeta3_all, y1_temp, y2_temp, c0_neginf, xbeta3, gamma, r3Uniq, r3Uniq_count, i, u3, mu03, sigSq03, a03, b03, tau3,rr);
                gsl_vector_set(y1_NA, i, 0);
                c_uniq_h3(r3, r3Uniq, r3Uniq_count, mu3_all, zeta3_all, mu3_vec, zeta3_vec, y1_NA, &u3);
            }else
            {
                gsl_vector_set(y1_NA, i, 0);
            }
        }
        
        gsl_vector_set(y1, i, y1_temp);
        gsl_vector_set(y2, i, y2_temp);
        
    }
    
    /*     */
    c_uniq_h3(r3, r3Uniq, r3Uniq_count, mu3_all, zeta3_all, mu3_vec, zeta3_vec, y1_NA, nClass_DP3);
    
    
    gsl_vector_free(xbeta1);
    gsl_vector_free(xbeta2);
    gsl_vector_free(xbeta3);
    return;
}






/* updating r1, mu1, eta1 using DPM of Normals */

void BAFT_DPscr_update_mu_zeta1(gsl_vector *c0,
                                gsl_vector *c0_neginf,
                                gsl_matrix *X1,
                                gsl_vector *y1,
                                gsl_vector *y2,
                                gsl_vector *y1_NA,
                                gsl_vector *beta1,
                                gsl_vector *gamma,
                                gsl_vector *r1,
                                gsl_vector *mu1_all,
                                gsl_vector *zeta1_all,
                                gsl_vector *r1Uniq,
                                gsl_vector *r1Uniq_count,
                                double tau1,
                                double a01,
                                double b01,
                                double mu01,
                                double sigSq01,
                                double beta01_prop_var,
                                double zeta1_prop_var,
                                int *accept_beta01,
                                int *accept_zeta1,
                                int *nClass_DP1,
                                gsl_rng *rr)
{
    int i, j, k, n_ir, r_ind, wh_ri, u_, LL_negInf;
    double b_mc, sum_prob, val, mu, zeta, eta, cond;
    double mean, sd, mean_prop, sd_prop;
    
    int u = *nClass_DP1;
    int n = y1 -> size;
    
    gsl_vector *xbeta1 = gsl_vector_calloc(n);
    gsl_blas_dgemv(CblasNoTrans, 1, X1, beta1, 0, xbeta1);
    
    gsl_vector *mu_temp = gsl_vector_calloc(n);
    gsl_vector *zeta_temp = gsl_vector_calloc(n);
    
    gsl_vector *mu_dummy = gsl_vector_calloc(n);
    gsl_vector *zeta_dummy = gsl_vector_calloc(n);
    gsl_vector *rUniq_dummy = gsl_vector_calloc(n);
    
    /*************************************************************/
    
    /* Step 1: updating latent classes (r) */
    
    /*************************************************************/
    
    /* */
    for(i = 0; i < n; i++)
    {
        eta = gsl_vector_get(xbeta1, i) + gsl_vector_get(gamma, i);
        LL_negInf = (int) gsl_vector_get(c0_neginf, i);
        
        /* identfy "u" unique values */
        c_uniq(r1, r1Uniq, r1Uniq_count, mu1_all, zeta1_all, &u);
        
        for(j = 0; j < u; j++)
        {
            if(gsl_vector_get(r1, i) == gsl_vector_get(r1Uniq, j)) wh_ri = j;
        }
        
        if(gsl_vector_get(r1Uniq_count, wh_ri) > 1) /* not singleton */
        {
            u_ = u;
        }else /* singleton */
        {
            u_ = u-1;
        }
        
        gsl_vector *prob = gsl_vector_calloc(u_+1);
        
        
        zeta = INFINITY;
        while(zeta == INFINITY || isnan(zeta))
        {
            zeta = rgamma(a01, 1/b01);
        }
        
        mu = NAN;
        while(mu == INFINITY || isnan(mu))
        {
            mu = rnorm(mu01, sqrt(sigSq01));
        }
        
        if(gsl_vector_get(r1Uniq_count, wh_ri) > 1) /* not singleton: proposing to create a new component */
        {
            for(j = 0; j < u_; j++)
            {
                n_ir = gsl_vector_get(r1Uniq_count, j);
                if(gsl_vector_get(r1, i) == gsl_vector_get(r1Uniq, j)) n_ir -= 1;
                
                mean = gsl_vector_get(mu1_all, j) + eta;
                sd = sqrt(1/gsl_vector_get(zeta1_all, j));
                
                if(gsl_vector_get(y1_NA, i) == 0)
                {
                    val = dnorm(gsl_vector_get(y1, i), mean, sd, 0);
                }else if(gsl_vector_get(y1_NA, i) == 1)
                {
                    val = pnorm(gsl_vector_get(y2, i), mean, sd, 0, 0);
                }
                
                if(LL_negInf == 0)
                {
                    cond = pnorm(gsl_vector_get(c0, i), mean, sd, 0, 0);
                    if(cond < exp(-700))
                    {
                        val = 0;
                    }else
                    {
                        val /= pnorm(gsl_vector_get(c0, i), mean, sd, 0, 0);
                    }
                }
                
                val *= (double) n_ir / (double)(n-1+tau1);
                gsl_vector_set(prob, j, val);
            }
            
            mean = mu + eta;
            sd = sqrt(1/zeta);
            
            if(gsl_vector_get(y1_NA, i) == 0)
            {
                val = dnorm(gsl_vector_get(y1, i), mean, sd, 0);
            }else if(gsl_vector_get(y1_NA, i) == 1)
            {
                val = pnorm(gsl_vector_get(y2, i), mean, sd, 0, 0);
            }
            if(LL_negInf == 0)
            {
                cond = pnorm(gsl_vector_get(c0, i), mean, sd, 0, 0);
                if(cond < exp(-700))
                {
                    val = 0;
                }else
                {
                    val /= pnorm(gsl_vector_get(c0, i), mean, sd, 0, 0);
                }
            }
            val *= (double) tau1 / (double)(n-1+tau1);
            
            gsl_vector_set(prob, u_, val);
            
            sum_prob = 0;
            for(k = 0; k < u_+1; k++)
            {
                sum_prob += gsl_vector_get(prob, k);
            }
            
            if(sum_prob > pow(10, -300))
            {
                b_mc = 1/sum_prob;
                gsl_vector_scale(prob, b_mc);
                
                r_ind = c_multinom_sample(rr, prob, u_+1);
                
                if(r_ind <= u_)
                {
                    gsl_vector_set(r1, i, gsl_vector_get(r1Uniq, r_ind-1));
                }else if(r_ind == u_+1)
                {
                    gsl_vector_set(r1, i, gsl_vector_max(r1Uniq)+1);
                    gsl_vector_set(r1Uniq, u_, gsl_vector_max(r1Uniq)+1);
                    gsl_vector_set(mu1_all, u_, mu);
                    gsl_vector_set(zeta1_all, u_, zeta);
                    u += 1;
                }
            }
            
        }else /* singletion */
        {
            gsl_vector_memcpy(mu_temp, mu1_all);
            gsl_vector_memcpy(zeta_temp, zeta1_all);
            
            gsl_vector_set(mu_temp, wh_ri, mu);
            gsl_vector_set(zeta_temp, wh_ri, zeta);
            
            for(j = 0; j < u_+1; j++)
            {
                mean = gsl_vector_get(mu_temp, j) + eta;
                sd = sqrt(1/gsl_vector_get(zeta_temp, j));
                
                n_ir = gsl_vector_get(r1Uniq_count, j);
                if(gsl_vector_get(r1, i) == gsl_vector_get(r1Uniq, j)) n_ir -= 1;
                
                if(gsl_vector_get(y1_NA, i) == 0)
                {
                    val = dnorm(gsl_vector_get(y1, i), mean, sd, 0);
                }else if(gsl_vector_get(y1_NA, i) == 1)
                {
                    val = pnorm(gsl_vector_get(y2, i), mean, sd, 0, 0);
                }
                
                if(LL_negInf == 0)
                {
                    cond = pnorm(gsl_vector_get(c0, i), mean, sd, 0, 0);
                    if(cond < exp(-700))
                    {
                        val = 0;
                    }else
                    {
                        val /= pnorm(gsl_vector_get(c0, i), mean, sd, 0, 0);
                    }
                }
                
                if(j == wh_ri)
                {
                    val *= (double) tau1 / (double)(n-1+tau1);
                }else
                {
                    val *= (double) n_ir / (double)(n-1+tau1);
                }
                gsl_vector_set(prob, j, val);
            }
            
            sum_prob = 0;
            for(k = 0; k < u_+1; k++)
            {
                sum_prob += gsl_vector_get(prob, k);
            }
            
            if(sum_prob > pow(10, -300))
            {
                b_mc = 1/sum_prob;
                gsl_vector_scale(prob, b_mc);
                
                r_ind = c_multinom_sample(rr, prob, u_+1);
                gsl_vector_set(r1, i, gsl_vector_get(r1Uniq, r_ind-1));
                
                if(r_ind-1 != wh_ri)
                {
                    if(wh_ri < u_)
                    {
                        gsl_vector_memcpy(rUniq_dummy, r1Uniq);
                        gsl_vector_memcpy(mu_dummy, mu1_all);
                        gsl_vector_memcpy(zeta_dummy, zeta1_all);
                        for(k=wh_ri;k<u_;k++)
                        {
                            gsl_vector_set(r1Uniq, k, gsl_vector_get(rUniq_dummy, k+1));
                            gsl_vector_set(mu_temp, k, gsl_vector_get(mu_dummy, k+1));
                            gsl_vector_set(zeta_temp, k, gsl_vector_get(zeta_dummy, k+1));
                        }
                    }
                    gsl_vector_set(r1Uniq, u_, 0);
                    gsl_vector_set(mu_temp, u_, 0);
                    gsl_vector_set(zeta_temp, u_, 0);
                    u -= 1;
                }
                gsl_vector_memcpy(mu1_all, mu_temp);
                gsl_vector_memcpy(zeta1_all, zeta_temp);
            }
        }
        
        gsl_vector_free(prob);
        
        for(j = u_+1; j < n; j++)
        {
            gsl_vector_set(r1Uniq, j, 0);
            gsl_vector_set(mu1_all, j, 0);
            gsl_vector_set(zeta1_all, j, 0);
        }
    }
    
    /*************************************************************/
    
    /* Step 2: update (beta0, sigSq) using the posterior distribution that is based on {epsilon_i :i∈{k:rk =r}}. */
    
    /*************************************************************/
    
    /* identfy "u" unique values */
    c_uniq(r1, r1Uniq, r1Uniq_count, mu1_all, zeta1_all, &u);
    *nClass_DP1 = u;
    
    double mu_prop, zeta_prop, mu_ini, zeta_ini, temp, temp_prop;
    double loglh, loglh_prop, logprior, logprior_prop, logprop, logprop_prop;
    double logR;
    int uu;
    
    for(j = 0; j < u; j++)
    {
        mu_ini = gsl_vector_get(mu1_all, j);
        zeta_ini = gsl_vector_get(zeta1_all, j);
        
        /* update mu */
        mu_prop = rnorm((double) mu_ini, (double) sqrt(beta01_prop_var));
        
        loglh=0;
        loglh_prop=0;
        for(i = 0; i < n; i++)
        {
            if(gsl_vector_get(r1, i) == gsl_vector_get(r1Uniq, j))
            {
                LL_negInf = (int) gsl_vector_get(c0_neginf, i);
                eta = gsl_vector_get(xbeta1, i) + gsl_vector_get(gamma, i);
                mean = mu_ini + eta;
                sd = sqrt(1/zeta_ini);
                mean_prop = mu_prop + eta;
                
                if(gsl_vector_get(y1_NA, i) == 0)
                {
                    temp = dnorm(gsl_vector_get(y1, i), mean, sd, 1);
                    temp_prop = dnorm(gsl_vector_get(y1, i), mean_prop, sd, 1);
                    
                }else if(gsl_vector_get(y1_NA, i) == 1)
                {
                    temp = pnorm(gsl_vector_get(y2, i), mean, sd, 0, 1);
                    temp_prop = pnorm(gsl_vector_get(y2, i), mean_prop, sd, 0, 1);
                }
                
                if(LL_negInf == 0)
                {
                    temp -= pnorm(gsl_vector_get(c0, i), mean, sd, 0, 1);
                    temp_prop -= pnorm(gsl_vector_get(c0, i), mean_prop, sd, 0, 1);
                }
                loglh += temp;
                loglh_prop += temp_prop;
            }
        }
        
        logprior = dnorm(mu_ini, mu01, sqrt(sigSq01), 1);
        logprior_prop = dnorm(mu_prop, mu01, sqrt(sigSq01), 1);
        
        logprop = dnorm(mu_ini, mu_prop, sqrt(beta01_prop_var), 1);
        logprop_prop = dnorm(mu_prop, mu_ini, sqrt(beta01_prop_var), 1);
        
        logR = loglh_prop - loglh + logprior_prop - logprior + logprop - logprop_prop;
        uu = log(runif(0, 1)) < logR;
        
        if(uu == 1)
        {
            gsl_vector_set(mu1_all, j, mu_prop);
            *accept_beta01 += 1;
        }
        
        /* update zeta */
        mu_ini = gsl_vector_get(mu1_all, j);
        
        zeta_prop = INFINITY;
        
        while(zeta_prop == INFINITY || isnan(zeta_prop))
        {
            zeta_prop = rgamma(pow(zeta_ini, 2)/zeta1_prop_var, zeta1_prop_var/zeta_ini);
        }
        
        loglh=0;
        loglh_prop=0;
        for(i = 0; i < n; i++)
        {
            if(gsl_vector_get(r1, i) == gsl_vector_get(r1Uniq, j))
            {
                LL_negInf = (int) gsl_vector_get(c0_neginf, i);
                eta = gsl_vector_get(xbeta1, i) + gsl_vector_get(gamma, i);
                mean = mu_ini + eta;
                sd = sqrt(1/zeta_ini);
                sd_prop = sqrt(1/zeta_prop);
                
                if(gsl_vector_get(y1_NA, i) == 0)
                {
                    temp = dnorm(gsl_vector_get(y1, i), mean, sd, 1);
                    temp_prop = dnorm(gsl_vector_get(y1, i), mean, sd_prop, 1);
                    
                }else if(gsl_vector_get(y1_NA, i) == 1)
                {
                    temp = pnorm(gsl_vector_get(y2, i), mean, sd, 0, 1);
                    temp_prop = pnorm(gsl_vector_get(y2, i), mean, sd_prop, 0, 1);
                }
                
                if(LL_negInf == 0)
                {
                    temp -= pnorm(gsl_vector_get(c0, i), mean, sd, 0, 1);
                    temp_prop -= pnorm(gsl_vector_get(c0, i), mean, sd_prop, 0, 1);
                }
                loglh += temp;
                loglh_prop += temp_prop;
            }
        }
        
        logprior = dgamma(zeta_ini, a01, 1/b01, 1);
        logprior_prop = dgamma(zeta_prop, a01, 1/b01, 1);
        
        logprop = dgamma(zeta_ini, pow(zeta_prop, 2)/zeta1_prop_var, zeta1_prop_var/zeta_prop, 1);
        logprop_prop = dgamma(zeta_prop, pow(zeta_ini, 2)/zeta1_prop_var, zeta1_prop_var/zeta_ini, 1);
        
        
        logR = loglh_prop - loglh + logprior_prop - logprior + logprop - logprop_prop;
        uu = log(runif(0, 1)) < logR;
        
        if(uu == 1)
        {
            gsl_vector_set(zeta1_all, j, zeta_prop);
            *accept_zeta1 += 1;
        }
    }
    
    gsl_vector_free(xbeta1);
    gsl_vector_free(mu_temp);
    gsl_vector_free(zeta_temp);
    gsl_vector_free(mu_dummy);
    gsl_vector_free(zeta_dummy);
    gsl_vector_free(rUniq_dummy);
    
    return;
    
}











/* updating r2, mu2, eta2 using DPM of Normals */

void BAFT_DPscr_update_mu_zeta2(gsl_vector *c0,
                                gsl_vector *c0_neginf,
                                gsl_matrix *X2,
                                gsl_vector *y1,
                                gsl_vector *y2,
                                gsl_vector *y1_NA,
                                gsl_vector *beta2,
                                gsl_vector *gamma,
                                gsl_vector *r2,
                                gsl_vector *mu2_all,
                                gsl_vector *zeta2_all,
                                gsl_vector *r2Uniq,
                                gsl_vector *r2Uniq_count,
                                double tau2,
                                double a02,
                                double b02,
                                double mu02,
                                double sigSq02,
                                double beta02_prop_var,
                                double zeta2_prop_var,
                                int *accept_beta02,
                                int *accept_zeta2,
                                int *nClass_DP2,
                                gsl_rng *rr)
{
    int i, j, k, n_ir, r_ind, wh_ri, u_, LL_negInf;
    double b_mc, sum_prob, val, mu, zeta, eta, cond;
    double mean, sd, mean_prop, sd_prop;
    
    int u = *nClass_DP2;
    int n = y2 -> size;
    
    gsl_vector *xbeta2 = gsl_vector_calloc(n);
    gsl_blas_dgemv(CblasNoTrans, 1, X2, beta2, 0, xbeta2);
    
    gsl_vector *mu_temp = gsl_vector_calloc(n);
    gsl_vector *zeta_temp = gsl_vector_calloc(n);
    
    gsl_vector *mu_dummy = gsl_vector_calloc(n);
    gsl_vector *zeta_dummy = gsl_vector_calloc(n);
    gsl_vector *rUniq_dummy = gsl_vector_calloc(n);
    
    /*************************************************************/
    
    /* Step 1: updating latent classes (r) */
    
    /*************************************************************/
    
    /* */
    for(i = 0; i < n; i++)
    {
        eta = gsl_vector_get(xbeta2, i) + gsl_vector_get(gamma, i);
        LL_negInf = (int) gsl_vector_get(c0_neginf, i);
        
        /* identfy "u" unique values */
        c_uniq(r2, r2Uniq, r2Uniq_count, mu2_all, zeta2_all, &u);
        
        for(j = 0; j < u; j++)
        {
            if(gsl_vector_get(r2, i) == gsl_vector_get(r2Uniq, j)) wh_ri = j;
        }
        
        if(gsl_vector_get(r2Uniq_count, wh_ri) > 1) /* not singleton */
        {
            u_ = u;
        }else /* singleton */
        {
            u_ = u-1;
        }
        
        gsl_vector *prob = gsl_vector_calloc(u_+1);
        
        zeta = INFINITY;
        
        while(zeta == INFINITY || isnan(zeta))
        {
            zeta = rgamma(a02, 1/b02);
        }
        
        mu = NAN;
        
        while(mu == INFINITY || isnan(mu))
        {
            mu = rnorm(mu02, sqrt(sigSq02));
        }
        
        
        if(gsl_vector_get(r2Uniq_count, wh_ri) > 1) /* not singleton: proposing to create a new component */
        {
            for(j = 0; j < u_; j++)
            {
                n_ir = gsl_vector_get(r2Uniq_count, j);
                if(gsl_vector_get(r2, i) == gsl_vector_get(r2Uniq, j)) n_ir -= 1;
                
                mean = gsl_vector_get(mu2_all, j) + eta;
                sd = sqrt(1/gsl_vector_get(zeta2_all, j));
                
                if(gsl_vector_get(y1_NA, i) == 0)
                {
                    val = pnorm(gsl_vector_get(y1, i), mean, sd, 0, 0);
                }else if(gsl_vector_get(y1_NA, i) == 1)
                {
                    val = dnorm(gsl_vector_get(y2, i), mean, sd, 0);
                }
                
                if(LL_negInf == 0)
                {
                    cond = pnorm(gsl_vector_get(c0, i), mean, sd, 0, 0);
                    if(cond < exp(-700))
                    {
                        val = 0;
                    }else
                    {
                        val /= pnorm(gsl_vector_get(c0, i), mean, sd, 0, 0);
                    }
                }
                
                val *= (double) n_ir / (double)(n-1+tau2);
                gsl_vector_set(prob, j, val);
            }
            
            mean = mu + eta;
            sd = sqrt(1/zeta);
            
            if(gsl_vector_get(y1_NA, i) == 0)
            {
                val = pnorm(gsl_vector_get(y1, i), mean, sd, 0, 0);
            }else if(gsl_vector_get(y1_NA, i) == 1)
            {
                val = dnorm(gsl_vector_get(y2, i), mean, sd, 0);
            }
            if(LL_negInf == 0)
            {
                cond = pnorm(gsl_vector_get(c0, i), mean, sd, 0, 0);
                if(cond < exp(-700))
                {
                    val = 0;
                }else
                {
                    val /= pnorm(gsl_vector_get(c0, i), mean, sd, 0, 0);
                }
            }
            val *= (double) tau2 / (double)(n-1+tau2);
            
            gsl_vector_set(prob, u_, val);
            
            sum_prob = 0;
            for(k = 0; k < u_+1; k++)
            {
                sum_prob += gsl_vector_get(prob, k);
            }
            
            if(sum_prob > pow(10, -300))
            {
                b_mc = 1/sum_prob;
                gsl_vector_scale(prob, b_mc);
                
                r_ind = c_multinom_sample(rr, prob, u_+1);
                
                if(r_ind <= u_)
                {
                    gsl_vector_set(r2, i, gsl_vector_get(r2Uniq, r_ind-1));
                }else if(r_ind == u_+1)
                {
                    gsl_vector_set(r2, i, gsl_vector_max(r2Uniq)+1);
                    gsl_vector_set(r2Uniq, u_, gsl_vector_max(r2Uniq)+1);
                    gsl_vector_set(mu2_all, u_, mu);
                    gsl_vector_set(zeta2_all, u_, zeta);
                    u += 1;
                }
            }
            
        }else /* singletion */
        {
            gsl_vector_memcpy(mu_temp, mu2_all);
            gsl_vector_memcpy(zeta_temp, zeta2_all);
            
            gsl_vector_set(mu_temp, wh_ri, mu);
            gsl_vector_set(zeta_temp, wh_ri, zeta);
            
            for(j = 0; j < u_+1; j++)
            {
                mean = gsl_vector_get(mu_temp, j) + eta;
                sd = sqrt(1/gsl_vector_get(zeta_temp, j));
                
                n_ir = gsl_vector_get(r2Uniq_count, j);
                if(gsl_vector_get(r2, i) == gsl_vector_get(r2Uniq, j)) n_ir -= 1;
                
                if(gsl_vector_get(y1_NA, i) == 0)
                {
                    val = pnorm(gsl_vector_get(y1, i), mean, sd, 0, 0);
                }else if(gsl_vector_get(y1_NA, i) == 1)
                {
                    val = dnorm(gsl_vector_get(y2, i), mean, sd, 0);
                }
                
                if(LL_negInf == 0)
                {
                    cond = pnorm(gsl_vector_get(c0, i), mean, sd, 0, 0);
                    if(cond < exp(-700))
                    {
                        val = 0;
                    }else
                    {
                        val /= pnorm(gsl_vector_get(c0, i), mean, sd, 0, 0);
                    }
                }
                
                if(j == wh_ri)
                {
                    val *= (double) tau2 / (double)(n-1+tau2);
                }else
                {
                    val *= (double) n_ir / (double)(n-1+tau2);
                }
                gsl_vector_set(prob, j, val);
            }
            
            sum_prob = 0;
            for(k = 0; k < u_+1; k++)
            {
                sum_prob += gsl_vector_get(prob, k);
            }
            
            if(sum_prob > pow(10, -300))
            {
                b_mc = 1/sum_prob;
                gsl_vector_scale(prob, b_mc);
                
                r_ind = c_multinom_sample(rr, prob, u_+1);
                gsl_vector_set(r2, i, gsl_vector_get(r2Uniq, r_ind-1));
                
                if(r_ind-1 != wh_ri)
                {
                    if(wh_ri < u_)
                    {
                        gsl_vector_memcpy(rUniq_dummy, r2Uniq);
                        gsl_vector_memcpy(mu_dummy, mu2_all);
                        gsl_vector_memcpy(zeta_dummy, zeta2_all);
                        for(k=wh_ri;k<u_;k++)
                        {
                            gsl_vector_set(r2Uniq, k, gsl_vector_get(rUniq_dummy, k+1));
                            gsl_vector_set(mu_temp, k, gsl_vector_get(mu_dummy, k+1));
                            gsl_vector_set(zeta_temp, k, gsl_vector_get(zeta_dummy, k+1));
                        }
                    }
                    gsl_vector_set(r2Uniq, u_, 0);
                    gsl_vector_set(mu_temp, u_, 0);
                    gsl_vector_set(zeta_temp, u_, 0);
                    u -= 1;
                }
                gsl_vector_memcpy(mu2_all, mu_temp);
                gsl_vector_memcpy(zeta2_all, zeta_temp);
            }
        }
        
        gsl_vector_free(prob);
        
        for(j = u_+1; j < n; j++)
        {
            gsl_vector_set(r2Uniq, j, 0);
            gsl_vector_set(mu2_all, j, 0);
            gsl_vector_set(zeta2_all, j, 0);
        }
    }
    
    /*************************************************************/
    
    /* Step 2: update (beta0, sigSq) using the posterior distribution that is based on {epsilon_i :i∈{k:rk =r}}. */
    
    /*************************************************************/
    
    /* identfy "u" unique values */
    c_uniq(r2, r2Uniq, r2Uniq_count, mu2_all, zeta2_all, &u);
    *nClass_DP2 = u;
    
    double mu_prop, zeta_prop, mu_ini, zeta_ini, temp, temp_prop;
    double loglh, loglh_prop, logprior, logprior_prop, logprop, logprop_prop;
    double logR;
    int uu;
    
    for(j = 0; j < u; j++)
    {
        mu_ini = gsl_vector_get(mu2_all, j);
        zeta_ini = gsl_vector_get(zeta2_all, j);
        
        /* update mu */
        mu_prop = rnorm((double) mu_ini, (double) sqrt(beta02_prop_var));
        
        loglh=0;
        loglh_prop=0;
        for(i = 0; i < n; i++)
        {
            if(gsl_vector_get(r2, i) == gsl_vector_get(r2Uniq, j))
            {
                LL_negInf = (int) gsl_vector_get(c0_neginf, i);
                eta = gsl_vector_get(xbeta2, i) + gsl_vector_get(gamma, i);
                mean = mu_ini + eta;
                sd = sqrt(1/zeta_ini);
                mean_prop = mu_prop + eta;
                
                if(gsl_vector_get(y1_NA, i) == 0)
                {
                    temp = pnorm(gsl_vector_get(y1, i), mean, sd, 0, 1);
                    temp_prop = pnorm(gsl_vector_get(y1, i), mean_prop, sd, 0, 1);
                    
                }else if(gsl_vector_get(y1_NA, i) == 1)
                {
                    temp = dnorm(gsl_vector_get(y2, i), mean, sd, 1);
                    temp_prop = dnorm(gsl_vector_get(y2, i), mean_prop, sd, 1);
                }
                
                if(LL_negInf == 0)
                {
                    temp -= pnorm(gsl_vector_get(c0, i), mean, sd, 0, 1);
                    temp_prop -= pnorm(gsl_vector_get(c0, i), mean_prop, sd, 0, 1);
                }
                loglh += temp;
                loglh_prop += temp_prop;
            }
        }
        
        logprior = dnorm(mu_ini, mu02, sqrt(sigSq02), 1);
        logprior_prop = dnorm(mu_prop, mu02, sqrt(sigSq02), 1);
        
        logprop = dnorm(mu_ini, mu_prop, sqrt(beta02_prop_var), 1);
        logprop_prop = dnorm(mu_prop, mu_ini, sqrt(beta02_prop_var), 1);
        
        logR = loglh_prop - loglh + logprior_prop - logprior + logprop - logprop_prop;
        uu = log(runif(0, 1)) < logR;
        
        if(uu == 1)
        {
            gsl_vector_set(mu2_all, j, mu_prop);
            *accept_beta02 += 1;
        }
        
        /* update zeta */
        mu_ini = gsl_vector_get(mu2_all, j);
        
        zeta_prop = INFINITY;
        while(zeta_prop == INFINITY || isnan(zeta_prop))
        {
            zeta_prop = rgamma(pow(zeta_ini, 2)/zeta2_prop_var, zeta2_prop_var/zeta_ini);
        }
        
        loglh=0;
        loglh_prop=0;
        for(i = 0; i < n; i++)
        {
            if(gsl_vector_get(r2, i) == gsl_vector_get(r2Uniq, j))
            {
                LL_negInf = (int) gsl_vector_get(c0_neginf, i);
                eta = gsl_vector_get(xbeta2, i) + gsl_vector_get(gamma, i);
                mean = mu_ini + eta;
                sd = sqrt(1/zeta_ini);
                sd_prop = sqrt(1/zeta_prop);
                
                if(gsl_vector_get(y1_NA, i) == 0)
                {
                    temp = pnorm(gsl_vector_get(y1, i), mean, sd, 0, 1);
                    temp_prop = pnorm(gsl_vector_get(y1, i), mean, sd_prop, 0, 1);
                    
                }else if(gsl_vector_get(y1_NA, i) == 1)
                {
                    temp = dnorm(gsl_vector_get(y2, i), mean, sd, 1);
                    temp_prop = dnorm(gsl_vector_get(y2, i), mean, sd_prop, 1);
                }
                
                if(LL_negInf == 0)
                {
                    temp -= pnorm(gsl_vector_get(c0, i), mean, sd, 0, 1);
                    temp_prop -= pnorm(gsl_vector_get(c0, i), mean, sd_prop, 0, 1);
                }
                loglh += temp;
                loglh_prop += temp_prop;
            }
        }
        
        logprior = dgamma(zeta_ini, a02, 1/b02, 1);
        logprior_prop = dgamma(zeta_prop, a02, 1/b02, 1);
        
        logprop = dgamma(zeta_ini, pow(zeta_prop, 2)/zeta2_prop_var, zeta2_prop_var/zeta_prop, 1);
        logprop_prop = dgamma(zeta_prop, pow(zeta_ini, 2)/zeta2_prop_var, zeta2_prop_var/zeta_ini, 1);
        
        logR = loglh_prop - loglh + logprior_prop - logprior + logprop - logprop_prop;
        uu = log(runif(0, 1)) < logR;
        
        if(uu == 1)
        {
            gsl_vector_set(zeta2_all, j, zeta_prop);
            *accept_zeta2 += 1;
        }
    }
    
    gsl_vector_free(xbeta2);
    gsl_vector_free(mu_temp);
    gsl_vector_free(zeta_temp);
    gsl_vector_free(mu_dummy);
    gsl_vector_free(zeta_dummy);
    gsl_vector_free(rUniq_dummy);
    
    return;
    
}






/* updating r3, mu3, eta3 using DPM of Normals */

void BAFT_DPscr_update_mu_zeta3(gsl_vector *c0,
                                gsl_vector *c0_neginf,
                                gsl_matrix *X3,
                                gsl_vector *y1,
                                gsl_vector *y2,
                                gsl_vector *y1_NA,
                                gsl_vector *beta3,
                                gsl_vector *gamma,
                                gsl_vector *r3,
                                gsl_vector *mu3_all,
                                gsl_vector *mu3_vec,
                                gsl_vector *zeta3_all,
                                gsl_vector *zeta3_vec,
                                gsl_vector *r3Uniq,
                                gsl_vector *r3Uniq_count,
                                double tau3,
                                double a03,
                                double b03,
                                double mu03,
                                double sigSq03,
                                double beta03_prop_var,
                                double zeta3_prop_var,
                                int *accept_beta03,
                                int *accept_zeta3,
                                int *nClass_DP3,
                                gsl_rng *rr)
{
    int i, j, k, n_ir, r_ind, wh_ri, u_, LL_negInf;
    double b_mc, sum_prob, val, mu, zeta, eta;
    double mean, sd, mean_prop, sd_prop, yStar;
    
    int u = *nClass_DP3;
    int n = y1 -> size;
    
    gsl_vector *xbeta3 = gsl_vector_calloc(n);
    gsl_blas_dgemv(CblasNoTrans, 1, X3, beta3, 0, xbeta3);
    
    gsl_vector *mu_temp = gsl_vector_calloc(n);
    gsl_vector *zeta_temp = gsl_vector_calloc(n);
    
    gsl_vector *mu_dummy = gsl_vector_calloc(n);
    gsl_vector *zeta_dummy = gsl_vector_calloc(n);
    gsl_vector *rUniq_dummy = gsl_vector_calloc(n);
    
    
    
    /*************************************************************/
    
    /* Step 1: updating latent classes (r) */
    
    /*************************************************************/
    
    /* */
    for(i = 0; i < n; i++)
    {
        if((gsl_vector_get(y1, i) < gsl_vector_get(y2, i)) && (gsl_vector_get(y1_NA, i) == 0))
        {
            yStar = gsl_vector_get(y2, i) + log(1-exp(gsl_vector_get(y1, i)-gsl_vector_get(y2, i)));
            
            eta = gsl_vector_get(xbeta3, i) + gsl_vector_get(gamma, i);
            LL_negInf = (int) gsl_vector_get(c0_neginf, i);
            
            /* identfy "u" unique values */
            c_uniq_h3(r3, r3Uniq, r3Uniq_count, mu3_all, zeta3_all, mu3_vec, zeta3_vec, y1_NA, &u);
            
            
            for(j = 0; j < u; j++)
            {
                if(gsl_vector_get(r3, i) == gsl_vector_get(r3Uniq, j)) wh_ri = j;
            }
            
            if(gsl_vector_get(r3Uniq_count, wh_ri) > 1) /* not singleton */
            {
                u_ = u;
            }else /* singleton */
            {
                u_ = u-1;
            }
            
            gsl_vector *prob = gsl_vector_calloc(u_+1);
            
            zeta = INFINITY;
            
            while(zeta == INFINITY || isnan(zeta))
            {
                zeta = rgamma(a03, 1/b03);
            }
            
            mu = NAN;
            
            while(mu == INFINITY || isnan(mu))
            {
                mu = rnorm(mu03, sqrt(sigSq03));
            }
            
            
            if(gsl_vector_get(r3Uniq_count, wh_ri) > 1) /* not singleton: proposing to create a new component */
            {
                for(j = 0; j < u_; j++)
                {
                    n_ir = gsl_vector_get(r3Uniq_count, j);
                    if(gsl_vector_get(r3, i) == gsl_vector_get(r3Uniq, j)) n_ir -= 1;
                    
                    mean = gsl_vector_get(mu3_all, j) + eta;
                    sd = sqrt(1/gsl_vector_get(zeta3_all, j));
                    
                    val = dnorm(yStar, mean, sd, 0);
                    val *= (double) n_ir / (double)(n-1+tau3);
                    gsl_vector_set(prob, j, val);
                }
                
                mean = mu + eta;
                sd = sqrt(1/zeta);
                
                val = dnorm(yStar, mean, sd, 0);
                val *= (double) tau3 / (double)(n-1+tau3);
                
                gsl_vector_set(prob, u_, val);
                
                sum_prob = 0;
                for(k = 0; k < u_+1; k++)
                {
                    sum_prob += gsl_vector_get(prob, k);
                }
                
                if(sum_prob > pow(10, -300))
                {
                    b_mc = 1/sum_prob;
                    gsl_vector_scale(prob, b_mc);
                    
                    r_ind = c_multinom_sample(rr, prob, u_+1);
                    
                    if(r_ind <= u_)
                    {
                        gsl_vector_set(r3, i, gsl_vector_get(r3Uniq, r_ind-1));
                        gsl_vector_set(mu3_vec, i, gsl_vector_get(mu3_all, r_ind-1));
                        gsl_vector_set(zeta3_vec, i, gsl_vector_get(zeta3_all, r_ind-1));
                    }else if(r_ind == u_+1)
                    {
                        gsl_vector_set(r3, i, gsl_vector_max(r3Uniq)+1);
                        gsl_vector_set(mu3_vec, i, mu);
                        gsl_vector_set(zeta3_vec, i, zeta);
                        gsl_vector_set(r3Uniq, u_, gsl_vector_max(r3Uniq)+1);
                        gsl_vector_set(mu3_all, u_, mu);
                        gsl_vector_set(zeta3_all, u_, zeta);
                        u += 1;
                    }
                }
            }else /* singletion */
            {
                gsl_vector_memcpy(mu_temp, mu3_all);
                gsl_vector_memcpy(zeta_temp, zeta3_all);
                
                gsl_vector_set(mu_temp, wh_ri, mu);
                gsl_vector_set(zeta_temp, wh_ri, zeta);
                
                for(j = 0; j < u_+1; j++)
                {
                    mean = gsl_vector_get(mu_temp, j) + eta;
                    sd = sqrt(1/gsl_vector_get(zeta_temp, j));
                    
                    n_ir = gsl_vector_get(r3Uniq_count, j);
                    if(gsl_vector_get(r3, i) == gsl_vector_get(r3Uniq, j)) n_ir -= 1;
                    
                    val = dnorm(yStar, mean, sd, 0);
                    
                    if(j == wh_ri)
                    {
                        val *= (double) tau3 / (double)(n-1+tau3);
                    }else
                    {
                        val *= (double) n_ir / (double)(n-1+tau3);
                    }
                    gsl_vector_set(prob, j, val);
                }
                
                sum_prob = 0;
                for(k = 0; k < u_+1; k++)
                {
                    sum_prob += gsl_vector_get(prob, k);
                }
                
                if(sum_prob > pow(10, -300))
                {
                    b_mc = 1/sum_prob;
                    gsl_vector_scale(prob, b_mc);
                    
                    r_ind = c_multinom_sample(rr, prob, u_+1);
                    gsl_vector_set(r3, i, gsl_vector_get(r3Uniq, r_ind-1));
                    
                    if(r_ind-1 != wh_ri)
                    {
                        if(wh_ri < u_)
                        {
                            gsl_vector_memcpy(rUniq_dummy, r3Uniq);
                            gsl_vector_memcpy(mu_dummy, mu3_all);
                            gsl_vector_memcpy(zeta_dummy, zeta3_all);
                            for(k=wh_ri;k<u_;k++)
                            {
                                gsl_vector_set(r3Uniq, k, gsl_vector_get(rUniq_dummy, k+1));
                                gsl_vector_set(mu_temp, k, gsl_vector_get(mu_dummy, k+1));
                                gsl_vector_set(zeta_temp, k, gsl_vector_get(zeta_dummy, k+1));
                            }
                        }
                        gsl_vector_set(mu3_vec, i, gsl_vector_get(mu_temp, r_ind-1));
                        gsl_vector_set(zeta3_vec, i, gsl_vector_get(zeta_temp, r_ind-1));
                        gsl_vector_set(r3Uniq, u_, 0);
                        gsl_vector_set(mu_temp, u_, 0);
                        gsl_vector_set(zeta_temp, u_, 0);
                        u -= 1;
                    }
                    gsl_vector_memcpy(mu3_all, mu_temp);
                    gsl_vector_memcpy(zeta3_all, zeta_temp);
                }
            }
            
            gsl_vector_free(prob);
            
            for(j = u_+1; j < n; j++)
            {
                gsl_vector_set(r3Uniq, j, 0);
                gsl_vector_set(mu3_all, j, 0);
                gsl_vector_set(zeta3_all, j, 0);
            }
        }else
        {
            gsl_vector_set(r3, i, 0);
        }
    }
    
    /*************************************************************/
    
    /* Step 2: update (beta0, sigSq) using the posterior distribution that is based on {epsilon_i :i∈{k:rk =r}}. */
    
    /*************************************************************/
    
    /* identfy "u" unique values */
    c_uniq_h3(r3, r3Uniq, r3Uniq_count, mu3_all, zeta3_all, mu3_vec, zeta3_vec, y1_NA, &u);
    *nClass_DP3 = u;
    
    double mu_prop, zeta_prop, mu_ini, zeta_ini, temp, temp_prop;
    double loglh, loglh_prop, logprior, logprior_prop, logprop, logprop_prop;
    double logR;
    int uu;
    
    for(j = 0; j < u; j++)
    {
        mu_ini = gsl_vector_get(mu3_all, j);
        zeta_ini = gsl_vector_get(zeta3_all, j);
        
        /* update mu */
        mu_prop = rnorm((double) mu_ini, (double) sqrt(beta03_prop_var));
        
        loglh=0;
        loglh_prop=0;
        for(i = 0; i < n; i++)
        {
            if((gsl_vector_get(y1, i) < gsl_vector_get(y2, i)) && (gsl_vector_get(y1_NA, i) == 0))
            {
                if(gsl_vector_get(r3, i) == gsl_vector_get(r3Uniq, j))
                {
                    yStar = gsl_vector_get(y2, i) + log(1-exp(gsl_vector_get(y1, i)-gsl_vector_get(y2, i)));
                    eta = gsl_vector_get(xbeta3, i) + gsl_vector_get(gamma, i);
                    mean = mu_ini + eta;
                    sd = sqrt(1/zeta_ini);
                    mean_prop = mu_prop + eta;
                    
                    temp = dnorm(yStar, mean, sd, 1);
                    temp_prop = dnorm(yStar, mean_prop, sd, 1);
                    
                    loglh += temp;
                    loglh_prop += temp_prop;
                }
            }
        }
        logprior = dnorm(mu_ini, mu03, sqrt(sigSq03), 1);
        logprior_prop = dnorm(mu_prop, mu03, sqrt(sigSq03), 1);
        
        logprop = dnorm(mu_ini, mu_prop, sqrt(beta03_prop_var), 1);
        logprop_prop = dnorm(mu_prop, mu_ini, sqrt(beta03_prop_var), 1);
        
        logR = loglh_prop - loglh + logprior_prop - logprior + logprop - logprop_prop;
        uu = log(runif(0, 1)) < logR;
        
        if(uu == 1)
        {
            for(i = 0; i < n; i++)
            {
                if(gsl_vector_get(r3, i) == gsl_vector_get(r3Uniq, j))
                {
                    gsl_vector_set(mu3_vec, i, mu_prop);
                }
            }
            gsl_vector_set(mu3_all, j, mu_prop);
            *accept_beta03 += 1;
        }
        
        /* update zeta */
        mu_ini = gsl_vector_get(mu3_all, j);
        
        zeta_prop = INFINITY;
        while(zeta_prop == INFINITY || isnan(zeta_prop))
        {
            zeta_prop = rgamma(pow(zeta_ini, 2)/zeta3_prop_var, zeta3_prop_var/zeta_ini);
        }
        
        
        loglh=0;
        loglh_prop=0;
        for(i = 0; i < n; i++)
        {
            if((gsl_vector_get(y1, i) < gsl_vector_get(y2, i)) && (gsl_vector_get(y1_NA, i) == 0))
            {
                if(gsl_vector_get(r3, i) == gsl_vector_get(r3Uniq, j))
                {
                    yStar = gsl_vector_get(y2, i) + log(1-exp(gsl_vector_get(y1, i)-gsl_vector_get(y2, i)));
                    eta = gsl_vector_get(xbeta3, i) + gsl_vector_get(gamma, i);
                    mean = mu_ini + eta;
                    sd = sqrt(1/zeta_ini);
                    sd_prop = sqrt(1/zeta_prop);
                    
                    temp = dnorm(yStar, mean, sd, 1);
                    temp_prop = dnorm(yStar, mean, sd_prop, 1);
                    
                    loglh += temp;
                    loglh_prop += temp_prop;
                }
            }
            
        }
        logprior = dgamma(zeta_ini, a03, 1/b03, 1);
        logprior_prop = dgamma(zeta_prop, a03, 1/b03, 1);
        
        logprop = dgamma(zeta_ini, pow(zeta_prop, 2)/zeta3_prop_var, zeta3_prop_var/zeta_prop, 1);
        logprop_prop = dgamma(zeta_prop, pow(zeta_ini, 2)/zeta3_prop_var, zeta3_prop_var/zeta_ini, 1);
        
        logR = loglh_prop - loglh + logprior_prop - logprior + logprop - logprop_prop;
        uu = log(runif(0, 1)) < logR;
        
        if(uu == 1)
        {
            for(i = 0; i < n; i++)
            {
                if(gsl_vector_get(r3, i) == gsl_vector_get(r3Uniq, j))
                {
                    gsl_vector_set(zeta3_vec, i, zeta_prop);
                }
            }
            gsl_vector_set(zeta3_all, j, zeta_prop);
            *accept_zeta3 += 1;
        }
    }
    
    gsl_vector_free(xbeta3);
    gsl_vector_free(mu_temp);
    gsl_vector_free(zeta_temp);
    gsl_vector_free(mu_dummy);
    gsl_vector_free(zeta_dummy);
    gsl_vector_free(rUniq_dummy);
    
    return;
    
}










/* Updating gamma */

void BAFT_DPscr_update_gamma(gsl_vector *y1_NA,
                             gsl_vector *c0,
                             gsl_vector *c0_neginf,
                             gsl_matrix *X1,
                             gsl_matrix *X2,
                             gsl_matrix *X3,
                             gsl_vector *y1,
                             gsl_vector *y2,
                             gsl_vector *beta1,
                             gsl_vector *beta2,
                             gsl_vector *beta3,
                             gsl_vector *gamma,
                             gsl_vector *r1,
                             gsl_vector *r2,
                             gsl_vector *r3,
                             gsl_vector *mu1_all,
                             gsl_vector *mu2_all,
                             gsl_vector *mu3_all,
                             gsl_vector *zeta1_all,
                             gsl_vector *zeta2_all,
                             gsl_vector *zeta3_all,
                             gsl_vector *r1Uniq,
                             gsl_vector *r2Uniq,
                             gsl_vector *r3Uniq,
                             gsl_vector *r1Uniq_count,
                             gsl_vector *r2Uniq_count,
                             gsl_vector *r3Uniq_count,
                             int *nClass_DP1,
                             int *nClass_DP2,
                             int *nClass_DP3,
                             double theta,
                             double gamma_prop_var,
                             gsl_vector *accept_gamma)
{
    int i, uu, k, ii1, ii2, ii3;
    double gamma_prop, gamma_ini;
    double eta1, eta1_prop, eta2, eta2_prop, eta3, eta3_prop;
    double loglh, loglh_prop, logR;
    double mean1, mean1_prop, sd1;
    double mean2, mean2_prop, sd2;
    double mean3, mean3_prop, sd3;
    double logprior, logprior_prop;
    double yStar;
    
    int n = X1 -> size1;
    int u1 = *nClass_DP1;
    int u2 = *nClass_DP2;
    int u3 = *nClass_DP3;
    
    gsl_vector *xbeta1 = gsl_vector_calloc(n);
    gsl_vector *xbeta2 = gsl_vector_calloc(n);
    gsl_vector *xbeta3 = gsl_vector_calloc(n);
    
    gsl_blas_dgemv(CblasNoTrans, 1, X1, beta1, 0, xbeta1);
    gsl_blas_dgemv(CblasNoTrans, 1, X2, beta2, 0, xbeta2);
    gsl_blas_dgemv(CblasNoTrans, 1, X3, beta3, 0, xbeta3);
    
    
    for(i=0;i<n;i++)
    {
        gamma_ini = gsl_vector_get(gamma, i);
        gamma_prop = rnorm(gamma_ini, sqrt(gamma_prop_var));
        
        eta1 = gsl_vector_get(xbeta1, i) + gamma_ini;
        eta1_prop = gsl_vector_get(xbeta1, i) + gamma_prop;
        
        eta2 = gsl_vector_get(xbeta2, i) + gamma_ini;
        eta2_prop = gsl_vector_get(xbeta2, i) + gamma_prop;
        
        for(k = 0; k < u1; k++)
        {
            if(gsl_vector_get(r1, i) == gsl_vector_get(r1Uniq, k)) ii1 = k;
        }
        
        for(k = 0; k < u2; k++)
        {
            if(gsl_vector_get(r2, i) == gsl_vector_get(r2Uniq, k)) ii2 = k;
        }
        
        mean1 = eta1 + gsl_vector_get(mu1_all, ii1);
        mean1_prop = eta1_prop + gsl_vector_get(mu1_all, ii1);
        sd1 = (double) pow(gsl_vector_get(zeta1_all, ii1), -0.5);
        
        mean2 = eta2 + gsl_vector_get(mu2_all, ii2);
        mean2_prop = eta2_prop + gsl_vector_get(mu2_all, ii2);
        sd2 = (double) pow(gsl_vector_get(zeta2_all, ii2), -0.5);
        
        if(gsl_vector_get(y1_NA, i) == 0)
        {
            if(gsl_vector_get(y1, i) < gsl_vector_get(y2, i))
            {
                eta3 = gsl_vector_get(xbeta3, i) + gamma_ini;
                eta3_prop = gsl_vector_get(xbeta3, i) + gamma_prop;
                
                for(k = 0; k < u3; k++)
                {
                    if(gsl_vector_get(r3, i) == gsl_vector_get(r3Uniq, k)) ii3 = k;
                }
                mean3 = eta3 + gsl_vector_get(mu3_all, ii3);
                mean3_prop = eta3_prop + gsl_vector_get(mu3_all, ii3);
                sd3 = (double) pow(gsl_vector_get(zeta3_all, ii3), -0.5);
            }
        }
        
        loglh = 0;
        loglh_prop = 0;
        if(gsl_vector_get(y1_NA, i) == 0)
        {
            loglh += dnorm(gsl_vector_get(y1, i), mean1, sd1, 1);
            loglh_prop += dnorm(gsl_vector_get(y1, i), mean1_prop, sd1, 1);
            
            loglh += pnorm(gsl_vector_get(y1, i), mean2, sd2, 0, 1);
            loglh_prop += pnorm(gsl_vector_get(y1, i), mean2_prop, sd2, 0, 1);
            
            if(gsl_vector_get(y1, i) < gsl_vector_get(y2, i))
            {
                yStar = gsl_vector_get(y2, i) + log(1-exp(gsl_vector_get(y1, i)-gsl_vector_get(y2, i)));
                
                loglh += dnorm(yStar, mean3, sd3, 1);
                loglh_prop += dnorm(yStar, mean3_prop, sd3, 1);
            }
            
        }else
        {
            loglh += pnorm(gsl_vector_get(y2, i), mean1, sd1, 0, 1);
            loglh_prop += pnorm(gsl_vector_get(y2, i), mean1_prop, sd1, 0, 1);
            
            loglh += dnorm(gsl_vector_get(y2, i), mean2, sd2, 1);
            loglh_prop += dnorm(gsl_vector_get(y2, i), mean2_prop, sd2, 1);
        }
        
        if(gsl_vector_get(c0_neginf, i) == 0)
        {
            loglh -= pnorm(gsl_vector_get(c0, i), mean1, sd1, 0, 1);
            loglh_prop -= pnorm(gsl_vector_get(c0, i), mean1_prop, sd1, 0, 1);
            
            loglh -= pnorm(gsl_vector_get(c0, i), mean2, sd2, 0, 1);
            loglh_prop -= pnorm(gsl_vector_get(c0, i), mean2_prop, sd2, 0, 1);
        }
        
        logprior = dnorm(gamma_ini, 0, sqrt(theta), 1);
        logprior_prop = dnorm(gamma_prop, 0, sqrt(theta), 1);
        
        logR = loglh_prop - loglh + logprior_prop - logprior;
        uu = log(runif(0, 1)) < logR;
        
        if(uu == 1)
        {
            gsl_vector_set(gamma, i, gamma_prop);
            gsl_vector_set(accept_gamma, i, gsl_vector_get(accept_gamma, i) + 1);
        }
    }
    
    gsl_vector_free(xbeta1);
    gsl_vector_free(xbeta2);
    gsl_vector_free(xbeta3);
    
    return;
}






/* Updating theta */

void BAFT_DPscr_update_theta(gsl_vector *gamma,
                             double *theta,
                             double a_theta,
                             double b_theta)
{
    int i;
    int n = gamma -> size;
    double a = a_theta + n/2;
    double b = b_theta;
    
    for(i=0;i<n;i++)
    {
        b += 0.5*pow(gsl_vector_get(gamma,i),2);
    }
    
    c_rigamma(theta, a, b);
}





/* updating precision parameter of DP prior: tau */

void BAFT_DPscr_update_tau(int *n,
                           double *tau,
                           double aTau,
                           double bTau,
                           int *nClass_DP)
{
    double eta, pEta, tau_shape, tau_rate, tau_scale, dist1Ind;
    
    eta = rbeta(*tau+1, *n);
    
    pEta = (aTau + (double) *nClass_DP - 1)/(*n * bTau - *n * log(eta) + aTau + (double) *nClass_DP - 1);
    
    dist1Ind = rbinom(1, pEta);
    
    tau_shape = aTau + (double) *nClass_DP - 1 + dist1Ind;
    tau_rate = bTau - log(eta);
    tau_scale = 1/tau_rate;
    
    *tau = rgamma(tau_shape, tau_scale);
    
    return;
}








