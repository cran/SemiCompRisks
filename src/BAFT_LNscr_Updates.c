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

#include "BAFT_LNscr.h"



/* Updating y1 and y2 */

void BAFT_LNscr_update_y(gsl_vector *y1L,
                         gsl_vector *y1U,
                         gsl_vector *y2L,
                         gsl_vector *y2U,
                         gsl_vector *y1L_neginf,
                         gsl_vector *y2L_neginf,
                         gsl_vector *y1U_posinf,
                         gsl_vector *y2U_posinf,
                         gsl_vector *y1_NA,
                         gsl_matrix *X1,
                         gsl_matrix *X2,
                         gsl_matrix *X3,
                         gsl_vector *y1,
                         gsl_vector *y2,
                         gsl_vector *beta1,
                         gsl_vector *beta2,
                         gsl_vector *beta3,
                         gsl_vector *gamma,
                         double beta01,
                         double beta02,
                         double beta03,
                         double sigSq1,
                         double sigSq2,
                         double sigSq3)
{
    double eta1, eta2, eta3, y1_temp, y2_temp, gap, CenNew, CenNew1, CenNew2;
    int i;
    int n = y1 -> size;
    
    gsl_vector *xbeta1 = gsl_vector_calloc(n);
    gsl_vector *xbeta2 = gsl_vector_calloc(n);
    gsl_vector *xbeta3 = gsl_vector_calloc(n);
    gsl_blas_dgemv(CblasNoTrans, 1, X1, beta1, 0, xbeta1);
    gsl_blas_dgemv(CblasNoTrans, 1, X2, beta2, 0, xbeta2);
    gsl_blas_dgemv(CblasNoTrans, 1, X3, beta3, 0, xbeta3);
    
    for(i=0;i<n;i++)
    {
        gap = 0;
        
        eta1 = beta01+gsl_vector_get(xbeta1, i)+gsl_vector_get(gamma, i);
        eta2 = beta02+gsl_vector_get(xbeta2, i)+gsl_vector_get(gamma, i);
        eta3 = beta03+gsl_vector_get(xbeta3, i)+gsl_vector_get(gamma, i);
        
        /*  both end points are right-censored */
        
        if(gsl_vector_get(y1U_posinf, i) == 1 && gsl_vector_get(y2U_posinf, i) == 1)
        {
            c_rtnorm(eta1, sqrt(sigSq1), gsl_vector_get(y1L, i), gsl_vector_get(y1U, i), gsl_vector_get(y1L_neginf, i), 1, &y1_temp);
            c_rtnorm(eta2, sqrt(sigSq2), gsl_vector_get(y2L, i), gsl_vector_get(y2U, i), gsl_vector_get(y2L_neginf, i), 1, &y2_temp);
            if(y1_temp < y2_temp)
            {
                gap = rnorm(eta3, sqrt(sigSq3));
                y2_temp = log(exp(gap) + exp(y1_temp));
                gsl_vector_set(y1_NA, i, 0);
            }else
            {
                y1_temp = y2_temp;
                gsl_vector_set(y1_NA, i, 1);
            }
        }
        
        /*  non-terminal event is right-censored and terminal event is observed (or interval-censored) */
        
        if(gsl_vector_get(y1U_posinf, i) == 1 && gsl_vector_get(y2U_posinf, i) == 0)
        {
            c_rtnorm(eta1, sqrt(sigSq1), gsl_vector_get(y1L, i), gsl_vector_get(y1U, i), gsl_vector_get(y1L_neginf, i), 1, &y1_temp);
            
            /*  i) terminal event is interval-censored */
            if(gsl_vector_get(y2L, i) < gsl_vector_get(y2U, i))
            {
                c_rtnorm(eta2, sqrt(sigSq2), gsl_vector_get(y2L, i), gsl_vector_get(y2U, i), gsl_vector_get(y2L_neginf, i), 0, &y2_temp);
                
                if(y1_temp < y2_temp)
                {
                    CenNew = log(exp(gsl_vector_get(y2U, i))-exp(y1_temp));
                    c_rtnorm(eta3, sqrt(sigSq3), -10000000000, CenNew, 1, 0, &gap);
                    y2_temp = log(exp(gap) + exp(y1_temp));
                    gsl_vector_set(y1_NA, i, 0);
                    
                }else
                {
                    y1_temp = y2_temp;
                    gsl_vector_set(y1_NA, i, 1);
                }/*  ii) terminal event is observed */
            }else if(gsl_vector_get(y2L, i) == gsl_vector_get(y2U, i))
            {
                y2_temp = gsl_vector_get(y2L, i);
                
                if(y1_temp < y2_temp)
                {
                    gsl_vector_set(y1_NA, i, 0);
                    
                }else
                {
                    y1_temp = y2_temp;
                    gsl_vector_set(y1_NA, i, 1);
                }
            }
        }
        
        /*  non-terminal event is observed (or interval-censored) and terminal event is right-censored */
        
        if(gsl_vector_get(y1U_posinf, i) == 0 && gsl_vector_get(y2U_posinf, i) == 1)
        {
            /*  i) non-terminal event is interval-censored */
            if(gsl_vector_get(y1L, i) < gsl_vector_get(y1U, i))
            {
                c_rtnorm(eta1, sqrt(sigSq1), gsl_vector_get(y1L, i), gsl_vector_get(y1U, i), 0, gsl_vector_get(y1L_neginf, i), &y1_temp);
                /*  ii) non-terminal event is observed */
            }else if(gsl_vector_get(y1L, i) == gsl_vector_get(y1U, i))
            {
                y1_temp = gsl_vector_get(y1L, i);
            }
            if(y1_temp < gsl_vector_get(y2L, i))
            {
                CenNew = log(exp(gsl_vector_get(y2L, i))-exp(y1_temp));
                c_rtnorm(eta3, sqrt(sigSq3), CenNew, 10000000000, 0, 1, &gap);
            }else
            {
                gap = rnorm(eta3, sqrt(sigSq3));
            }
            
            y2_temp = log(exp(gap) + exp(y1_temp));
            gsl_vector_set(y1_NA, i, 0);
        }
        
        /*  both events are observed (or interval-censored)*/
        if(gsl_vector_get(y1U_posinf, i) == 0 && gsl_vector_get(y2U_posinf, i) == 0)
        {
            /*  non-terminal event is interval-censored */
            if(gsl_vector_get(y1L, i) < gsl_vector_get(y1U, i))
            {
                c_rtnorm(eta1, sqrt(sigSq1), gsl_vector_get(y1L, i), gsl_vector_get(y1U, i), gsl_vector_get(y1L_neginf, i), 0, &y1_temp);
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
                    CenNew1 = log(exp(gsl_vector_get(y2L, i))-exp(y1_temp));
                    CenNew2 = log(exp(gsl_vector_get(y2U, i))-exp(y1_temp));
                    c_rtnorm(eta3, sqrt(sigSq3), CenNew1, CenNew2, 0, 0, &gap);
                }else if(y1_temp < gsl_vector_get(y2U, i) && y1_temp > gsl_vector_get(y2L, i))
                {
                    CenNew2 = log(exp(gsl_vector_get(y2U, i))-exp(y1_temp));
                    c_rtnorm(eta3, sqrt(sigSq3), -1000000000, CenNew2, 1, 0, &gap);
                }
                y2_temp = log(exp(gap) + exp(y1_temp));
                
                /*  terminal event is observed */
            }else if(gsl_vector_get(y2L, i) == gsl_vector_get(y2U, i))
            {
                y2_temp = gsl_vector_get(y2L, i);
                gap = log(exp(y2_temp) - exp(y1_temp));
            }
            gsl_vector_set(y1_NA, i, 0);
        }
        gsl_vector_set(y1, i, y1_temp);
        gsl_vector_set(y2, i, y2_temp);
    }
    
    gsl_vector_free(xbeta1);
    gsl_vector_free(xbeta2);
    gsl_vector_free(xbeta3);
    return;
}




/* Updating beta1 */

void BAFT_LNscr_update_beta1(gsl_vector *y1_NA,
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
                             double beta01,
                             double beta02,
                             double beta03,
                             double sigSq1,
                             double sigSq2,
                             double sigSq3,
                             double beta1_prop_var,
                             gsl_vector *accept_beta1)
{
    int i, j, u;
    double loglh, loglh_prop, logR, temp;
    
    int n = X1 -> size1;
    int p1 = X1 -> size2;
    
    gsl_vector *beta1_prop = gsl_vector_calloc(p1);
    
    for(j=0;j<p1;j++)
    {
        loglh = 0;
        loglh_prop = 0;
        
        gsl_vector_memcpy(beta1_prop, beta1);
        gsl_vector_set(beta1_prop, j, rnorm(gsl_vector_get(beta1, j), sqrt(beta1_prop_var)));
        
        for(i=0;i<n;i++)
        {
            if(gsl_vector_get(y1_NA, i) == 0)
            {
                log_Jpdf_Upper_BAFT_LN(i, gsl_vector_get(y1, i), gsl_vector_get(y2, i), gsl_vector_get(c0, i), c0_neginf, X1, X2, X3, beta1, beta2, beta3, gamma, beta01, beta02, beta03, sigSq1, sigSq2, sigSq3, &temp);
                loglh += temp;
                
                log_Jpdf_Upper_BAFT_LN(i, gsl_vector_get(y1, i), gsl_vector_get(y2, i), gsl_vector_get(c0, i), c0_neginf, X1, X2, X3, beta1_prop, beta2, beta3, gamma, beta01, beta02, beta03, sigSq1, sigSq2, sigSq3, &temp);
                loglh_prop += temp;
            }else
            {
                log_Jpdf_Lower_BAFT_LN(i, gsl_vector_get(y2, i), gsl_vector_get(c0, i), c0_neginf, X1, X2, beta1, beta2, gamma, beta01, beta02, sigSq1, sigSq2, &temp);
                loglh += temp;
                log_Jpdf_Lower_BAFT_LN(i, gsl_vector_get(y2, i), gsl_vector_get(c0, i), c0_neginf, X1, X2, beta1_prop, beta2, gamma, beta01, beta02, sigSq1, sigSq2, &temp);
                loglh_prop += temp;
            }
            
        }
        
        logR = loglh_prop - loglh;
        u = log(runif(0, 1)) < logR;
        
        if(u == 1)
        {
            gsl_vector_memcpy(beta1, beta1_prop);
            gsl_vector_set(accept_beta1, j, gsl_vector_get(accept_beta1, j) + 1);
        }
    }
    
    gsl_vector_free(beta1_prop);
    return;
}







/* Updating beta2 */

void BAFT_LNscr_update_beta2(gsl_vector *y1_NA,
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
                             double beta01,
                             double beta02,
                             double beta03,
                             double sigSq1,
                             double sigSq2,
                             double sigSq3,
                             double beta2_prop_var,
                             gsl_vector *accept_beta2)
{
    int i, j, u;
    double loglh, loglh_prop, logR, temp;
    
    int n = X2 -> size1;
    int p2 = X2 -> size2;
    
    gsl_vector *beta2_prop = gsl_vector_calloc(p2);
    
    for(j=0;j<p2;j++)
    {
        loglh = 0;
        loglh_prop = 0;
        
        gsl_vector_memcpy(beta2_prop, beta2);
        gsl_vector_set(beta2_prop, j, rnorm(gsl_vector_get(beta2, j), sqrt(beta2_prop_var)));
        
        for(i=0;i<n;i++)
        {
            if(gsl_vector_get(y1_NA, i) == 0)
            {
                log_Jpdf_Upper_BAFT_LN(i, gsl_vector_get(y1, i), gsl_vector_get(y2, i), gsl_vector_get(c0, i), c0_neginf, X1, X2, X3, beta1, beta2, beta3, gamma, beta01, beta02, beta03, sigSq1, sigSq2, sigSq3, &temp);
                loglh += temp;
                
                log_Jpdf_Upper_BAFT_LN(i, gsl_vector_get(y1, i), gsl_vector_get(y2, i), gsl_vector_get(c0, i), c0_neginf, X1, X2, X3, beta1, beta2_prop, beta3, gamma, beta01, beta02, beta03, sigSq1, sigSq2, sigSq3, &temp);
                loglh_prop += temp;
            }else
            {
                log_Jpdf_Lower_BAFT_LN(i, gsl_vector_get(y2, i), gsl_vector_get(c0, i), c0_neginf, X1, X2, beta1, beta2, gamma, beta01, beta02, sigSq1, sigSq2, &temp);
                loglh += temp;
                
                log_Jpdf_Lower_BAFT_LN(i, gsl_vector_get(y2, i), gsl_vector_get(c0, i), c0_neginf, X1, X2, beta1, beta2_prop, gamma, beta01, beta02, sigSq1, sigSq2, &temp);
                loglh_prop += temp;
            }
        }
        
        logR = loglh_prop - loglh;
        u = log(runif(0, 1)) < logR;
        
        if(u == 1)
        {
            gsl_vector_memcpy(beta2, beta2_prop);
            gsl_vector_set(accept_beta2, j, gsl_vector_get(accept_beta2, j) + 1);
        }
    }
    
    gsl_vector_free(beta2_prop);
    return;
}






/* Updating beta3 */

void BAFT_LNscr_update_beta3(gsl_vector *y1_NA,
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
                             double beta01,
                             double beta02,
                             double beta03,
                             double sigSq1,
                             double sigSq2,
                             double sigSq3,
                             double beta3_prop_var,
                             gsl_vector *accept_beta3)
{
    int i, j, u;
    double loglh, loglh_prop, logR, temp;
    
    int n = X3 -> size1;
    int p3 = X3 -> size2;
    
    gsl_vector *beta3_prop = gsl_vector_calloc(p3);
    
    for(j=0;j<p3;j++)
    {
        loglh = 0;
        loglh_prop = 0;
        
        gsl_vector_memcpy(beta3_prop, beta3);
        gsl_vector_set(beta3_prop, j, rnorm(gsl_vector_get(beta3, j), sqrt(beta3_prop_var)));
        
        for(i=0;i<n;i++)
        {
            if(gsl_vector_get(y1_NA, i) == 0)
            {
                log_Jpdf_Upper_BAFT_LN(i, gsl_vector_get(y1, i), gsl_vector_get(y2, i), gsl_vector_get(c0, i), c0_neginf, X1, X2, X3, beta1, beta2, beta3, gamma, beta01, beta02, beta03, sigSq1, sigSq2, sigSq3, &temp);
                loglh += temp;
                
                log_Jpdf_Upper_BAFT_LN(i, gsl_vector_get(y1, i), gsl_vector_get(y2, i), gsl_vector_get(c0, i), c0_neginf, X1, X2, X3, beta1, beta2, beta3_prop, gamma, beta01, beta02, beta03, sigSq1, sigSq2, sigSq3, &temp);
                loglh_prop += temp;
            }
        }
        
        logR = loglh_prop - loglh;
        u = log(runif(0, 1)) < logR;
        
        if(u == 1)
        {
            gsl_vector_memcpy(beta3, beta3_prop);
            gsl_vector_set(accept_beta3, j, gsl_vector_get(accept_beta3, j) + 1);
        }
    }
    
    gsl_vector_free(beta3_prop);
    return;
}




/* Updating beta01 */

void BAFT_LNscr_update_beta01(gsl_vector *y1_NA,
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
                              double *beta01,
                              double beta02,
                              double beta03,
                              double sigSq1,
                              double sigSq2,
                              double sigSq3,
                              double beta01_prop_var,
                              int *accept_beta01)
{
    int i, u;
    double loglh, loglh_prop, logR, temp, logprior, logprior_prop;
    
    int n = X1 -> size1;
    
    loglh = 0;
    loglh_prop = 0;
    
    double beta01_prop = rnorm(*beta01, sqrt(beta01_prop_var));
    
    for(i=0;i<n;i++)
    {
        if(gsl_vector_get(y1_NA, i) == 0)
        {
            log_Jpdf_Upper_BAFT_LN(i, gsl_vector_get(y1, i), gsl_vector_get(y2, i), gsl_vector_get(c0, i), c0_neginf, X1, X2, X3, beta1, beta2, beta3, gamma, *beta01, beta02, beta03, sigSq1, sigSq2, sigSq3, &temp);
            loglh += temp;
            
            log_Jpdf_Upper_BAFT_LN(i, gsl_vector_get(y1, i), gsl_vector_get(y2, i), gsl_vector_get(c0, i), c0_neginf, X1, X2, X3, beta1, beta2, beta3, gamma, beta01_prop, beta02, beta03, sigSq1, sigSq2, sigSq3, &temp);
            loglh_prop += temp;
        }else
        {
            log_Jpdf_Lower_BAFT_LN(i, gsl_vector_get(y2, i), gsl_vector_get(c0, i), c0_neginf, X1, X2, beta1, beta2, gamma, *beta01, beta02, sigSq1, sigSq2, &temp);
            loglh += temp;
            
            log_Jpdf_Lower_BAFT_LN(i, gsl_vector_get(y2, i), gsl_vector_get(c0, i), c0_neginf, X1, X2, beta1, beta2, gamma, beta01_prop, beta02, sigSq1, sigSq2, &temp);
            loglh_prop += temp;
        }
        
    }
    
    logprior = dnorm(*beta01, 0, pow(10,6)*sqrt(sigSq1), 1);
    logprior_prop = dnorm(beta01_prop, 0, pow(10,6)*sqrt(sigSq1), 1);;
    logR = loglh_prop - loglh + logprior_prop - logprior;
    
    u = log(runif(0, 1)) < logR;
    
    if(u == 1)
    {
        *beta01 = beta01_prop;
        *accept_beta01 += 1;
    }
    
    return;
}






/* Updating beta02 */

void BAFT_LNscr_update_beta02(gsl_vector *y1_NA,
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
                              double beta01,
                              double *beta02,
                              double beta03,
                              double sigSq1,
                              double sigSq2,
                              double sigSq3,
                              double beta02_prop_var,
                              int *accept_beta02)
{
    int i, u;
    double loglh, loglh_prop, logR, temp, logprior, logprior_prop;
    
    int n = X1 -> size1;
    
    loglh = 0;
    loglh_prop = 0;
    
    double beta02_prop = rnorm(*beta02, sqrt(beta02_prop_var));
    
    for(i=0;i<n;i++)
    {
        if(gsl_vector_get(y1_NA, i) == 0)
        {
            log_Jpdf_Upper_BAFT_LN(i, gsl_vector_get(y1, i), gsl_vector_get(y2, i), gsl_vector_get(c0, i), c0_neginf, X1, X2, X3, beta1, beta2, beta3, gamma, beta01, *beta02, beta03, sigSq1, sigSq2, sigSq3, &temp);
            loglh += temp;
            
            log_Jpdf_Upper_BAFT_LN(i, gsl_vector_get(y1, i), gsl_vector_get(y2, i), gsl_vector_get(c0, i), c0_neginf, X1, X2, X3, beta1, beta2, beta3, gamma, beta01, beta02_prop, beta03, sigSq1, sigSq2, sigSq3, &temp);
            loglh_prop += temp;
        }else
        {
            log_Jpdf_Lower_BAFT_LN(i, gsl_vector_get(y2, i), gsl_vector_get(c0, i), c0_neginf, X1, X2, beta1, beta2, gamma, beta01, *beta02, sigSq1, sigSq2, &temp);
            loglh += temp;
            
            log_Jpdf_Lower_BAFT_LN(i, gsl_vector_get(y2, i), gsl_vector_get(c0, i), c0_neginf, X1, X2, beta1, beta2, gamma, beta01, beta02_prop, sigSq1, sigSq2, &temp);
            loglh_prop += temp;
        }
    }
    
    logprior = dnorm(*beta02, 0, pow(10,6)*sqrt(sigSq2), 1);
    logprior_prop = dnorm(beta02_prop, 0, pow(10,6)*sqrt(sigSq2), 1);;
    logR = loglh_prop - loglh + logprior_prop - logprior;
    
    u = log(runif(0, 1)) < logR;
    
    if(u == 1)
    {
        *beta02 = beta02_prop;
        *accept_beta02 += 1;
    }
    
    return;
}




/* Updating beta03 */

void BAFT_LNscr_update_beta03(gsl_vector *y1_NA,
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
                              double beta01,
                              double beta02,
                              double *beta03,
                              double sigSq1,
                              double sigSq2,
                              double sigSq3,
                              double beta03_prop_var,
                              int *accept_beta03)
{
    int i, u;
    double loglh, loglh_prop, logR, temp, logprior, logprior_prop;
    
    int n = X1 -> size1;
    
    loglh = 0;
    loglh_prop = 0;
    
    double beta03_prop = rnorm(*beta03, sqrt(beta03_prop_var));
    
    for(i=0;i<n;i++)
    {
        if(gsl_vector_get(y1_NA, i) == 0)
        {
            log_Jpdf_Upper_BAFT_LN(i, gsl_vector_get(y1, i), gsl_vector_get(y2, i), gsl_vector_get(c0, i), c0_neginf, X1, X2, X3, beta1, beta2, beta3, gamma, beta01, beta02, *beta03, sigSq1, sigSq2, sigSq3, &temp);
            loglh += temp;
            
            log_Jpdf_Upper_BAFT_LN(i, gsl_vector_get(y1, i), gsl_vector_get(y2, i), gsl_vector_get(c0, i), c0_neginf, X1, X2, X3, beta1, beta2, beta3, gamma, beta01, beta02, beta03_prop, sigSq1, sigSq2, sigSq3, &temp);
            loglh_prop += temp;
        }
    }
    
    logprior = dnorm(*beta03, 0, pow(10,6)*sqrt(sigSq3), 1);
    logprior_prop = dnorm(beta03_prop, 0, pow(10,6)*sqrt(sigSq3), 1);;
    logR = loglh_prop - loglh + logprior_prop - logprior;
    
    u = log(runif(0, 1)) < logR;
    
    if(u == 1)
    {
        *beta03 = beta03_prop;
        *accept_beta03 += 1;
    }
    
    return;
}





/* Updating sigSq1 */

void BAFT_LNscr_update_sigSq1(gsl_vector *y1_NA,
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
                              double beta01,
                              double beta02,
                              double beta03,
                              double *sigSq1,
                              double sigSq2,
                              double sigSq3,
                              double a_sigSq1,
                              double b_sigSq1,
                              double sigSq1_prop_var,
                              int *accept_sigSq1)
{
    int i, u;
    double loglh, loglh_prop, logR, temp, xeta1_prop, sigSq1_prop;
    double logprior, logprior_prop, logprior1, logprior1_prop;
    int n = X1 -> size1;
    
    loglh = 0;
    loglh_prop = 0;
    
    xeta1_prop = rnorm(log(*sigSq1), sqrt(sigSq1_prop_var));
    sigSq1_prop = exp(xeta1_prop);
    
    
    for(i=0;i<n;i++)
    {
        if(gsl_vector_get(y1_NA, i) == 0)
        {
            log_Jpdf_Upper_BAFT_LN(i, gsl_vector_get(y1, i), gsl_vector_get(y2, i), gsl_vector_get(c0, i), c0_neginf, X1, X2, X3, beta1, beta2, beta3, gamma, beta01, beta02, beta03, *sigSq1, sigSq2, sigSq3, &temp);
            loglh += temp;
            
            log_Jpdf_Upper_BAFT_LN(i, gsl_vector_get(y1, i), gsl_vector_get(y2, i), gsl_vector_get(c0, i), c0_neginf, X1, X2, X3, beta1, beta2, beta3, gamma, beta01, beta02, beta03, sigSq1_prop, sigSq2, sigSq3, &temp);
            loglh_prop += temp;
        }else
        {
            log_Jpdf_Lower_BAFT_LN(i, gsl_vector_get(y2, i), gsl_vector_get(c0, i), c0_neginf, X1, X2, beta1, beta2, gamma, beta01, beta02, *sigSq1, sigSq2, &temp);
            loglh += temp;
            
            log_Jpdf_Lower_BAFT_LN(i, gsl_vector_get(y2, i), gsl_vector_get(c0, i), c0_neginf, X1, X2, beta1, beta2, gamma, beta01, beta02, sigSq1_prop, sigSq2, &temp);
            loglh_prop += temp;
        }
    }
    
    logprior = (-a_sigSq1-1)*log(*sigSq1)-b_sigSq1 /(*sigSq1);
    logprior_prop = (-a_sigSq1-1)*log(sigSq1_prop)-b_sigSq1/sigSq1_prop;
    
    logprior1 = dnorm(beta01, 0, pow(10,6)*sqrt(*sigSq1), 1);
    logprior1_prop = dnorm(beta01, 0, pow(10,6)*sqrt(sigSq1_prop), 1);
    
    logR = loglh_prop - loglh + logprior_prop - logprior + logprior1_prop - logprior1 + log(*sigSq1) - xeta1_prop;
    u = log(runif(0, 1)) < logR;
    
    if(u == 1)
    {
        *sigSq1 = sigSq1_prop;
        *accept_sigSq1 += 1;
    }
    
    return;
}








/* Updating sigSq2 */

void BAFT_LNscr_update_sigSq2(gsl_vector *y1_NA,
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
                              double beta01,
                              double beta02,
                              double beta03,
                              double sigSq1,
                              double *sigSq2,
                              double sigSq3,
                              double a_sigSq2,
                              double b_sigSq2,
                              double sigSq2_prop_var,
                              int *accept_sigSq2)
{
    int i, u;
    double loglh, loglh_prop, logR, temp, xeta2_prop, sigSq2_prop;
    double logprior, logprior_prop, logprior1, logprior1_prop;
    int n = X1 -> size1;
    
    loglh = 0;
    loglh_prop = 0;
    
    xeta2_prop = rnorm(log(*sigSq2), sqrt(sigSq2_prop_var));
    sigSq2_prop = exp(xeta2_prop);
    
    
    for(i=0;i<n;i++)
    {
        if(gsl_vector_get(y1_NA, i) == 0)
        {
            log_Jpdf_Upper_BAFT_LN(i, gsl_vector_get(y1, i), gsl_vector_get(y2, i), gsl_vector_get(c0, i), c0_neginf, X1, X2, X3, beta1, beta2, beta3, gamma, beta01, beta02, beta03, sigSq1, *sigSq2, sigSq3, &temp);
            loglh += temp;
            
            log_Jpdf_Upper_BAFT_LN(i, gsl_vector_get(y1, i), gsl_vector_get(y2, i), gsl_vector_get(c0, i), c0_neginf, X1, X2, X3, beta1, beta2, beta3, gamma, beta01, beta02, beta03, sigSq1, sigSq2_prop, sigSq3, &temp);
            loglh_prop += temp;
        }else
        {
            log_Jpdf_Lower_BAFT_LN(i, gsl_vector_get(y2, i), gsl_vector_get(c0, i), c0_neginf, X1, X2, beta1, beta2, gamma, beta01, beta02, sigSq1, *sigSq2, &temp);
            loglh += temp;
            
            log_Jpdf_Lower_BAFT_LN(i, gsl_vector_get(y2, i), gsl_vector_get(c0, i), c0_neginf, X1, X2, beta1, beta2, gamma, beta01, beta02, sigSq1, sigSq2_prop, &temp);
            loglh_prop += temp;
        }
    }
    
    logprior = (-a_sigSq2-1)*log(*sigSq2)-b_sigSq2 /(*sigSq2);
    logprior_prop = (-a_sigSq2-1)*log(sigSq2_prop)-b_sigSq2/sigSq2_prop;
    
    logprior1 = dnorm(beta02, 0, pow(10,6)*sqrt(*sigSq2), 1);
    logprior1_prop = dnorm(beta02, 0, pow(10,6)*sqrt(sigSq2_prop), 1);
    
    logR = loglh_prop - loglh + logprior_prop - logprior + logprior1_prop - logprior1 + log(*sigSq2) - xeta2_prop;
    u = log(runif(0, 1)) < logR;
    
    if(u == 1)
    {
        *sigSq2 = sigSq2_prop;
        *accept_sigSq2 += 1;
    }
    
    return;
}











/* Updating sigSq3 */

void BAFT_LNscr_update_sigSq3(gsl_vector *y1_NA,
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
                              double beta01,
                              double beta02,
                              double beta03,
                              double sigSq1,
                              double sigSq2,
                              double *sigSq3,
                              double a_sigSq3,
                              double b_sigSq3,
                              double sigSq3_prop_var,
                              int *accept_sigSq3)
{
    int i, u;
    double loglh, loglh_prop, logR, temp, xeta3_prop, sigSq3_prop;
    double logprior, logprior_prop, logprior1, logprior1_prop;
    int n = X1 -> size1;
    
    loglh = 0;
    loglh_prop = 0;
    
    xeta3_prop = rnorm(log(*sigSq3), sqrt(sigSq3_prop_var));
    sigSq3_prop = exp(xeta3_prop);
    
    
    for(i=0;i<n;i++)
    {
        if(gsl_vector_get(y1_NA, i) == 0)
        {
            log_Jpdf_Upper_BAFT_LN(i, gsl_vector_get(y1, i), gsl_vector_get(y2, i), gsl_vector_get(c0, i), c0_neginf, X1, X2, X3, beta1, beta2, beta3, gamma, beta01, beta02, beta03, sigSq1, sigSq2, *sigSq3, &temp);
            loglh += temp;
            
            log_Jpdf_Upper_BAFT_LN(i, gsl_vector_get(y1, i), gsl_vector_get(y2, i), gsl_vector_get(c0, i), c0_neginf, X1, X2, X3, beta1, beta2, beta3, gamma, beta01, beta02, beta03, sigSq1, sigSq2, sigSq3_prop, &temp);
            loglh_prop += temp;
        }
    }
    
    logprior = (-a_sigSq3-1)*log(*sigSq3)-b_sigSq3 /(*sigSq3);
    logprior_prop = (-a_sigSq3-1)*log(sigSq3_prop)-b_sigSq3/sigSq3_prop;
    
    logprior1 = dnorm(beta03, 0, pow(10,6)*sqrt(*sigSq3), 1);
    logprior1_prop = dnorm(beta03, 0, pow(10,6)*sqrt(sigSq3_prop), 1);
    
    logR = loglh_prop - loglh + logprior_prop - logprior +  logprior1_prop - logprior1 + log(*sigSq3) - xeta3_prop;
    u = log(runif(0, 1)) < logR;
    
    if(u == 1)
    {
        *sigSq3 = sigSq3_prop;
        *accept_sigSq3 += 1;
    }
    
    return;
}











/* Updating gamma */

void BAFT_LNscr_update_gamma(gsl_vector *y1_NA,
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
                             double beta01,
                             double beta02,
                             double beta03,
                             double sigSq1,
                             double sigSq2,
                             double sigSq3,
                             double theta,
                             double gamma_prop_var,
                             gsl_vector *accept_gamma)
{
    int i, u;
    double loglh, loglh_prop, logR, temp;
    double logprior, logprior_prop;
    
    int n = X1 -> size1;
    
    gsl_vector *gamma_prop = gsl_vector_calloc(n);
    
    for(i=0;i<n;i++)
    {
        loglh = 0;
        loglh_prop = 0;
        
        gsl_vector_memcpy(gamma_prop, gamma);
        gsl_vector_set(gamma_prop, i, rnorm(gsl_vector_get(gamma, i), sqrt(gamma_prop_var)));
        
        if(gsl_vector_get(y1_NA, i) == 0)
        {
            log_Jpdf_Upper_BAFT_LN(i, gsl_vector_get(y1, i), gsl_vector_get(y2, i), gsl_vector_get(c0, i), c0_neginf, X1, X2, X3, beta1, beta2, beta3, gamma, beta01, beta02, beta03, sigSq1, sigSq2, sigSq3, &temp);
            loglh = temp;
            
            log_Jpdf_Upper_BAFT_LN(i, gsl_vector_get(y1, i), gsl_vector_get(y2, i), gsl_vector_get(c0, i), c0_neginf, X1, X2, X3, beta1, beta2, beta3, gamma_prop, beta01, beta02, beta03, sigSq1, sigSq2, sigSq3, &temp);
            loglh_prop = temp;
        }else
        {
            log_Jpdf_Lower_BAFT_LN(i, gsl_vector_get(y2, i), gsl_vector_get(c0, i), c0_neginf, X1, X2, beta1, beta2, gamma, beta01, beta02, sigSq1, sigSq2, &temp);
            loglh = temp;
            
            log_Jpdf_Lower_BAFT_LN(i, gsl_vector_get(y2, i), gsl_vector_get(c0, i), c0_neginf, X1, X2, beta1, beta2, gamma_prop, beta01, beta02, sigSq1, sigSq2, &temp);
            loglh_prop = temp;
        }
        
        logprior = dnorm(gsl_vector_get(gamma, i), 0, sqrt(theta), 1);
        logprior_prop = dnorm(gsl_vector_get(gamma_prop, i), 0, sqrt(theta), 1);
        
        logR = loglh_prop - loglh + logprior_prop - logprior;
        u = log(runif(0, 1)) < logR;
                
        if(u == 1)
        {
            gsl_vector_memcpy(gamma, gamma_prop);
            gsl_vector_set(accept_gamma, i, gsl_vector_get(accept_gamma, i) + 1);
        }
    }
    
    gsl_vector_free(gamma_prop);
    return;
}






/* Updating theta */

void BAFT_LNscr_update_theta(gsl_vector *gamma,
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





