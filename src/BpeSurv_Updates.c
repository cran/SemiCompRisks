#include <stdio.h>
#include <math.h>

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_sort_vector.h"
#include "gsl/gsl_heapsort.h"

#include "R.h"
#include "Rmath.h"

#include "BpeSurv.h"






/* updating regression parameter: beta */

/**/

void BpeSur_updateRP2(gsl_vector *beta,
               gsl_vector *xbeta,
               gsl_vector *accept_beta,
               gsl_vector *lambda,
               gsl_vector *s,
               gsl_vector *survTime,
               gsl_vector *survEvent,
               gsl_matrix *survCov,
               int J)
{
    double D1, D2, logLH;
    double D1_prop, D2_prop, logLH_prop;
    double beta_prop_me, beta_prop_var, temp_prop;
    double beta_prop_me_prop, beta_prop_var_prop;
    double logProp_IniToProp, logProp_PropToIni;
    double logR, Del;
    int u, m, i, j;
    
    int p = beta -> size;
    int n = survTime -> size;
    
    gsl_vector *beta_prop = gsl_vector_calloc(p);
    
    m = (int) runif(0, p);
    
    /* m = 3; */
    
    logLH = 0; D1 = 0; D2 = 0;
    logLH_prop = 0; D1_prop = 0; D2_prop = 0;
    
    for(i = 0; i < n; i++)
    {
        if(gsl_vector_get(survEvent, i) == 1)
        {
            logLH += gsl_vector_get(xbeta, i);
            D1 += gsl_matrix_get(survCov, i, m);
        }

        for(j = 0; j < J+1; j++)
        {
            if(j > 0)
            {
                Del = c_max(0, (c_min(gsl_vector_get(s, j), gsl_vector_get(survTime, i)) - gsl_vector_get(s, j-1)));
            }
            if(j == 0)
            {
                Del = c_max(0, c_min(gsl_vector_get(s, j), gsl_vector_get(survTime, i)) - 0);
            }
            if(Del > 0)
            {
                logLH   += - Del*exp(gsl_vector_get(lambda, j))*exp(gsl_vector_get(xbeta, i));
                D1      += - Del*exp(gsl_vector_get(lambda, j))*exp(gsl_vector_get(xbeta, i))*gsl_matrix_get(survCov, i, m);
                D2      += - Del*exp(gsl_vector_get(lambda, j))*exp(gsl_vector_get(xbeta, i))*pow(gsl_matrix_get(survCov, i, m), 2);
            }
        }
    }
    

        
    
    beta_prop_me    = gsl_vector_get(beta, m) - D1/D2;
    beta_prop_var   = - pow(2.4, 2)/D2;
    
    temp_prop = rnorm(beta_prop_me, sqrt(beta_prop_var));
    
    gsl_vector_memcpy(beta_prop, beta);
    gsl_vector_set(beta_prop, m, temp_prop);
    
    gsl_vector *xbeta_prop = gsl_vector_calloc(n);
    
    gsl_blas_dgemv(CblasNoTrans, 1, survCov, beta_prop, 0, xbeta_prop);
    
    for(i = 0; i < n; i++)
    {
        if(gsl_vector_get(survEvent, i) == 1)
        {
            logLH_prop += gsl_vector_get(xbeta_prop, i);
            D1_prop += gsl_matrix_get(survCov, i, m);
        }
        
        for(j = 0; j < J+1; j++)
        {
            if(j > 0)
            {
                Del = c_max(0, (c_min(gsl_vector_get(s, j), gsl_vector_get(survTime, i)) - gsl_vector_get(s, j-1)));
            }
            if(j == 0)
            {
                Del = c_max(0, c_min(gsl_vector_get(s, j), gsl_vector_get(survTime, i)) - 0);
            }
            if(Del > 0)
            {
                logLH_prop   += - Del*exp(gsl_vector_get(lambda, j))*exp(gsl_vector_get(xbeta_prop, i));
                D1_prop      += - Del*exp(gsl_vector_get(lambda, j))*exp(gsl_vector_get(xbeta_prop, i))*gsl_matrix_get(survCov, i, m);
                D2_prop      += - Del*exp(gsl_vector_get(lambda, j))*exp(gsl_vector_get(xbeta_prop, i))*pow(gsl_matrix_get(survCov, i, m), 2);
            }
        }
    }
    
    beta_prop_me_prop   = temp_prop - D1_prop/D2_prop;
    beta_prop_var_prop  = - pow(2.4, 2)/D2_prop;
    
    logProp_IniToProp = dnorm(temp_prop, beta_prop_me, sqrt(beta_prop_var), 1);
    logProp_PropToIni = dnorm(gsl_vector_get(beta, m), beta_prop_me_prop, sqrt(beta_prop_var_prop), 1);
    
    logR = logLH_prop - logLH + logProp_PropToIni - logProp_IniToProp;
    
    u = log(runif(0, 1)) < logR;
    
    if(u == 1)
    {
        gsl_vector_set(beta, m, temp_prop);
        gsl_vector_swap(xbeta, xbeta_prop);
        gsl_vector_set(accept_beta, m, (gsl_vector_get(accept_beta, m) + u));
    }
    
    gsl_vector_free(beta_prop);
    gsl_vector_free(xbeta_prop);
    
      
    return;
}








/* updating regression parameter: beta */

/**/

void BpeSur_updateRP1(gsl_vector *beta,
              gsl_vector *xbeta,
              gsl_vector *accept_beta,              
              gsl_vector *lambda,
              gsl_vector *survTime,
              gsl_vector *survEvent,
              gsl_matrix *survCov,
              gsl_matrix *ind_r,
              gsl_matrix *ind_d,
              gsl_matrix *Delta,
              int J)
{
    double D1, D2, logLH;
    double D1_prop, D2_prop, logLH_prop;
    double beta_prop_me, beta_prop_var, temp_prop;
    double beta_prop_me_prop, beta_prop_var_prop;
    double logProp_IniToProp, logProp_PropToIni;
    double logR;
    int u, m, i, j;
    
    int p = beta -> size;
    int n = survTime -> size;
    
    gsl_vector *beta_prop = gsl_vector_calloc(p);
    
    m = (int) runif(0, p);
    
    logLH = 0; D1 = 0; D2 = 0;
    logLH_prop = 0; D1_prop = 0; D2_prop = 0;
    
    for(j = 0; j < J+1; j++)
    {
        for(i = 0; i < n; i++)
        {
            if(gsl_matrix_get(ind_r, i, j) == 1)
            {
                logLH   += -exp(gsl_vector_get(lambda, j))*gsl_matrix_get(Delta, i, j)*exp(gsl_vector_get(xbeta, i));
                D1      += -exp(gsl_vector_get(lambda, j))*gsl_matrix_get(Delta, i, j)*exp(gsl_vector_get(xbeta, i))*gsl_matrix_get(survCov, i, m);
                D2      += -exp(gsl_vector_get(lambda, j))*gsl_matrix_get(Delta, i, j)*exp(gsl_vector_get(xbeta, i))* pow(gsl_matrix_get(survCov, i, m), 2);
            }
            if(gsl_matrix_get(ind_d, i, j) == 1)
            {
                logLH   += gsl_vector_get(xbeta, i);
                D1      += gsl_matrix_get(survCov, i, m);
            }
        }
    }
    
    
    beta_prop_me    = gsl_vector_get(beta, m) - D1/D2;
    beta_prop_var   = - pow(2.4, 2)/D2;
    
    temp_prop = rnorm(beta_prop_me, sqrt(beta_prop_var));
    
    gsl_vector_memcpy(beta_prop, beta);
    gsl_vector_set(beta_prop, m, temp_prop);
    
    gsl_vector *xbeta_prop = gsl_vector_calloc(n);
    
    gsl_blas_dgemv(CblasNoTrans, 1, survCov, beta_prop, 0, xbeta_prop);
    
    for(j = 0; j < J+1; j++)
    {
        for(i = 0; i < n; i++)
        {
            if(gsl_matrix_get(ind_r, i, j) == 1)
            {
                logLH_prop   += -exp(gsl_vector_get(lambda, j))*gsl_matrix_get(Delta, i, j)*exp(gsl_vector_get(xbeta_prop, i));
                D1_prop      += -exp(gsl_vector_get(lambda, j))*gsl_matrix_get(Delta, i, j)*exp(gsl_vector_get(xbeta_prop, i))*gsl_matrix_get(survCov, i, m);
                D2_prop      += -exp(gsl_vector_get(lambda, j))*gsl_matrix_get(Delta, i, j)*exp(gsl_vector_get(xbeta_prop, i))* pow(gsl_matrix_get(survCov, i, m), 2);
            }
            if(gsl_matrix_get(ind_d, i, j) == 1)
            {
                logLH_prop   += gsl_vector_get(xbeta_prop, i);
                D1_prop      += gsl_matrix_get(survCov, i, m);
            }
        }
    }
    
    beta_prop_me_prop   = temp_prop - D1_prop/D2_prop;
    beta_prop_var_prop  = - pow(2.4, 2)/D2_prop;
    
    logProp_IniToProp = dnorm(temp_prop, beta_prop_me, sqrt(beta_prop_var), 1);
    logProp_PropToIni = dnorm(gsl_vector_get(beta, m), beta_prop_me_prop, sqrt(beta_prop_var_prop), 1);
    
    logR = logLH_prop - logLH + logProp_PropToIni - logProp_IniToProp;
    
    u = log(runif(0, 1)) < logR;
    
    if(u == 1)
    {
        gsl_vector_set(beta, m, temp_prop);
        gsl_vector_swap(xbeta, xbeta_prop);
        gsl_vector_set(accept_beta, m, (gsl_vector_get(accept_beta, m) + u));
    }
    
    gsl_vector_free(beta_prop);
    gsl_vector_free(xbeta_prop);
    
    
    
    return;
}


























/* updating log-baseline hazard function parameter: lambda */



void BpeSur_updateBH2(gsl_vector *lambda,
               gsl_vector *s,
               gsl_vector *xbeta,
               gsl_vector *survTime,
               gsl_vector *survEvent,
               gsl_matrix *Sigma_lam,
               gsl_matrix *invSigma_lam,
               gsl_matrix *W,
               gsl_matrix *Q,
               double mu_lam,
               double sigSq_lam,
               int J)
{
    double D1, D2, logLH, Del, inc;
    double D1_prop, D2_prop, logLH_prop;
    double lambda_prop_me, lambda_prop_var, temp_prop;
    double lambda_prop_me_prop, lambda_prop_var_prop;
    double logPrior, logPrior_prop;
    double logProp_IniToProp, logProp_PropToIni;
    double logR;
    int u, i, j;
    
    double nu_lam, nu_lam_prop;
    
    int n = xbeta -> size;
    
    j = (int) runif(0, J+1);

    
    if(J+1 > 1)
    {
        if(j == 0) nu_lam = mu_lam + gsl_matrix_get(W, 0, 1) * (gsl_vector_get(lambda, 1) - mu_lam);
        if(j == J) nu_lam = mu_lam + gsl_matrix_get(W, J, J-1) * (gsl_vector_get(lambda, J-1) - mu_lam);
        if(j != 0 && j !=J) nu_lam = mu_lam + gsl_matrix_get(W, j, j-1) * (gsl_vector_get(lambda, j-1) - mu_lam) + gsl_matrix_get(W, j, j+1) * (gsl_vector_get(lambda, j+1) - mu_lam);
    }
    
    if(J+1 == 1)
    {
        nu_lam = mu_lam;
    }
    
    logLH = 0; D1 = 0; D2 = 0;
    logLH_prop = 0; D1_prop = 0; D2_prop = 0;
  
  
    
    
    for(i = 0; i < n; i++)
    {
        
        if(gsl_vector_get(survEvent, i) == 1)
        {
            if(j == 0 && gsl_vector_get(survTime, i) <= gsl_vector_get(s, 0))
            {
                logLH += gsl_vector_get(lambda, j);
                D1 += 1;
            }
            if(j != 0 && gsl_vector_get(survTime, i) > gsl_vector_get(s, j-1) && gsl_vector_get(survTime, i) <= gsl_vector_get(s, j))
            {
                logLH += gsl_vector_get(lambda, j);
                D1 += 1;                
            }
            
        }

        if(j > 0)
        {
            Del = c_max(0, (c_min(gsl_vector_get(s, j), gsl_vector_get(survTime, i)) - gsl_vector_get(s, j-1)));
        }
        if(j == 0)
        {
            Del = c_max(0, c_min(gsl_vector_get(s, j), gsl_vector_get(survTime, i)) - 0);
        }
        if(Del > 0)
        {
            inc = - Del*exp(gsl_vector_get(lambda, j))*exp(gsl_vector_get(xbeta, i));
            logLH   += inc;
            D1      += inc;
            D2      += inc;
        }
    }
    

    D1      += -1/(sigSq_lam * gsl_matrix_get(Q, j, j))*(gsl_vector_get(lambda, j)-nu_lam);
    D2      += -1/(sigSq_lam * gsl_matrix_get(Q, j, j));
    
    
    lambda_prop_me    = gsl_vector_get(lambda, j) - D1/D2;
    lambda_prop_var   = - pow(2.4, 2)/D2;
    
    temp_prop = rnorm(lambda_prop_me, sqrt(lambda_prop_var));
    
    gsl_vector *lambda_prop = gsl_vector_calloc(J+1);
    
    gsl_vector_view lambda_sub = gsl_vector_subvector(lambda, 0, J+1);
    
    gsl_vector_memcpy(lambda_prop, &lambda_sub.vector);
    gsl_vector_set(lambda_prop, j, temp_prop);
    
    if(J+1 > 1)
    {
        if(j == 0) nu_lam_prop = mu_lam + gsl_matrix_get(W, 0, 1) * (gsl_vector_get(lambda_prop, 1) - mu_lam);
        if(j == J) nu_lam_prop = mu_lam + gsl_matrix_get(W, J, J-1) * (gsl_vector_get(lambda_prop, J-1) - mu_lam);
        if(j != 0 && j !=J) nu_lam_prop = mu_lam + gsl_matrix_get(W, j, j-1) * (gsl_vector_get(lambda_prop, j-1) - mu_lam) + gsl_matrix_get(W, j, j+1) * (gsl_vector_get(lambda_prop, j+1) - mu_lam);
    }
    
    if(J+1 == 1)
    {
        nu_lam_prop = mu_lam;
    }
    


    
    for(i = 0; i < n; i++)
    {
        
        if(gsl_vector_get(survEvent, i) == 1)
        {
            if(j == 0 && gsl_vector_get(survTime, i) <= gsl_vector_get(s, 0))
            {
                logLH_prop += gsl_vector_get(lambda_prop, j);
                D1_prop += 1;
            }
            if(j != 0 && gsl_vector_get(survTime, i) > gsl_vector_get(s, j-1) && gsl_vector_get(survTime, i) <= gsl_vector_get(s, j))
            {
                logLH_prop += gsl_vector_get(lambda_prop, j);
                D1_prop += 1;
            }
            
        }
        
        if(j > 0)
        {
            Del = c_max(0, (c_min(gsl_vector_get(s, j), gsl_vector_get(survTime, i)) - gsl_vector_get(s, j-1)));
        }
        if(j == 0)
        {
            Del = c_max(0, c_min(gsl_vector_get(s, j), gsl_vector_get(survTime, i)) - 0);
        }
        if(Del > 0)
        {
            inc = - Del*exp(gsl_vector_get(lambda_prop, j))*exp(gsl_vector_get(xbeta, i));
            logLH_prop   += inc;
            D1_prop      += inc;
            D2_prop      += inc;
        }
    }
    
    
    D1_prop += -1/(sigSq_lam * gsl_matrix_get(Q, j, j))*(gsl_vector_get(lambda_prop, j)-nu_lam);
    D2_prop += -1/(sigSq_lam * gsl_matrix_get(Q, j, j));
  

    
    
    lambda_prop_me_prop    = gsl_vector_get(lambda_prop, j) - D1_prop/D2_prop;
    lambda_prop_var_prop   = - pow(2.4, 2)/D2_prop;
    
    gsl_matrix_view invS_sub = gsl_matrix_submatrix(invSigma_lam, 0, 0, J+1, J+1);
    
    if(J+1 > 1)
    {
        c_dmvnorm(&lambda_sub.vector, mu_lam, sqrt(sigSq_lam), &invS_sub.matrix, &logPrior);
        c_dmvnorm(lambda_prop, mu_lam, sqrt(sigSq_lam), &invS_sub.matrix, &logPrior_prop);
    }
    if(J+1 == 1)
    {
        logPrior        = dnorm(gsl_vector_get(lambda, j), mu_lam, sqrt(sigSq_lam*gsl_matrix_get(Sigma_lam, 0, 0)), 1);
        logPrior_prop   = dnorm(temp_prop, mu_lam, sqrt(sigSq_lam*gsl_matrix_get(Sigma_lam, 0, 0)), 1);
    }
    
    logProp_IniToProp = dnorm(temp_prop, lambda_prop_me, sqrt(lambda_prop_var), 1);
    logProp_PropToIni = dnorm(gsl_vector_get(lambda, j), lambda_prop_me_prop, sqrt(lambda_prop_var_prop), 1);
    
    logR = logLH_prop - logLH + logPrior_prop - logPrior +  logProp_PropToIni - logProp_IniToProp;
    
    
    u = log(runif(0, 1)) < logR;
  
    if(u == 1) gsl_vector_set(lambda, j, temp_prop);
    
    gsl_vector_free(lambda_prop);
 
    
      
    
    return;
}












/* updating log-baseline hazard function parameter: lambda */



void BpeSur_updateBH1(gsl_vector *lambda,
              gsl_vector *xbeta,
              gsl_matrix *ind_r,
              gsl_matrix *Delta,
              gsl_vector *nEvent,
              gsl_matrix *Sigma_lam,
              gsl_matrix *invSigma_lam,
              gsl_matrix *W,
              gsl_matrix *Q,
              double mu_lam,
              double sigSq_lam,
              int J)
{
    double D1, D2, logLH;
    double D1_prop, D2_prop, logLH_prop;
    double lambda_prop_me, lambda_prop_var, temp_prop;
    double lambda_prop_me_prop, lambda_prop_var_prop;
    double logPrior, logPrior_prop;
    double logProp_IniToProp, logProp_PropToIni;
    double logR;
    int u, i, j;
    
    double nu_lam, nu_lam_prop;
    
    int n = xbeta -> size;
    
    j = (int) runif(0, J+1);
    
    
    if(J+1 > 1)
    {
        if(j == 0) nu_lam = mu_lam + gsl_matrix_get(W, 0, 1) * (gsl_vector_get(lambda, 1) - mu_lam);
        if(j == J) nu_lam = mu_lam + gsl_matrix_get(W, J, J-1) * (gsl_vector_get(lambda, J-1) - mu_lam);
        if(j != 0 && j !=J) nu_lam = mu_lam + gsl_matrix_get(W, j, j-1) * (gsl_vector_get(lambda, j-1) - mu_lam) + gsl_matrix_get(W, j, j+1) * (gsl_vector_get(lambda, j+1) - mu_lam);
    }
    
    if(J+1 == 1)
    {
        nu_lam = mu_lam;
    }
    
    logLH = 0; D1 = 0; D2 = 0;
    logLH_prop = 0; D1_prop = 0; D2_prop = 0;
    
    
    for(i = 0; i < n; i++)
    {
        if(gsl_matrix_get(ind_r, i, j) == 1)
        {
            logLH += -exp(gsl_vector_get(lambda, j))*gsl_matrix_get(Delta, i, j)*exp(gsl_vector_get(xbeta, i));
            D1 += -exp(gsl_vector_get(lambda, j))*gsl_matrix_get(Delta, i, j)*exp(gsl_vector_get(xbeta, i));
            D2 += -exp(gsl_vector_get(lambda, j))*gsl_matrix_get(Delta, i, j)*exp(gsl_vector_get(xbeta, i));
        }
    }
    
    logLH   += gsl_vector_get(lambda, j) * gsl_vector_get(nEvent, j);
    D1      += gsl_vector_get(nEvent, j)-1/(sigSq_lam * gsl_matrix_get(Q, j, j))*(gsl_vector_get(lambda, j)-nu_lam);
    D2      += -1/(sigSq_lam * gsl_matrix_get(Q, j, j));
    
    
    lambda_prop_me    = gsl_vector_get(lambda, j) - D1/D2;
    lambda_prop_var   = - pow(2.4, 2)/D2;
    
    temp_prop = rnorm(lambda_prop_me, sqrt(lambda_prop_var));
    
    
   
    gsl_vector *lambda_prop = gsl_vector_calloc(J+1);
    
    gsl_vector_view lambda_sub = gsl_vector_subvector(lambda, 0, J+1);
    
    gsl_vector_memcpy(lambda_prop, &lambda_sub.vector);
    gsl_vector_set(lambda_prop, j, temp_prop);
    
    if(J+1 > 1)
    {
        if(j == 0) nu_lam_prop = mu_lam + gsl_matrix_get(W, 0, 1) * (gsl_vector_get(lambda_prop, 1) - mu_lam);
        if(j == J) nu_lam_prop = mu_lam + gsl_matrix_get(W, J, J-1) * (gsl_vector_get(lambda_prop, J-1) - mu_lam);
        if(j != 0 && j !=J) nu_lam_prop = mu_lam + gsl_matrix_get(W, j, j-1) * (gsl_vector_get(lambda_prop, j-1) - mu_lam) + gsl_matrix_get(W, j, j+1) * (gsl_vector_get(lambda_prop, j+1) - mu_lam);
    }
    
    if(J+1 == 1)
    {
        nu_lam_prop = mu_lam;
    }
    
    
    for(i = 0; i < n; i++)
    {
        if(gsl_matrix_get(ind_r, i, j) == 1)
        {
            logLH_prop += -exp(gsl_vector_get(lambda_prop, j))*gsl_matrix_get(Delta, i, j)*exp(gsl_vector_get(xbeta, i));
            D1_prop += -exp(gsl_vector_get(lambda_prop, j))*gsl_matrix_get(Delta, i, j)*exp(gsl_vector_get(xbeta, i));
            D2_prop += -exp(gsl_vector_get(lambda_prop, j))*gsl_matrix_get(Delta, i, j)*exp(gsl_vector_get(xbeta, i));
        }
    }
    
    logLH_prop   += gsl_vector_get(lambda_prop, j) * gsl_vector_get(nEvent, j);
    D1_prop      += gsl_vector_get(nEvent, j)-1/(sigSq_lam * gsl_matrix_get(Q, j, j))*(gsl_vector_get(lambda_prop, j)-nu_lam);
    D2_prop      += -1/(sigSq_lam * gsl_matrix_get(Q, j, j));
    
    lambda_prop_me_prop    = gsl_vector_get(lambda_prop, j) - D1_prop/D2_prop;
    lambda_prop_var_prop   = - pow(2.4, 2)/D2_prop;
    
    gsl_matrix_view invS_sub = gsl_matrix_submatrix(invSigma_lam, 0, 0, J+1, J+1);
    
    if(J+1 > 1)
    {
        c_dmvnorm(&lambda_sub.vector, mu_lam, sqrt(sigSq_lam), &invS_sub.matrix, &logPrior);
        c_dmvnorm(lambda_prop, mu_lam, sqrt(sigSq_lam), &invS_sub.matrix, &logPrior_prop);
    }
    if(J+1 == 1)
    {
        logPrior        = dnorm(gsl_vector_get(lambda, j), mu_lam, sqrt(sigSq_lam*gsl_matrix_get(Sigma_lam, 0, 0)), 1);
        logPrior_prop   = dnorm(temp_prop, mu_lam, sqrt(sigSq_lam*gsl_matrix_get(Sigma_lam, 0, 0)), 1);
    }
    
    logProp_IniToProp = dnorm(temp_prop, lambda_prop_me, sqrt(lambda_prop_var), 1);
    logProp_PropToIni = dnorm(gsl_vector_get(lambda, j), lambda_prop_me_prop, sqrt(lambda_prop_var_prop), 1);
    
    logR = logLH_prop - logLH + logPrior_prop - logPrior +  logProp_PropToIni - logProp_IniToProp;
    
    u = log(runif(0, 1)) < logR;

    
    if(u == 1) gsl_vector_set(lambda, j, temp_prop);

    gsl_vector_free(lambda_prop);
    
      
    return;
}





















/* Updating second stage survival components: mu_lam and sigSq_lam */

void BpeSur_updateSP(double *mu_lam,
              double *sigSq_lam,
              gsl_vector *lambda,
              gsl_matrix *Sigma_lam,
              gsl_matrix *invSigma_lam,
              double a,
              double b,
              int J)
{
    double num, den, sigSH, sigRT, sigSC, tau, mu_lam_mean, mu_lam_var;
    
    gsl_vector *ones = gsl_vector_calloc(J+1);
    gsl_vector_set_all(ones, 1);
    
    gsl_matrix_view invSlam_sub = gsl_matrix_submatrix(invSigma_lam, 0, 0, J+1, J+1);
    gsl_vector_view lam_sub     = gsl_vector_subvector(lambda, 0, J+1);
    
    c_quadform_vMu(ones, &invSlam_sub.matrix, &lam_sub.vector, &num);
    c_quadform_vMv(ones, &invSlam_sub.matrix, &den);
    
    mu_lam_mean = num/den;
    mu_lam_var = *sigSq_lam/den;
    
    *mu_lam = rnorm(mu_lam_mean, sqrt(mu_lam_var));
    
    gsl_vector *diff = gsl_vector_calloc(J+1);
    gsl_vector_set_all(diff, *mu_lam);
    gsl_vector_sub(diff, &lam_sub.vector);
    c_quadform_vMv(diff, &invSlam_sub.matrix, &sigRT);
    sigRT /= 2;
    sigRT += b;
    sigSC = 1/sigRT;
    sigSH = a + (double) (J+1)/2;
    tau = rgamma(sigSH, sigSC);
    *sigSq_lam = 1/tau;
    
    gsl_vector_free(ones);
    gsl_vector_free(diff);
    
    return;
}










/* Updating the number of splits and their positions: J and s (Birth move) */


void BpeSur_updateBI2(gsl_vector *s,
               int *J,
               int *accept_BI,
               gsl_vector *survTime,
               gsl_vector *survEvent,
               gsl_vector *xbeta,
               gsl_matrix *Sigma_lam,
               gsl_matrix *invSigma_lam,
               gsl_matrix *W,
               gsl_matrix *Q,
               gsl_vector *lambda,
               gsl_vector *s_propBI,
               int num_s_propBI,
               double delPert,
               int alpha,
               double c_lam,
               double mu_lam,
               double sigSq_lam,
               double s_max)
{
    int count, num_s_propBI_fin, skip;
    int star_inx, j_old, J_new, i, j, u;
    double s_star, Upert, newLam1, newLam2;
    double logLH, logLH_prop, Del;
    double logPrior, logPrior_prop, logPriorR, logPropR;
    double logJacob, logR;
    
    int n = xbeta -> size;
    
    count = 0;
    
    gsl_vector *interInx = gsl_vector_calloc(num_s_propBI);
    
    for(i = 0; i < num_s_propBI; i++)
    {
        for(j = 0; j < *J+1; j++)
        {
            if(gsl_vector_get(s_propBI, i) == gsl_vector_get(s, j))
            {
                count += 1;
                gsl_vector_set(interInx, count-1, i);
            }
        }
    }
    
    gsl_vector *s_propBI_fin = gsl_vector_calloc(num_s_propBI-count);
    
    num_s_propBI_fin = s_propBI_fin -> size;
    skip = 0;
    
    if(count > 0)
    {
        for(i = 0; i < num_s_propBI; i++)
        {
            if(i != gsl_vector_get(interInx, skip))
            {
                gsl_vector_set(s_propBI_fin, i-skip, gsl_vector_get(s_propBI, i));
            }
            if(i == gsl_vector_get(interInx, skip)) skip += 1;
        }
    }
    if(count == 0) gsl_vector_memcpy(s_propBI_fin, s_propBI);
    
    
    star_inx = (int) runif(0, num_s_propBI_fin);
    s_star = gsl_vector_get(s_propBI_fin, star_inx);
        
    j_old = -1;
    i = 0;
    
    while(j_old < 0)
    {
        if(gsl_vector_get(s, i) >= s_star) j_old = i;
        else i += 1;
    }
    
    gsl_vector *s_new = gsl_vector_calloc(*J+2);
    for(i = 0; i < *J+1; i++)
    {
        gsl_vector_set(s_new, i, gsl_vector_get(s, i));
    }
    gsl_vector_set(s_new, *J+1, s_star);
    gsl_sort_vector(s_new);
    
    J_new = *J+1;
    
    Upert = runif(0.5 - delPert, 0.5 + delPert);
    
    if(j_old != 0)
    {
        newLam1 = gsl_vector_get(lambda, j_old) - (gsl_vector_get(s, j_old) - s_star)/(gsl_vector_get(s, j_old) - gsl_vector_get(s, j_old-1)) * log((1-Upert)/Upert);
        newLam2 = gsl_vector_get(lambda, j_old) + (s_star - gsl_vector_get(s, j_old-1))/(gsl_vector_get(s, j_old) - gsl_vector_get(s, j_old-1)) * log((1-Upert)/Upert);
    }
    
    if(j_old == 0)
    {
        newLam1 = gsl_vector_get(lambda, j_old) - (gsl_vector_get(s, j_old) - s_star)/(gsl_vector_get(s, j_old) - 0) * log((1-Upert)/Upert);
        newLam2 = gsl_vector_get(lambda, j_old) + (s_star - 0)/(gsl_vector_get(s, j_old) - 0) * log((1-Upert)/Upert);
    }
    
    gsl_vector *lambda_new = gsl_vector_calloc(*J+2);
    
    skip = 0;
    for(i = 0; i < *J+2; i++)
    {
        if(i == j_old)
        {
            gsl_vector_set(lambda_new, i, newLam1);
        }
        else if(i == j_old + 1)
        {
            skip += 1;
            gsl_vector_set(lambda_new, i, newLam2);
        }
        else gsl_vector_set(lambda_new, i, gsl_vector_get(lambda, i-skip));
    }
    
    
    gsl_matrix *Sigma_lam_new       = gsl_matrix_calloc(J_new+1, J_new+1);
    gsl_matrix *invSigma_lam_new    = gsl_matrix_calloc(J_new+1, J_new+1);
    gsl_matrix *W_new               = gsl_matrix_calloc(J_new+1, J_new+1);
    gsl_matrix *Q_new               = gsl_matrix_calloc(J_new+1, J_new+1);
    
    cal_Sigma(Sigma_lam_new, invSigma_lam_new, W_new, Q_new, s_new, c_lam, J_new);
    
    logLH = 0; logLH_prop = 0;
    logPrior = 0; logPrior_prop = 0; logPropR = 0;
    

    
    for(i = 0; i < n; i++)
    {
        if(gsl_vector_get(survEvent, i) == 1)
        {
            if(j_old == 0 && gsl_vector_get(survTime, i) <= gsl_vector_get(s, 0))
            {
                logLH += gsl_vector_get(lambda, j_old);
            }
            if(j_old != 0 && gsl_vector_get(survTime, i) > gsl_vector_get(s, j_old-1) && gsl_vector_get(survTime, i) <= gsl_vector_get(s, j_old))
            {
                logLH += gsl_vector_get(lambda, j_old);                
            }
        }
        if(j_old > 0)
        {
            Del = c_max(0, (c_min(gsl_vector_get(s, j_old), gsl_vector_get(survTime, i)) - gsl_vector_get(s, j_old-1)));
        }
        if(j_old == 0)
        {
            Del = c_max(0, c_min(gsl_vector_get(s, j_old), gsl_vector_get(survTime, i)) - 0);
        }
        if(Del > 0)
        {
            logLH   += - Del*exp(gsl_vector_get(lambda, j_old))*exp(gsl_vector_get(xbeta, i));
        }
    }
 
    
    
    for(j = j_old; j < j_old + 2; j++)
    {
        for(i = 0; i < n; i++)
        {
            if(gsl_vector_get(survEvent, i) == 1)
            {
                if(j == 0 && gsl_vector_get(survTime, i) <= gsl_vector_get(s_new, 0))
                {
                    logLH_prop += gsl_vector_get(lambda_new, j);
                }
                if(j != 0 && gsl_vector_get(survTime, i) > gsl_vector_get(s_new, j-1) && gsl_vector_get(survTime, i) <= gsl_vector_get(s_new, j))
                {
                    logLH_prop += gsl_vector_get(lambda_new, j);
                }
            }
            if(j > 0)
            {
                Del = c_max(0, (c_min(gsl_vector_get(s_new, j), gsl_vector_get(survTime, i)) - gsl_vector_get(s_new, j-1)));
            }
            if(j == 0)
            {
                Del = c_max(0, c_min(gsl_vector_get(s_new, j), gsl_vector_get(survTime, i)) - 0);
            }
            if(Del > 0)
            {
                logLH_prop   += - Del*exp(gsl_vector_get(lambda_new, j))*exp(gsl_vector_get(xbeta, i));
            }
        }
    }
    
    
    gsl_vector_view lambda_sub  = gsl_vector_subvector(lambda, 0, *J+1);
    gsl_matrix_view invS_sub    = gsl_matrix_submatrix(invSigma_lam, 0, 0, *J+1, *J+1);
    
    if(*J+1 != 1)
    {
        c_dmvnorm(&lambda_sub.vector, mu_lam, sqrt(sigSq_lam), &invS_sub.matrix, &logPrior);
        c_dmvnorm(lambda_new, mu_lam, sqrt(sigSq_lam), invSigma_lam_new, &logPrior_prop);
        
        if(j_old != 0)
        {
            logPrior_prop += log( (2*(*J) + 3)*(2*(*J) + 2) * pow(gsl_vector_get(s, *J), -2) * (s_star - gsl_vector_get(s, j_old-1)) * (gsl_vector_get(s, j_old) - s_star)/(gsl_vector_get(s, j_old) - gsl_vector_get(s, j_old-1)) );
        }
        if(j_old == 0)
        {
            logPrior_prop += log( (2*(*J) + 3)*(2*(*J) + 2) * pow(gsl_vector_get(s, *J), -2) * (s_star - 0) * (gsl_vector_get(s, j_old) - s_star)/(gsl_vector_get(s, j_old) - 0) );
        }
    }
    
    if(*J+1 == 1)
    {
        logPrior = dnorm(gsl_vector_get(lambda, 0), mu_lam, sqrt(sigSq_lam*gsl_matrix_get(Sigma_lam, 0, 0)), 1);
        c_dmvnorm(lambda_new, mu_lam, sqrt(sigSq_lam), invSigma_lam_new, &logPrior_prop);
        
        logPrior_prop += log( (2*(*J) + 3)*(2*(*J) + 2) * pow(gsl_vector_get(s, *J), -2) * (s_star - 0) * (gsl_vector_get(s, j_old) - s_star)/(gsl_vector_get(s, j_old) - 0) );
    }
    
    
    logPriorR = log((double) alpha/((*J) + 1)) + logPrior_prop - logPrior;
    
    logPropR = log(num_s_propBI_fin/alpha) - dunif(Upert, 0.5-delPert, 0.5+delPert, 1);
    
    logJacob = log(1/(1-Upert)/Upert);
    
    logR = logLH_prop - logLH + logPriorR + logPropR + logJacob;
    
    u = log(runif(0, 1)) < logR;
    
    
    if(u == 1)
    {
        gsl_matrix_view Sigma_lam_save      = gsl_matrix_submatrix(Sigma_lam, 0, 0, J_new+1, J_new+1);
        gsl_matrix_view invSigma_lam_save   = gsl_matrix_submatrix(invSigma_lam, 0, 0, J_new+1, J_new+1);
        gsl_matrix_view W_save              = gsl_matrix_submatrix(W, 0, 0, J_new+1, J_new+1);
        gsl_matrix_view Q_save              = gsl_matrix_submatrix(Q, 0, 0, J_new+1, J_new+1);
        gsl_vector_view s_save              = gsl_vector_subvector(s, 0, J_new+1);
        gsl_vector_view lambda_save         = gsl_vector_subvector(lambda, 0, J_new+1);
        
        *accept_BI += 1;
        *J = J_new;
        
        gsl_matrix_memcpy(&Sigma_lam_save.matrix, Sigma_lam_new);
        gsl_matrix_memcpy(&invSigma_lam_save.matrix, invSigma_lam_new);
        gsl_matrix_memcpy(&W_save.matrix, W_new);
        gsl_matrix_memcpy(&Q_save.matrix, Q_new);
        gsl_vector_memcpy(&s_save.vector, s_new);
        gsl_vector_memcpy(&lambda_save.vector, lambda_new);
    }
    
    gsl_vector_free(interInx);
    gsl_vector_free(s_propBI_fin);
    gsl_vector_free(s_new);
    gsl_vector_free(lambda_new);
    gsl_matrix_free(Sigma_lam_new);
    gsl_matrix_free(invSigma_lam_new);
    gsl_matrix_free(W_new);
    gsl_matrix_free(Q_new);
    
        
    return;
}















/* Updating the number of splits and their positions: J and s (Birth move) */


void BpeSur_updateBI1(gsl_vector *s,
              int *J,
              int *accept_BI,
              gsl_vector *survTime,
              gsl_vector *survEvent,
              gsl_vector *case0yleq,
              gsl_vector *case0ygeq,
              gsl_vector *case1yleq,
              gsl_vector *case1ygeq,
              gsl_vector *xbeta,
              gsl_matrix *ind_r,
              gsl_matrix *ind_d,
              gsl_vector *nEvent,
              gsl_matrix *Delta,
              gsl_matrix *Sigma_lam,
              gsl_matrix *invSigma_lam,
              gsl_matrix *W,
              gsl_matrix *Q,
              gsl_vector *lambda,
              gsl_vector *s_propBI,
              int num_s_propBI,
              double delPert,
              int alpha,
              double c_lam,
              double mu_lam,
              double sigSq_lam,              
              double s_max)
{
    int count, num_s_propBI_fin, skip;
    int star_inx, j_old, J_new, i, j, u;
    double s_star, Upert, newLam1, newLam2;
    double logLH, logLH_prop;
    double logPrior, logPrior_prop, logPriorR, logPropR;
    double logJacob, logR;
    
    int n = xbeta -> size;
    
    count = 0;
    
    gsl_vector *interInx = gsl_vector_calloc(num_s_propBI);
    
    for(i = 0; i < num_s_propBI; i++)
    {
        for(j = 0; j < *J+1; j++)
        {
            if(gsl_vector_get(s_propBI, i) == gsl_vector_get(s, j))
            {
                count += 1;
                gsl_vector_set(interInx, count-1, i);
            }
        }
    }
    
    gsl_vector *s_propBI_fin = gsl_vector_calloc(num_s_propBI-count);
    
    num_s_propBI_fin = s_propBI_fin -> size;
    skip = 0;
    
    if(count > 0)
    {
        for(i = 0; i < num_s_propBI; i++)
        {
            if(i != gsl_vector_get(interInx, skip))
            {
                gsl_vector_set(s_propBI_fin, i-skip, gsl_vector_get(s_propBI, i));
            }
            if(i == gsl_vector_get(interInx, skip)) skip += 1;
        }
    }
    if(count == 0) gsl_vector_memcpy(s_propBI_fin, s_propBI);
    
    star_inx = (int) runif(0, num_s_propBI_fin);
    s_star = gsl_vector_get(s_propBI_fin, star_inx);
    
    j_old = -1;
    i = 0;
    
    while(j_old < 0)
    {
        if(gsl_vector_get(s, i) >= s_star) j_old = i;
        else i += 1;
    }
    
    gsl_vector *s_new = gsl_vector_calloc(*J+2);
    for(i = 0; i < *J+1; i++)
    {
        gsl_vector_set(s_new, i, gsl_vector_get(s, i));
    }
    gsl_vector_set(s_new, *J+1, s_star);
    gsl_sort_vector(s_new);
    
    J_new = *J+1;
    
    Upert = runif(0.5 - delPert, 0.5 + delPert);
    
    if(j_old != 0)
    {
        newLam1 = gsl_vector_get(lambda, j_old) - (gsl_vector_get(s, j_old) - s_star)/(gsl_vector_get(s, j_old) - gsl_vector_get(s, j_old-1)) * log((1-Upert)/Upert);
        newLam2 = gsl_vector_get(lambda, j_old) + (s_star - gsl_vector_get(s, j_old-1))/(gsl_vector_get(s, j_old) - gsl_vector_get(s, j_old-1)) * log((1-Upert)/Upert);
    }
    
    if(j_old == 0)
    {
        newLam1 = gsl_vector_get(lambda, j_old) - (gsl_vector_get(s, j_old) - s_star)/(gsl_vector_get(s, j_old) - 0) * log((1-Upert)/Upert);
        newLam2 = gsl_vector_get(lambda, j_old) + (s_star - 0)/(gsl_vector_get(s, j_old) - 0) * log((1-Upert)/Upert);
    }
    
    gsl_vector *lambda_new = gsl_vector_calloc(*J+2);
    
    skip = 0;
    for(i = 0; i < *J+2; i++)
    {
        if(i == j_old)
        {
            gsl_vector_set(lambda_new, i, newLam1);
        }
        else if(i == j_old + 1)
        {
            skip += 1;
            gsl_vector_set(lambda_new, i, newLam2);
        }
        else gsl_vector_set(lambda_new, i, gsl_vector_get(lambda, i-skip));
    }
    
    gsl_matrix *ind_d_new   = gsl_matrix_calloc(n, J_new+1);
    gsl_matrix *ind_r_new   = gsl_matrix_calloc(n, J_new+1);
    gsl_vector *nEvent_new  = gsl_vector_calloc(J_new+1);
    
    
    set_Ind(ind_d_new, ind_r_new, nEvent_new, s_new, survTime, survEvent, case0yleq, case0ygeq, case1yleq, case1ygeq, s_max, J_new);
    
    
    gsl_matrix *Delta_new = gsl_matrix_alloc(n, J_new+1);
    cal_Delta(Delta_new, survTime, s_new, J_new);
    
    gsl_matrix *Sigma_lam_new       = gsl_matrix_calloc(J_new+1, J_new+1);
    gsl_matrix *invSigma_lam_new    = gsl_matrix_calloc(J_new+1, J_new+1);
    gsl_matrix *W_new               = gsl_matrix_calloc(J_new+1, J_new+1);
    gsl_matrix *Q_new               = gsl_matrix_calloc(J_new+1, J_new+1);
    
    cal_Sigma(Sigma_lam_new, invSigma_lam_new, W_new, Q_new, s_new, c_lam, J_new);
    
    logLH = 0; logLH_prop = 0;
    logPrior = 0; logPrior_prop = 0; logPropR = 0;
    
    
    for(i = 0; i < n; i++)
    {
        if(gsl_matrix_get(ind_r, i, j_old) == 1)
        {
            logLH   += -exp(gsl_vector_get(lambda, j_old))*gsl_matrix_get(Delta, i, j_old)*exp(gsl_vector_get(xbeta, i));
        }
        if(gsl_matrix_get(ind_d, i, j_old) == 1)
        {
            logLH   += gsl_vector_get(lambda, j_old);
        }
    }
    
    
    
    
    for(j = j_old; j < j_old + 2; j++)
    {
        for(i = 0; i < n; i++)
        {
            if(gsl_matrix_get(ind_r_new, i, j) == 1)
            {
                logLH_prop   += -exp(gsl_vector_get(lambda_new, j))*gsl_matrix_get(Delta_new, i, j)*exp(gsl_vector_get(xbeta, i));
            }
            if(gsl_matrix_get(ind_d_new, i, j) == 1)
            {
                logLH_prop   += gsl_vector_get(lambda_new, j);
            }
        }
    }
    
    
    gsl_vector_view lambda_sub  = gsl_vector_subvector(lambda, 0, *J+1);
    gsl_matrix_view invS_sub    = gsl_matrix_submatrix(invSigma_lam, 0, 0, *J+1, *J+1);
    
    if(*J+1 != 1)
    {
        c_dmvnorm(&lambda_sub.vector, mu_lam, sqrt(sigSq_lam), &invS_sub.matrix, &logPrior);
        c_dmvnorm(lambda_new, mu_lam, sqrt(sigSq_lam), invSigma_lam_new, &logPrior_prop);
        
        if(j_old != 0)
        {
            logPrior_prop += log( (2*(*J) + 3)*(2*(*J) + 2) * pow(gsl_vector_get(s, *J), -2) * (s_star - gsl_vector_get(s, j_old-1)) * (gsl_vector_get(s, j_old) - s_star)/(gsl_vector_get(s, j_old) - gsl_vector_get(s, j_old-1)) );
        }
        if(j_old == 0)
        {
            logPrior_prop += log( (2*(*J) + 3)*(2*(*J) + 2) * pow(gsl_vector_get(s, *J), -2) * (s_star - 0) * (gsl_vector_get(s, j_old) - s_star)/(gsl_vector_get(s, j_old) - 0) );
        }
    }
    
    if(*J+1 == 1)
    {
        logPrior = dnorm(gsl_vector_get(lambda, 0), mu_lam, sqrt(sigSq_lam*gsl_matrix_get(Sigma_lam, 0, 0)), 1);
        c_dmvnorm(lambda_new, mu_lam, sqrt(sigSq_lam), invSigma_lam_new, &logPrior_prop);
        
        logPrior_prop += log( (2*(*J) + 3)*(2*(*J) + 2) * pow(gsl_vector_get(s, *J), -2) * (s_star - 0) * (gsl_vector_get(s, j_old) - s_star)/(gsl_vector_get(s, j_old) - 0) );
    }
    
    
    logPriorR = log((double) alpha/((*J) + 1)) + logPrior_prop - logPrior;
     
    logPropR = log(num_s_propBI_fin/alpha) - dunif(Upert, 0.5-delPert, 0.5+delPert, 1);
    
    logJacob = log(1/(1-Upert)/Upert);
    
    logR = logLH_prop - logLH + logPriorR + logPropR + logJacob;
    
    u = log(runif(0, 1)) < logR;
    
    
    if(u == 1)
    {
        gsl_matrix_view ind_r_save          = gsl_matrix_submatrix(ind_r, 0, 0, n, J_new+1);
        gsl_matrix_view ind_d_save          = gsl_matrix_submatrix(ind_d, 0, 0, n, J_new+1);
        gsl_vector_view nEvent_save         = gsl_vector_subvector(nEvent, 0, J_new+1);
        gsl_matrix_view Delta_save          = gsl_matrix_submatrix(Delta, 0, 0, n, J_new+1);
        gsl_matrix_view Sigma_lam_save      = gsl_matrix_submatrix(Sigma_lam, 0, 0, J_new+1, J_new+1);
        gsl_matrix_view invSigma_lam_save   = gsl_matrix_submatrix(invSigma_lam, 0, 0, J_new+1, J_new+1);
        gsl_matrix_view W_save              = gsl_matrix_submatrix(W, 0, 0, J_new+1, J_new+1);
        gsl_matrix_view Q_save              = gsl_matrix_submatrix(Q, 0, 0, J_new+1, J_new+1);
        gsl_vector_view s_save              = gsl_vector_subvector(s, 0, J_new+1);
        gsl_vector_view lambda_save         = gsl_vector_subvector(lambda, 0, J_new+1);
        
        *accept_BI += 1;
        *J = J_new;
        
        gsl_matrix_memcpy(&ind_r_save.matrix, ind_r_new);
        gsl_matrix_memcpy(&ind_d_save.matrix, ind_d_new);
        gsl_vector_memcpy(&nEvent_save.vector, nEvent_new);
        gsl_matrix_memcpy(&Delta_save.matrix, Delta_new);
        gsl_matrix_memcpy(&Sigma_lam_save.matrix, Sigma_lam_new);
        gsl_matrix_memcpy(&invSigma_lam_save.matrix, invSigma_lam_new);
        gsl_matrix_memcpy(&W_save.matrix, W_new);
        gsl_matrix_memcpy(&Q_save.matrix, Q_new);
        gsl_vector_memcpy(&s_save.vector, s_new);
        gsl_vector_memcpy(&lambda_save.vector, lambda_new);
    }

    gsl_vector_free(interInx);
    gsl_vector_free(s_propBI_fin);
    gsl_vector_free(s_new);
    gsl_vector_free(lambda_new);
    gsl_matrix_free(ind_d_new);
    gsl_matrix_free(ind_r_new);
    gsl_vector_free(nEvent_new);
    gsl_matrix_free(Delta_new);
    gsl_matrix_free(Sigma_lam_new);
    gsl_matrix_free(invSigma_lam_new);
    gsl_matrix_free(W_new);
    gsl_matrix_free(Q_new);
    
    return;
}













/* Updating the number of splits and their positions: J and s (Death move) */


void BpeSur_updateDI2(gsl_vector *s,
               int *J,
               int *accept_DI,
               gsl_vector *survTime,
               gsl_vector *survEvent,
               gsl_vector *xbeta,
               gsl_matrix *Sigma_lam,
               gsl_matrix *invSigma_lam,
               gsl_matrix *W,
               gsl_matrix *Q,
               gsl_vector *lambda,
               gsl_vector *s_propBI,
               int num_s_propBI,
               double delPert,
               int alpha,
               double c_lam,
               double mu_lam,
               double sigSq_lam,
               double s_max,
               int J_max)
{
    
    
    int skip, i, j;
    int j_old, J_new, u;
    double s_star, Upert, newLam;
    double logLH, logLH_prop, Del;
    double logPrior, logPrior_prop, logPriorR, logPropR;
    double logJacob, logR;
    
    int n = xbeta -> size;
    
    j_old = (int) runif(0, *J);
    
    gsl_vector *s_new = gsl_vector_calloc(*J);
    
    skip = 0;
    for(i = 0; i < *J+1; i++)
    {
        if(i == j_old) skip += 1;
        else gsl_vector_set(s_new, i-skip, gsl_vector_get(s, i));
    }
    
    J_new = *J-1;
    
    Upert = 1/(exp(gsl_vector_get(lambda, j_old+1) - gsl_vector_get(lambda, j_old)) + 1);
    
    if(j_old != 0)
    {
        newLam = ((gsl_vector_get(s, j_old) - gsl_vector_get(s, j_old-1)) * gsl_vector_get(lambda, j_old) + (gsl_vector_get(s, j_old+1) - gsl_vector_get(s, j_old)) * gsl_vector_get(lambda, j_old+1)) / (gsl_vector_get(s, j_old+1) - gsl_vector_get(s, j_old-1));
    }
    
    if(j_old == 0)
    {
        newLam = ((gsl_vector_get(s, j_old) - 0) * gsl_vector_get(lambda, j_old) + (gsl_vector_get(s, j_old+1) - gsl_vector_get(s, j_old)) * gsl_vector_get(lambda, j_old+1)) / (gsl_vector_get(s, j_old+1) - 0);
    }
    
    gsl_vector *lambda_new = gsl_vector_calloc(J_new+1);
    
    skip = 0;
    for(i = 0; i < J_new+1; i++)
    {
        if(i == j_old){
            gsl_vector_set(lambda_new, i, newLam);
            skip += 1;
        }
        else gsl_vector_set(lambda_new, i, gsl_vector_get(lambda, i+skip));
    }
    
    gsl_matrix *Sigma_lam_new       = gsl_matrix_calloc(J_new+1, J_new+1);
    gsl_matrix *invSigma_lam_new    = gsl_matrix_calloc(J_new+1, J_new+1);
    gsl_matrix *W_new               = gsl_matrix_calloc(J_new+1, J_new+1);
    gsl_matrix *Q_new               = gsl_matrix_calloc(J_new+1, J_new+1);
    
    cal_Sigma(Sigma_lam_new, invSigma_lam_new, W_new, Q_new, s_new, c_lam, J_new);
    
    logLH = 0; logLH_prop = 0;
    logPrior = 0; logPrior_prop = 0; logPropR = 0;
    
    
    for(j = j_old; j < j_old + 2; j++)
    {
        for(i = 0; i < n; i++)
        {
            if(gsl_vector_get(survEvent, i) == 1)
            {
                if(j == 0 && gsl_vector_get(survTime, i) <= gsl_vector_get(s, 0))
                {
                    logLH += gsl_vector_get(lambda, j);
                }
                if(j != 0 && gsl_vector_get(survTime, i) > gsl_vector_get(s, j-1) && gsl_vector_get(survTime, i) <= gsl_vector_get(s, j))
                {
                    logLH += gsl_vector_get(lambda, j);
                }
            }
            if(j > 0)
            {
                Del = c_max(0, (c_min(gsl_vector_get(s, j), gsl_vector_get(survTime, i)) - gsl_vector_get(s, j-1)));
            }
            if(j == 0)
            {
                Del = c_max(0, c_min(gsl_vector_get(s, j), gsl_vector_get(survTime, i)) - 0);
            }
            if(Del > 0)
            {
                logLH   += - Del*exp(gsl_vector_get(lambda, j))*exp(gsl_vector_get(xbeta, i));
            }
        }
    }
    
       
    
    for(i = 0; i < n; i++)
    {
        if(gsl_vector_get(survEvent, i) == 1)
        {
            if(j_old == 0 && gsl_vector_get(survTime, i) <= gsl_vector_get(s_new, 0))
            {
                logLH_prop += gsl_vector_get(lambda_new, j_old);
            }
            if(j_old != 0 && gsl_vector_get(survTime, i) > gsl_vector_get(s_new, j_old-1) && gsl_vector_get(survTime, i) <= gsl_vector_get(s_new, j_old))
            {
                logLH_prop += gsl_vector_get(lambda_new, j_old);
            }
        }
        if(j_old > 0)
        {
            Del = c_max(0, (c_min(gsl_vector_get(s_new, j_old), gsl_vector_get(survTime, i)) - gsl_vector_get(s_new, j_old-1)));
        }
        if(j_old == 0)
        {
            Del = c_max(0, c_min(gsl_vector_get(s_new, j_old), gsl_vector_get(survTime, i)) - 0);
        }
        if(Del > 0)
        {
            logLH_prop   += - Del*exp(gsl_vector_get(lambda_new, j_old))*exp(gsl_vector_get(xbeta, i));
        }
    }
    
    
    gsl_vector_view lambda_sub  = gsl_vector_subvector(lambda, 0, *J+1);
    gsl_matrix_view invS_sub    = gsl_matrix_submatrix(invSigma_lam, 0, 0, *J+1, *J+1);
    
    
    if(*J+1 != 2)
    {
        c_dmvnorm(&lambda_sub.vector, mu_lam, sqrt(sigSq_lam), &invS_sub.matrix, &logPrior);
        c_dmvnorm(lambda_new, mu_lam, sqrt(sigSq_lam), invSigma_lam_new, &logPrior_prop);
        
        if(j_old != 0)
        {
            logPrior_prop += log(( (double) 1/(2*(*J) + 1)/(2*(*J)))*pow(gsl_vector_get(s, *J), 2)*(gsl_vector_get(s, j_old+1) - gsl_vector_get(s, j_old-1))/(gsl_vector_get(s, j_old) - gsl_vector_get(s, j_old-1))/(gsl_vector_get(s, j_old+1) - gsl_vector_get(s, j_old)));
        }
        if(j_old == 0)
        {
            logPrior_prop += log(( (double) 1/(2*(*J)+1)/(2*(*J)))*pow(gsl_vector_get(s, *J), 2)*(gsl_vector_get(s, j_old+1) - 0)/(gsl_vector_get(s, j_old) - 0)/(gsl_vector_get(s, j_old+1) - gsl_vector_get(s, j_old)));
        }
    }
    
    
    if(*J+1 == 2)
    {
        logPrior_prop = dnorm(gsl_vector_get(lambda_new, 0), mu_lam, sqrt(sigSq_lam*gsl_matrix_get(Sigma_lam_new, 0, 0)), 1);
        c_dmvnorm(lambda, mu_lam, sqrt(sigSq_lam), invSigma_lam, &logPrior);
        
        logPrior_prop += log(( (double) 1/(2*(*J)+1)/(2*(*J)))*pow(gsl_vector_get(s, *J), 2)*(gsl_vector_get(s, j_old+1) - 0)/(gsl_vector_get(s, j_old) - 0)/(gsl_vector_get(s, j_old+1) - gsl_vector_get(s, j_old)));
    }
    
    
    logPriorR = log((double) *J/alpha) + logPrior_prop - logPrior;
    
    logPropR = log((double) alpha/num_s_propBI) - dunif(Upert, 0.5-delPert, 0.5+delPert, 1);
    
    logJacob = log((1-Upert)*Upert);
    
    logR = logLH_prop - logLH + logPriorR + logPropR + logJacob;
    
    u = log(runif(0, 1)) < logR;
    
    
    if(u == 1)
    {
 
        gsl_matrix_view Sigma_lam_save      = gsl_matrix_submatrix(Sigma_lam, 0, 0, J_new+1, J_new+1);
        gsl_matrix_view invSigma_lam_save   = gsl_matrix_submatrix(invSigma_lam, 0, 0, J_new+1, J_new+1);
        gsl_matrix_view W_save              = gsl_matrix_submatrix(W, 0, 0, J_new+1, J_new+1);
        gsl_matrix_view Q_save              = gsl_matrix_submatrix(Q, 0, 0, J_new+1, J_new+1);
        gsl_vector_view s_save              = gsl_vector_subvector(s, 0, J_new+1);
        gsl_vector_view lambda_save         = gsl_vector_subvector(lambda, 0, J_new+1);
        
        gsl_matrix_memcpy(&Sigma_lam_save.matrix, Sigma_lam_new);
        gsl_matrix_memcpy(&invSigma_lam_save.matrix, invSigma_lam_new);
        gsl_matrix_memcpy(&W_save.matrix, W_new);
        gsl_matrix_memcpy(&Q_save.matrix, Q_new);
        gsl_vector_memcpy(&s_save.vector, s_new);
        gsl_vector_memcpy(&lambda_save.vector, lambda_new);
        

        gsl_vector *zeroVec_J = gsl_vector_calloc(J_max+1);
        

        gsl_matrix_set_col(Sigma_lam, *J, zeroVec_J);
        gsl_matrix_set_row(Sigma_lam, *J, zeroVec_J);
        gsl_matrix_set_col(invSigma_lam, *J, zeroVec_J);
        gsl_matrix_set_row(invSigma_lam, *J, zeroVec_J);
        gsl_matrix_set_col(W, *J, zeroVec_J);
        gsl_matrix_set_row(W, *J, zeroVec_J);
        gsl_matrix_set_col(Q, *J, zeroVec_J);
        gsl_matrix_set_row(Q, *J, zeroVec_J);
        gsl_vector_set(s, *J, 0);
        gsl_vector_set(lambda, *J, 0);
        
        *accept_DI += 1;
        *J = J_new;
        
        gsl_vector_free(zeroVec_J);
        
    }
    
    gsl_vector_free(s_new);
    gsl_vector_free(lambda_new);
    gsl_matrix_free(Sigma_lam_new);
    gsl_matrix_free(invSigma_lam_new);
    gsl_matrix_free(W_new);
    gsl_matrix_free(Q_new);
    
    
    return;
}























































/* Updating the number of splits and their positions: J and s (Death move) */ 


void BpeSur_updateDI1(gsl_vector *s,
              int *J,
              int *accept_DI,
              gsl_vector *survTime,
              gsl_vector *survEvent,
              gsl_vector *case0yleq,
              gsl_vector *case0ygeq,
              gsl_vector *case1yleq,
              gsl_vector *case1ygeq,
              gsl_vector *xbeta,
              gsl_matrix *ind_r,
              gsl_matrix *ind_d,
              gsl_vector *nEvent,
              gsl_matrix *Delta,
              gsl_matrix *Sigma_lam,
              gsl_matrix *invSigma_lam,
              gsl_matrix *W,
              gsl_matrix *Q,
              gsl_vector *lambda,
              gsl_vector *s_propBI,
              int num_s_propBI,
              double delPert,
              int alpha,
              double c_lam,
              double mu_lam,
              double sigSq_lam,
              double s_max,
              int J_max)
{
    
    
    int skip, i, j;
    int j_old, J_new, u;
    double s_star, Upert, newLam;
    double logLH, logLH_prop;
    double logPrior, logPrior_prop, logPriorR, logPropR;
    double logJacob, logR;
    
    int n = xbeta -> size;
    
    j_old = (int) runif(0, *J);
    
    gsl_vector *s_new = gsl_vector_calloc(*J);
    
    skip = 0;
    for(i = 0; i < *J+1; i++)
    {
        if(i == j_old) skip += 1;
        else gsl_vector_set(s_new, i-skip, gsl_vector_get(s, i));
    }
    
    J_new = *J-1;
    
    Upert = 1/(exp(gsl_vector_get(lambda, j_old+1) - gsl_vector_get(lambda, j_old)) + 1);
    
    if(j_old != 0)
    {
        newLam = ((gsl_vector_get(s, j_old) - gsl_vector_get(s, j_old-1)) * gsl_vector_get(lambda, j_old) + (gsl_vector_get(s, j_old+1) - gsl_vector_get(s, j_old)) * gsl_vector_get(lambda, j_old+1)) / (gsl_vector_get(s, j_old+1) - gsl_vector_get(s, j_old-1));
    }
    
    if(j_old == 0)
    {
        newLam = ((gsl_vector_get(s, j_old) - 0) * gsl_vector_get(lambda, j_old) + (gsl_vector_get(s, j_old+1) - gsl_vector_get(s, j_old)) * gsl_vector_get(lambda, j_old+1)) / (gsl_vector_get(s, j_old+1) - 0);
    }
    
    gsl_vector *lambda_new = gsl_vector_calloc(J_new+1);
    
    skip = 0;
    for(i = 0; i < J_new+1; i++)
    {
        if(i == j_old){
            gsl_vector_set(lambda_new, i, newLam);
            skip += 1;
        }
        else gsl_vector_set(lambda_new, i, gsl_vector_get(lambda, i+skip));
    }
    
    gsl_matrix *ind_d_new   = gsl_matrix_calloc(n, J_new+1);
    gsl_matrix *ind_r_new   = gsl_matrix_calloc(n, J_new+1);
    gsl_vector *nEvent_new  = gsl_vector_calloc(J_new+1);
    
    set_Ind(ind_d_new, ind_r_new, nEvent_new, s_new, survTime, survEvent, case0yleq, case0ygeq, case1yleq, case1ygeq, s_max, J_new);
    
    gsl_matrix *Delta_new = gsl_matrix_alloc(n, J_new+1);
    cal_Delta(Delta_new, survTime, s_new, J_new);
    
    gsl_matrix *Sigma_lam_new       = gsl_matrix_calloc(J_new+1, J_new+1);
    gsl_matrix *invSigma_lam_new    = gsl_matrix_calloc(J_new+1, J_new+1);
    gsl_matrix *W_new               = gsl_matrix_calloc(J_new+1, J_new+1);
    gsl_matrix *Q_new               = gsl_matrix_calloc(J_new+1, J_new+1);
    
    cal_Sigma(Sigma_lam_new, invSigma_lam_new, W_new, Q_new, s_new, c_lam, J_new);
    
    logLH = 0; logLH_prop = 0;
    logPrior = 0; logPrior_prop = 0; logPropR = 0;
    
    
    for(j = j_old; j < j_old + 2; j++)
    {
        for(i = 0; i < n; i++)
        {
            if(gsl_matrix_get(ind_r, i, j) == 1)
            {
                logLH   += -exp(gsl_vector_get(lambda, j))*gsl_matrix_get(Delta, i, j)*exp(gsl_vector_get(xbeta, i));
            }
            if(gsl_matrix_get(ind_d, i, j) == 1)
            {
                logLH   += gsl_vector_get(lambda, j);
            }
        }
    }
    
    
    for(i = 0; i < n; i++)
    {
        if(gsl_matrix_get(ind_r_new, i, j_old) == 1)
        {
            logLH_prop   += -exp(gsl_vector_get(lambda_new, j_old))*gsl_matrix_get(Delta_new, i, j_old)*exp(gsl_vector_get(xbeta, i));
        }
        if(gsl_matrix_get(ind_d_new, i, j_old) == 1)
        {
            logLH_prop   += gsl_vector_get(lambda_new, j_old);
        }
    }
    
    gsl_vector_view lambda_sub  = gsl_vector_subvector(lambda, 0, *J+1);
    gsl_matrix_view invS_sub    = gsl_matrix_submatrix(invSigma_lam, 0, 0, *J+1, *J+1);
    
    
    if(*J+1 != 2)
    {
        c_dmvnorm(&lambda_sub.vector, mu_lam, sqrt(sigSq_lam), &invS_sub.matrix, &logPrior);
        c_dmvnorm(lambda_new, mu_lam, sqrt(sigSq_lam), invSigma_lam_new, &logPrior_prop);
        
        if(j_old != 0)
        {
            logPrior_prop += log(( (double) 1/(2*(*J) + 1)/(2*(*J)))*pow(gsl_vector_get(s, *J), 2)*(gsl_vector_get(s, j_old+1) - gsl_vector_get(s, j_old-1))/(gsl_vector_get(s, j_old) - gsl_vector_get(s, j_old-1))/(gsl_vector_get(s, j_old+1) - gsl_vector_get(s, j_old)));
        }
        if(j_old == 0)
        {
            logPrior_prop += log(( (double) 1/(2*(*J)+1)/(2*(*J)))*pow(gsl_vector_get(s, *J), 2)*(gsl_vector_get(s, j_old+1) - 0)/(gsl_vector_get(s, j_old) - 0)/(gsl_vector_get(s, j_old+1) - gsl_vector_get(s, j_old)));
        }
    }
    
    
    if(*J+1 == 2)
    {
        logPrior_prop = dnorm(gsl_vector_get(lambda_new, 0), mu_lam, sqrt(sigSq_lam*gsl_matrix_get(Sigma_lam_new, 0, 0)), 1);
        c_dmvnorm(lambda, mu_lam, sqrt(sigSq_lam), invSigma_lam, &logPrior);
        
        logPrior_prop += log(( (double) 1/(2*(*J)+1)/(2*(*J)))*pow(gsl_vector_get(s, *J), 2)*(gsl_vector_get(s, j_old+1) - 0)/(gsl_vector_get(s, j_old) - 0)/(gsl_vector_get(s, j_old+1) - gsl_vector_get(s, j_old)));
    }
    
    
    logPriorR = log((double) *J/alpha) + logPrior_prop - logPrior;
    
    logPropR = log((double) alpha/num_s_propBI) - dunif(Upert, 0.5-delPert, 0.5+delPert, 1);
    
    logJacob = log((1-Upert)*Upert);
    
    logR = logLH_prop - logLH + logPriorR + logPropR + logJacob;
    
    u = log(runif(0, 1)) < logR;
    
  
    
    if(u == 1)
    {
        gsl_matrix_view ind_r_save          = gsl_matrix_submatrix(ind_r, 0, 0, n, J_new+1);
        gsl_matrix_view ind_d_save          = gsl_matrix_submatrix(ind_d, 0, 0, n, J_new+1);
        gsl_vector_view nEvent_save         = gsl_vector_subvector(nEvent, 0, J_new+1);
        gsl_matrix_view Delta_save          = gsl_matrix_submatrix(Delta, 0, 0, n, J_new+1);
        gsl_matrix_view Sigma_lam_save      = gsl_matrix_submatrix(Sigma_lam, 0, 0, J_new+1, J_new+1);
        gsl_matrix_view invSigma_lam_save   = gsl_matrix_submatrix(invSigma_lam, 0, 0, J_new+1, J_new+1);
        gsl_matrix_view W_save              = gsl_matrix_submatrix(W, 0, 0, J_new+1, J_new+1);
        gsl_matrix_view Q_save              = gsl_matrix_submatrix(Q, 0, 0, J_new+1, J_new+1);
        gsl_vector_view s_save              = gsl_vector_subvector(s, 0, J_new+1);
        gsl_vector_view lambda_save         = gsl_vector_subvector(lambda, 0, J_new+1);

        gsl_matrix_memcpy(&ind_r_save.matrix, ind_r_new);
        gsl_matrix_memcpy(&ind_d_save.matrix, ind_d_new);
        gsl_vector_memcpy(&nEvent_save.vector, nEvent_new);
        gsl_matrix_memcpy(&Delta_save.matrix, Delta_new);
        gsl_matrix_memcpy(&Sigma_lam_save.matrix, Sigma_lam_new);
        gsl_matrix_memcpy(&invSigma_lam_save.matrix, invSigma_lam_new);
        gsl_matrix_memcpy(&W_save.matrix, W_new);
        gsl_matrix_memcpy(&Q_save.matrix, Q_new);
        gsl_vector_memcpy(&s_save.vector, s_new);
        gsl_vector_memcpy(&lambda_save.vector, lambda_new);
        
        gsl_vector *zeroVec_n = gsl_vector_calloc(n);
        gsl_vector *zeroVec_J = gsl_vector_calloc(J_max+1);
        
        gsl_matrix_set_col(ind_r, *J, zeroVec_n);
        gsl_matrix_set_col(ind_d, *J, zeroVec_n);
        gsl_vector_set(nEvent, *J, 0);
        gsl_matrix_set_col(Delta, *J, zeroVec_n);
        gsl_matrix_set_col(Sigma_lam, *J, zeroVec_J);
        gsl_matrix_set_row(Sigma_lam, *J, zeroVec_J);
        gsl_matrix_set_col(invSigma_lam, *J, zeroVec_J);
        gsl_matrix_set_row(invSigma_lam, *J, zeroVec_J);
        gsl_matrix_set_col(W, *J, zeroVec_J);
        gsl_matrix_set_row(W, *J, zeroVec_J);
        gsl_matrix_set_col(Q, *J, zeroVec_J);
        gsl_matrix_set_row(Q, *J, zeroVec_J);
        gsl_vector_set(s, *J, 0);
        gsl_vector_set(lambda, *J, 0);
        
        *accept_DI += 1;
        *J = J_new;
        
        gsl_vector_free(zeroVec_n);
        gsl_vector_free(zeroVec_J);
        
    }
    
    gsl_vector_free(s_new);
    gsl_vector_free(lambda_new);
    gsl_matrix_free(ind_d_new);
    gsl_matrix_free(ind_r_new);
    gsl_vector_free(nEvent_new);
    gsl_matrix_free(Delta_new);
    gsl_matrix_free(Sigma_lam_new);
    gsl_matrix_free(invSigma_lam_new);
    gsl_matrix_free(W_new);
    gsl_matrix_free(Q_new);
    
    
    return;
}











