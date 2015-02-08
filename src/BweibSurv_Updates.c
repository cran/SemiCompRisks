#include <stdio.h>
#include <math.h>

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"

#include "R.h"
#include "Rmath.h"

#include "BweibSurv.h"






/* updating regression parameter: beta */

/**/

void BweibSurv_updateRP(gsl_vector *beta,
              double *alpha,
              double *kappa,
              gsl_vector *survTime,
              gsl_vector *survEvent,
              gsl_matrix *survCov,
              gsl_vector *accept_beta)
{
    double LP, D1, D2, logLH;
    double LP_prop, D1_prop, D2_prop, logLH_prop;
    double beta_prop_me, beta_prop_var, temp_prop;
    double beta_prop_me_prop, beta_prop_var_prop;
    double logProp_IniToProp, logProp_PropToIni;
    double logR;
    int u;
    
    int n = survTime -> size;
    int p = survCov -> size2;
    int i, j;
    
    gsl_vector *beta_prop = gsl_vector_calloc(p);

    
    j = (int) runif(0, p);
    
    logLH = 0; D1 = 0; D2 = 0;
    logLH_prop = 0; D1_prop = 0; D2_prop = 0;
    
    for(i = 0; i < n; i++)
    {
        gsl_vector_view Xi = gsl_matrix_row(survCov, i);
        gsl_blas_ddot(&Xi.vector, beta, &LP);
        if(gsl_vector_get(survEvent, i) == 1)
        {
            logLH   += LP;
            D1      += gsl_matrix_get(survCov, i, j);
        }
        logLH   += -(*kappa) * pow(gsl_vector_get(survTime, i), *alpha) * exp(LP);
        D1      += -(*kappa) * pow(gsl_vector_get(survTime, i), *alpha) * exp(LP) * gsl_matrix_get(survCov, i, j);
        D2      += -(*kappa) * pow(gsl_vector_get(survTime, i), *alpha) * exp(LP) * pow(gsl_matrix_get(survCov, i, j), 2);
    }
    
    beta_prop_me    = gsl_vector_get(beta, j) - D1/D2;
    beta_prop_var   = - pow(2.4, 2)/D2;
    
    temp_prop = rnorm(beta_prop_me, sqrt(beta_prop_var));
    
    gsl_vector_memcpy(beta_prop, beta);
    gsl_vector_set(beta_prop, j, temp_prop);
    
    for(i = 0; i < n; i++)
    {
        gsl_vector_view Xi = gsl_matrix_row(survCov, i);
        gsl_blas_ddot(&Xi.vector, beta_prop, &LP_prop);
        if(gsl_vector_get(survEvent, i) == 1)
        {
            logLH_prop   += LP_prop;
            D1_prop      += gsl_matrix_get(survCov, i, j);
        }
        logLH_prop   += -(*kappa) * pow(gsl_vector_get(survTime, i), *alpha) * exp(LP_prop);
        D1_prop      += -(*kappa) * pow(gsl_vector_get(survTime, i), *alpha) * exp(LP_prop) * gsl_matrix_get(survCov, i, j);
        D2_prop      += -(*kappa) * pow(gsl_vector_get(survTime, i), *alpha) * exp(LP_prop) * pow(gsl_matrix_get(survCov, i, j), 2);
    }
    
    beta_prop_me_prop   = temp_prop - D1_prop/D2_prop;
    beta_prop_var_prop  = - pow(2.4, 2)/D2_prop;
    
    logProp_IniToProp = dnorm(temp_prop, beta_prop_me, sqrt(beta_prop_var), 1);
    logProp_PropToIni = dnorm(gsl_vector_get(beta, j), beta_prop_me_prop, sqrt(beta_prop_var_prop), 1);
    
    logR = logLH_prop - logLH + logProp_PropToIni - logProp_IniToProp;
    
    u = log(runif(0, 1)) < logR;
    
    if(u == 1)
    {
        gsl_vector_set(beta, j, temp_prop);
        gsl_vector_set(accept_beta, j, (gsl_vector_get(accept_beta, j) + u));
    }
    
    gsl_vector_free(beta_prop);
    
    return;
}







/* updating shape parameter: alpha */

/* use the random walk proposal */

/**/

void BweibSurv_updateSC1(gsl_vector *beta,
               double *alpha,
               double *kappa,
               gsl_vector *survTime,
               gsl_vector *survEvent,
               gsl_matrix *survCov,
               double mhProp_alpha_var,
               double a,
               double b,
               int *accept_alpha)
{
    double LP, logLH, logLH_prop;
    double temp_prop;
    double logPrior, logPrior_prop;
    double logProp_IniToProp, logProp_PropToIni;
    double logR;
    int u;
    int n = survTime -> size;
    int i;
    
    logLH = 0; logLH_prop = 0;
    temp_prop = rgamma(pow(*alpha, 2)/mhProp_alpha_var, mhProp_alpha_var/(*alpha));
    
    for(i = 0; i < n; i++)
    {
        gsl_vector_view Xi = gsl_matrix_row(survCov, i);
        gsl_blas_ddot(&Xi.vector, beta, &LP);
        if(gsl_vector_get(survEvent, i) == 1)
        {
            logLH   += log(*alpha) + (*alpha - 1)*log(gsl_vector_get(survTime, i));
            logLH_prop   += log(temp_prop) + (temp_prop - 1)*log(gsl_vector_get(survTime, i));
        }
        logLH += - (*kappa) * pow(gsl_vector_get(survTime, i), *alpha) * exp(LP);
        logLH_prop += - (*kappa) * pow(gsl_vector_get(survTime, i), temp_prop) * exp(LP);
    }
    
     logPrior        = dgamma(*alpha, a, 1/b, 1);
     logPrior_prop   = dgamma(temp_prop, a, 1/b, 1);
     
     logProp_PropToIni = dgamma(*alpha, pow(temp_prop, 2)/mhProp_alpha_var, mhProp_alpha_var/(temp_prop), 1);
     logProp_IniToProp = dgamma(temp_prop, pow(*alpha, 2)/mhProp_alpha_var, mhProp_alpha_var/(*alpha), 1);
     
     logR = logLH_prop - logLH + logPrior_prop - logPrior + logProp_PropToIni - logProp_IniToProp;

    u = log(runif(0, 1)) <logR;
    
    if(u == 1)
    {
        *alpha = temp_prop;
        *accept_alpha += u;
    }
    
    return;
    
}











/* updating shape parameter: alpha */

/* use the prior as proposal */

/**/

void BweibSurv_updateSC2(gsl_vector *beta,
              double *alpha,
              double *kappa,
              gsl_vector *survTime,
              gsl_vector *survEvent,
              gsl_matrix *survCov,
              double mhProp_alpha_var,
              double a,
              double b,
              int *accept_alpha)
{
    double LP, logLH, logLH_prop;
    double temp_prop;
    double logR;
    int u;
    int n = survTime -> size;
    int i;
    
    logLH = 0; logLH_prop = 0;
    temp_prop = rgamma(a, 1/b);
    
    for(i = 0; i < n; i++)
    {
        gsl_vector_view Xi = gsl_matrix_row(survCov, i);
        gsl_blas_ddot(&Xi.vector, beta, &LP);
        if(gsl_vector_get(survEvent, i) == 1)
        {
            logLH   += log(*alpha) + (*alpha - 1)*log(gsl_vector_get(survTime, i));
            logLH_prop   += log(temp_prop) + (temp_prop - 1)*log(gsl_vector_get(survTime, i));
        }
        logLH += - (*kappa) * pow(gsl_vector_get(survTime, i), *alpha) * exp(LP);
        logLH_prop += - (*kappa) * pow(gsl_vector_get(survTime, i), temp_prop) * exp(LP);
    }
    
    
    logR = logLH_prop - logLH;
    
    u = log(runif(0, 1)) <logR;
    
    if(u == 1)
    {
        *alpha = temp_prop;
        *accept_alpha += u;
    }
    
    return;
    
}






/* updating shape parameter: kappa */

/**/

void BweibSurv_updateSH(gsl_vector *beta,
              double *alpha,
              double *kappa,
              gsl_vector *survTime,
              gsl_vector *survEvent,
              gsl_matrix *survCov,
              double c,
              double d)
{
    int n = survTime -> size;
    int i;
    double LP, Kappa_shape, Kappa_rate, Kappa_scale;
    
    gsl_vector *ones = gsl_vector_calloc(n);
    gsl_vector_set_all(ones, 1);
    
    gsl_blas_ddot(ones, survEvent, &Kappa_shape);
    Kappa_shape += c;
    
    gsl_vector_free(ones);
    
    Kappa_rate = 0;
    
    for(i = 0; i< n; i++)
    {
        gsl_vector_view Xi = gsl_matrix_row(survCov, i);
        gsl_blas_ddot(&Xi.vector, beta ,&LP);
        
        Kappa_rate += exp(LP) * pow(gsl_vector_get(survTime, i), *alpha);
    }
    Kappa_rate += d;
    Kappa_scale = 1/Kappa_rate;
    
    *kappa = rgamma(Kappa_shape, Kappa_scale);
    return;
}


