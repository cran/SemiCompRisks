#include <stdio.h>
#include <math.h>

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"

#include "R.h"
#include "Rmath.h"

#include "BweibScrSM.h"








/* updating regression parameter: beta1 */

/**/

void BweibScrSM_updateRP1(gsl_vector *beta1,
              double *alpha1,
              double *kappa1,
              gsl_vector *gamma,               
              gsl_vector *survTime1,
              gsl_vector *survEvent1,
              gsl_matrix *survCov1,
              gsl_vector *accept_beta1)
{
    double LP, D1, D2, logLH;
    double LP_prop, D1_prop, D2_prop, logLH_prop;
    double beta_prop_me, beta_prop_var, temp_prop;
    double beta_prop_me_prop, beta_prop_var_prop;
    double logProp_IniToProp, logProp_PropToIni;
    double logR;
    int u;
    
    int n = survTime1 -> size;
    int p = survCov1 -> size2;
    int i, j;
    
    gsl_vector *beta_prop = gsl_vector_calloc(p);
    
    
    j = (int) runif(0, p);
    
    logLH = 0; D1 = 0; D2 = 0;
    logLH_prop = 0; D1_prop = 0; D2_prop = 0;
    
    for(i = 0; i < n; i++)
    {
        gsl_vector_view Xi = gsl_matrix_row(survCov1, i);
        gsl_blas_ddot(&Xi.vector, beta1, &LP);
        if(gsl_vector_get(survEvent1, i) == 1)
        {
            logLH   += LP;
            D1      += gsl_matrix_get(survCov1, i, j);
        }
        logLH   += -gsl_vector_get(gamma, i) * (*kappa1) * pow(gsl_vector_get(survTime1, i), *alpha1) * exp(LP);
        D1      += -gsl_vector_get(gamma, i) * (*kappa1) * pow(gsl_vector_get(survTime1, i), *alpha1) * exp(LP) * gsl_matrix_get(survCov1, i, j);
        D2      += -gsl_vector_get(gamma, i) * (*kappa1) * pow(gsl_vector_get(survTime1, i), *alpha1) * exp(LP) * pow(gsl_matrix_get(survCov1, i, j), 2);
    }
    
    beta_prop_me    = gsl_vector_get(beta1, j) - D1/D2;
    beta_prop_var   = - pow(2.4, 2)/D2;
    
    temp_prop = rnorm(beta_prop_me, sqrt(beta_prop_var));
    
    gsl_vector_memcpy(beta_prop, beta1);
    gsl_vector_set(beta_prop, j, temp_prop);
    
    for(i = 0; i < n; i++)
    {
        gsl_vector_view Xi = gsl_matrix_row(survCov1, i);
        gsl_blas_ddot(&Xi.vector, beta_prop, &LP_prop);
        if(gsl_vector_get(survEvent1, i) == 1)
        {
            logLH_prop   += LP_prop;
            D1_prop      += gsl_matrix_get(survCov1, i, j);
        }
        logLH_prop   += -gsl_vector_get(gamma, i) * (*kappa1) * pow(gsl_vector_get(survTime1, i), *alpha1) * exp(LP_prop);
        D1_prop      += -gsl_vector_get(gamma, i) * (*kappa1) * pow(gsl_vector_get(survTime1, i), *alpha1) * exp(LP_prop) * gsl_matrix_get(survCov1, i, j);
        D2_prop      += -gsl_vector_get(gamma, i) * (*kappa1) * pow(gsl_vector_get(survTime1, i), *alpha1) * exp(LP_prop) * pow(gsl_matrix_get(survCov1, i, j), 2);
    }
    
    beta_prop_me_prop   = temp_prop - D1_prop/D2_prop;
    beta_prop_var_prop  = - pow(2.4, 2)/D2_prop;
    
    logProp_IniToProp = dnorm(temp_prop, beta_prop_me, sqrt(beta_prop_var), 1);
    logProp_PropToIni = dnorm(gsl_vector_get(beta1, j), beta_prop_me_prop, sqrt(beta_prop_var_prop), 1);
    
    logR = logLH_prop - logLH + logProp_PropToIni - logProp_IniToProp;
    
    u = log(runif(0, 1)) < logR;
    
    if(u == 1)
    {
        gsl_vector_set(beta1, j, temp_prop);
        gsl_vector_set(accept_beta1, j, (gsl_vector_get(accept_beta1, j) + u));
    }
    
    gsl_vector_free(beta_prop);
    
    return;
}










/* updating regression parameter: beta2 */

/**/

void BweibScrSM_updateRP2(gsl_vector *beta2,
               double *alpha2,
               double *kappa2,
               gsl_vector *gamma,
               gsl_vector *survTime1,
               gsl_vector *case01,
               gsl_matrix *survCov2,
               gsl_vector *accept_beta2)
{
    double LP, D1, D2, logLH;
    double LP_prop, D1_prop, D2_prop, logLH_prop;
    double beta_prop_me, beta_prop_var, temp_prop;
    double beta_prop_me_prop, beta_prop_var_prop;
    double logProp_IniToProp, logProp_PropToIni;
    double logR;
    int u;
    
    int n = survTime1 -> size;
    int p = survCov2 -> size2;
    int i, j;
    
    gsl_vector *beta_prop = gsl_vector_calloc(p);
    
    
    j = (int) runif(0, p);
    
    logLH = 0; D1 = 0; D2 = 0;
    logLH_prop = 0; D1_prop = 0; D2_prop = 0;
    
    for(i = 0; i < n; i++)
    {
        gsl_vector_view Xi = gsl_matrix_row(survCov2, i);
        gsl_blas_ddot(&Xi.vector, beta2, &LP);
        if(gsl_vector_get(case01, i) == 1)
        {
            logLH   += LP;
            D1      += gsl_matrix_get(survCov2, i, j);
        }
        logLH   += -gsl_vector_get(gamma, i) * (*kappa2) * pow(gsl_vector_get(survTime1, i), *alpha2) * exp(LP);
        D1      += -gsl_vector_get(gamma, i) * (*kappa2) * pow(gsl_vector_get(survTime1, i), *alpha2) * exp(LP) * gsl_matrix_get(survCov2, i, j);
        D2      += -gsl_vector_get(gamma, i) * (*kappa2) * pow(gsl_vector_get(survTime1, i), *alpha2) * exp(LP) * pow(gsl_matrix_get(survCov2, i, j), 2);
    }
    
    beta_prop_me    = gsl_vector_get(beta2, j) - D1/D2;
    beta_prop_var   = - pow(2.4, 2)/D2;
    
    temp_prop = rnorm(beta_prop_me, sqrt(beta_prop_var));
    
    gsl_vector_memcpy(beta_prop, beta2);
    gsl_vector_set(beta_prop, j, temp_prop);
    
    for(i = 0; i < n; i++)
    {
        gsl_vector_view Xi = gsl_matrix_row(survCov2, i);
        gsl_blas_ddot(&Xi.vector, beta_prop, &LP_prop);
        if(gsl_vector_get(case01, i) == 1)
        {
            logLH_prop   += LP_prop;
            D1_prop      += gsl_matrix_get(survCov2, i, j);
        }
        logLH_prop   += -gsl_vector_get(gamma, i) * (*kappa2) * pow(gsl_vector_get(survTime1, i), *alpha2) * exp(LP_prop);
        D1_prop      += -gsl_vector_get(gamma, i) * (*kappa2) * pow(gsl_vector_get(survTime1, i), *alpha2) * exp(LP_prop) * gsl_matrix_get(survCov2, i, j);
        D2_prop      += -gsl_vector_get(gamma, i) * (*kappa2) * pow(gsl_vector_get(survTime1, i), *alpha2) * exp(LP_prop) * pow(gsl_matrix_get(survCov2, i, j), 2);
    }
    
    beta_prop_me_prop   = temp_prop - D1_prop/D2_prop;
    beta_prop_var_prop  = - pow(2.4, 2)/D2_prop;
    
    logProp_IniToProp = dnorm(temp_prop, beta_prop_me, sqrt(beta_prop_var), 1);
    logProp_PropToIni = dnorm(gsl_vector_get(beta2, j), beta_prop_me_prop, sqrt(beta_prop_var_prop), 1);
    
    logR = logLH_prop - logLH + logProp_PropToIni - logProp_IniToProp;
    
    u = log(runif(0, 1)) < logR;
    
    if(u == 1)
    {
        gsl_vector_set(beta2, j, temp_prop);
        gsl_vector_set(accept_beta2, j, (gsl_vector_get(accept_beta2, j) + u));
    }
    
    gsl_vector_free(beta_prop);
    
    return;
}








/* updating regression parameter: beta3 */

/**/

void BweibScrSM_updateRP3(gsl_vector *beta3,
               double *alpha3,
               double *kappa3,
               gsl_vector *gamma,
               gsl_vector *yStar,
               gsl_vector *case11,
               gsl_matrix *survCov3,
               gsl_vector *accept_beta3)
{
    double LP, D1, D2, logLH;
    double LP_prop, D1_prop, D2_prop, logLH_prop;
    double beta_prop_me, beta_prop_var, temp_prop;
    double beta_prop_me_prop, beta_prop_var_prop;
    double logProp_IniToProp, logProp_PropToIni;
    double logR;
    int u;
    
    int n = yStar -> size;
    int p = survCov3 -> size2;
    int i, j;
    
    gsl_vector *beta_prop = gsl_vector_calloc(p);
    
    
    j = (int) runif(0, p);
    
    logLH = 0; D1 = 0; D2 = 0;
    logLH_prop = 0; D1_prop = 0; D2_prop = 0;
    
    for(i = 0; i < n; i++)
    {
        gsl_vector_view Xi = gsl_matrix_row(survCov3, i);
        gsl_blas_ddot(&Xi.vector, beta3, &LP);
        if(gsl_vector_get(case11, i) == 1)
        {
            logLH   += LP;
            D1      += gsl_matrix_get(survCov3, i, j);
        }
        logLH   += -gsl_vector_get(gamma, i) * (*kappa3) * pow(gsl_vector_get(yStar, i), *alpha3) * exp(LP);
        D1      += -gsl_vector_get(gamma, i) * (*kappa3) * pow(gsl_vector_get(yStar, i), *alpha3) * exp(LP) * gsl_matrix_get(survCov3, i, j);
        D2      += -gsl_vector_get(gamma, i) * (*kappa3) * pow(gsl_vector_get(yStar, i), *alpha3) * exp(LP) * pow(gsl_matrix_get(survCov3, i, j), 2);
    }
    
    beta_prop_me    = gsl_vector_get(beta3, j) - D1/D2;
    beta_prop_var   = - pow(2.4, 2)/D2;
    
    temp_prop = rnorm(beta_prop_me, sqrt(beta_prop_var));
    
    gsl_vector_memcpy(beta_prop, beta3);
    gsl_vector_set(beta_prop, j, temp_prop);
    
    for(i = 0; i < n; i++)
    {
        gsl_vector_view Xi = gsl_matrix_row(survCov3, i);
        gsl_blas_ddot(&Xi.vector, beta_prop, &LP_prop);
        if(gsl_vector_get(case11, i) == 1)
        {
            logLH_prop   += LP_prop;
            D1_prop      += gsl_matrix_get(survCov3, i, j);
        }
        logLH_prop   += -gsl_vector_get(gamma, i) * (*kappa3) * pow(gsl_vector_get(yStar, i), *alpha3) * exp(LP_prop);
        D1_prop      += -gsl_vector_get(gamma, i) * (*kappa3) * pow(gsl_vector_get(yStar, i), *alpha3) * exp(LP_prop) * gsl_matrix_get(survCov3, i, j);
        D2_prop      += -gsl_vector_get(gamma, i) * (*kappa3) * pow(gsl_vector_get(yStar, i), *alpha3) * exp(LP_prop) * pow(gsl_matrix_get(survCov3, i, j), 2);
    }
    
    beta_prop_me_prop   = temp_prop - D1_prop/D2_prop;
    beta_prop_var_prop  = - pow(2.4, 2)/D2_prop;
    
    logProp_IniToProp = dnorm(temp_prop, beta_prop_me, sqrt(beta_prop_var), 1);
    logProp_PropToIni = dnorm(gsl_vector_get(beta3, j), beta_prop_me_prop, sqrt(beta_prop_var_prop), 1);
    
    logR = logLH_prop - logLH + logProp_PropToIni - logProp_IniToProp;
    
    u = log(runif(0, 1)) < logR;
    
    if(u == 1)
    {
        gsl_vector_set(beta3, j, temp_prop);
        gsl_vector_set(accept_beta3, j, (gsl_vector_get(accept_beta3, j) + u));
    }
    
    gsl_vector_free(beta_prop);
    
    return;
}













/* updating shape parameter: alpha1 */

/* use the random walk proposal */

/**/

void BweibScrSM_updateSC1(gsl_vector *beta1,
               double *alpha1,
               double *kappa1,
               gsl_vector *gamma,
               gsl_vector *survTime1,
               gsl_vector *survEvent1,
               gsl_matrix *survCov1,
               double mhProp_alpha1_var,
               double a1,
               double b1,
               int *accept_alpha1)
{
    double LP, logLH, logLH_prop;
    double temp_prop;
    double logPrior, logPrior_prop;
    double logProp_IniToProp, logProp_PropToIni;
    double logR;
    int u;
    int n = survTime1 -> size;
    int i;
    
    logLH = 0; logLH_prop = 0;
    temp_prop = rgamma(pow(*alpha1, 2)/mhProp_alpha1_var, mhProp_alpha1_var/(*alpha1));
    
    for(i = 0; i < n; i++)
    {
        gsl_vector_view Xi = gsl_matrix_row(survCov1, i);
        gsl_blas_ddot(&Xi.vector, beta1, &LP);
        if(gsl_vector_get(survEvent1, i) == 1)
        {
            logLH       += log(*alpha1) + (*alpha1 - 1)*log(gsl_vector_get(survTime1, i));
            logLH_prop  += log(temp_prop) + (temp_prop - 1)*log(gsl_vector_get(survTime1, i));
        }
        logLH       += - gsl_vector_get(gamma, i) * (*kappa1) * pow(gsl_vector_get(survTime1, i), *alpha1) * exp(LP);
        logLH_prop  += - gsl_vector_get(gamma, i) * (*kappa1) * pow(gsl_vector_get(survTime1, i), temp_prop) * exp(LP);
    }
    
    logPrior        = dgamma(*alpha1, a1, 1/b1, 1);
    logPrior_prop   = dgamma(temp_prop, a1, 1/b1, 1);
    
    logProp_PropToIni = dgamma(*alpha1, pow(temp_prop, 2)/mhProp_alpha1_var, mhProp_alpha1_var/(temp_prop), 1);
    logProp_IniToProp = dgamma(temp_prop, pow(*alpha1, 2)/mhProp_alpha1_var, mhProp_alpha1_var/(*alpha1), 1);
    
    logR = logLH_prop - logLH + logPrior_prop - logPrior + logProp_PropToIni - logProp_IniToProp;
    
    u = log(runif(0, 1)) <logR;
    
    if(u == 1)
    {
        *alpha1 = temp_prop;
        *accept_alpha1 += u;
    }
    
    return;
    
}













/* updating shape parameter: alpha2 */

/* use the random walk proposal */

/**/

void BweibScrSM_updateSC2(gsl_vector *beta2,
               double *alpha2,
               double *kappa2,
               gsl_vector *gamma,
               gsl_vector *survTime1,
               gsl_vector *survTime2,
               gsl_vector *case01,
               gsl_matrix *survCov2,
               double mhProp_alpha2_var,
               double a2,
               double b2,
               int *accept_alpha2)
{
    double LP, logLH, logLH_prop;
    double temp_prop;
    double logPrior, logPrior_prop;
    double logProp_IniToProp, logProp_PropToIni;
    double logR;
    int u;
    int n = survTime1 -> size;
    int i;
    
    logLH = 0; logLH_prop = 0;
    temp_prop = rgamma(pow(*alpha2, 2)/mhProp_alpha2_var, mhProp_alpha2_var/(*alpha2));
    
    for(i = 0; i < n; i++)
    {
        gsl_vector_view Xi = gsl_matrix_row(survCov2, i);
        gsl_blas_ddot(&Xi.vector, beta2, &LP);
        if(gsl_vector_get(case01, i) == 1)
        {
            logLH       += log(*alpha2) + (*alpha2 - 1)*log(gsl_vector_get(survTime2, i));
            logLH_prop  += log(temp_prop) + (temp_prop - 1)*log(gsl_vector_get(survTime2, i));
        }
        logLH       += - gsl_vector_get(gamma, i) * (*kappa2) * pow(gsl_vector_get(survTime1, i), *alpha2) * exp(LP);
        logLH_prop  += - gsl_vector_get(gamma, i) * (*kappa2) * pow(gsl_vector_get(survTime1, i), temp_prop) * exp(LP);
    }
    
    logPrior        = dgamma(*alpha2, a2, 1/b2, 1);
    logPrior_prop   = dgamma(temp_prop, a2, 1/b2, 1);
    
    logProp_PropToIni = dgamma(*alpha2, pow(temp_prop, 2)/mhProp_alpha2_var, mhProp_alpha2_var/(temp_prop), 1);
    logProp_IniToProp = dgamma(temp_prop, pow(*alpha2, 2)/mhProp_alpha2_var, mhProp_alpha2_var/(*alpha2), 1);
    
    logR = logLH_prop - logLH + logPrior_prop - logPrior + logProp_PropToIni - logProp_IniToProp;
    
    u = log(runif(0, 1)) <logR;
    
    if(u == 1)
    {
        *alpha2 = temp_prop;
        *accept_alpha2 += u;
    }
    
    return;
    
}











/* updating shape parameter: alpha3 */

/* use the random walk proposal */

/**/

void BweibScrSM_updateSC3(gsl_vector *beta3,
               double *alpha3,
               double *kappa3,
               gsl_vector *gamma,
               gsl_vector *yStar,
               gsl_vector *case11,
               gsl_matrix *survCov3,
               double mhProp_alpha3_var,
               double a3,
               double b3,
               int *accept_alpha3)
{
    double LP, logLH, logLH_prop;
    double temp_prop;
    double logPrior, logPrior_prop;
    double logProp_IniToProp, logProp_PropToIni;
    double logR;
    int u;
    int n = yStar -> size;
    int i;
    
    logLH = 0; logLH_prop = 0;
    temp_prop = rgamma(pow(*alpha3, 2)/mhProp_alpha3_var, mhProp_alpha3_var/(*alpha3));
    
    for(i = 0; i < n; i++)
    {
        gsl_vector_view Xi = gsl_matrix_row(survCov3, i);
        gsl_blas_ddot(&Xi.vector, beta3, &LP);
        if(gsl_vector_get(case11, i) == 1 && gsl_vector_get(yStar, i) != 0)
        {
            logLH       += log(*alpha3) + (*alpha3 - 1)*log(gsl_vector_get(yStar, i));
            logLH_prop  += log(temp_prop) + (temp_prop - 1)*log(gsl_vector_get(yStar, i));
        }
        logLH       += - gsl_vector_get(gamma, i) * (*kappa3) * pow(gsl_vector_get(yStar, i), *alpha3) * exp(LP);
        logLH_prop  += - gsl_vector_get(gamma, i) * (*kappa3) * pow(gsl_vector_get(yStar, i), temp_prop) * exp(LP);
    }
    
    logPrior        = dgamma(*alpha3, a3, 1/b3, 1);
    logPrior_prop   = dgamma(temp_prop, a3, 1/b3, 1);
    
    logProp_PropToIni = dgamma(*alpha3, pow(temp_prop, 2)/mhProp_alpha3_var, mhProp_alpha3_var/(temp_prop), 1);
    logProp_IniToProp = dgamma(temp_prop, pow(*alpha3, 2)/mhProp_alpha3_var, mhProp_alpha3_var/(*alpha3), 1);
    
    logR = logLH_prop - logLH + logPrior_prop - logPrior + logProp_PropToIni - logProp_IniToProp;
    
    u = log(runif(0, 1)) <logR;
    
    if(u == 1)
    {
        *alpha3 = temp_prop;
        *accept_alpha3 += u;
    }
    
    return;
    
}








/* updating shape parameter: kappa1 */

/**/

void BweibScrSM_updateSH1(gsl_vector *beta1,
              double *alpha1,
              double *kappa1,
              gsl_vector *gamma,
              gsl_vector *survTime1,
              gsl_vector *survEvent1,
              gsl_matrix *survCov1,
              double c1,
              double d1)
{
    int n = survTime1 -> size;
    int i;
    double LP, Kappa1_shape, Kappa1_rate, Kappa1_scale;
    
    gsl_vector *ones = gsl_vector_calloc(n);
    gsl_vector_set_all(ones, 1);
    
    gsl_blas_ddot(ones, survEvent1, &Kappa1_shape);
    Kappa1_shape += c1;
    
    gsl_vector_free(ones);
    
    Kappa1_rate = 0;
    
    for(i = 0; i< n; i++)
    {
        gsl_vector_view Xi = gsl_matrix_row(survCov1, i);
        gsl_blas_ddot(&Xi.vector, beta1 ,&LP);
        
        Kappa1_rate += gsl_vector_get(gamma, i) * pow(gsl_vector_get(survTime1, i), *alpha1) * exp(LP);
    }
    Kappa1_rate += d1;
    Kappa1_scale = 1/Kappa1_rate;
    
    *kappa1 = rgamma(Kappa1_shape, Kappa1_scale);
    return;
}







/* updating shape parameter: kappa2 */

/**/

void BweibScrSM_updateSH2(gsl_vector *beta2,
               double *alpha2,
               double *kappa2,
               gsl_vector *gamma,
               gsl_vector *survTime1,
               gsl_vector *case01,
               gsl_matrix *survCov2,
               double c2,
               double d2)
{
    int n = survTime1 -> size;
    int i;
    double LP, Kappa2_shape, Kappa2_rate, Kappa2_scale;
    
    gsl_vector *ones = gsl_vector_calloc(n);
    gsl_vector_set_all(ones, 1);
    
    gsl_blas_ddot(ones, case01, &Kappa2_shape);
    Kappa2_shape += c2;
    
    gsl_vector_free(ones);
    
    Kappa2_rate = 0;
    
    for(i = 0; i< n; i++)
    {
        gsl_vector_view Xi = gsl_matrix_row(survCov2, i);
        gsl_blas_ddot(&Xi.vector, beta2 ,&LP);
        
        Kappa2_rate += gsl_vector_get(gamma, i) * pow(gsl_vector_get(survTime1, i), *alpha2) * exp(LP);
    }
    Kappa2_rate += d2;
    Kappa2_scale = 1/Kappa2_rate;
    
    *kappa2 = rgamma(Kappa2_shape, Kappa2_scale);
    return;
}



/* updating shape parameter: kappa3 */

/**/

void BweibScrSM_updateSH3(gsl_vector *beta3,
               double *alpha3,
               double *kappa3,
               gsl_vector *gamma,
               gsl_vector *yStar,
               gsl_vector *case11,
               gsl_matrix *survCov3,
               double c3,
               double d3)
{
    int n = yStar -> size;
    int i;
    double LP, Kappa3_shape, Kappa3_rate, Kappa3_scale;
    
    gsl_vector *ones = gsl_vector_calloc(n);
    gsl_vector_set_all(ones, 1);
    
    gsl_blas_ddot(ones, case11, &Kappa3_shape);
    Kappa3_shape += c3;
    
    gsl_vector_free(ones);
    
    Kappa3_rate = 0;
    
    for(i = 0; i< n; i++)
    {
        gsl_vector_view Xi = gsl_matrix_row(survCov3, i);
        gsl_blas_ddot(&Xi.vector, beta3 ,&LP);
        
        Kappa3_rate += gsl_vector_get(gamma, i) * pow(gsl_vector_get(yStar, i), *alpha3) * exp(LP);
    }
    Kappa3_rate += d3;
    Kappa3_scale = 1/Kappa3_rate;
    
    *kappa3 = rgamma(Kappa3_shape, Kappa3_scale);
    return;
}













/* updating variance parameter: theta */

/* use the random walk proposal */

/**/

void BweibScrSM_updateDP(gsl_vector *gamma,
              double *theta,
              double mhProp_theta_var,
              double psi,
              double omega,
              int *accept_theta)
{
    double logPost, logPost_prop;
    double logProp_IniToProp, logProp_PropToIni;
    double logR;
    int n = gamma -> size;
    int u;
    int i;
    double xi = 1/(*theta);
    double temp_prop;
    
    logPost = 0; logPost_prop = 0;
    temp_prop = rgamma(pow(xi, 2)/mhProp_theta_var, mhProp_theta_var/(xi));
    
    for(i = 0; i < n; i++)
    {
        logPost      += xi * (log(gsl_vector_get(gamma, i)) - gsl_vector_get(gamma, i));
        logPost_prop += temp_prop * (log(gsl_vector_get(gamma, i)) - gsl_vector_get(gamma, i));
    }
    
    logPost         += (n * xi + psi - 1)*log(xi) - xi*omega - n*lgamma(xi);
    logPost_prop    += (n * temp_prop + psi - 1)*log(temp_prop) - temp_prop*omega - n*lgamma(temp_prop);
    
    logProp_PropToIni = dgamma(xi, pow(temp_prop, 2)/mhProp_theta_var, mhProp_theta_var/(temp_prop), 1);
    logProp_IniToProp = dgamma(temp_prop, pow(xi, 2)/mhProp_theta_var, mhProp_theta_var/(xi), 1);
    
    logR = logPost_prop - logPost + logProp_PropToIni - logProp_IniToProp;
    
    u = log(runif(0, 1)) <logR;
    
    if(u == 1)
    {
        *theta = 1/temp_prop;
        *accept_theta += u;
    }
    
    return;
}












/* updating frailty parameter: gamma */

/**/


 
void BweibScrSM_updateFP(gsl_vector *beta1,
              gsl_vector *beta2,
              gsl_vector *beta3,
              double alpha1,
              double alpha2,
              double alpha3,
              double kappa1,
              double kappa2,
              double kappa3,
              double theta,
              gsl_vector *gamma,
              gsl_vector *survTime1,
              gsl_vector *yStar,
              gsl_vector *survEvent1,
              gsl_vector *survEvent2,
              gsl_matrix *survCov1,
              gsl_matrix *survCov2,
              gsl_matrix *survCov3)
{
    int n = survTime1 -> size;
    int i;
    double gamma_shape, gamma_rate, gamma_scale;

    for(i = 0; i< n; i++)
    {
        gamma_shape = gsl_vector_get(survEvent1, i) + gsl_vector_get(survEvent2, i) + 1/theta;
        gamma_rate  = BweibScrSM_wFunc(i, beta1, beta2, beta3, alpha1, alpha2, alpha3, kappa1, kappa2, kappa3, survTime1, yStar, survCov1, survCov2, survCov3) + 1/theta;
        
        gamma_scale = 1/gamma_rate;
        gsl_vector_set(gamma, i, rgamma(gamma_shape, gamma_scale));
        
    }
    return;
}











