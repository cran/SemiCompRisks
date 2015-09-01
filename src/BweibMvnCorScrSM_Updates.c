#include <stdio.h>
#include <math.h>

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"


#include "R.h"
#include "Rmath.h"

#include "BweibMvnCorScrSM.h"






/* updating regression parameter: beta1 */

/**/

void BweibMvnCorScrSM_updateRP1(gsl_vector *beta1,
                           double *alpha1,
                           double *kappa1,
                           gsl_vector *gamma,
                           gsl_vector *V1,
                           gsl_vector *survTime1,
                           gsl_vector *survEvent1,
                           gsl_vector *cluster,
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
    int i, j, jj;
    
    gsl_vector *beta_prop = gsl_vector_calloc(p);
    
    
    j = (int) runif(0, p);
    
    
    /*
    j = 7;
    */
    
    logLH = 0; D1 = 0; D2 = 0;
    logLH_prop = 0; D1_prop = 0; D2_prop = 0;
    

    
    for(i = 0; i < n; i++)
    {
        gsl_vector_view Xi = gsl_matrix_row(survCov1, i);
        gsl_blas_ddot(&Xi.vector, beta1, &LP);
        jj = (int) gsl_vector_get(cluster, i) - 1;
        if(gsl_vector_get(survEvent1, i) == 1)
        {
            logLH   += LP;
            D1      += gsl_matrix_get(survCov1, i, j);
        }
        logLH   += -gsl_vector_get(gamma, i) * (*kappa1) * pow(gsl_vector_get(survTime1, i), *alpha1) * exp(LP + gsl_vector_get(V1, jj));
        D1      += -gsl_vector_get(gamma, i) * (*kappa1) * pow(gsl_vector_get(survTime1, i), *alpha1) * exp(LP + gsl_vector_get(V1, jj)) * gsl_matrix_get(survCov1, i, j);
        D2      += -gsl_vector_get(gamma, i) * (*kappa1) * pow(gsl_vector_get(survTime1, i), *alpha1) * exp(LP + gsl_vector_get(V1, jj)) * pow(gsl_matrix_get(survCov1, i, j), 2);
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
        jj = (int) gsl_vector_get(cluster, i) - 1;
        if(gsl_vector_get(survEvent1, i) == 1)
        {
            logLH_prop   += LP_prop;
            D1_prop      += gsl_matrix_get(survCov1, i, j);
        }
        logLH_prop   += -gsl_vector_get(gamma, i) * (*kappa1) * pow(gsl_vector_get(survTime1, i), *alpha1) * exp(LP_prop + gsl_vector_get(V1, jj));
        D1_prop      += -gsl_vector_get(gamma, i) * (*kappa1) * pow(gsl_vector_get(survTime1, i), *alpha1) * exp(LP_prop + gsl_vector_get(V1, jj)) * gsl_matrix_get(survCov1, i, j);
        D2_prop      += -gsl_vector_get(gamma, i) * (*kappa1) * pow(gsl_vector_get(survTime1, i), *alpha1) * exp(LP_prop + gsl_vector_get(V1, jj)) * pow(gsl_matrix_get(survCov1, i, j), 2);
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

void BweibMvnCorScrSM_updateRP2(gsl_vector *beta2,
                        double *alpha2,
                        double *kappa2,
                        double *nu2,
                        gsl_vector *gamma,
                        gsl_vector *V2,
                        gsl_vector *survTime1,
                        gsl_vector *case01,
                        gsl_vector *cluster,                        
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
    int i, j, jj;
    
    gsl_vector *beta_prop = gsl_vector_calloc(p);
    
    
    j = (int) runif(0, p);
    
    logLH = 0; D1 = 0; D2 = 0;
    logLH_prop = 0; D1_prop = 0; D2_prop = 0;
    
    for(i = 0; i < n; i++)
    {
        gsl_vector_view Xi = gsl_matrix_row(survCov2, i);
        gsl_blas_ddot(&Xi.vector, beta2, &LP);
        jj = (int) gsl_vector_get(cluster, i) - 1;        
        if(gsl_vector_get(case01, i) == 1)
        {
            logLH   += LP;
            D1      += gsl_matrix_get(survCov2, i, j);
        }
        logLH   += - pow(gsl_vector_get(gamma, i), *nu2) * (*kappa2) * pow(gsl_vector_get(survTime1, i), *alpha2) * exp(LP + gsl_vector_get(V2, jj));
        D1      += - pow(gsl_vector_get(gamma, i), *nu2) * (*kappa2) * pow(gsl_vector_get(survTime1, i), *alpha2) * exp(LP + gsl_vector_get(V2, jj)) * gsl_matrix_get(survCov2, i, j);
        D2      += - pow(gsl_vector_get(gamma, i), *nu2) * (*kappa2) * pow(gsl_vector_get(survTime1, i), *alpha2) * exp(LP + gsl_vector_get(V2, jj)) * pow(gsl_matrix_get(survCov2, i, j), 2);
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
        jj = (int) gsl_vector_get(cluster, i) - 1;        
        if(gsl_vector_get(case01, i) == 1)
        {
            logLH_prop   += LP_prop;
            D1_prop      += gsl_matrix_get(survCov2, i, j);
        }
        logLH_prop   += - pow(gsl_vector_get(gamma, i), *nu2) * (*kappa2) * pow(gsl_vector_get(survTime1, i), *alpha2) * exp(LP_prop + gsl_vector_get(V2, jj));
        D1_prop      += - pow(gsl_vector_get(gamma, i), *nu2) * (*kappa2) * pow(gsl_vector_get(survTime1, i), *alpha2) * exp(LP_prop + gsl_vector_get(V2, jj)) * gsl_matrix_get(survCov2, i, j);
        D2_prop      += - pow(gsl_vector_get(gamma, i), *nu2) * (*kappa2) * pow(gsl_vector_get(survTime1, i), *alpha2) * exp(LP_prop + gsl_vector_get(V2, jj)) * pow(gsl_matrix_get(survCov2, i, j), 2);
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

void BweibMvnCorScrSM_updateRP3(gsl_vector *beta3,
                        double *alpha3,
                        double *kappa3,
                        double *nu3,
                        gsl_vector *gamma,
                        gsl_vector *V3,
                        gsl_vector *yStar,
                        gsl_vector *case11,
                        gsl_vector *cluster,
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
    int i, j, jj;
    
    gsl_vector *beta_prop = gsl_vector_calloc(p);
    
    
    j = (int) runif(0, p);
    
    logLH = 0; D1 = 0; D2 = 0;
    logLH_prop = 0; D1_prop = 0; D2_prop = 0;
    
    for(i = 0; i < n; i++)
    {
        gsl_vector_view Xi = gsl_matrix_row(survCov3, i);
        gsl_blas_ddot(&Xi.vector, beta3, &LP);
        jj = (int) gsl_vector_get(cluster, i) - 1;         
        if(gsl_vector_get(case11, i) == 1)
        {
            logLH   += LP;
            D1      += gsl_matrix_get(survCov3, i, j);
        }
        logLH   += - pow(gsl_vector_get(gamma, i), *nu3) * (*kappa3) * pow(gsl_vector_get(yStar, i), *alpha3) * exp(LP + gsl_vector_get(V3, jj));
        D1      += - pow(gsl_vector_get(gamma, i), *nu3) * (*kappa3) * pow(gsl_vector_get(yStar, i), *alpha3) * exp(LP + gsl_vector_get(V3, jj)) * gsl_matrix_get(survCov3, i, j);
        D2      += - pow(gsl_vector_get(gamma, i), *nu3) * (*kappa3) * pow(gsl_vector_get(yStar, i), *alpha3) * exp(LP + gsl_vector_get(V3, jj)) * pow(gsl_matrix_get(survCov3, i, j), 2);
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
        jj = (int) gsl_vector_get(cluster, i) - 1;         
        if(gsl_vector_get(case11, i) == 1)
        {
            logLH_prop   += LP_prop;
            D1_prop      += gsl_matrix_get(survCov3, i, j);
        }
        logLH_prop   += - pow(gsl_vector_get(gamma, i), *nu3) * (*kappa3) * pow(gsl_vector_get(yStar, i), *alpha3) * exp(LP_prop + gsl_vector_get(V3, jj));
        D1_prop      += - pow(gsl_vector_get(gamma, i), *nu3) * (*kappa3) * pow(gsl_vector_get(yStar, i), *alpha3) * exp(LP_prop + gsl_vector_get(V3, jj)) * gsl_matrix_get(survCov3, i, j);
        D2_prop      += - pow(gsl_vector_get(gamma, i), *nu3) * (*kappa3) * pow(gsl_vector_get(yStar, i), *alpha3) * exp(LP_prop + gsl_vector_get(V3, jj)) * pow(gsl_matrix_get(survCov3, i, j), 2);
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

void BweibMvnCorScrSM_updateSC1(gsl_vector *beta1,
                        double *alpha1,
                        double *kappa1,
                        gsl_vector *gamma,
                        gsl_vector *V1,
                        gsl_vector *survTime1,
                        gsl_vector *survEvent1,
                        gsl_vector *cluster,
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
    int i, jj;
    
    logLH = 0; logLH_prop = 0;
    temp_prop = rgamma(pow(*alpha1, 2)/mhProp_alpha1_var, mhProp_alpha1_var/(*alpha1));
    
    for(i = 0; i < n; i++)
    {
        gsl_vector_view Xi = gsl_matrix_row(survCov1, i);
        gsl_blas_ddot(&Xi.vector, beta1, &LP);
        jj = (int) gsl_vector_get(cluster, i) - 1;        
        if(gsl_vector_get(survEvent1, i) == 1)
        {
            logLH       += log(*alpha1) + (*alpha1 - 1)*log(gsl_vector_get(survTime1, i));
            logLH_prop  += log(temp_prop) + (temp_prop - 1)*log(gsl_vector_get(survTime1, i));
        }
        logLH       += - gsl_vector_get(gamma, i) * (*kappa1) * pow(gsl_vector_get(survTime1, i), *alpha1) * exp(LP + gsl_vector_get(V1, jj));
        logLH_prop  += - gsl_vector_get(gamma, i) * (*kappa1) * pow(gsl_vector_get(survTime1, i), temp_prop) * exp(LP + gsl_vector_get(V1, jj));
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








/* updating shape parameter: alpha1 */

/* use the random walk proposal with log transformation*/

/**/

void BweibMvnCorScrSM_updateSC1_rw2(gsl_vector *beta1,
                        double *alpha1,
                        double *kappa1,
                        gsl_vector *gamma,
                        gsl_vector *V1,
                        gsl_vector *survTime1,
                        gsl_vector *survEvent1,
                        gsl_vector *cluster,
                        gsl_matrix *survCov1,
                        double mhProp_alpha1_var,
                        double a1,
                        double b1,
                        int *accept_alpha1)
{
    double z1, LP, logLH, logLH_prop;
    double temp_prop;
    double logPrior, logPrior_prop;
    double logProp_IniToProp, logProp_PropToIni;
    double logR;
    int u;
    int n = survTime1 -> size;
    int i, jj;
    
    z1 = log(*alpha1);
    
    logLH = 0; logLH_prop = 0;
    temp_prop = rnorm(z1, sqrt(mhProp_alpha1_var));
    
    for(i = 0; i < n; i++)
    {
        gsl_vector_view Xi = gsl_matrix_row(survCov1, i);
        gsl_blas_ddot(&Xi.vector, beta1, &LP);
        jj = (int) gsl_vector_get(cluster, i) - 1;
        if(gsl_vector_get(survEvent1, i) == 1)
        {
            logLH       += log(*alpha1) + (*alpha1 - 1)*log(gsl_vector_get(survTime1, i));
            logLH_prop  += temp_prop + (exp(temp_prop) - 1)*log(gsl_vector_get(survTime1, i));
        }
        logLH       += - gsl_vector_get(gamma, i) * (*kappa1) * pow(gsl_vector_get(survTime1, i), *alpha1) * exp(LP + gsl_vector_get(V1, jj));
        logLH_prop  += - gsl_vector_get(gamma, i) * (*kappa1) * pow(gsl_vector_get(survTime1, i), exp(temp_prop)) * exp(LP + gsl_vector_get(V1, jj));
    }
    
    logPrior        = a1*z1 - b1* (*alpha1);
    logPrior_prop   = a1*temp_prop - b1* (exp(temp_prop));
    
    logProp_PropToIni = dnorm(z1, temp_prop, sqrt(mhProp_alpha1_var), 1);
    logProp_IniToProp = dnorm(temp_prop, z1, sqrt(mhProp_alpha1_var), 1);
    
    logR = logLH_prop - logLH + logPrior_prop - logPrior + logProp_PropToIni - logProp_IniToProp;
    

    
    u = log(runif(0, 1)) <logR;
    
    
    if(u == 1)
    {
        *alpha1 = exp(temp_prop);
        *accept_alpha1 += u;
    }
    
    return;
    
}




/* updating shape parameter: alpha1 */

/* use AMCMC */

/**/

void BweibMvnCorScrSM_updateSC1_amc(gsl_vector *beta1,
                        double *alpha1,
                        double *kappa1,
                        gsl_vector *gamma,
                        gsl_vector *V1,
                        gsl_vector *survTime1,
                        gsl_vector *survEvent1,
                        gsl_vector *cluster,
                        gsl_matrix *survCov1,
                        double mhProp_alpha1_var,
                        double a1,
                        double b1,
                        int *accept_alpha1)
{
    double LP, logLH, logLH_prop, D1, D2, D1_prop, D2_prop;
    double temp_prop;
    double alpha_prop_me, alpha_prop_var;
    double alpha_prop_me_prop, alpha_prop_var_prop;
    double logPrior, logPrior_prop;
    double logProp_IniToProp, logProp_PropToIni;
    double logR;
    int u;
    int n = survTime1 -> size;
    int i, jj;
    
    logLH = 0; logLH_prop = 0;
    D1 = 0; D2 = 0;
    D1_prop = 0; D2_prop = 0;
    
    
    for(i = 0; i < n; i++)
    {
        gsl_vector_view Xi = gsl_matrix_row(survCov1, i);
        gsl_blas_ddot(&Xi.vector, beta1, &LP);
        jj = (int) gsl_vector_get(cluster, i) - 1;
        if(gsl_vector_get(survEvent1, i) == 1)
        {
            logLH       += log(*alpha1) + (*alpha1 - 1)*log(gsl_vector_get(survTime1, i));
            D1 += 1/(*alpha1) + log(gsl_vector_get(survTime1, i));
            D2 += -1/pow(*alpha1, 2);
        }
        logLH   += - gsl_vector_get(gamma, i) * (*kappa1) * pow(gsl_vector_get(survTime1, i), *alpha1) * exp(LP + gsl_vector_get(V1, jj));
        D1      += - gsl_vector_get(gamma, i) * (*kappa1) * pow(gsl_vector_get(survTime1, i), *alpha1) * exp(LP + gsl_vector_get(V1, jj))*log(*alpha1);
        D2      += - gsl_vector_get(gamma, i) * (*kappa1) * pow(gsl_vector_get(survTime1, i), *alpha1) * exp(LP + gsl_vector_get(V1, jj))*(pow(gsl_vector_get(survTime1, i), *alpha1)*pow(log(*alpha1), 2) + pow(gsl_vector_get(survTime1, i), *alpha1)/(*alpha1));
    }
    D1 += (a1 - 1)/(*alpha1) - b1;
    D2 += -(a1 - 1)/(pow(*alpha1, 2));
    
    alpha_prop_me = *alpha1 - D1/D2;
    alpha_prop_var = fabs(-pow(2.4, 2)/D2);
    
    
    temp_prop = rgamma(pow(alpha_prop_me, 2)/alpha_prop_var, alpha_prop_var/alpha_prop_me);
    
    for(i = 0; i < n; i++)
    {
        gsl_vector_view Xi = gsl_matrix_row(survCov1, i);
        gsl_blas_ddot(&Xi.vector, beta1, &LP);
        jj = (int) gsl_vector_get(cluster, i) - 1;
        if(gsl_vector_get(survEvent1, i) == 1)
        {
            logLH_prop       += log(temp_prop) + (temp_prop - 1)*log(gsl_vector_get(survTime1, i));
            D1_prop += 1/(temp_prop) + log(gsl_vector_get(survTime1, i));
            D2_prop += -1/pow(temp_prop, 2);
        }
        logLH_prop   += - gsl_vector_get(gamma, i) * (*kappa1) * pow(gsl_vector_get(survTime1, i), temp_prop) * exp(LP + gsl_vector_get(V1, jj));
        D1_prop      += - gsl_vector_get(gamma, i) * (*kappa1) * pow(gsl_vector_get(survTime1, i), temp_prop) * exp(LP + gsl_vector_get(V1, jj))*log(temp_prop);
        D2_prop      += - gsl_vector_get(gamma, i) * (*kappa1) * pow(gsl_vector_get(survTime1, i), temp_prop) * exp(LP + gsl_vector_get(V1, jj))*(pow(gsl_vector_get(survTime1, i), temp_prop)*pow(log(temp_prop), 2) + pow(gsl_vector_get(survTime1, i), temp_prop)/(temp_prop));
    }
    D1_prop += (a1 - 1)/(temp_prop) - b1;
    D2_prop += -(a1 - 1)/(pow(temp_prop, 2));
    
    alpha_prop_me_prop = temp_prop - D1_prop/D2_prop;
    alpha_prop_var_prop = fabs(-pow(2.4, 2)/D2_prop);
    
    logPrior        = dgamma(*alpha1, a1, 1/b1, 1);
    logPrior_prop   = dgamma(temp_prop, a1, 1/b1, 1);
    
    logProp_PropToIni = dgamma(*alpha1, pow(alpha_prop_me_prop, 2)/alpha_prop_var_prop, alpha_prop_var_prop/(alpha_prop_me_prop), 1);
    logProp_IniToProp = dgamma(temp_prop, pow(alpha_prop_me, 2)/alpha_prop_var, alpha_prop_var/(alpha_prop_me), 1);
    
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

void BweibMvnCorScrSM_updateSC2(gsl_vector *beta2,
                        double *alpha2,
                        double *kappa2,
                        double *nu2,
                        gsl_vector *gamma,
                        gsl_vector *V2,
                        gsl_vector *survTime1,
                        gsl_vector *survTime2,
                        gsl_vector *case01,
                        gsl_vector *cluster,
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
    int i, jj;
    
    logLH = 0; logLH_prop = 0;
    temp_prop = rgamma(pow(*alpha2, 2)/mhProp_alpha2_var, mhProp_alpha2_var/(*alpha2));
    
    for(i = 0; i < n; i++)
    {
        gsl_vector_view Xi = gsl_matrix_row(survCov2, i);
        gsl_blas_ddot(&Xi.vector, beta2, &LP);
        jj = (int) gsl_vector_get(cluster, i) - 1;         
        if(gsl_vector_get(case01, i) == 1)
        {
            logLH       += log(*alpha2) + (*alpha2 - 1)*log(gsl_vector_get(survTime2, i));
            logLH_prop  += log(temp_prop) + (temp_prop - 1)*log(gsl_vector_get(survTime2, i));
        }
        logLH       += - pow(gsl_vector_get(gamma, i), *nu2) * (*kappa2) * pow(gsl_vector_get(survTime1, i), *alpha2) * exp(LP + gsl_vector_get(V2, jj));
        logLH_prop  += - pow(gsl_vector_get(gamma, i), *nu2) * (*kappa2) * pow(gsl_vector_get(survTime1, i), temp_prop) * exp(LP + gsl_vector_get(V2, jj));
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








/* updating shape parameter: alpha2 */

/* use the random walk proposal with log transformation*/

/**/

void BweibMvnCorScrSM_updateSC2_rw2(gsl_vector *beta2,
                        double *alpha2,
                        double *kappa2,
                        double *nu2,
                        gsl_vector *gamma,
                        gsl_vector *V2,
                        gsl_vector *survTime1,
                        gsl_vector *survTime2,
                        gsl_vector *case01,
                        gsl_vector *cluster,
                        gsl_matrix *survCov2,
                        double mhProp_alpha2_var,
                        double a2,
                        double b2,
                        int *accept_alpha2)
{
    double z2, LP, logLH, logLH_prop;
    double temp_prop;
    double logPrior, logPrior_prop;
    double logProp_IniToProp, logProp_PropToIni;
    double logR;
    int u;
    int n = survTime1 -> size;
    int i, jj;
    
    z2 = log(*alpha2);
    
    logLH = 0; logLH_prop = 0;
    temp_prop = rnorm(z2, sqrt(mhProp_alpha2_var));
    
    for(i = 0; i < n; i++)
    {
        gsl_vector_view Xi = gsl_matrix_row(survCov2, i);
        gsl_blas_ddot(&Xi.vector, beta2, &LP);
        jj = (int) gsl_vector_get(cluster, i) - 1;
        if(gsl_vector_get(case01, i) == 1)
        {
            logLH       += log(*alpha2) + (*alpha2 - 1)*log(gsl_vector_get(survTime2, i));
            logLH_prop  += temp_prop + (exp(temp_prop) - 1)*log(gsl_vector_get(survTime2, i));
        }
        logLH       += - pow(gsl_vector_get(gamma, i), *nu2) * (*kappa2) * pow(gsl_vector_get(survTime1, i), *alpha2) * exp(LP + gsl_vector_get(V2, jj));
        logLH_prop  += - pow(gsl_vector_get(gamma, i), *nu2) * (*kappa2) * pow(gsl_vector_get(survTime1, i), exp(temp_prop)) * exp(LP + gsl_vector_get(V2, jj));
    }
    
    logPrior        = a2*z2 - b2* (*alpha2);
    logPrior_prop   = a2*temp_prop - b2* (exp(temp_prop));
    
    logProp_PropToIni = dnorm(z2, temp_prop, sqrt(mhProp_alpha2_var), 1);
    logProp_IniToProp = dnorm(temp_prop, z2, sqrt(mhProp_alpha2_var), 1);
    
    logR = logLH_prop - logLH + logPrior_prop - logPrior + logProp_PropToIni - logProp_IniToProp;
    
    
    
    
    u = log(runif(0, 1)) <logR;
    
    if(u == 1)
    {
        *alpha2 = exp(temp_prop);
        *accept_alpha2 += u;
    }
    
    return;
    
}










/* updating shape parameter: alpha3 */

/* use the random walk proposal */

/**/

void BweibMvnCorScrSM_updateSC3(gsl_vector *beta3,
                        double *alpha3,
                        double *kappa3,
                        double *nu3,
                        gsl_vector *gamma,
                        gsl_vector *V3,
                        gsl_vector *survTime1,
                        gsl_vector *survTime2,
                        gsl_vector *case11,
                        gsl_vector *cluster,
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
    int n = survTime1 -> size;
    int i, jj;
    
    logLH = 0; logLH_prop = 0;
    temp_prop = rgamma(pow(*alpha3, 2)/mhProp_alpha3_var, mhProp_alpha3_var/(*alpha3));
    
    for(i = 0; i < n; i++)
    {
        gsl_vector_view Xi = gsl_matrix_row(survCov3, i);
        gsl_blas_ddot(&Xi.vector, beta3, &LP);
        jj = (int) gsl_vector_get(cluster, i) - 1;          
        if(gsl_vector_get(case11, i) == 1)
        {
            logLH       += log(*alpha3) + (*alpha3 - 1)*log(gsl_vector_get(survTime2, i));
            logLH_prop  += log(temp_prop) + (temp_prop - 1)*log(gsl_vector_get(survTime2, i));
        }
        logLH       += - pow(gsl_vector_get(gamma, i), *nu3) * (*kappa3) * (pow(gsl_vector_get(survTime2, i), *alpha3) - pow(gsl_vector_get(survTime1, i), *alpha3)) * exp(LP + gsl_vector_get(V3, jj));
        logLH_prop  += - pow(gsl_vector_get(gamma, i), *nu3) * (*kappa3) * (pow(gsl_vector_get(survTime2, i), temp_prop) - pow(gsl_vector_get(survTime1, i), temp_prop)) * exp(LP + gsl_vector_get(V3, jj));
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







/* updating shape parameter: alpha3 */

/* use the random walk proposal with log transformation*/

/**/

void BweibMvnCorScrSM_updateSC3_rw2(gsl_vector *beta3,
                        double *alpha3,
                        double *kappa3,
                        double *nu3,
                        gsl_vector *gamma,
                        gsl_vector *V3,
                        gsl_vector *yStar,
                        gsl_vector *case11,
                        gsl_vector *cluster,
                        gsl_matrix *survCov3,
                        double mhProp_alpha3_var,
                        double a3,
                        double b3,
                        int *accept_alpha3)
{
    double z3, LP, logLH, logLH_prop;
    double temp_prop;
    double logPrior, logPrior_prop;
    double logProp_IniToProp, logProp_PropToIni;
    double logR;
    int u;
    int n = yStar -> size;
    int i, jj;
    
    z3 = log(*alpha3);
    
    logLH = 0; logLH_prop = 0;
    temp_prop = rnorm(z3, sqrt(mhProp_alpha3_var));
    
    for(i = 0; i < n; i++)
    {
        gsl_vector_view Xi = gsl_matrix_row(survCov3, i);
        gsl_blas_ddot(&Xi.vector, beta3, &LP);
        jj = (int) gsl_vector_get(cluster, i) - 1;
        if(gsl_vector_get(case11, i) == 1 && gsl_vector_get(yStar, i) != 0)
        {
            logLH       += log(*alpha3) + (*alpha3 - 1)*log(gsl_vector_get(yStar, i));
            logLH_prop  += temp_prop + (exp(temp_prop) - 1)*log(gsl_vector_get(yStar, i));
        }
        logLH       += - pow(gsl_vector_get(gamma, i), *nu3) * (*kappa3) * pow(gsl_vector_get(yStar, i), *alpha3) * exp(LP + gsl_vector_get(V3, jj));
        logLH_prop  += - pow(gsl_vector_get(gamma, i), *nu3) * (*kappa3) * pow(gsl_vector_get(yStar, i), exp(temp_prop)) * exp(LP + gsl_vector_get(V3, jj));
    }
    
    logPrior        = a3*z3 - b3* (*alpha3);
    logPrior_prop   = a3*temp_prop - b3* (exp(temp_prop));
    
    logProp_PropToIni = dnorm(z3, temp_prop, sqrt(mhProp_alpha3_var), 1);
    logProp_IniToProp = dnorm(temp_prop, z3, sqrt(mhProp_alpha3_var), 1);
    
    logR = logLH_prop - logLH + logPrior_prop - logPrior + logProp_PropToIni - logProp_IniToProp;
    
  
    
    u = log(runif(0, 1)) <logR;
   
    
    if(u == 1)
    {
        *alpha3 = exp(temp_prop);
        *accept_alpha3 += u;
    }
    
    return;
    
}





/* updating shape parameter: kappa1 */

/**/

void BweibMvnCorScrSM_updateSH1(gsl_vector *beta1,
                        double *alpha1,
                        double *kappa1,
                        gsl_vector *gamma,
                        gsl_vector *V1,
                        gsl_vector *survTime1,
                        gsl_vector *survEvent1,
                        gsl_vector *cluster,
                        gsl_matrix *survCov1,
                        double c1,
                        double d1)
{
    int n = survTime1 -> size;
    int i, jj;
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
        jj = (int) gsl_vector_get(cluster, i) - 1;
        
        Kappa1_rate += gsl_vector_get(gamma, i) * pow(gsl_vector_get(survTime1, i), *alpha1) * exp(LP + gsl_vector_get(V1, jj));
    }
    Kappa1_rate += d1;
    Kappa1_scale = 1/Kappa1_rate;
    
    *kappa1 = rgamma(Kappa1_shape, Kappa1_scale);
    return;
}







/* updating shape parameter: kappa2 */

/**/

void BweibMvnCorScrSM_updateSH2(gsl_vector *beta2,
                        double *alpha2,
                        double *kappa2,
                        double *nu2,
                        gsl_vector *gamma,
                        gsl_vector *V2,
                        gsl_vector *survTime1,
                        gsl_vector *case01,
                        gsl_vector *cluster,
                        gsl_matrix *survCov2,
                        double c2,
                        double d2)
{
    int n = survTime1 -> size;
    int i, jj;
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
        jj = (int) gsl_vector_get(cluster, i) - 1;        
        
        Kappa2_rate += pow(gsl_vector_get(gamma, i), *nu2) * pow(gsl_vector_get(survTime1, i), *alpha2) * exp(LP + gsl_vector_get(V2, jj));
    }
    Kappa2_rate += d2;
    Kappa2_scale = 1/Kappa2_rate;
    
    *kappa2 = rgamma(Kappa2_shape, Kappa2_scale);
    return;
}



/* updating shape parameter: kappa3 */

/**/

void BweibMvnCorScrSM_updateSH3(gsl_vector *beta3,
                        double *alpha3,
                        double *kappa3,
                        double *nu3,
                        gsl_vector *gamma,
                        gsl_vector *V3,
                        gsl_vector *yStar,
                        gsl_vector *case11,
                        gsl_vector *cluster,
                        gsl_matrix *survCov3,
                        double c3,
                        double d3)
{
    int n = yStar -> size;
    int i, jj;
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
        jj = (int) gsl_vector_get(cluster, i) - 1;        
        
        Kappa3_rate += pow(gsl_vector_get(gamma, i), *nu3) * pow(gsl_vector_get(yStar, i), *alpha3) * exp(LP + gsl_vector_get(V3, jj));
    }
    Kappa3_rate += d3;
    Kappa3_scale = 1/Kappa3_rate;
    
    *kappa3 = rgamma(Kappa3_shape, Kappa3_scale);
    return;
}













/* updating variance parameter: theta */

/* use the random walk proposal */

/**/

void BweibMvnCorScrSM_updateDP(gsl_vector *gamma,
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


 
void BweibMvnCorScrSM_updateFP_Gibbs(gsl_vector *beta1,
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
                             gsl_vector *V1,
                             gsl_vector *V2,
                             gsl_vector *V3,
                             gsl_vector *survTime1,
                             gsl_vector *yStar,
                             gsl_vector *survEvent1,
                             gsl_vector *survEvent2,
                             gsl_vector *cluster,                             
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
        gamma_rate  = BweibMvnCorScrSM_wFunc_old(i, beta1, beta2, beta3, alpha1, alpha2, alpha3, kappa1, kappa2, kappa3, V1, V2, V3, survTime1, yStar, cluster, survCov1, survCov2, survCov3) + 1/theta;
        
        gamma_scale = 1/gamma_rate;
        gsl_vector_set(gamma, i, rgamma(gamma_shape, gamma_scale));
        
    }
    return;
}







/* updating nu2 and nu3 */




void BweibMvnCorScrSM_updateMP(gsl_vector *beta2,
                       gsl_vector *beta3,
                       double *alpha2,
                       double *alpha3,
                       double *kappa2,
                       double *kappa3,
                       double *nu2,
                       double *nu3,
                       gsl_vector *gamma,
                       gsl_vector *V2,
                       gsl_vector *V3,
                       gsl_vector *survTime1,
                       gsl_vector *survTime2,
                       gsl_vector *case01,
                       gsl_vector *case11,
                       gsl_vector *cluster,
                       gsl_matrix *survCov2,
                       gsl_matrix *survCov3,
                       int *accept_nu2,
                       int *accept_nu3)
{
    double LP2, LP3, D1, D2, logLH;
    double D1_prop, D2_prop, logLH_prop;
    double nu2_prop_me, nu2_prop_var, nu2_prop;
    double nu3_prop_me, nu3_prop_var, nu3_prop;
    double nu2_prop_me_prop, nu2_prop_var_prop;
    double nu3_prop_me_prop, nu3_prop_var_prop;
    double logProp_IniToProp, logProp_PropToIni;
    double logR;
    int i, jj, u;
    
    int n = survTime1 -> size;


    
    
    /* update nu2 */
    
    logLH = 0; D1 = 0; D2 = 0;
    logLH_prop = 0; D1_prop = 0; D2_prop = 0;
    
    for(i = 0; i < n; i++)
    {
        gsl_vector_view Xi2 = gsl_matrix_row(survCov2, i);
        gsl_blas_ddot(&Xi2.vector, beta2, &LP2);
        jj = (int) gsl_vector_get(cluster, i) - 1;
        if(gsl_vector_get(case01, i) == 1)
        {
            logLH   += *nu2 * log(gsl_vector_get(gamma, i));
            D1      += log(gsl_vector_get(gamma, i));
        }
        logLH   += - pow(gsl_vector_get(gamma, i), *nu2) * (*kappa2) * pow(gsl_vector_get(survTime1, i), *alpha2) * exp(LP2 + gsl_vector_get(V2, jj));
        D1      += - pow(gsl_vector_get(gamma, i), *nu2) * (*kappa2) * pow(gsl_vector_get(survTime1, i), *alpha2) * exp(LP2 + gsl_vector_get(V2, jj)) * log(gsl_vector_get(gamma, i));
        D2      += - pow(gsl_vector_get(gamma, i), *nu2) * (*kappa2) * pow(gsl_vector_get(survTime1, i), *alpha2) * exp(LP2 + gsl_vector_get(V2, jj)) * pow(log(gsl_vector_get(gamma, i)), 2);
    }
    
    nu2_prop_me    = *nu2 - D1/D2;
    nu2_prop_var   = - pow(2.4, 2)/D2;
    
    nu2_prop = rnorm(nu2_prop_me, sqrt(nu2_prop_var));
    
    for(i = 0; i < n; i++)
    {
        gsl_vector_view Xi2 = gsl_matrix_row(survCov2, i);
        gsl_blas_ddot(&Xi2.vector, beta2, &LP2);
        jj = (int) gsl_vector_get(cluster, i) - 1;
        if(gsl_vector_get(case01, i) == 1)
        {
            logLH_prop   += nu2_prop * log(gsl_vector_get(gamma, i));
            D1_prop      += log(gsl_vector_get(gamma, i));
        }
        logLH_prop   += - pow(gsl_vector_get(gamma, i), nu2_prop) * (*kappa2) * pow(gsl_vector_get(survTime1, i), *alpha2) * exp(LP2 + gsl_vector_get(V2, jj));
        D1_prop      += - pow(gsl_vector_get(gamma, i), nu2_prop) * (*kappa2) * pow(gsl_vector_get(survTime1, i), *alpha2) * exp(LP2 + gsl_vector_get(V2, jj)) * log(gsl_vector_get(gamma, i));
        D2_prop      += - pow(gsl_vector_get(gamma, i), nu2_prop) * (*kappa2) * pow(gsl_vector_get(survTime1, i), *alpha2) * exp(LP2 + gsl_vector_get(V2, jj)) * pow(log(gsl_vector_get(gamma, i)), 2);
    }
    
    nu2_prop_me_prop   = nu2_prop - D1_prop/D2_prop;
    nu2_prop_var_prop  = - pow(2.4, 2)/D2_prop;
    
    logProp_IniToProp = dnorm(nu2_prop, nu2_prop_me, sqrt(nu2_prop_var), 1);
    logProp_PropToIni = dnorm(*nu2, nu2_prop_me_prop, sqrt(nu2_prop_var_prop), 1);
    
    logR = logLH_prop - logLH + logProp_PropToIni - logProp_IniToProp;
    
    u = log(runif(0, 1)) < logR;
    
    if(u == 1)
    {
        *nu2 = nu2_prop;
        *accept_nu2 += u;
    }
       
    /* update nu3 */
    
    logLH = 0; D1 = 0; D2 = 0;
    logLH_prop = 0; D1_prop = 0; D2_prop = 0;
    
    for(i = 0; i < n; i++)
    {
        gsl_vector_view Xi3 = gsl_matrix_row(survCov3, i);
        gsl_blas_ddot(&Xi3.vector, beta3, &LP3);
        jj = (int) gsl_vector_get(cluster, i) - 1;
        if(gsl_vector_get(case11, i) == 1)
        {
            logLH   += *nu3 * log(gsl_vector_get(gamma, i));
            D1      += log(gsl_vector_get(gamma, i));
        }
        logLH   += - pow(gsl_vector_get(gamma, i), *nu3) * (*kappa3) * (pow(gsl_vector_get(survTime2, i), *alpha3) - pow(gsl_vector_get(survTime1, i), *alpha3)) * exp(LP3 + gsl_vector_get(V3, jj));
        D1      += - pow(gsl_vector_get(gamma, i), *nu3) * (*kappa3) * (pow(gsl_vector_get(survTime2, i), *alpha3) - pow(gsl_vector_get(survTime1, i), *alpha3)) * exp(LP3 + gsl_vector_get(V3, jj)) * log(gsl_vector_get(gamma, i));
        D2      += - pow(gsl_vector_get(gamma, i), *nu3) * (*kappa3) * (pow(gsl_vector_get(survTime2, i), *alpha3) - pow(gsl_vector_get(survTime1, i), *alpha3)) * exp(LP3 + gsl_vector_get(V3, jj)) * pow(log(gsl_vector_get(gamma, i)), 2);
    }
    
    nu3_prop_me    = *nu3 - D1/D2;
    nu3_prop_var   = - pow(2.4, 2)/D2;
    
    nu3_prop = rnorm(nu3_prop_me, sqrt(nu3_prop_var));
    
    for(i = 0; i < n; i++)
    {
        gsl_vector_view Xi3 = gsl_matrix_row(survCov3, i);
        gsl_blas_ddot(&Xi3.vector, beta3, &LP3);
        jj = (int) gsl_vector_get(cluster, i) - 1;
        if(gsl_vector_get(case11, i) == 1)
        {
            logLH_prop   += nu3_prop * log(gsl_vector_get(gamma, i));
            D1_prop      += log(gsl_vector_get(gamma, i));
        }
        logLH_prop   += - pow(gsl_vector_get(gamma, i), nu3_prop) * (*kappa3) * (pow(gsl_vector_get(survTime2, i), *alpha3) - pow(gsl_vector_get(survTime1, i), *alpha3)) * exp(LP3 + gsl_vector_get(V3, jj));
        D1_prop      += - pow(gsl_vector_get(gamma, i), nu3_prop) * (*kappa3) * (pow(gsl_vector_get(survTime2, i), *alpha3) - pow(gsl_vector_get(survTime1, i), *alpha3)) * exp(LP3 + gsl_vector_get(V3, jj)) * log(gsl_vector_get(gamma, i));
        D2_prop      += - pow(gsl_vector_get(gamma, i), nu3_prop) * (*kappa3) * (pow(gsl_vector_get(survTime2, i), *alpha3) - pow(gsl_vector_get(survTime1, i), *alpha3)) * exp(LP3 + gsl_vector_get(V3, jj)) * pow(log(gsl_vector_get(gamma, i)), 2);
    }
    
    nu3_prop_me_prop   = nu3_prop - D1_prop/D2_prop;
    nu3_prop_var_prop  = - pow(2.4, 2)/D2_prop;
    
    logProp_IniToProp = dnorm(nu3_prop, nu3_prop_me, sqrt(nu3_prop_var), 1);
    logProp_PropToIni = dnorm(*nu3, nu3_prop_me_prop, sqrt(nu3_prop_var_prop), 1);
    
    logR = logLH_prop - logLH + logProp_PropToIni - logProp_IniToProp;
    
    u = log(runif(0, 1)) < logR;
    
       
    if(u == 1)
    {
        *nu3 = nu3_prop;
        *accept_nu3 += u;
    }

    return;
    
}
















/* updating frailty parameter: gamma */

/**/



void BweibMvnCorScrSM_updateFP(gsl_vector *beta1,
                       gsl_vector *beta2,
                       gsl_vector *beta3,
                       double alpha1,
                       double alpha2,
                       double alpha3,
                       double kappa1,
                       double kappa2,
                       double kappa3,
                       double nu2,
                       double nu3,
                       double theta,
                       gsl_vector *gamma,
                       gsl_vector *V1,
                       gsl_vector *V2,
                       gsl_vector *V3,
                       gsl_vector *survTime1,
                       gsl_vector *survTime2,
                       gsl_vector *survEvent1,
                       gsl_vector *survEvent2,
                       gsl_vector *cluster,
                       gsl_matrix *survCov1,
                       gsl_matrix *survCov2,
                       gsl_matrix *survCov3,
                       gsl_vector *accept_gamma,
                       double mhProp_gamma_var,
                       int numGamUpdate,
                       gsl_vector *mhGam_chk,
                       int *ChgProp)
{
    double logLH, logLH_prop;
    double gam, temp_prop, del1, del2;
    double logPrior, logPrior_prop;
    double logProp_IniToProp, logProp_PropToIni;
    double logR;
    int u, i;
    int n = survTime1 -> size;
    
    for(i = 0; i < n; i++)
    {
        if(gsl_vector_get(mhGam_chk, i) == 0)
        {
            gam = gsl_vector_get(gamma, i);
            del1 = gsl_vector_get(survEvent1, i);
            del2 = gsl_vector_get(survEvent2, i);
            logLH = 0; logLH_prop = 0;
            
            temp_prop = rgamma(pow(gam, 2)/mhProp_gamma_var, mhProp_gamma_var/(gam));
            
            if(temp_prop == 0){
                gsl_vector_set(mhGam_chk, i, 1);
                *ChgProp = 1;
            }
            else
            {
                logLH       += (del1 + nu2 * del2 + (nu3-nu2) * del1 * del2 + 1/theta - 1) * log(gam);
                logLH_prop  += (del1 + nu2 * del2 + (nu3-nu2) * del1 * del2 + 1/theta - 1) * log(temp_prop);
                
                logLH += -BweibMvnCorScrSM_wFunc(i, beta1, beta2, beta3, alpha1, alpha2, alpha3, kappa1, kappa2, kappa3, nu2, nu3, gam, V1, V2, V3, survTime1, survTime2, cluster, survCov1, survCov2, survCov3);
                logLH_prop += -BweibMvnCorScrSM_wFunc(i, beta1, beta2, beta3, alpha1, alpha2, alpha3, kappa1, kappa2, kappa3, nu2, nu3, temp_prop, V1, V2, V3, survTime1, survTime2, cluster, survCov1, survCov2, survCov3);
                
                logLH       += 1/theta * gam;
                logLH_prop  += 1/theta * temp_prop;
                
                logPrior        = dgamma(gam, 1/theta, theta, 1);
                logPrior_prop   = dgamma(temp_prop, 1/theta, theta, 1);
                
                logProp_PropToIni = dgamma(gam, pow(temp_prop, 2)/mhProp_gamma_var, mhProp_gamma_var/(temp_prop), 1);
                logProp_IniToProp = dgamma(temp_prop, pow(gam, 2)/mhProp_gamma_var, mhProp_gamma_var/(gam), 1);
                
                logR = logLH_prop - logLH + logPrior_prop - logPrior + logProp_PropToIni - logProp_IniToProp;
                
                u = log(runif(0, 1)) <logR;

                
                if(u == 1)
                {
                    gsl_vector_set(gamma, i, temp_prop);
                    gsl_vector_set(accept_gamma, i, (gsl_vector_get(accept_gamma, i) + u));
                }
            }
        }
        
        if(gsl_vector_get(mhGam_chk, i) == 1)
        {
            gam = gsl_vector_get(gamma, i);
            del1 = gsl_vector_get(survEvent1, i);
            del2 = gsl_vector_get(survEvent2, i);
            logLH = 0; logLH_prop = 0;
            
   
            
            temp_prop = rgamma(1, 1);
            
            
            logLH       += (del1 + nu2 * del2 + (nu3-nu2) * del1 * del2 + 1/theta - 1) * log(gam);
            logLH_prop  += (del1 + nu2 * del2 + (nu3-nu2) * del1 * del2 + 1/theta - 1) * log(temp_prop);
            
            logLH += -BweibMvnCorScrSM_wFunc(i, beta1, beta2, beta3, alpha1, alpha2, alpha3, kappa1, kappa2, kappa3, nu2, nu3, gam, V1, V2, V3, survTime1, survTime2, cluster, survCov1, survCov2, survCov3);
            logLH_prop += -BweibMvnCorScrSM_wFunc(i, beta1, beta2, beta3, alpha1, alpha2, alpha3, kappa1, kappa2, kappa3, nu2, nu3, temp_prop, V1, V2, V3, survTime1, survTime2, cluster, survCov1, survCov2, survCov3);
            
            logLH       += 1/theta * gam;
            logLH_prop  += 1/theta * temp_prop;
            
            logPrior        = dgamma(gam, 1/theta, theta, 1);
            logPrior_prop   = dgamma(temp_prop, 1/theta, theta, 1);
            
            logProp_PropToIni = dgamma(gam, 1, 1, 1);
            logProp_IniToProp = dgamma(temp_prop, 1, 1, 1);
            
            logR = logLH_prop - logLH + logPrior_prop - logPrior + logProp_PropToIni - logProp_IniToProp;
            
            u = log(runif(0, 1)) <logR;
            
 
            
            
            if(u == 1)
            {
                gsl_vector_set(gamma, i, temp_prop);
                gsl_vector_set(accept_gamma, i, (gsl_vector_get(accept_gamma, i) + u));
            }
        }
    }
    
   
    

    return;
}














/* updating cluster-specific random effect using RW-MH: V1 */

/**/



void BweibMvnCorScrSM_updateCP1(gsl_vector *beta1,
                        double alpha1,
                        double kappa1,
                        gsl_vector *gamma,
                        gsl_vector *V1,
                        gsl_vector *V2,
                        gsl_vector *V3,
                        gsl_matrix *Sigma_V,
                        gsl_vector *survTime1,
                        gsl_vector *survEvent1,
                        gsl_vector *cluster,
                        gsl_matrix *survCov1,
                        gsl_vector *n_j,
                        gsl_vector *accept_V1,
                        double mhProp_V1_var)
{
    double LP, logLH, logLH_prop;
    double logPrior, logPrior_prop, temp_prop;
    double logProp_IniToProp, logProp_PropToIni;
    double logR;
    int u, i, j;

    int J = V1  -> size;
    
    int startInx = 0;
    int endInx = 0;
    
    gsl_vector *Vj = gsl_vector_calloc(3);
    gsl_vector *Vj_prop = gsl_vector_calloc(3);
    gsl_vector *zeroVec = gsl_vector_calloc(3);
    
    gsl_matrix *invSigma_V = gsl_matrix_calloc(3,3);
    matrixInv(Sigma_V, invSigma_V);
    
    for(j = 0; j < J; j++)
    {
        logLH = 0; logLH_prop = 0;
        
        gsl_vector_set(Vj, 0, gsl_vector_get(V1, j));
        gsl_vector_set(Vj, 1, gsl_vector_get(V2, j));
        gsl_vector_set(Vj, 2, gsl_vector_get(V3, j));
        gsl_vector_memcpy(Vj_prop, Vj);
        
        temp_prop = rnorm(gsl_vector_get(V1, j), sqrt(mhProp_V1_var));
        gsl_vector_set(Vj_prop, 0, temp_prop);
        
        endInx += (int) gsl_vector_get(n_j, j);
        
        for(i = startInx; i < endInx; i++)
        {
            gsl_vector_view Xi = gsl_matrix_row(survCov1, i);
            gsl_blas_ddot(&Xi.vector, beta1, &LP);
            
            if(gsl_vector_get(survEvent1, i) == 1)
            {
                logLH += gsl_vector_get(V1, j);
                logLH_prop += temp_prop;
            }
            logLH +=  -gsl_vector_get(gamma, i) * (kappa1) * pow(gsl_vector_get(survTime1, i), alpha1) * exp(LP + gsl_vector_get(V1, j));
            logLH_prop +=  -gsl_vector_get(gamma, i) * (kappa1) * pow(gsl_vector_get(survTime1, i), alpha1) * exp(LP + temp_prop);
        }
        startInx = endInx;
        
        c_dmvnorm2(Vj, zeroVec, 1, invSigma_V, &logPrior);
        c_dmvnorm2(Vj_prop, zeroVec, 1, invSigma_V, &logPrior_prop);
        
        logProp_PropToIni = dnorm(gsl_vector_get(V1, j), temp_prop, sqrt(mhProp_V1_var), 1);
        logProp_IniToProp = dnorm(temp_prop, gsl_vector_get(V1, j), sqrt(mhProp_V1_var), 1);
        
        logR = logLH_prop - logLH + logPrior_prop - logPrior + logProp_PropToIni - logProp_IniToProp;
        
        u = log(runif(0, 1)) <logR;
        
        if(u == 1)
        {
            gsl_vector_set(V1, j, temp_prop);
            gsl_vector_set(accept_V1, j, (gsl_vector_get(accept_V1, j) + u));
        }
    }
    
    gsl_vector_free(Vj);
    gsl_vector_free(Vj_prop);
    gsl_vector_free(zeroVec);
    gsl_matrix_free(invSigma_V);
    
    return;
}








/* updating cluster-specific random effect using AMCMC: V1 */

/**/



void BweibMvnCorScrSM_updateCP1_amcmc(gsl_vector *beta1,
                                    double alpha1,
                                    double kappa1,
                                    gsl_vector *gamma,
                                    gsl_vector *V1,
                                    gsl_vector *V2,
                                    gsl_vector *V3,
                                    gsl_matrix *Sigma_V,
                                    gsl_vector *survTime1,
                                    gsl_vector *survEvent1,
                                    gsl_vector *cluster,
                                    gsl_matrix *survCov1,
                                    gsl_vector *n_j,
                                    gsl_vector *accept_V1,
                                    double mhProp_V1_var)
{
    double LP, logLH, logLH_prop;
    double D1, D2, D1_prop, D2_prop;
    double logPrior, logPrior_prop, temp_prop;
    double logProp_IniToProp, logProp_PropToIni;
    double logR, term, v_prop_me, v_prop_var, v_prop_me_prop, v_prop_var_prop;
    int u, i, j;

    int J = V1  -> size;
   
    
    int startInx = 0;
    int endInx = 0;
    
    gsl_vector *Vj = gsl_vector_calloc(3);
    gsl_vector *Vj_prop = gsl_vector_calloc(3);
    gsl_vector *zeroVec = gsl_vector_calloc(3);
    
    gsl_matrix *invSigma_V = gsl_matrix_calloc(3,3);
    matrixInv(Sigma_V, invSigma_V);
    
    for(j = 0; j < J; j++)
    {
        logLH = 0; logLH_prop = 0;
        D1 = 0; D2 = 0;
        D1_prop = 0; D2_prop = 0;
        
        gsl_vector_set(Vj, 0, gsl_vector_get(V1, j));
        gsl_vector_set(Vj, 1, gsl_vector_get(V2, j));
        gsl_vector_set(Vj, 2, gsl_vector_get(V3, j));
        gsl_vector_memcpy(Vj_prop, Vj);
        
        endInx += (int) gsl_vector_get(n_j, j);
        
        for(i = startInx; i < endInx; i++)
        {
            gsl_vector_view Xi = gsl_matrix_row(survCov1, i);
            gsl_blas_ddot(&Xi.vector, beta1, &LP);
            
            if(gsl_vector_get(survEvent1, i) == 1)
            {
                logLH += gsl_vector_get(V1, j);
                D1 += 1;
            }
            term = -gsl_vector_get(gamma, i) * (kappa1) * pow(gsl_vector_get(survTime1, i), alpha1) * exp(LP + gsl_vector_get(V1, j));
            logLH +=  term;
            D1 += term;
            D2 += term;
        }
  
        D1 += -gsl_matrix_get(invSigma_V, 0, 0)*gsl_vector_get(Vj, 0)-gsl_matrix_get(invSigma_V, 1, 0)*gsl_vector_get(Vj, 1) - gsl_matrix_get(invSigma_V, 2, 0)*gsl_vector_get(Vj, 2);
        D2 += -gsl_matrix_get(invSigma_V, 0, 0);
        
        
        v_prop_me    = gsl_vector_get(V1, j) - D1/D2;
        v_prop_var   = - pow(2.4, 2)/D2;
        
        temp_prop = rnorm(v_prop_me, sqrt(v_prop_var));
        gsl_vector_set(Vj_prop, 0, temp_prop);
        
        for(i = startInx; i < endInx; i++)
        {
            gsl_vector_view Xi = gsl_matrix_row(survCov1, i);
            gsl_blas_ddot(&Xi.vector, beta1, &LP);
            
            if(gsl_vector_get(survEvent1, i) == 1)
            {
                logLH_prop += temp_prop;
                D1_prop += 1;
            }
            term = -gsl_vector_get(gamma, i) * (kappa1) * pow(gsl_vector_get(survTime1, i), alpha1) * exp(LP + temp_prop);
            logLH_prop +=  term;
            D1_prop += term;
            D2_prop += term;
        }

        D1_prop += -gsl_matrix_get(invSigma_V, 0, 0)*temp_prop-gsl_matrix_get(invSigma_V, 1, 0)*gsl_vector_get(Vj, 1) - gsl_matrix_get(invSigma_V, 2, 0)*gsl_vector_get(Vj, 2);
        D2_prop += -gsl_matrix_get(invSigma_V, 0, 0);
        
        v_prop_me_prop   = temp_prop - D1_prop/D2_prop;
        v_prop_var_prop  = - pow(2.4, 2)/D2_prop;
        
        startInx = endInx;
        

        
        c_dmvnorm2(Vj, zeroVec, 1, invSigma_V, &logPrior);
        c_dmvnorm2(Vj_prop, zeroVec, 1, invSigma_V, &logPrior_prop);
        
        logProp_PropToIni = dnorm(gsl_vector_get(V1, j), v_prop_me_prop, sqrt(v_prop_var_prop), 1);
        logProp_IniToProp = dnorm(temp_prop, v_prop_me, sqrt(v_prop_var), 1);
        
        logR = logLH_prop - logLH + logPrior_prop - logPrior + logProp_PropToIni - logProp_IniToProp;
        
        u = log(runif(0, 1)) <logR;
        
        if(u == 1)
        {
            gsl_vector_set(V1, j, temp_prop);
            gsl_vector_set(accept_V1, j, (gsl_vector_get(accept_V1, j) + u));
        }
        
       
    }
    
    gsl_vector_free(Vj);
    gsl_vector_free(Vj_prop);
    gsl_vector_free(zeroVec);
    gsl_matrix_free(invSigma_V);
    
    return;
}






/* updating cluster-specific random effect using RW-MH: V2 */

/**/



void BweibMvnCorScrSM_updateCP2(gsl_vector *beta2,
                        double alpha2,
                        double kappa2,
                        double nu2,
                        gsl_vector *gamma,
                        gsl_vector *V1,
                        gsl_vector *V2,
                        gsl_vector *V3,
                        gsl_matrix *Sigma_V,
                        gsl_vector *survTime1,
                        gsl_vector *case01,
                        gsl_vector *cluster,
                        gsl_matrix *survCov2,
                        gsl_vector *n_j,
                        gsl_vector *accept_V2,
                        double mhProp_V2_var)
{
    double LP, logLH, logLH_prop;
    double logPrior, logPrior_prop, temp_prop;
    double logProp_IniToProp, logProp_PropToIni;
    double logR;
    int u, i, j;

    int J = V2  -> size;
    
    int startInx = 0;
    int endInx = 0;
    
    gsl_vector *Vj = gsl_vector_calloc(3);
    gsl_vector *Vj_prop = gsl_vector_calloc(3);
    gsl_vector *zeroVec = gsl_vector_calloc(3);
    
    gsl_matrix *invSigma_V = gsl_matrix_calloc(3,3);
    matrixInv(Sigma_V, invSigma_V);
    
    for(j = 0; j < J; j++)
    {
        logLH = 0; logLH_prop = 0;
        
        gsl_vector_set(Vj, 0, gsl_vector_get(V1, j));
        gsl_vector_set(Vj, 1, gsl_vector_get(V2, j));
        gsl_vector_set(Vj, 2, gsl_vector_get(V3, j));
        gsl_vector_memcpy(Vj_prop, Vj);
        
        temp_prop = rnorm(gsl_vector_get(V2, j), sqrt(mhProp_V2_var));
        gsl_vector_set(Vj_prop, 1, temp_prop);
        
        endInx += (int) gsl_vector_get(n_j, j);
        
        for(i = startInx; i < endInx; i++)
        {
            gsl_vector_view Xi = gsl_matrix_row(survCov2, i);
            gsl_blas_ddot(&Xi.vector, beta2, &LP);
            
            if(gsl_vector_get(case01, i) == 1)
            {
                logLH += gsl_vector_get(V2, j);
                logLH_prop += temp_prop;
            }
            logLH +=  -pow(gsl_vector_get(gamma, i), nu2) * (kappa2) * pow(gsl_vector_get(survTime1, i), alpha2) * exp(LP + gsl_vector_get(V2, j));
            logLH_prop +=  -pow(gsl_vector_get(gamma, i), nu2) * (kappa2) * pow(gsl_vector_get(survTime1, i), alpha2) * exp(LP + temp_prop);
        }
        startInx = endInx;
        
        c_dmvnorm2(Vj, zeroVec, 1, invSigma_V, &logPrior);
        c_dmvnorm2(Vj_prop, zeroVec, 1, invSigma_V, &logPrior_prop);
        
        logProp_PropToIni = dnorm(gsl_vector_get(V2, j), temp_prop, sqrt(mhProp_V2_var), 1);
        logProp_IniToProp = dnorm(temp_prop, gsl_vector_get(V2, j), sqrt(mhProp_V2_var), 1);
        
        logR = logLH_prop - logLH + logPrior_prop - logPrior + logProp_PropToIni - logProp_IniToProp;
        
        u = log(runif(0, 1)) <logR;
        
        if(u == 1)
        {
            gsl_vector_set(V2, j, temp_prop);
            gsl_vector_set(accept_V2, j, (gsl_vector_get(accept_V2, j) + u));
        }
    }
    
    gsl_vector_free(Vj);
    gsl_vector_free(Vj_prop);
    gsl_vector_free(zeroVec);
    gsl_matrix_free(invSigma_V);
    
    return;
}










/* updating cluster-specific random effect using AMCMC: V2 */

/**/



void BweibMvnCorScrSM_updateCP2_amcmc(gsl_vector *beta2,
                        double alpha2,
                        double kappa2,
                        double nu2,
                        gsl_vector *gamma,
                        gsl_vector *V1,
                        gsl_vector *V2,
                        gsl_vector *V3,
                        gsl_matrix *Sigma_V,
                        gsl_vector *survTime1,
                        gsl_vector *case01,
                        gsl_vector *cluster,
                        gsl_matrix *survCov2,
                        gsl_vector *n_j,
                        gsl_vector *accept_V2,
                        double mhProp_V2_var)
{
    double LP, logLH, logLH_prop;
    double D1, D2, D1_prop, D2_prop;    
    double logPrior, logPrior_prop, temp_prop;
    double logProp_IniToProp, logProp_PropToIni;
    double logR, term, v_prop_me, v_prop_var, v_prop_me_prop, v_prop_var_prop;
    int u, i, j;

    int J = V2  -> size;    
    
    int startInx = 0;
    int endInx = 0;
    
    gsl_vector *Vj = gsl_vector_calloc(3);
    gsl_vector *Vj_prop = gsl_vector_calloc(3);
    gsl_vector *zeroVec = gsl_vector_calloc(3);
    
    gsl_matrix *invSigma_V = gsl_matrix_calloc(3,3);
    matrixInv(Sigma_V, invSigma_V);
    
    for(j = 0; j < J; j++)
    {
        logLH = 0; logLH_prop = 0;
        D1 = 0; D2 = 0;
        D1_prop = 0; D2_prop = 0;
        
        gsl_vector_set(Vj, 0, gsl_vector_get(V1, j));
        gsl_vector_set(Vj, 1, gsl_vector_get(V2, j));
        gsl_vector_set(Vj, 2, gsl_vector_get(V3, j));
        gsl_vector_memcpy(Vj_prop, Vj);
        
        endInx += (int) gsl_vector_get(n_j, j);
        
        for(i = startInx; i < endInx; i++)
        {
            gsl_vector_view Xi = gsl_matrix_row(survCov2, i);
            gsl_blas_ddot(&Xi.vector, beta2, &LP);
            
            if(gsl_vector_get(case01, i) == 1)
            {
                logLH += gsl_vector_get(V2, j);
                D1 += 1;
            }
            term = -pow(gsl_vector_get(gamma, i), nu2) * (kappa2) * pow(gsl_vector_get(survTime1, i), alpha2) * exp(LP + gsl_vector_get(V2, j));
            logLH +=  term;
            D1 += term;
            D2 += term;
        }
        
        D1 += -gsl_matrix_get(invSigma_V, 1, 1)*gsl_vector_get(Vj, 1)-gsl_matrix_get(invSigma_V, 0, 1)*gsl_vector_get(Vj, 0) - gsl_matrix_get(invSigma_V, 2, 1)*gsl_vector_get(Vj, 2);
        D2 += -gsl_matrix_get(invSigma_V, 1, 1);
        
        v_prop_me    = gsl_vector_get(V2, j) - D1/D2;
        v_prop_var   = - pow(2.4, 2)/D2;
        
        temp_prop = rnorm(v_prop_me, sqrt(v_prop_var));
        gsl_vector_set(Vj_prop, 1, temp_prop);
        
        for(i = startInx; i < endInx; i++)
        {
            gsl_vector_view Xi = gsl_matrix_row(survCov2, i);
            gsl_blas_ddot(&Xi.vector, beta2, &LP);
            
            if(gsl_vector_get(case01, i) == 1)
            {
                logLH_prop += temp_prop;
                D1_prop += 1;
            }
            term = -pow(gsl_vector_get(gamma, i), nu2) * (kappa2) * pow(gsl_vector_get(survTime1, i), alpha2) * exp(LP + temp_prop);
            logLH_prop +=  term;
            D1_prop += term;
            D2_prop += term;
        }
        
        D1_prop += -gsl_matrix_get(invSigma_V, 1, 1)*temp_prop - gsl_matrix_get(invSigma_V, 0, 1)*gsl_vector_get(Vj, 0) - gsl_matrix_get(invSigma_V, 2, 1)*gsl_vector_get(Vj, 2);
        D2_prop += -gsl_matrix_get(invSigma_V, 1, 1);
        
        v_prop_me_prop   = temp_prop - D1_prop/D2_prop;
        v_prop_var_prop  = - pow(2.4, 2)/D2_prop;
         
        startInx = endInx;
        
        c_dmvnorm2(Vj, zeroVec, 1, invSigma_V, &logPrior);
        c_dmvnorm2(Vj_prop, zeroVec, 1, invSigma_V, &logPrior_prop);
        
        logProp_PropToIni = dnorm(gsl_vector_get(V2, j), v_prop_me_prop, sqrt(v_prop_var_prop), 1);
        logProp_IniToProp = dnorm(temp_prop, v_prop_me, sqrt(v_prop_var), 1);
        
        logR = logLH_prop - logLH + logPrior_prop - logPrior + logProp_PropToIni - logProp_IniToProp;
        
        u = log(runif(0, 1)) <logR;
        
        if(u == 1)
        {
            gsl_vector_set(V2, j, temp_prop);
            gsl_vector_set(accept_V2, j, (gsl_vector_get(accept_V2, j) + u));
        }
    }
    
    gsl_vector_free(Vj);
    gsl_vector_free(Vj_prop);
    gsl_vector_free(zeroVec);
    gsl_matrix_free(invSigma_V);
    
    return;
}













/* updating cluster-specific random effect using RW-MH: V3 */

/**/



void BweibMvnCorScrSM_updateCP3(gsl_vector *beta3,
                        double alpha3,
                        double kappa3,
                        double nu3,
                        gsl_vector *gamma,
                        gsl_vector *V1,
                        gsl_vector *V2,
                        gsl_vector *V3,
                        gsl_matrix *Sigma_V,
                        gsl_vector *survTime1,
                        gsl_vector *survTime2,
                        gsl_vector *case11,
                        gsl_vector *cluster,
                        gsl_matrix *survCov3,
                        gsl_vector *n_j,
                        gsl_vector *accept_V3,
                        double mhProp_V3_var)
{
    double LP, logLH, logLH_prop;
    double logPrior, logPrior_prop, temp_prop;
    double logProp_IniToProp, logProp_PropToIni;
    double logR;
    int u, i, j;

    int J = V3  -> size;
    
    int startInx = 0;
    int endInx = 0;
    
    gsl_vector *Vj = gsl_vector_calloc(3);
    gsl_vector *Vj_prop = gsl_vector_calloc(3);
    gsl_vector *zeroVec = gsl_vector_calloc(3);
    
    gsl_matrix *invSigma_V = gsl_matrix_calloc(3,3);
    matrixInv(Sigma_V, invSigma_V);
    
    for(j = 0; j < J; j++)
    {
        logLH = 0; logLH_prop = 0;
        
        gsl_vector_set(Vj, 0, gsl_vector_get(V1, j));
        gsl_vector_set(Vj, 1, gsl_vector_get(V2, j));
        gsl_vector_set(Vj, 2, gsl_vector_get(V3, j));
        gsl_vector_memcpy(Vj_prop, Vj);
        
        temp_prop = rnorm(gsl_vector_get(V3, j), sqrt(mhProp_V3_var));
        gsl_vector_set(Vj_prop, 2, temp_prop);
        
        endInx += (int) gsl_vector_get(n_j, j);
        
        for(i = startInx; i < endInx; i++)
        {
            gsl_vector_view Xi = gsl_matrix_row(survCov3, i);
            gsl_blas_ddot(&Xi.vector, beta3, &LP);
            
            if(gsl_vector_get(case11, i) == 1)
            {
                logLH += gsl_vector_get(V3, j);
                logLH_prop += temp_prop;
            }
            logLH += -pow(gsl_vector_get(gamma, i), nu3) * (kappa3) * (pow(gsl_vector_get(survTime2, i), alpha3)-pow(gsl_vector_get(survTime1, i), alpha3)) * exp(LP + gsl_vector_get(V3, j));
            logLH_prop += -pow(gsl_vector_get(gamma, i), nu3) * (kappa3) * (pow(gsl_vector_get(survTime2, i), alpha3)-pow(gsl_vector_get(survTime1, i), alpha3)) * exp(LP + temp_prop);
        }
        startInx = endInx;
        
        c_dmvnorm2(Vj, zeroVec, 1, invSigma_V, &logPrior);
        c_dmvnorm2(Vj_prop, zeroVec, 1, invSigma_V, &logPrior_prop);
        
        logProp_PropToIni = dnorm(gsl_vector_get(V3, j), temp_prop, sqrt(mhProp_V3_var), 1);
        logProp_IniToProp = dnorm(temp_prop, gsl_vector_get(V3, j), sqrt(mhProp_V3_var), 1);
        
        logR = logLH_prop - logLH + logPrior_prop - logPrior + logProp_PropToIni - logProp_IniToProp;
        
        u = log(runif(0, 1)) <logR;
        
        if(u == 1)
        {
            gsl_vector_set(V3, j, temp_prop);
            gsl_vector_set(accept_V3, j, (gsl_vector_get(accept_V3, j) + u));
        }
    }
    
    gsl_vector_free(Vj);
    gsl_vector_free(Vj_prop);
    gsl_vector_free(zeroVec);
    gsl_matrix_free(invSigma_V);
    
    return;
}







/* updating cluster-specific random effect using AMCMC: V3 */

/**/



void BweibMvnCorScrSM_updateCP3_amcmc(gsl_vector *beta3,
                        double alpha3,
                        double kappa3,
                        double nu3,
                        gsl_vector *gamma,
                        gsl_vector *V1,
                        gsl_vector *V2,
                        gsl_vector *V3,
                        gsl_matrix *Sigma_V,
                        gsl_vector *yStar,
                        gsl_vector *case11,
                        gsl_vector *cluster,
                        gsl_matrix *survCov3,
                        gsl_vector *n_j,
                        gsl_vector *accept_V3,
                        double mhProp_V3_var)
{
    double LP, logLH, logLH_prop;
    double D1, D2, D1_prop, D2_prop;        
    double logPrior, logPrior_prop, temp_prop;
    double logProp_IniToProp, logProp_PropToIni;
    double logR, term, v_prop_me, v_prop_var, v_prop_me_prop, v_prop_var_prop;
    int u, i, j;

    int J = V3  -> size;
    
    int startInx = 0;
    int endInx = 0;
    
    gsl_vector *Vj = gsl_vector_calloc(3);
    gsl_vector *Vj_prop = gsl_vector_calloc(3);
    gsl_vector *zeroVec = gsl_vector_calloc(3);
    
    gsl_matrix *invSigma_V = gsl_matrix_calloc(3,3);
    matrixInv(Sigma_V, invSigma_V);
    
    for(j = 0; j < J; j++)
    {
        logLH = 0; logLH_prop = 0;
        D1 = 0; D2 = 0;
        D1_prop = 0; D2_prop = 0;
        
        gsl_vector_set(Vj, 0, gsl_vector_get(V1, j));
        gsl_vector_set(Vj, 1, gsl_vector_get(V2, j));
        gsl_vector_set(Vj, 2, gsl_vector_get(V3, j));
        gsl_vector_memcpy(Vj_prop, Vj);
        
        endInx += (int) gsl_vector_get(n_j, j);
        
        for(i = startInx; i < endInx; i++)
        {
            gsl_vector_view Xi = gsl_matrix_row(survCov3, i);
            gsl_blas_ddot(&Xi.vector, beta3, &LP);
            
            if(gsl_vector_get(case11, i) == 1)
            {
                logLH += gsl_vector_get(V3, j);
                D1 += 1;
            }
            term = -pow(gsl_vector_get(gamma, i), nu3) * (kappa3) * pow(gsl_vector_get(yStar, i), alpha3) * exp(LP + gsl_vector_get(V3, j));
            logLH +=  term;
            D1 += term;
            D2 += term;
        }
        
        D1 += -gsl_matrix_get(invSigma_V, 2, 2)*gsl_vector_get(Vj, 2)-gsl_matrix_get(invSigma_V, 0, 2)*gsl_vector_get(Vj, 0) - gsl_matrix_get(invSigma_V, 1, 2)*gsl_vector_get(Vj, 1);
        D2 += -gsl_matrix_get(invSigma_V, 2, 2);
        
        v_prop_me    = gsl_vector_get(V3, j) - D1/D2;
        v_prop_var   = - pow(2.4, 2)/D2;
        
        temp_prop = rnorm(v_prop_me, sqrt(v_prop_var));
        gsl_vector_set(Vj_prop, 2, temp_prop);
        
        for(i = startInx; i < endInx; i++)
        {
            gsl_vector_view Xi = gsl_matrix_row(survCov3, i);
            gsl_blas_ddot(&Xi.vector, beta3, &LP);
            
            if(gsl_vector_get(case11, i) == 1)
            {
                logLH_prop += temp_prop;
                D1_prop += 1;
            }
            term = -pow(gsl_vector_get(gamma, i), nu3) * (kappa3) * pow(gsl_vector_get(yStar, i), alpha3) * exp(LP + temp_prop);
            logLH_prop +=  term;
            D1_prop += term;
            D2_prop += term;
        }

        D1_prop += -gsl_matrix_get(invSigma_V, 2, 2)*temp_prop - gsl_matrix_get(invSigma_V, 0, 2)*gsl_vector_get(Vj, 0) - gsl_matrix_get(invSigma_V, 1, 2)*gsl_vector_get(Vj, 1);
        D2_prop += -gsl_matrix_get(invSigma_V, 2, 2);
        
        v_prop_me_prop   = temp_prop - D1_prop/D2_prop;
        v_prop_var_prop  = - pow(2.4, 2)/D2_prop;
        
        startInx = endInx;
        
        c_dmvnorm2(Vj, zeroVec, 1, invSigma_V, &logPrior);
        c_dmvnorm2(Vj_prop, zeroVec, 1, invSigma_V, &logPrior_prop);
        
        logProp_PropToIni = dnorm(gsl_vector_get(V3, j), v_prop_me_prop, sqrt(v_prop_var_prop), 1);
        logProp_IniToProp = dnorm(temp_prop, v_prop_me, sqrt(v_prop_var), 1);
        
        logR = logLH_prop - logLH + logPrior_prop - logPrior + logProp_PropToIni - logProp_IniToProp;
        
        u = log(runif(0, 1)) <logR;
        
        if(u == 1)
        {
            gsl_vector_set(V3, j, temp_prop);
            gsl_vector_set(accept_V3, j, (gsl_vector_get(accept_V3, j) + u));
        }
    }
    
    gsl_vector_free(Vj);
    gsl_vector_free(Vj_prop);
    gsl_vector_free(zeroVec);
    gsl_matrix_free(invSigma_V);
    
    return;
}












/* updating variance covariance matrix of cluster-specific random effect: Sigma_V */

/**/



void BweibMvnCorScrSM_updateVP(gsl_vector *V1,
                             gsl_vector *V2,
                             gsl_vector *V3,
                             gsl_matrix *Sigma_V,
                             double rho_v,
                             gsl_matrix *Psi_v)
{
    int j, df;
    int J = V1 -> size;
    
    gsl_vector *Vj = gsl_vector_calloc(3);
    gsl_matrix *VV = gsl_matrix_calloc(3,3);
    gsl_matrix *Sum = gsl_matrix_calloc(3,3);
    gsl_matrix_memcpy(Sum, Psi_v);
    
    for(j = 0; j < J; j++)
    {
        gsl_vector_set(Vj, 0, gsl_vector_get(V1, j));
        gsl_vector_set(Vj, 1, gsl_vector_get(V2, j));
        gsl_vector_set(Vj, 2, gsl_vector_get(V3, j));
        
        gsl_blas_dger(1, Vj, Vj, VV);
        gsl_matrix_add(Sum, VV);
        gsl_matrix_set_zero(VV);
    }
    
    df = (int) rho_v + J;
    
        
    
    c_riwishart(df, Sum, Sigma_V);
    
    gsl_vector_free(Vj);
    gsl_matrix_free(VV);
    gsl_matrix_free(Sum);
    
    return;
}

