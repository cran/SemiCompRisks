#include <stdio.h>
#include <math.h>

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"

#include "R.h"
#include "Rmath.h"

#include "BweibCorSurv.h"






/* updating regression parameter: beta */

/**/

void BweibCorSurv_updateRP(gsl_vector *beta,
                           double *alpha,
                           double *kappa,
                           gsl_vector *V,
                           gsl_vector *survTime,
                           gsl_vector *survEvent,
                           gsl_vector *cluster,
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
    int i, j, jj;
    
    gsl_vector *beta_prop = gsl_vector_calloc(p);
        
    j = (int) runif(0, p);
    
    logLH = 0; D1 = 0; D2 = 0;
    logLH_prop = 0; D1_prop = 0; D2_prop = 0;
    
    for(i = 0; i < n; i++)
    {
        gsl_vector_view Xi = gsl_matrix_row(survCov, i);
        gsl_blas_ddot(&Xi.vector, beta, &LP);
        jj = (int) gsl_vector_get(cluster, i) - 1;
        if(gsl_vector_get(survEvent, i) == 1)
        {
            logLH   += LP;
            D1      += gsl_matrix_get(survCov, i, j);
        }
        logLH   += -(*kappa) * pow(gsl_vector_get(survTime, i), *alpha) * exp(LP + gsl_vector_get(V, jj));
        D1      += -(*kappa) * pow(gsl_vector_get(survTime, i), *alpha) * exp(LP + gsl_vector_get(V, jj)) * gsl_matrix_get(survCov, i, j);
        D2      += -(*kappa) * pow(gsl_vector_get(survTime, i), *alpha) * exp(LP + gsl_vector_get(V, jj)) * pow(gsl_matrix_get(survCov, i, j), 2);
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
        jj = (int) gsl_vector_get(cluster, i) - 1;
        if(gsl_vector_get(survEvent, i) == 1)        
        {
            logLH_prop   += LP_prop;
            D1_prop      += gsl_matrix_get(survCov, i, j);
        }
        logLH_prop   += -(*kappa) * pow(gsl_vector_get(survTime, i), *alpha) * exp(LP_prop + gsl_vector_get(V, jj));
        D1_prop      += -(*kappa) * pow(gsl_vector_get(survTime, i), *alpha) * exp(LP_prop + gsl_vector_get(V, jj)) * gsl_matrix_get(survCov, i, j);
        D2_prop      += -(*kappa) * pow(gsl_vector_get(survTime, i), *alpha) * exp(LP_prop + gsl_vector_get(V, jj)) * pow(gsl_matrix_get(survCov, i, j), 2);
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

/* use the random walk proposal with log transformation*/

/**/

void BweibSurv_updateSH_rw2(gsl_vector *beta,
                        double *alpha,
                        double *kappa,
                        gsl_vector *V,
                        gsl_vector *survTime,
                        gsl_vector *survEvent,
                        gsl_vector *cluster,
                        gsl_matrix *survCov,
                        double mhProp_alpha_var,
                        double a,
                        double b,
                        int *accept_alpha)
{
    double z1, LP, logLH, logLH_prop;
    double temp_prop;
    double logPrior, logPrior_prop;
    double logProp_IniToProp, logProp_PropToIni;
    double logR;
    int u;
    int n = survTime -> size;
    int i, jj;
    
    z1 = log(*alpha);
    
    logLH = 0; logLH_prop = 0;
    temp_prop = rnorm(z1, sqrt(mhProp_alpha_var));
    
    for(i = 0; i < n; i++)
    {
        gsl_vector_view Xi = gsl_matrix_row(survCov, i);
        gsl_blas_ddot(&Xi.vector, beta, &LP);
        jj = (int) gsl_vector_get(cluster, i) - 1;
        if(gsl_vector_get(survEvent, i) == 1)
        {
            logLH       += log(*alpha) + (*alpha - 1)*log(gsl_vector_get(survTime, i));
            logLH_prop  += temp_prop + (exp(temp_prop) - 1)*log(gsl_vector_get(survTime, i));
        }
        logLH       += -(*kappa) * pow(gsl_vector_get(survTime, i), *alpha) * exp(LP + gsl_vector_get(V, jj));
        logLH_prop  += -(*kappa) * pow(gsl_vector_get(survTime, i), exp(temp_prop)) * exp(LP + gsl_vector_get(V, jj));
    }
    
    logPrior        = a*z1 - b* (*alpha);
    logPrior_prop   = a*temp_prop - b* (exp(z1));
    
    logProp_PropToIni = dnorm(z1, temp_prop, sqrt(mhProp_alpha_var), 1);
    logProp_IniToProp = dnorm(temp_prop, z1, sqrt(mhProp_alpha_var), 1);
    
    logR = logLH_prop - logLH + logPrior_prop - logPrior + logProp_PropToIni - logProp_IniToProp;
    
    
    u = log(runif(0, 1)) <logR;
       
    if(u == 1)
    {
        *alpha = exp(temp_prop);
        *accept_alpha += u;
    }
    
    return;
    
}




/* updating scale parameter: kappa */

/**/

void BweibSurv_updateSC(gsl_vector *beta,
                        double *alpha,
                        double *kappa,
                        gsl_vector *V,
                        gsl_vector *survTime,
                        gsl_vector *survEvent,
                        gsl_vector *cluster,
                        gsl_matrix *survCov,
                        double c,
                        double d)
{
    int n = survTime -> size;
    int i, jj;
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
        jj = (int) gsl_vector_get(cluster, i) - 1;
        
        Kappa_rate += pow(gsl_vector_get(survTime, i), *alpha) * exp(LP + gsl_vector_get(V, jj));
    }
    Kappa_rate += d;
    Kappa_scale = 1/Kappa_rate;
    
    *kappa = rgamma(Kappa_shape, Kappa_scale);
    return;
}






/* updating cluster-specific random effect using RW-MH: V */

/**/



void BweibSurv_updateCP(gsl_vector *beta,
                        double alpha,
                        double kappa,
                        gsl_vector *V,
                        double zeta,
                        gsl_vector *survTime,
                        gsl_vector *survEvent,
                        gsl_vector *cluster,
                        gsl_matrix *survCov,
                        gsl_vector *n_j,
                        gsl_vector *accept_V,
                        double mhProp_V_var)
{
    double LP, logLH, logLH_prop;
    double logPrior, logPrior_prop, temp_prop;
    double logProp_IniToProp, logProp_PropToIni;
    double logR;
    int u, i, j;

    int J = V  -> size;
    
    int startInx = 0;
    int endInx = 0;
    
    for(j = 0; j < J; j++)
    {
        logLH = 0; logLH_prop = 0;
        
        temp_prop = rnorm(gsl_vector_get(V, j), sqrt(mhProp_V_var));
        
        endInx += (int) gsl_vector_get(n_j, j);
        
        for(i = startInx; i < endInx; i++)
        {
            gsl_vector_view Xi = gsl_matrix_row(survCov, i);
            gsl_blas_ddot(&Xi.vector, beta, &LP);
            
            if(gsl_vector_get(survEvent, i) == 1)
            {
                logLH += gsl_vector_get(V, j);
                logLH_prop += temp_prop;
            }
            logLH +=  -(kappa) * pow(gsl_vector_get(survTime, i), alpha) * exp(LP + gsl_vector_get(V, j));
            logLH_prop +=  -(kappa) * pow(gsl_vector_get(survTime, i), alpha) * exp(LP + temp_prop);
        }
        startInx = endInx;
        
        logPrior        = -zeta*pow(gsl_vector_get(V, j), 2)/2;
        logPrior_prop   = -zeta*pow(temp_prop, 2)/2;
        
        logProp_PropToIni = dnorm(gsl_vector_get(V, j), temp_prop, sqrt(mhProp_V_var), 1);
        logProp_IniToProp = dnorm(temp_prop, gsl_vector_get(V, j), sqrt(mhProp_V_var), 1);
        
        logR = logLH_prop - logLH + logPrior_prop - logPrior + logProp_PropToIni - logProp_IniToProp;
        
        u = log(runif(0, 1)) <logR;
        
        if(u == 1)
        {
            gsl_vector_set(V, j, temp_prop);
            gsl_vector_set(accept_V, j, (gsl_vector_get(accept_V, j) + u));
        }
    }
    
    return;
}










/* updating variance covariance matrix of cluster-specific random effect: Sigma_V */

/**/



void BweibSurv_updateVP(gsl_vector *V,
                       double *zeta,
                       double rho1,
                       double rho2)
{
    int j;
    int J = V -> size;
    double zeta_shape, zeta_rate, zeta_scale;
    
    zeta_shape = rho1 + (double) J/2;
    
    zeta_rate = 0;
    
    for(j = 0; j < J; j++)
    {
        zeta_rate += pow(gsl_vector_get(V, j), 2);
    }
    zeta_rate /= 2;
    zeta_rate += rho2;
    zeta_scale = 1/zeta_rate;
    
    *zeta = rgamma(zeta_shape, zeta_scale);
    
   
    return;
}






