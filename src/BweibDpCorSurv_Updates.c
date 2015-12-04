#include <stdio.h>
#include <math.h>

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_sf.h"

#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

#include "R.h"
#include "Rmath.h"

#include "BweibDpCorSurv.h"






/* updating regression parameter: beta */

/**/

void BweibDpCorSurv_updateRP(gsl_vector *beta,
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

void BweibDpCorSurv_updateSH_rw2(gsl_vector *beta,
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

void BweibDpCorSurv_updateSC(gsl_vector *beta,
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

















/* updating cluster-specific random effects
 : the prior is used for the proposal density */

void BweibDpCorSurv_updateCP(gsl_vector *beta,
                             double alpha,
                             double kappa,
                           gsl_vector *V,
                           gsl_vector *survTime,
                           gsl_vector *survEvent,
                           gsl_vector *cluster,
                           gsl_matrix *survCov,
                           gsl_vector *n_j,
                           gsl_vector *mu_all,
                           gsl_vector *zeta_all,
                           gsl_vector *c,
                           gsl_vector *accept_V,
                             double mhProp_V_var,
                           double mu0,
                           double zeta0,
                           double a0,
                           double b0,
                           double tau,
                           int *nClass_DP,
                           gsl_rng *rr)
{
    int i, j, jj, k, u, n_jc, c_ind;
    double prob2, b_mc, sum_prob, val, mu, zeta;
    

    int J = V -> size;
    
    
    gsl_vector *cUniq = gsl_vector_calloc(J);
    gsl_vector *cUniq_count = gsl_vector_calloc(J);
    gsl_vector *cTemp = gsl_vector_calloc(J);
    
    gsl_vector *prob1 = gsl_vector_calloc(J+1);
    
    double Vbar, Vsum, muA, zetaA, aA, bA, tempSum;

    
    
    gsl_vector_set_zero(mu_all);
    gsl_vector_set_zero(zeta_all);
    
    
    
    /*************************************************************/
    
    /* Step 1: updating latent classes (c) */
    
    /*************************************************************/
    
    
    for(jj = 0; jj < J; jj++)
    {
        
        /* identfy "u" unique values */
        
        gsl_vector_memcpy(cTemp, c);
        
        u = 1;
        
        for(i = 0; i < J; i++)
        {
            if(i == 0)
            {
                gsl_vector_set(cUniq, u-1, gsl_vector_get(cTemp, i));
                
                for(j = i; j < J; j++)
                {
                    if(gsl_vector_get(cTemp, j) == gsl_vector_get(cUniq, u-1))
                    {
                        gsl_vector_set(cUniq_count, u-1, gsl_vector_get(cUniq_count, u-1)+1);
                        gsl_vector_set(cTemp, j, 0);
                    }
                }
            }
            
            if(i != 0 && gsl_vector_get(cTemp, i) != 0)
            {
                u += 1;
                gsl_vector_set(cUniq, u-1, gsl_vector_get(cTemp, i));
                
                for(j = i; j < J; j++)
                {
                    if(gsl_vector_get(cTemp, j) == gsl_vector_get(cUniq, u-1))
                    {
                        gsl_vector_set(cUniq_count, u-1, gsl_vector_get(cUniq_count, u-1)+1);
                        gsl_vector_set(cTemp, j, 0);
                    }
                }
                
            }
        }

        
        /* calculating probabilities for each of "u" latent class */
        
        
        
        gsl_vector_set_zero(prob1);
        
        for(j = 0; j < u; j++)
        {
            
            n_jc = gsl_vector_get(cUniq_count, j);
            
            if(gsl_vector_get(c, jj) == gsl_vector_get(cUniq, j)) n_jc -= 1;
            
            zetaA = zeta0 + gsl_vector_get(cUniq_count, j);
            
            aA = a0 + gsl_vector_get(cUniq_count, j)/2;
            
            Vsum = 0;
            
            for(k = 0; k < J; k++)
            {
                if(gsl_vector_get(c, k) == gsl_vector_get(cUniq, j) && k != jj)
                {
                    Vsum += gsl_vector_get(V, k);
                }
            }
            Vbar = Vsum;
            
            if(n_jc != 0)
            {
                Vbar = Vsum/n_jc;
            }
            
            muA = (mu0*zeta0 + Vsum) / (zeta0 + gsl_vector_get(cUniq_count, j));
                        
            bA = b0;
            bA += (zeta0*gsl_vector_get(cUniq_count, j)*pow(Vbar - mu0, 2))/ (2*(zeta0 + gsl_vector_get(cUniq_count, j)));
            
            for(k = 0; k < J; k++)
            {
                if(gsl_vector_get(c, k) == gsl_vector_get(cUniq, j))
                {
                    tempSum = gsl_vector_get(V, k) - Vbar;
                    bA += pow(tempSum, 2)/2;
                }
            }
            
            
            val = (double) n_jc / (double)(J-1+tau) * Qfunc_univ(gsl_vector_get(V, jj), muA, zetaA, aA, bA);
            
            gsl_vector_set(prob1, j, val);
            
        }
        
        
        prob2 = tau/(double)(J-1+tau) * Qfunc_univ(gsl_vector_get(V, jj), mu0, zeta0, a0, b0);
        
        sum_prob = 0;
        for(j = 0; j < u; j++) sum_prob += gsl_vector_get(prob1, j);
        sum_prob += prob2;
        b_mc = 1/sum_prob;
        
        gsl_vector_scale(prob1, b_mc);
        prob2 *= b_mc;
        
        gsl_vector_set(prob1, u, prob2);
        
            
        /* sample c based on the probabilities */
        
        c_ind = c_multinom_sample(rr, prob1, u+1);
        
        
        if(c_ind <= u)
        {
            gsl_vector_set(c, jj, gsl_vector_get(cUniq, c_ind-1));
        }
        if(c_ind > u)
        {
            gsl_vector_set(c, jj, gsl_vector_max(cUniq)+1);
        }
                
        
        /* initializing vectors and matrices */
        
        gsl_vector_set_zero(cUniq);
        gsl_vector_set_zero(cUniq_count);
        
    }
    
  
    
    /*************************************************************/
    
    /* Step 2: update (mu, zeta) using the posterior distribution that is based on {Vj :j∈{k:ck =c}}. */
    
    /*************************************************************/
    
    
    
    /* identfy "u" unique values */
    
    gsl_vector_memcpy(cTemp, c);
    
    u = 1;
    
    for(i = 0; i < J; i++)
    {
        if(i == 0)
        {
            gsl_vector_set(cUniq, u-1, gsl_vector_get(cTemp, i));
            
            for(j = i; j < J; j++)
            {
                if(gsl_vector_get(cTemp, j) == gsl_vector_get(cUniq, u-1))
                {
                    gsl_vector_set(cUniq_count, u-1, gsl_vector_get(cUniq_count, u-1)+1);
                    gsl_vector_set(cTemp, j, 0);
                }
            }
        }
        
        if(i != 0 && gsl_vector_get(cTemp, i) != 0)
        {
            u += 1;
            gsl_vector_set(cUniq, u-1, gsl_vector_get(cTemp, i));
            
            for(j = i; j < J; j++)
            {
                if(gsl_vector_get(cTemp, j) == gsl_vector_get(cUniq, u-1))
                {
                    gsl_vector_set(cUniq_count, u-1, gsl_vector_get(cUniq_count, u-1)+1);
                    gsl_vector_set(cTemp, j, 0);
                }
            }
            
        }
    }
    
      
    *nClass_DP = u;
    
    
    
    for(j = 0; j < u; j++)
    {
        
        n_jc = gsl_vector_get(cUniq_count, j);        
        
        zetaA = zeta0 + gsl_vector_get(cUniq_count, j);
        
        aA = a0 + gsl_vector_get(cUniq_count, j)/2;

        Vsum = 0;
        
        for(k = 0; k < J; k++)
        {
            if(gsl_vector_get(c, k) == gsl_vector_get(cUniq, j))
            {
                Vsum += gsl_vector_get(V, k);
            }
        }
        Vbar = Vsum;
        
        if(n_jc != 0)
        {
            Vbar = Vsum/n_jc;
        }
        
        muA = (mu0*zeta0 + Vsum) / zetaA;
        
        
        bA = b0;
        bA += (zeta0*gsl_vector_get(cUniq_count, j)*pow(Vbar - mu0, 2))/ (2*(zeta0 + gsl_vector_get(cUniq_count, j)));
        
        for(k = 0; k < J; k++)
        {
            if(gsl_vector_get(c, k) == gsl_vector_get(cUniq, j))
            {
                tempSum = gsl_vector_get(V, k) - Vbar;
                bA += pow(tempSum, 2)/2;
            }
        }
        
        zeta = rgamma(aA, 1/bA);
        mu = rnorm(muA, sqrt(1/(zetaA*zeta)));
        
        gsl_vector_set(mu_all, j, mu);
        gsl_vector_set(zeta_all, j, zeta);
        
    }
    
    
    
    /*************************************************************/
    
    /* Step 3: updating the Vj using MH algorithm */
    
    /*************************************************************/
    
    
    double LP, logLH, logLH_prop;
    double logPrior, logPrior_prop, temp_prop;
    double logProp_IniToProp, logProp_PropToIni;
    double logR, mu_temp, zeta_temp;
    int uu;
    
    int startInx = 0;
    int endInx = 0;
    
    for(j = 0; j < J; j++)
    {
        for(i = 0; i < u; i++)
        {
            if(gsl_vector_get(c, j) == gsl_vector_get(cUniq, i)) jj = i;
        }
        
        mu_temp = gsl_vector_get(mu_all, jj);
        zeta_temp = gsl_vector_get(zeta_all, jj);
        
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

        } /* the end of the loop with i */
               
        startInx = endInx;
        
        logPrior = -(double) zeta_temp/2*pow(gsl_vector_get(V, j), 2);
        logPrior_prop = -(double) zeta_temp/2*pow(temp_prop, 2);
        
        logProp_PropToIni = dnorm(gsl_vector_get(V, j), temp_prop, sqrt(mhProp_V_var), 1);
        logProp_IniToProp = dnorm(temp_prop, gsl_vector_get(V, j), sqrt(mhProp_V_var), 1);;
        
        logR = logLH_prop - logLH + logPrior_prop - logPrior + logProp_PropToIni - logProp_IniToProp;
        
        uu = log(runif(0, 1)) < logR;
                
        if(uu == 1)
        {
            gsl_vector_set(V, j, temp_prop);
            gsl_vector_set(accept_V, j, (gsl_vector_get(accept_V, j) + uu));
        }
        

        
    }/* the end of the loop with j */
    
    
    
    gsl_vector_free(cUniq);
    gsl_vector_free(cUniq_count);
    gsl_vector_free(cTemp);
    gsl_vector_free(prob1);
    
    
    return;
    
    
    
    
}













/* updating precision parameter of DP prior: tau */

void BweibDpCorSurv_updatePP(int *n,
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





