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

#include "BAFT_DPsurv.h"




/* updating mu and zeta using DPM of Normals */

void BAFT_DPsurv_update_mu_zeta(gsl_vector *yL,
                                gsl_vector *yU,
                                gsl_vector *yU_posinf,
                                gsl_vector *c0,
                                gsl_vector *c0_neginf,
                                gsl_matrix *X,
                                gsl_vector *y,
                                gsl_vector *beta,
                                gsl_vector *r,
                                gsl_vector *mu_all,
                                gsl_vector *zeta_all,
                                gsl_vector *rUniq,
                                gsl_vector *rUniq_count,
                                double tau,
                                double a0,
                                double b0,
                                double mu0,
                                double sigSq0,
                                double beta0_prop_var,
                                double zeta_prop_var,
                                int *accept_beta0,
                                int *accept_sigSq,
                                int *nClass_DP,
                                gsl_rng *rr)
{
    int i, j, k, n_ir, r_ind, wh_ri, u_, LL_negInf;
    double b_mc, sum_prob, val, mu, zeta, eta, cond;
    double mean, sd, mean_prop, sd_prop;
    
    int u = *nClass_DP;
    int n = y -> size;
    
    gsl_vector *xbeta = gsl_vector_calloc(n);
    gsl_blas_dgemv(CblasNoTrans, 1, X, beta, 0, xbeta);
    
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
        
        eta = gsl_vector_get(xbeta, i);
        LL_negInf = (int) gsl_vector_get(c0_neginf, i);
        
        /* identfy "u" unique values */
        c_uniq(r, rUniq, rUniq_count, mu_all, zeta_all, &u);
        
        for(j = 0; j < u; j++)
        {
            if(gsl_vector_get(r, i) == gsl_vector_get(rUniq, j)) wh_ri = j;
        }
        
        if(gsl_vector_get(rUniq_count, wh_ri) > 1) /* not singleton */
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
            zeta = rgamma(a0, 1/b0);
        }
        
        mu = NAN;
        while(mu == INFINITY || isnan(mu))
        {
            mu = rnorm(mu0, sqrt(sigSq0));
        }
        
        if(gsl_vector_get(rUniq_count, wh_ri) > 1) /* not singleton: proposing to create a new component */
        {
            for(j = 0; j < u_; j++)
            {
                n_ir = gsl_vector_get(rUniq_count, j);
                if(gsl_vector_get(r, i) == gsl_vector_get(rUniq, j)) n_ir -= 1;
                
                mean = gsl_vector_get(mu_all, j) + eta;
                sd = sqrt(1/gsl_vector_get(zeta_all, j));
                
                val = dnorm(gsl_vector_get(y, i), mean, sd, 0);
                
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
                
                val *= (double) n_ir / (double)(n-1+tau);
                gsl_vector_set(prob, j, val);
            }
            
            mean = mu + eta;
            sd = sqrt(1/zeta);
            
            val = dnorm(gsl_vector_get(y, i), mean, sd, 0);
            
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
            
            val *= (double) tau  / (double)(n-1+tau);
            
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
                    gsl_vector_set(r, i, gsl_vector_get(rUniq, r_ind-1));
                }else if(r_ind == u_+1)
                {
                    gsl_vector_set(r, i, gsl_vector_max(rUniq)+1);
                    gsl_vector_set(rUniq, u_, gsl_vector_max(rUniq)+1);
                    gsl_vector_set(mu_all, u_, mu);
                    gsl_vector_set(zeta_all, u_, zeta);
                    u += 1;
                }
            }
            
        }else /* singletion */
        {
            gsl_vector_memcpy(mu_temp, mu_all);
            gsl_vector_memcpy(zeta_temp, zeta_all);
            
            gsl_vector_set(mu_temp, wh_ri, mu);
            gsl_vector_set(zeta_temp, wh_ri, zeta);
            
            for(j = 0; j < u_+1; j++)
            {
                mean = gsl_vector_get(mu_temp, j) + eta;
                sd = sqrt(1/gsl_vector_get(zeta_temp, j));
                
                n_ir = gsl_vector_get(rUniq_count, j);
                if(gsl_vector_get(r, i) == gsl_vector_get(rUniq, j)) n_ir -= 1;
                
                val = dnorm(gsl_vector_get(y, i), mean, sd, 0);
                
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
                    val *= (double) tau / (double)(n-1+tau);
                }else
                {
                    val *= (double) n_ir / (double)(n-1+tau);
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
                
                gsl_vector_set(r, i, gsl_vector_get(rUniq, r_ind-1));
                
                if(r_ind-1 != wh_ri)
                {
                    if(wh_ri < u_)
                    {
                        gsl_vector_memcpy(rUniq_dummy, rUniq);
                        gsl_vector_memcpy(mu_dummy, mu_temp);
                        gsl_vector_memcpy(zeta_dummy, zeta_temp);
                        for(k=wh_ri;k<u_;k++)
                        {
                            gsl_vector_set(rUniq, k, gsl_vector_get(rUniq_dummy, k+1));
                            gsl_vector_set(mu_temp, k, gsl_vector_get(mu_dummy, k+1));
                            gsl_vector_set(zeta_temp, k, gsl_vector_get(zeta_dummy, k+1));
                        }
                    }
                    gsl_vector_set(rUniq, u_, 0);
                    gsl_vector_set(mu_temp, u_, 0);
                    gsl_vector_set(zeta_temp, u_, 0);
                    u -= 1;
                }
                
                gsl_vector_memcpy(mu_all, mu_temp);
                gsl_vector_memcpy(zeta_all, zeta_temp);
            }
        }
        
        gsl_vector_free(prob);
        
        for(j = u_+1; j < n; j++)
        {
            gsl_vector_set(rUniq, j, 0);
            gsl_vector_set(mu_all, j, 0);
            gsl_vector_set(zeta_all, j, 0);
        }
        
    }
    
    
    
    
    /*************************************************************/
    
    /* Step 2: update (beta0, sigSq) using the posterior distribution that is based on {epsilon_i :i∈{k:rk =r}}. */
    
    /*************************************************************/
    
    /* identfy "u" unique values */
    c_uniq(r, rUniq, rUniq_count, mu_all, zeta_all, &u);
    *nClass_DP = u;
    
    double mu_prop, zeta_prop, mu_ini, zeta_ini, temp, temp_prop;
    double loglh, loglh_prop, logprior, logprior_prop, logprop, logprop_prop;
    double logR;
    int uu;
    
    for(k = 0; k < 1; k++)
    {
        for(j = 0; j < u; j++)
        {
            mu_ini = gsl_vector_get(mu_all, j);
            zeta_ini = gsl_vector_get(zeta_all, j);
            
            /* update mu */
            mu_prop = rnorm((double) mu_ini, (double) sqrt(beta0_prop_var));
            
            loglh=0;
            loglh_prop=0;
            for(i = 0; i < n; i++)
            {
                if(gsl_vector_get(r, i) == gsl_vector_get(rUniq, j))
                {
                    LL_negInf = (int) gsl_vector_get(c0_neginf, i);
                    eta = gsl_vector_get(xbeta, i);
                    mean = mu_ini + eta;
                    mean_prop = mu_prop + eta;
                    sd = sqrt(1/zeta_ini);
                    
                    temp = dnorm(gsl_vector_get(y, i), mean, sd, 1);
                    temp_prop = dnorm(gsl_vector_get(y, i), mean_prop, sd, 1);
                    
                    if(LL_negInf == 0)
                    {
                        temp -= pnorm(gsl_vector_get(c0, i), mean, sd, 0, 1);
                        temp_prop -= pnorm(gsl_vector_get(c0, i), mean_prop, sd, 0, 1);
                    }
                    
                    loglh += temp;
                    loglh_prop += temp_prop;
                }
            }
            
            logprior = dnorm(mu_ini, mu0, sqrt(sigSq0), 1);
            logprior_prop = dnorm(mu_prop, mu0, sqrt(sigSq0), 1);
            
            logprop = dnorm(mu_ini, mu_prop, sqrt(beta0_prop_var), 1);
            logprop_prop = dnorm(mu_prop, mu_ini, sqrt(beta0_prop_var), 1);
            
            logR = loglh_prop - loglh + logprop_prop - logprior + logprop - logprop_prop;
            uu = log(runif(0, 1)) < logR;
            if(uu == 1)
            {
                gsl_vector_set(mu_all, j, mu_prop);
                *accept_beta0 += 1;
            }
            
            /* update zeta */
            mu_ini = gsl_vector_get(mu_all, j);
            
            zeta_prop = INFINITY;
            
            while(zeta_prop == INFINITY || isnan(zeta_prop))
            {
                zeta_prop = rgamma(pow(zeta_ini, 2)/zeta_prop_var, zeta_prop_var/zeta_ini);
            }
            
            loglh=0;
            loglh_prop=0;
            for(i = 0; i < n; i++)
            {
                if(gsl_vector_get(r, i) == gsl_vector_get(rUniq, j))
                {
                    LL_negInf = (int) gsl_vector_get(c0_neginf, i);
                    eta = gsl_vector_get(xbeta, i);
                    mean = mu_ini + eta;
                    sd = sqrt(1/zeta_ini);
                    sd_prop = sqrt(1/zeta_prop);
                    
                    temp = dnorm(gsl_vector_get(y, i), mean, sd, 1);
                    temp_prop = dnorm(gsl_vector_get(y, i), mean, sd_prop, 1);
                    
                    if(LL_negInf == 0)
                    {
                        temp -= pnorm(gsl_vector_get(c0, i), mean, sd, 0, 1);
                        temp_prop -= pnorm(gsl_vector_get(c0, i), mean, sd_prop, 0, 1);
                    }
                    
                    loglh += temp;
                    loglh_prop += temp_prop;
                }
            }
            logprior = dgamma(zeta_ini, a0, 1/b0, 1);
            logprior_prop = dgamma(zeta_prop, a0, 1/b0, 1);
            
            logprop += dgamma(zeta_ini, pow(zeta_prop, 2)/zeta_prop_var, zeta_prop_var/zeta_prop, 1);
            logprop_prop += dgamma(zeta_prop, pow(zeta_ini, 2)/zeta_prop_var, zeta_prop_var/zeta_ini, 1);
            
            logR = loglh_prop - loglh + logprior_prop - logprior + logprop - logprop_prop;
            uu = log(runif(0, 1)) < logR;
            if(uu == 1)
            {
                gsl_vector_set(zeta_all, j, zeta_prop);
                *accept_sigSq += 1;
            }
            
        }
    }
    
    gsl_vector_free(mu_dummy);
    gsl_vector_free(zeta_dummy);
    gsl_vector_free(rUniq_dummy);
    gsl_vector_free(xbeta);
    gsl_vector_free(mu_temp);
    gsl_vector_free(zeta_temp);
    
    return;
    
}






/* Updating y */

void BAFT_DPsurv_update_y(gsl_vector *yL,
                          gsl_vector *yU,
                          gsl_vector *yL_neginf,
                          gsl_vector *yU_posinf,
                          gsl_vector *c0,
                          gsl_matrix *X,
                          gsl_vector *y,
                          gsl_vector *beta,
                          gsl_vector *r,
                          gsl_vector *mu_all,
                          gsl_vector *zeta_all,
                          gsl_vector *rUniq,
                          gsl_vector *rUniq_count,
                          int *nClass_DP,
                          gsl_rng *rr)
{
    double eta, sample, mixProb_den, mu_k, sigSq_k;
    int i, k, selInd;
    int n = y -> size;
    
    int u = *nClass_DP;
    
    gsl_vector *mixProb = gsl_vector_calloc(u);
    gsl_vector *xbeta = gsl_vector_calloc(n);
    gsl_blas_dgemv(CblasNoTrans, 1, X, beta, 0, xbeta);
    
    for(i=0;i<n;i++)
    {
        eta = gsl_vector_get(xbeta, i);
        
        /* interval censored or right censored */
        if(gsl_vector_get(yU, i) != gsl_vector_get(yL, i))
        {
            mixProb_den = 0;
            for(k=0;k<u;k++)
            {
                mu_k = (double) gsl_vector_get(mu_all, k);
                sigSq_k = (double) pow(gsl_vector_get(zeta_all, k), -1);
                
                if(gsl_vector_get(yL_neginf, i) == 0 && gsl_vector_get(yU_posinf, i) == 0)
                {
                    gsl_vector_set(mixProb, k, gsl_vector_get(rUniq_count, k) * (pnorm(gsl_vector_get(yU, i), eta + mu_k, sqrt(sigSq_k), 1, 0) - pnorm(gsl_vector_get(yL, i), eta + mu_k, sqrt(sigSq_k), 1, 0)));
                }else if(gsl_vector_get(yL_neginf, i) == 0 && gsl_vector_get(yU_posinf, i) == 1)
                {
                    gsl_vector_set(mixProb, k, gsl_vector_get(rUniq_count, k) * pnorm(gsl_vector_get(yL, i), eta + mu_k, sqrt(sigSq_k), 0, 0));
                }else if(gsl_vector_get(yL_neginf, i) == 1 && gsl_vector_get(yU_posinf, i) == 0)
                {
                    gsl_vector_set(mixProb, k, gsl_vector_get(rUniq_count, k) * pnorm(gsl_vector_get(yU, i), eta + mu_k, sqrt(sigSq_k), 1, 0));
                }else if(gsl_vector_get(yL_neginf, i) == 1 && gsl_vector_get(yU_posinf, i) == 1)
                {
                    gsl_vector_set(mixProb, k, gsl_vector_get(rUniq_count, k));
                }
                
                mixProb_den += gsl_vector_get(mixProb, k);
            }
            
            if(mixProb_den == 0)
            {
                for(k=0; k<u; k++)
                {
                    gsl_vector_set(mixProb, k, (double) 1/u);
                }
            }else
            {
                gsl_vector_scale(mixProb, 1/mixProb_den);
            }
            
            selInd = c_multinom_sample(rr, mixProb, u);
            
            c_rtnorm(eta + gsl_vector_get(mu_all, selInd-1), sqrt(1/gsl_vector_get(zeta_all, selInd-1)), gsl_vector_get(yL, i), gsl_vector_get(yU, i), gsl_vector_get(yL_neginf, i), gsl_vector_get(yU_posinf, i), &sample);
            gsl_vector_set(y, i, sample);
            /*  observed */
            
        }else if(gsl_vector_get(yU, i) == gsl_vector_get(yL, i))
        {
            gsl_vector_set(y, i, gsl_vector_get(yU, i));
        }
    }
    
    gsl_vector_free(xbeta);
    gsl_vector_free(mixProb);
    return;
}







/* Updating beta */

void BAFT_DPsurv_update_beta(gsl_vector *yL,
                             gsl_vector *yU,
                             gsl_vector *yU_posinf,
                             gsl_vector *c0,
                             gsl_vector *c0_neginf,
                             gsl_matrix *X,
                             gsl_vector *y,
                             gsl_vector *beta,
                             gsl_vector *r,
                             gsl_vector *mu_all,
                             gsl_vector *zeta_all,
                             gsl_vector *rUniq,
                             gsl_vector *rUniq_count,
                             int *nClass_DP,
                             double beta_prop_var,
                             gsl_vector *accept_beta)
{
    int i, j, uu;
    double eta, eta_prop, loglh, loglh_prop, logR;
    double temp, temp_prop;
    
    int n = X -> size1;
    int p = X -> size2;
    
    int u = *nClass_DP;
    
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
        eta = gsl_vector_get(xbeta, i);
        eta_prop = gsl_vector_get(xbeta_prop, i);
        
        if(gsl_vector_get(c0_neginf, i) == 0)
        {
            temp = fmixTN(gsl_vector_get(y, i), gsl_vector_get(c0, i), 1000000, 0, 1, mu_all, zeta_all, rUniq_count, u, eta);
            temp_prop = fmixTN(gsl_vector_get(y, i), gsl_vector_get(c0, i), 1000000, 0, 1, mu_all, zeta_all, rUniq_count, u, eta_prop);
        }else
        {
            temp = fmixTN(gsl_vector_get(y, i), -1000000, 1000000, 1, 1, mu_all, zeta_all, rUniq_count, u, eta);
            temp_prop = fmixTN(gsl_vector_get(y, i), -1000000, 1000000, 1, 1, mu_all, zeta_all, rUniq_count, u, eta_prop);
        }
        loglh += log(temp);
        loglh_prop += log(temp_prop);
    }
    
    logR = loglh_prop - loglh;
    uu = log(runif(0, 1)) < logR;
    if(uu == 1)
    {
        gsl_vector_memcpy(beta, beta_prop);
        gsl_vector_set(accept_beta, j, gsl_vector_get(accept_beta, j) + 1);
    }
    
    gsl_vector_free(beta_prop);
    gsl_vector_free(xbeta);
    gsl_vector_free(xbeta_prop);
    
    return;
}



/* updating precision parameter of DP prior: tau */

void BAFT_DPsurv_update_tau(int *n,
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












