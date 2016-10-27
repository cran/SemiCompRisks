/*
 TO COMPILE USE THE CODE:
 
 R CMD SHLIB BAFT_DPuni.c BAFT_DPuni_Updates.c BAFT_DPuni_Utilities.c -lgsl -lgslcblas
 
 */

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
#include "BAFT_DPsurv.h"


/* */
void BAFT_DPunimcmc(double Ymat[],
                    double yLInf[],
                    double yUInf[],
                    double c0Inf[],
                    double Xmat[],
                    int *n,
                    int *p,
                    double hyperP[],
                    double mcmcP[],
                    double startValues[],
                    int *numReps,
                    int *thin,
                    double *burninPerc,
                    double samples_y[],
                    double samples_beta[],
                    double samples_r[],
                    double samples_mu[],
                    double samples_sigSq[],
                    double samples_tau[],
                    double samples_misc[])
{
    GetRNGstate();
    time_t now;
    int i, j, M;
    
    const gsl_rng_type * TT;
    gsl_rng * rr;
    
    gsl_rng_env_setup();
    
    TT = gsl_rng_default;
    rr = gsl_rng_alloc(TT);
    
    /* Data */
    
    gsl_vector *yL = gsl_vector_calloc(*n);
    gsl_vector *yL_neginf = gsl_vector_calloc(*n);
    gsl_vector *yU = gsl_vector_calloc(*n);
    gsl_vector *yU_posinf = gsl_vector_calloc(*n);
    gsl_vector *c0 = gsl_vector_calloc(*n);
    gsl_vector *c0_neginf = gsl_vector_calloc(*n);
    
    int nP;
    
    if(*p > 0) nP = *p;
    if(*p == 0) nP = 1;
    
    gsl_matrix *X = gsl_matrix_calloc(*n, (nP));
    
    
    for(i = 0; i < *n; i++)
    {
        gsl_vector_set(yL, i, Ymat[(0 * *n) + i]);
        gsl_vector_set(yU, i, Ymat[(1 * *n) + i]);
        gsl_vector_set(c0, i, Ymat[(2 * *n) + i]);
        gsl_vector_set(yL_neginf, i, yLInf[i]);
        gsl_vector_set(yU_posinf, i, yUInf[i]);
        gsl_vector_set(c0_neginf, i, c0Inf[i]);
        if(*p >0)
        {
            for(j = 0; j < *p; j++)
            {
                gsl_matrix_set(X, i, j, Xmat[(j* *n) + i]);
            }
        }
    }
    
    /* Hyperparameters */
    
    double a0 = hyperP[0];
    double b0 = hyperP[1];
    double aTau = hyperP[2];
    double bTau = hyperP[3];
    double mu0 = hyperP[4];
    double sigSq0  = hyperP[5];
    
    double beta_prop_var = mcmcP[0];
    double beta0_prop_var = mcmcP[1];
    double zeta_prop_var = mcmcP[2];
    
    /* Starting values */
    
    gsl_vector *y = gsl_vector_calloc(*n);
    gsl_vector *beta = gsl_vector_calloc(nP);
    gsl_vector *r = gsl_vector_calloc(*n);
    double tau = startValues[*n+*p+*n];
    
    for(i = 0; i < *n; i++)
    {
        gsl_vector_set(y, i, startValues[i]);
    }
    if(*p > 0)
    {
        for(i = 0; i < *p; i++) gsl_vector_set(beta, i, startValues[*n+i]);
    }
    for(i = 0; i < *n; i++)
    {
        gsl_vector_set(r, i, startValues[*n+*p+i]);
    }
    
    
    /* Variables required for storage of samples */
    
    int StoreInx;
    gsl_vector *accept_beta = gsl_vector_calloc(nP);
    int accept_beta0 = 0;
    int accept_sigSq = 0;
    
    gsl_vector *mu_all = gsl_vector_calloc(*n);
    gsl_vector *zeta_all = gsl_vector_calloc(*n);
    
    gsl_vector *rUniq = gsl_vector_calloc(*n);
    gsl_vector *rUniq_count = gsl_vector_calloc(*n);
    
    int nClass_DP;
    
    c_uniq1(r, rUniq, rUniq_count, &nClass_DP);
    
    for(i = 0; i < nClass_DP; i++)
    {
        gsl_vector_set(mu_all, i, startValues[*n+*p+*n+1+i]);
        gsl_vector_set(zeta_all, i, startValues[*n+*p+*n+1+*n+i]);
    }
    
    
    /* selecting a move */
    /* move: 1= DPM, 2= tau, 3= y, 4= beta */
    
    double pDPM = 0.2;
    double pTau = 0.2;
    double pY = 0.1;
    double pBeta = 1 - (pDPM + pTau + pY);
    int move;
    gsl_vector *moveProb = gsl_vector_calloc(4);
    gsl_vector_set(moveProb, 0, pDPM);
    gsl_vector_set(moveProb, 1, pTau);
    gsl_vector_set(moveProb, 2, pY);
    gsl_vector_set(moveProb, 3, pBeta);
    
    BAFT_DPsurv_update_y(yL, yU, yL_neginf, yU_posinf, c0, X, y, beta, r, mu_all, zeta_all, rUniq, rUniq_count, &nClass_DP, rr);
    
    
    for(M = 0; M < *numReps; M++)
    {
        move = c_multinom_sample(rr, moveProb, 4);
        
        if(move == 1)
        {
            /* update mu and zeta         */
            BAFT_DPsurv_update_mu_zeta(yL, yU, yU_posinf, c0, c0_neginf, X, y, beta, r, mu_all, zeta_all, rUniq, rUniq_count, tau, a0, b0, mu0, sigSq0, beta0_prop_var, zeta_prop_var, &accept_beta0, &accept_sigSq, &nClass_DP, rr);
        }
        
        
        if(move == 2)
        {
            /* update tau */
            BAFT_DPsurv_update_tau(n, &tau, aTau, bTau, &nClass_DP);
        }
        
        
        
        if(move == 3)
        {
            /* Updating y  */
            BAFT_DPsurv_update_y(yL, yU, yL_neginf, yU_posinf, c0, X, y, beta, r, mu_all, zeta_all, rUniq, rUniq_count, &nClass_DP, rr);
        }
        
        
        
        if(move == 4)
        {
            /* Updating beta */
            if(*p > 0)
            {
                BAFT_DPsurv_update_beta(yL, yU, yU_posinf, c0, c0_neginf, X, y, beta, r, mu_all, zeta_all, rUniq, rUniq_count, &nClass_DP, beta_prop_var, accept_beta);
            }
        }
        
        
        /* Storing posterior samples */
        
        if( ( (M+1) % *thin ) == 0 && (M+1) > (*numReps * *burninPerc))
        {
            StoreInx = (M+1)/(*thin)- (*numReps * *burninPerc)/(*thin);
            
            for(i = 0; i < *n; i++) samples_y[(StoreInx - 1) * (*n) + i] = gsl_vector_get(y, i);
            
            if(*p >0)
            {
                for(j = 0; j < *p; j++) samples_beta[(StoreInx - 1) * (*p) + j] = gsl_vector_get(beta, j);
            }
            samples_tau[StoreInx - 1] = tau;
            
            for(j = 0; j < *n; j++) samples_r[(StoreInx - 1) * (*n) + j] = gsl_vector_get(r, j);
            
            for(j = 0; j < *n; j++) samples_mu[(StoreInx - 1) * (*n) + j] = gsl_vector_get(mu_all, j);
            
            for(j = 0; j < *n; j++) samples_sigSq[(StoreInx - 1) * (*n) + j] = 1/gsl_vector_get(zeta_all, j);
        }
        
        if(M == (*numReps - 1))
        {
            for(i = 0; i < *p; i++)
            {
                samples_misc[i] = (int) gsl_vector_get(accept_beta, i);
            }
            
            samples_misc[*p] = (int) accept_beta0;
            samples_misc[*p+1] = (int) accept_sigSq;
        }
        
        if( ( (M+1) % 1000 ) == 0)
        {
            time(&now);
            Rprintf("iteration: %d: %s\n", M+1, ctime(&now));
            R_FlushConsole();
            R_ProcessEvents();
        }
    }
    PutRNGstate();
    return;
}





















