/*
 TO COMPILE USE THE CODE:
 
 R CMD SHLIB BAFTuni.c BAFTuni_Updates.c BAFTuni_Utilities.c -lgsl -lgslcblas
 
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
#include "BAFT_LNsurv.h"


/* */
void BAFTunimcmc(double Ymat[],
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
                 double samples_beta0[],
                 double samples_sigSq[],
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
    
    double a_sigSq = hyperP[0];
    double b_sigSq = hyperP[1];
    
    double beta_prop_var = mcmcP[0];
    double beta0_prop_var = mcmcP[1];
    double sigSq_prop_var = mcmcP[2];
    
    /* Starting values */
    
    gsl_vector *y = gsl_vector_calloc(*n);
    gsl_vector *beta = gsl_vector_calloc(nP);
    double beta0 = startValues[*n+*p];
    double sigSq = startValues[*n+*p+1];
    
    for(i = 0; i < *n; i++)
    {
        gsl_vector_set(y, i, startValues[i]);
    }
    if(*p > 0)
    {
        for(i = 0; i < *p; i++) gsl_vector_set(beta, i, startValues[*n+i]);
    }
    
    
    
    /* Variables required for storage of samples */
    
    int StoreInx;
    gsl_vector *accept_beta = gsl_vector_calloc(nP);
    int accept_beta0 = 0;
    int accept_sigSq = 0;
    
    /* selecting a move */
    /* move: 1= y, 2= beta, 3= beta0, 4= sigSq */
    
    double pY = 0.1;
    double pBeta = 0.3;
    double pBeta0 = 0.3;
    double pSigSq = 1 - (pY + pBeta + pBeta0);
    int move;
    gsl_vector *moveProb = gsl_vector_calloc(4);
    gsl_vector_set(moveProb, 0, pY);
    gsl_vector_set(moveProb, 1, pBeta);
    gsl_vector_set(moveProb, 2, pBeta0);
    gsl_vector_set(moveProb, 3, pSigSq);
    
    
    BAFT_LNsurv_update_y(yL, yU, yU_posinf, c0, X, y, beta, beta0, sigSq);
    
    
    for(M = 0; M < *numReps; M++)
    {
        move = c_multinom_sample(rr, moveProb, 4);
        
        if(move == 1)
        {
            /* Updating y    */
            BAFT_LNsurv_update_y(yL, yU, yU_posinf, c0, X, y, beta, beta0, sigSq);
        }
        
        if(move == 2)
        {
            /* Updating beta          */
            BAFT_LNsurv_update_beta(yL, yU, yU_posinf, c0, c0_neginf, X, y, beta, beta0, sigSq, beta_prop_var, accept_beta);
        }
        
        if(move == 3)
        {
            /* Updating beta0 */
            BAFT_LNsurv_update_beta0(yL, yU, yU_posinf, c0, c0_neginf, X, y, beta, &beta0, sigSq, beta0_prop_var, &accept_beta0);
        }
        
        if(move == 4)
        {
            /* Updating sigSq         */
            BAFT_LNsurv_update_sigSq(yL, yU, yU_posinf, c0, c0_neginf, X, y, beta, beta0, &sigSq, a_sigSq, b_sigSq, sigSq_prop_var, &accept_sigSq);
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
            samples_beta0[StoreInx - 1] = beta0;
            samples_sigSq[StoreInx - 1] = sigSq;
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
        
        if( ( (M+1) % 10000 ) == 0)
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





















