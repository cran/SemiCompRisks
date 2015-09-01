/*
 TO COMPILE USE THE CODE:
 
 R CMD SHLIB BweibDpCorSurv.c BweibDpCorSurv_Updates.c BweibDpCorSurv_Utilities.c -lgsl -lgslcblas
 
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

#include "BweibDpCorSurv.h"






/* */
void BweibDpCorSurvmcmc(double survData[],
                  int *n,
                  int *p,
                  int *J,
                  double nj[],
                  double hyperParams[],
                  double mcmcParams[],
                  double startValues[],
                  int *numReps,
                  int *thin,
                  double *burninPerc,
                  double samples_beta[],
                  double samples_alpha[],
                  double samples_kappa[],
                  double samples_V[],
                  double samples_c[],
                        double samples_mu[],
                        double samples_zeta[],
                        double samples_tau[],
                        double samples_misc[],
                        double moveVec[])
{
    GetRNGstate();
    
    time_t now;       
    
    int i, j, MM;

    
    const gsl_rng_type * TT;
    gsl_rng * rr;
    
    gsl_rng_env_setup();
    
    TT = gsl_rng_default;
    rr = gsl_rng_alloc(TT);
    
    
    /* Survival Data */
    
    gsl_vector *survTime    = gsl_vector_alloc(*n);
    gsl_vector *survEvent   = gsl_vector_alloc(*n);
    gsl_vector *cluster      = gsl_vector_alloc(*n);
    for(i = 0; i < *n; i++)
    {
        gsl_vector_set(survTime, i, survData[(0 * *n) + i]);
        gsl_vector_set(survEvent, i, survData[(1* *n) + i]);
        gsl_vector_set(cluster, i, survData[(2* *n) + i]);
    }

    int nP;
    
    if(*p > 0) nP = *p;
    if(*p == 0) nP = 1;
    
    gsl_matrix *survCov     = gsl_matrix_calloc(*n, nP);
    
    if(*p >0)
    {
        for(i = 0; i < *n; i++)
        {
            for(j = 0; j < *(p); j++)
            {
                gsl_matrix_set(survCov, i, j, survData[((3+j)* *n) + i]);
            }
        }
    }
        
    gsl_vector *n_j = gsl_vector_calloc(*J);
    
    for(j = 0; j < *J; j++)
    {
        gsl_vector_set(n_j, j, nj[j]);
    }
    

    /* Hyperparameters */
    
    double a       = hyperParams[0];
    double b       = hyperParams[1];
    double c_kappa = hyperParams[2];
    double d       = hyperParams[3];
    double mu0     = hyperParams[4];
    double zeta0   = hyperParams[5];
    double a0       = hyperParams[6];
    double b0       = hyperParams[7];
    double aTau    = hyperParams[8];
    double bTau    = hyperParams[9];
    
    

    /* varialbes for M-H step */
    
    double mhProp_alpha_var = mcmcParams[0];
    double mhProp_V_var     = mcmcParams[1];
    
    
    
    /* Starting values */
    
    gsl_vector *beta = gsl_vector_calloc(nP);
    
    if(*p > 0)
    {
        for(j = 0; j < *p; j++) gsl_vector_set(beta, j, startValues[j]);
    }
    
    double alpha = startValues[*p];
    double kappa = startValues[*p + 1];
        
    gsl_vector *V = gsl_vector_calloc(*J);

    for(j = 0; j < *J; j++)
    {
        gsl_vector_set(V, j, startValues[*p + 2 + j]);
    }
    
    gsl_vector *c = gsl_vector_calloc(*J);
    
    for(i = 0; i < *J; i++)
    {
        gsl_vector_set(c, i, startValues[*p + 2 + *J + i]);
    }
    
    double tau = startValues[*p + 2 + *J + *J];
    
    
 

    /* Variables required for storage of samples */
    
    int StoreInx;
    
    gsl_vector *accept_beta = gsl_vector_calloc(nP);
    gsl_vector *accept_V    = gsl_vector_calloc(*J);

    int accept_alpha = 0;
    
    
    gsl_vector *mu_all = gsl_vector_calloc(*J);
    gsl_vector *zeta_all = gsl_vector_calloc(*J);
    
    int nClass_DP;

    
    /* Compute probabilities for various types of moves */
    
    double pRP, pSH, pSC, pCP, choice;
    int move, numUpdate;
    
    numUpdate = 3;
    if(*p > 0) numUpdate += 1;
    
    /*     */
    pCP = (double) 0.3;

    
    double probSub = (1 - pCP)/(numUpdate-1);
    
    pRP = (*p > 0) ? probSub : 0;
    pSC = probSub;
    pSH  = 1-(pRP + pSC + pCP);
    

   
    for(MM = 0; MM < *numReps; MM++)
    {
        /* selecting a move */
        /* move: 1=RP, 2=SH, 3=SC, 4=CP */
        
        choice  = runif(0, 1);
        move    = 1;
        if(choice > pRP) move = 2;
        if(choice > pRP + pSH) move = 3;
        if(choice > pRP + pSH + pSC) move = 4;
        
        moveVec[MM] = (double) move;

        
        /* updating regression parameter: beta */
        
        if(move == 1)
        {
            BweibDpCorSurv_updateRP(beta, &alpha, &kappa, V, survTime, survEvent, cluster, survCov, accept_beta);
        }
        
        /* updating shape parameter: alpha */
        
        if(move == 2)
        {
            BweibDpCorSurv_updateSH_rw2(beta, &alpha, &kappa, V, survTime, survEvent, cluster, survCov, mhProp_alpha_var, a, b, &accept_alpha);
        }

       

        /* updating scale parameter: kappa */
        
        if(move == 3)
        {
            BweibDpCorSurv_updateSC(beta, &alpha, &kappa, V, survTime, survEvent, cluster, survCov, c_kappa, d);
        }

  

        
        
        /* updating cluster-specific random effect: V */
    
        if(move == 4)
        {
            BweibDpCorSurv_updateCP(beta, alpha, kappa, V, survTime, survEvent, cluster, survCov, n_j, mu_all, zeta_all, c, accept_V, mhProp_V_var, mu0, zeta0, a0, b0, tau, &nClass_DP, rr);
            
            BweibDpCorSurv_updatePP(n, &tau, aTau, bTau, &nClass_DP);
 
        }
        
        

        /*        */
        
        
        /* Storing posterior samples */
        
        
        if( ( (MM+1) % *thin ) == 0 && (MM+1) > (*numReps * *burninPerc))
        {
            StoreInx = (MM+1)/(*thin)- (*numReps * *burninPerc)/(*thin);
 
            samples_alpha[StoreInx - 1] = alpha;
            samples_kappa[StoreInx - 1] = kappa;

            if(*p >0)
            {
                for(j = 0; j < *p; j++) samples_beta[(StoreInx - 1) * (*p) + j] = gsl_vector_get(beta, j);
            }

            for(j = 0; j < *J; j++) samples_V[(StoreInx - 1) * (*J) + j] = gsl_vector_get(V, j);
            
            for(j = 0; j < *J; j++) samples_c[(StoreInx - 1) * (*J) + j] = gsl_vector_get(c, j);
            
            for(j = 0; j < *J; j++) samples_mu[(StoreInx - 1) * (*J) + j] = gsl_vector_get(mu_all, j);
            
            for(j = 0; j < *J; j++) samples_zeta[(StoreInx - 1) * (*J) + j] = gsl_vector_get(zeta_all, j);
            
            samples_tau[StoreInx - 1] = tau;


            if(MM == (*numReps - 1))
            {
                            /*                             */
                if(*p >0)
                {
                    for(j = 0; j < *p; j++) samples_misc[j] = (int) gsl_vector_get(accept_beta, j);
                }

                /*                 */
                samples_misc[*p]    = accept_alpha;
                /*    */
                for(i = 0; i < *J; i++) samples_misc[*p + 1 + i] = (int) gsl_vector_get(accept_V, i);

            }


        }
        
        
        if( ( (MM+1) % 10000 ) == 0)
        {
            time(&now);
            
            Rprintf("iteration: %d: %s\n", MM+1, ctime(&now));
            
            
            R_FlushConsole();
            R_ProcessEvents();
            
            
        }
        
        
        
        
        
    }

    
    PutRNGstate();
    return;
    
    
}





















