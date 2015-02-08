/*
 TO COMPILE USE THE CODE:

 R CMD SHLIB BpeSurv.c BpeSurv_Updates.c BpeSurv_Utilities.c -lgsl -lgslcblas
 
 */

#include <stdio.h>
#include <math.h>
#include <time.h>

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_sort_vector.h"
#include "gsl/gsl_heapsort.h"

#include "R.h"
#include "Rmath.h"

#include "BpeSurv.h"


/* */
void BpeSurvmcmc(double survData[],
                    int *n,
                    int *p,
                    double hyperParams[],
                    double startValues[],
                    double mcmcParams[],
                    int *numReps,
                    int *thin,
                    double *burninPerc,
                    double samples_beta[],
                    double samples_mu_lam[],
                    double samples_sigSq_lam[],
                    double samples_J[],
                    double samples_s[],
                    double samples_misc[],
                    double lambda_fin[])
{
    GetRNGstate();
    
    time_t now;    
    
    int i, j, M;
    
    /* Survival Data */
    
    gsl_vector *survTime    = gsl_vector_alloc(*n);
    gsl_vector *survEvent   = gsl_vector_alloc(*n);
    for(i = 0; i < *n; i++)
    {
        gsl_vector_set(survTime, i, survData[(0 * *n) + i]);
        gsl_vector_set(survEvent, i, survData[(1* *n) + i]);
    }
    
    int nP;
    
    if(*p > 0) nP = *p;
    if(*p == 0) nP = 1;
    
    gsl_matrix *survCov = gsl_matrix_calloc(*n, nP);
    
    
    if(*p > 0)
    {
        for(i = 0; i < *n; i++)
        {
            for(j = 0; j < *p; j++) gsl_matrix_set(survCov, i, j, survData[((2+j)* *n) + i]);
        }
    }
    
    /* Hyperparameters */
    
    double a        = hyperParams[0];
    double b        = hyperParams[1];
    double alpha    = hyperParams[2];
    double c_lam    = hyperParams[3];
    
    
    
    
    /* varialbes for birth and death moves */
    
    double C                = mcmcParams[0];
    double delPert          = mcmcParams[1];
    int num_s_propBI        = mcmcParams[2];
    int J_max               = mcmcParams[3];
    double s_max            = mcmcParams[4];    
    
    gsl_vector *s_propBI    = gsl_vector_calloc(num_s_propBI);
    for(j = 0; j < num_s_propBI; j++) gsl_vector_set(s_propBI, j, mcmcParams[6+j]);
   
    
    /* time points where lambda values are stored  */

    int nTime_lambda = mcmcParams[5];
    
    gsl_vector *time_lambda = gsl_vector_calloc(nTime_lambda);
    
    for(i = 0; i < nTime_lambda; i++)
    {
        gsl_vector_set(time_lambda, i, mcmcParams[num_s_propBI + 6 + i]);
    }
                       
  
    /* Starting values */
    
    gsl_vector *beta = gsl_vector_calloc(nP);
    if(*p > 0)
    {
        for(j = 0; j < *p; j++) gsl_vector_set(beta, j, startValues[j]);
    }
    
    int J               = startValues[*p];
    double mu_lam       = startValues[*p+1];
    double sigSq_lam    = startValues[*p+2];
    
    gsl_vector *lambda  = gsl_vector_calloc(J_max+1);
    for(j = 0; j < (J+1); j++) gsl_vector_set(lambda, j, startValues[*p+3+j]);
        
    gsl_vector *s       = gsl_vector_calloc(J_max+1);
    for(j = 0; j < (J+1); j++) gsl_vector_set(s, j, startValues[*p+J+4+j]);
    
    gsl_vector *xbeta = gsl_vector_calloc(*n);
    gsl_blas_dgemv(CblasNoTrans, 1, survCov, beta, 0, xbeta);


    
    /* Calculating Sigma_lam (from W and Q) */

    
    gsl_matrix *Sigma_lam       = gsl_matrix_calloc(J_max+1, J_max+1);
    gsl_matrix *invSigma_lam    = gsl_matrix_calloc(J_max+1, J_max+1);
    gsl_matrix *W               = gsl_matrix_calloc(J_max+1, J_max+1);
    gsl_matrix *Q               = gsl_matrix_calloc(J_max+1, J_max+1);
 

    cal_Sigma(Sigma_lam, invSigma_lam, W, Q, s, c_lam, J);
    

    
    /* Variables required for storage of samples */
    
    int StoreInx;
    
    gsl_vector *accept_beta = gsl_vector_calloc(nP);
    
    int accept_BI = 0;
    int accept_DI = 0;
    
    
    /* Compute probabilities for various types of moves */
    
    double pRP, pBH, pSP, pBI, pDI, choice;
    int move, numUpdate;
    
    gsl_vector *pB          = gsl_vector_calloc(J_max);
    gsl_vector *pD          = gsl_vector_calloc(J_max);
    gsl_vector *rho_lam_vec = gsl_vector_calloc(J_max);
    double  rho_lam;
    
    for(j = 0; j < J_max; j++)
    {
        gsl_vector_set(pB, j, c_min(1, alpha/(j+1 + 1)));
        gsl_vector_set(pD, j, c_min(1, (j+1)/alpha));
    }
    
    for(j = 0; j < J_max; j++) gsl_vector_set(rho_lam_vec, j, C/(gsl_vector_get(pB, j) + gsl_vector_get(pD, j)));
    
    rho_lam = gsl_vector_min(rho_lam_vec);
    
    numUpdate = 4;
    if(*p > 0) numUpdate += 1;
    
    
    for(M = 0; M < *numReps; M++)
    {
        
        if(J < J_max)
        {
            pBI = rho_lam * c_min(1, alpha/(J+1));
            pDI = rho_lam * c_min(1, J/alpha);
        }
        if(J >= J_max)
        {
            pBI = 0;
            pDI = rho_lam * 2;
        }
        
        pRP = (*p > 0) ? (double) (1-pBI-pDI)/(numUpdate-2) : 0;
        pBH = (double) (1-pBI-pDI)/(numUpdate-2);
        pSP = 1-(pBI+pDI+pRP+pBH);
        
    
        /* selecting a move */
        /* move: 1=RP, 2=BH, 3=SP, 4=BI, 5=DI */
        
        choice  = runif(0, 1);
        move    = 1;
        if(choice > pRP) move = 2;
        if(choice > pRP + pBH) move = 3;
        if(choice > pRP + pBH + pSP) move = 4;
        if(choice > pRP + pBH + pSP + pBI) move = 5; 
        
        
        /* updating regression parameter: beta */
        
        

        
        if(move == 1)
        {
            BpeSur_updateRP2(beta, xbeta, accept_beta, lambda, s, survTime, survEvent, survCov, J);

        }
        
        
        
        

        /* updating log-baseline hazard function parameter: lambda */

    
        
        if(move == 2)
        {
            BpeSur_updateBH2(lambda, s, xbeta, survTime, survEvent, Sigma_lam, invSigma_lam, W, Q, mu_lam, sigSq_lam, J);
        }
    

		/* Updating second stage survival components: mu_lam and sigSq_lam */
        
      

        if(move == 3)
        {
            BpeSur_updateSP(&mu_lam, &sigSq_lam, lambda, Sigma_lam, invSigma_lam, a, b, J);
        }

        
        /* Updating the number of splits and their positions: J and s (Birth move) */        
        


        if(move == 4)
        {
            BpeSur_updateBI2(s, &J, &accept_BI, survTime, survEvent, xbeta, Sigma_lam, invSigma_lam, W, Q, lambda, s_propBI, num_s_propBI, delPert, alpha, c_lam, mu_lam, sigSq_lam, s_max);
        }
        

        
        /* Updating the number of splits and their positions: J and s (Death move) */           
        


        if(move == 5)
        {
            BpeSur_updateDI2(s, &J, &accept_DI, survTime, survEvent, xbeta, Sigma_lam, invSigma_lam, W, Q, lambda, s_propBI, num_s_propBI, delPert, alpha, c_lam, mu_lam, sigSq_lam, s_max, J_max);
        }
        
        
        
        
        
        /* Storing posterior samples */
        

        if( ( (M+1) % *thin ) == 0 && (M+1) > (*numReps * *burninPerc))
        {
            StoreInx = (M+1)/(*thin)- (*numReps * *burninPerc)/(*thin);
            
            if(*p >0)
            {
                for(j = 0; j < *p; j++)
                {
                    samples_beta[(StoreInx - 1) * (*p) + j] = gsl_vector_get(beta, j);
                }
            }

            j = 0;
            
            for(i = 0; i < nTime_lambda; i++)
            {
                j = 0;
                while(gsl_vector_get(time_lambda, i) > gsl_vector_get(s, j))
                {
                    j += 1;
                }
                lambda_fin[(StoreInx - 1) * (nTime_lambda) + i] = gsl_vector_get(lambda, j);

            }

            
            samples_mu_lam[StoreInx - 1] = mu_lam;
            samples_sigSq_lam[StoreInx - 1] = sigSq_lam;
            samples_J[StoreInx - 1] = J;
            for(j = 0; j < J+1; j++)
            {
                samples_s[(StoreInx - 1) * (J_max+1) + j] = gsl_vector_get(s, j);
            }
  
        }

        
        if(M == (*numReps - 1))
        {
            if(*p > 0)
            {
                for(j = 0; j < *p; j++)
                {
                    samples_misc[j] = (int) gsl_vector_get(accept_beta, j);
                }
            }
            samples_misc[*p] = accept_BI;
            samples_misc[*p+1] = accept_DI;
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






