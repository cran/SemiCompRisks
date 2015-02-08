/*
 TO COMPILE USE THE CODE:
 
 R CMD SHLIB BpeMvnCorSurv.c BpeMvnCorSurv_Updates.c BpeMvnCorSurv_Utilities.c -lgsl -lgslcblas
 
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

#include "BpeMvnCorSurv.h"







/* */
void BpeMvnCorSurvmcmc(double survData[],
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
                       double samples_mu_lam[],
                       double samples_sigSq_lam[],
                       double samples_K[],
                       double samples_s[],
                       double samples_V[],
                       double samples_zeta[],
                       double samples_misc[],
                       double lambda_fin[],
                       double dev[])
{
    GetRNGstate();
    
    time_t now;        
    
    int i, j, MM, lastChgProp;
    int ChgProp = 0;
    
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
    double alpha   = hyperParams[2];
    double c_lam   = hyperParams[3];
    double rho1    = hyperParams[4];
    double rho2    = hyperParams[5];
    
    
    /* varialbes for birth and death moves */
    
    double C                = mcmcParams[0];
    double delPert          = mcmcParams[1];
    int num_s_propBI        = mcmcParams[2];
    int K_max               = mcmcParams[3];
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
    
    double mhProp_V_var     = mcmcParams[num_s_propBI + 6 + nTime_lambda];    

    
    
    /* Starting values */
    
    gsl_vector *beta = gsl_vector_calloc(nP);
    
    if(*p > 0)
    {
        for(j = 0; j < *p; j++) gsl_vector_set(beta, j, startValues[j]);
    }
    
    int K               = startValues[*p];
    double mu_lam       = startValues[*p+1];
    double sigSq_lam    = startValues[*p+2];
    
    gsl_vector *lambda  = gsl_vector_calloc(K_max+1);
    for(j = 0; j < (K+1); j++) gsl_vector_set(lambda, j, startValues[*p+3+j]);
    
    gsl_vector *s       = gsl_vector_calloc(K_max+1);
    for(j = 0; j < (K+1); j++) gsl_vector_set(s, j, startValues[*p+K+4+j]);
        
    gsl_vector *V = gsl_vector_calloc(*J);

    for(j = 0; j < *J; j++)
    {
        gsl_vector_set(V, j, startValues[*p+K+4+K+1+j]);
    }
    
    double zeta = startValues[*p+K+4+K+1+*J];
    
    
    
    gsl_vector *xbeta = gsl_vector_calloc(*n);
    gsl_blas_dgemv(CblasNoTrans, 1, survCov, beta, 0, xbeta);
    

    

    /* Calculating Sigma_lam (from W and Q) */
    
    
    gsl_matrix *Sigma_lam       = gsl_matrix_calloc(K_max+1, K_max+1);
    gsl_matrix *invSigma_lam    = gsl_matrix_calloc(K_max+1, K_max+1);
    gsl_matrix *W               = gsl_matrix_calloc(K_max+1, K_max+1);
    gsl_matrix *Q               = gsl_matrix_calloc(K_max+1, K_max+1);
    
    
    cal_Sigma(Sigma_lam, invSigma_lam, W, Q, s, c_lam, K);
    
    
    
    /* Variables required for storage of samples */
    
    int StoreInx;
    
    gsl_vector *accept_beta = gsl_vector_calloc(nP);
    gsl_vector *accept_V    = gsl_vector_calloc(*J);

    int accept_BI = 0;
    int accept_DI = 0;
    
    double logLH;        

    
    /* Compute probabilities for various types of moves */
    
    double pRP, pBH, pSP, pBI, pDI, pCP, pVP, choice;
    int move, numUpdate;
    
    gsl_vector *pB          = gsl_vector_calloc(K_max);
    gsl_vector *pD          = gsl_vector_calloc(K_max);
    gsl_vector *rho_lam_vec = gsl_vector_calloc(K_max);
    double  rho_lam;
    
    for(j = 0; j < K_max; j++)
    {
        gsl_vector_set(pB, j, c_min(1, alpha/(j+1 + 1)));
        gsl_vector_set(pD, j, c_min(1, (j+1)/alpha));
    }
    
    for(j = 0; j < K_max; j++) gsl_vector_set(rho_lam_vec, j, C/(gsl_vector_get(pB, j) + gsl_vector_get(pD, j)));
    
    rho_lam = gsl_vector_min(rho_lam_vec);
    
    numUpdate = 6;
    if(*p > 0) numUpdate += 1;
    
    
   
    
    for(MM = 0; MM < *numReps; MM++)
    {
        
        if(K < K_max)
        {
            pBI = rho_lam * c_min(1, alpha/(K+1));
            pDI = rho_lam * c_min(1, K/alpha);
        }
        if(K >= K_max)
        {
            pBI = 0;
            pDI = rho_lam * 2;
        }
        

        
        pRP = (*p > 0) ? (double) (1-pBI-pDI)/(numUpdate-2) : 0;
        pBH = (double) (1-pBI-pDI)/(numUpdate-2);
        pSP = (double) (1-pBI-pDI)/(numUpdate-2);
        pCP = (double) (1-pBI-pDI)/(numUpdate-2);
        pVP = 1-(pBI+pDI+pRP+pBH+pSP+pCP);
        
        
        /* selecting a move */
        /* move: 1=RP, 2=BH, 3=SP, 4=BI, 5=DI, 6 = CP, 7=VP */
        
        choice  = runif(0, 1);
        move    = 1;
        if(choice > pRP) move = 2;
        if(choice > pRP + pBH) move = 3;
        if(choice > pRP + pBH + pSP) move = 4;
        if(choice > pRP + pBH + pSP + pBI) move = 5;
        if(choice > pRP + pBH + pSP + pBI + pDI) move = 6;
        if(choice > pRP + pBH + pSP + pBI + pDI+ pCP) move = 7;
        
                
        /* updating regression parameter: beta
        
        move = 100;*/
        
        if(move == 1)
        {
            BpeMvnCorSurv_updateRP(beta, xbeta, lambda, s, K, V, survTime, survEvent, cluster, survCov, accept_beta);
        }

        
        
        /* updating log-baseline hazard function parameter: lambda   
        
        move = 200;*/
        
        if(move == 2)
        {
            BpeMvnCorSurv_updateBH(lambda, s, xbeta, V, survTime, survEvent, cluster, Sigma_lam, invSigma_lam, W, Q, mu_lam, sigSq_lam, K);
                 
        }

        
        
        /* updating second stage survival components: mu_lam and sigSq_lam
        
         move = 300;*/
        
        
        if(move == 3)
        {
            BpeMvnCorSurv_updateSP(&mu_lam, &sigSq_lam, lambda, Sigma_lam, invSigma_lam, a, b, K);
        }

        
        
        /* Updating the number of splits and their positions: K and s (Birth move) 
         
         move = 400;*/
        
        if(move == 4)
        {
            BpeMvnCorSurv_updateBI(s, &K, &accept_BI, survTime, survEvent, xbeta, V, cluster, Sigma_lam, invSigma_lam, W, Q, lambda, s_propBI, num_s_propBI, delPert, alpha, c_lam, mu_lam, sigSq_lam, s_max);
        }
        
        
        /* Updating the number of splits and their positions: K and s (Death move)  
         
         move = 500;*/
        
        if(move == 5)
        {
            BpeMvnCorSurv_updateDI(s, &K, &accept_DI, survTime, survEvent, xbeta, V, cluster, Sigma_lam, invSigma_lam, W, Q, lambda, s_propBI, num_s_propBI, delPert, alpha, c_lam, mu_lam, sigSq_lam, s_max, K_max);
        }
        
        

        
        /* updating cluster-specific random effect: V
         
         move = 600;*/
        
        if(move == 6)
        {
            BpeMvnCorSurv_updateCP(beta, V, lambda, s, K, zeta, survTime, survEvent, cluster, survCov, n_j, accept_V, mhProp_V_var);
            
        }
        
        /* updating precition parameter of the random effect distribution: zeta 
        
        move = 700;*/
        
        if(move == 7)
        {
            BpeMvnCorSurv_updateVP(V, &zeta, rho1, rho2);
        }
         
        

        /*        */
        
        
        /* Storing posterior samples         */
        
        
        if( ( (MM+1) % *thin ) == 0 && (MM+1) > (*numReps * *burninPerc))
        {
            StoreInx = (MM+1)/(*thin)- (*numReps * *burninPerc)/(*thin);

            if(*p >0)
            {
                for(j = 0; j < *p; j++) samples_beta[(StoreInx - 1) * (*p) + j] = gsl_vector_get(beta, j);
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
            samples_K[StoreInx - 1] = K;
            for(j = 0; j < K+1; j++)
            {
                samples_s[(StoreInx - 1) * (K_max+1) + j] = gsl_vector_get(s, j);
            }

            for(j = 0; j < *J; j++) samples_V[(StoreInx - 1) * (*J) + j] = gsl_vector_get(V, j);
            
            samples_zeta[StoreInx - 1] = zeta;
            
            BpeMvnCorSurv_logLH(beta, xbeta, lambda, s, V, survTime, survEvent, survCov, cluster, K, &logLH);
            
            dev[StoreInx - 1] = -2*logLH;     
            
        }


            if(MM == (*numReps - 1))
            {
                if(*p >0)
                {
                    for(j = 0; j < *p; j++) samples_misc[j] = (int) gsl_vector_get(accept_beta, j);
                }

                samples_misc[*p] = accept_BI;
                samples_misc[*p+1] = accept_DI;

                for(i = 0; i < *J; i++) samples_misc[*p+2 + i] = (int) gsl_vector_get(accept_V, i);

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





















