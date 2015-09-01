/*
 TO COMPILE USE THE CODE:
 
 R CMD SHLIB BweibCorSurv.c BweibCorSurv_Updates.c BweibCorSurv_Utilities.c -lgsl -lgslcblas
 
 */

#include <stdio.h>
#include <math.h>
#include <time.h>

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_sort_vector.h"


#include "R.h"
#include "Rmath.h"

#include "BweibCorSurv.h"





/* */
void BweibCorSurvmcmc(double survData[],
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
                  double samples_zeta[],
                      double samples_misc[],
                      double moveVec[])
{
    GetRNGstate();
    
    time_t now;    
    
    int i, j, MM;

    
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
    double c       = hyperParams[2];
    double d       = hyperParams[3];
    double rho1    = hyperParams[4];
    double rho2    = hyperParams[5];
    
    

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
    
    double zeta = startValues[*p + 2 + *J];
    
    

    /* Variables required for storage of samples */
    
    int StoreInx;
    
    gsl_vector *accept_beta = gsl_vector_calloc(nP);
    gsl_vector *accept_V    = gsl_vector_calloc(*J);

    int accept_alpha = 0;

    
    /* Compute probabilities for various types of moves */
    
    double pRP, pSH, pSC, pCP, pVP, choice;
    int move, numUpdate;
    
    numUpdate = 4;
    if(*p > 0) numUpdate += 1;
    
    
    pSH = (double) 0.3;
    pCP = (double) 0.3;
    
    double probSub = (1 - pSH - pCP)/(numUpdate-2);
    
    
    pRP = (*p > 0) ? probSub : 0;
    pSC = probSub;
    pVP  = 1-(pRP + pSH + pSC + pCP);
    

   
    for(MM = 0; MM < *numReps; MM++)
    {
        /* selecting a move */
        /* move: 1=RP, 2=SH, 3=SC, 4=CP, 5=VP*/
        
        choice  = runif(0, 1);
        move    = 1;
        if(choice > pRP) move = 2;
        if(choice > pRP + pSH) move = 3;
        if(choice > pRP + pSH + pSC) move = 4;
        if(choice > pRP + pSH + pSC + pCP) move = 5;
        
        moveVec[MM] = (double) move;

        
        /* updating regression parameter: beta
        
        move = 100;*/
        
        if(move == 1)
        {
            BweibCorSurv_updateRP(beta, &alpha, &kappa, V, survTime, survEvent, cluster, survCov, accept_beta);
        }
        
        /* updating shape parameter: alpha
        
        move = 40;                */  
        
        if(move == 2)
        {
            BweibSurv_updateSH_rw2(beta, &alpha, &kappa, V, survTime, survEvent, cluster, survCov, mhProp_alpha_var, a, b, &accept_alpha);
        }

       

        /* updating scale parameter: kappa
  
        move = 70;     */
        
        if(move == 3)
        {
            BweibSurv_updateSC(beta, &alpha, &kappa, V, survTime, survEvent, cluster, survCov, c, d);
        }

  

        
        
        /* updating cluster-specific random effect: V 
         
         move = 130;*/
    
        if(move == 4)
        {
            BweibSurv_updateCP(beta, alpha, kappa, V, zeta, survTime, survEvent, cluster, survCov, n_j, accept_V, mhProp_V_var);
 
        }
        
        
        /* updating precition parameter of the random effect distribution: zeta
        
        move = 160;*/
        
        if(move == 5)
        {
            BweibSurv_updateVP(V, &zeta, rho1, rho2);
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
            
            samples_zeta[StoreInx - 1] = zeta;


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





















