/*
 TO COMPILE USE THE CODE:
 
 R CMD SHLIB BweibSurv.c BweibSurv_Updates.c BweibSurv_Utilities.c -lgsl -lgslcblas
 
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

#include "BweibSurv.h"


/* */
void BweibSurvmcmc(double survData[],
                    int *n,
                    int *p,
                    double hyperParams[],
                    double *mcmcParams,
                    double startValues[],
                    int *numReps,
                    int *thin,
                    double *burninPerc,
                    double samples_beta[],
                    double samples_alpha[],
                    double samples_kappa[],
                    double samples_misc[])
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
    
    double a = hyperParams[0];
    double b = hyperParams[1];
    double c = hyperParams[2];
    double d = hyperParams[3];

    /* varialbes for M-H step */
    
    double mhProp_alpha_var = *mcmcParams;
    
    
    /* Starting values */
    
    gsl_vector *beta = gsl_vector_calloc(nP);
    if(*p > 0)
    {
        for(j = 0; j < *p; j++) gsl_vector_set(beta, j, startValues[j]);
    }

    double alpha = startValues[*p];
    double kappa = startValues[*p+1];
    
    
    /* Variables required for storage of samples */
    
    int StoreInx;
    
    gsl_vector *accept_beta = gsl_vector_calloc(nP);
    
    int accept_alpha = 0;
    

    for(M = 0; M < *numReps; M++)
    {
        /* updating regression parameter: beta */
        if(*p > 0)
        {
            BweibSurv_updateRP(beta, &alpha, &kappa, survTime, survEvent, survCov, accept_beta);
        }

        /* updating scale parameter: alpha */
        
        BweibSurv_updateSC1(beta, &alpha, &kappa, survTime, survEvent, survCov, mhProp_alpha_var, a, b, &accept_alpha);
        
        
        /* updating shape parameter: kappa */
        
        BweibSurv_updateSH(beta, &alpha, &kappa, survTime, survEvent, survCov, c, d);
            
        
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
            
            samples_alpha[StoreInx - 1] = alpha;
            samples_kappa[StoreInx - 1] = kappa;
        }
        
        if(M == (*numReps - 1)){
            if(*p >0)
            {
                for(j = 0; j < *p; j++)
                {
                    samples_misc[j] = (int) gsl_vector_get(accept_beta, j);
                }
            }
            samples_misc[*p] = accept_alpha;
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






