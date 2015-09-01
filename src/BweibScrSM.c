/*
 TO COMPILE USE THE CODE:
 
 R CMD SHLIB BweibScr.c BweibScr_Updates.c BweibScr_Utilities.c -lgsl -lgslcblas
 
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

#include "BweibScrSM.h"





/* */
void BweibScrSMmcmc(double survData[],
                  int *n,
                  int *p1,
                  int *p2,
                  int *p3,
                  double hyperParams[],
                  double mcmcParams[],
                  double startValues[],
                  int *numReps,
                  int *thin,
                  double *burninPerc,
                  int *nGam_save,
                  double samples_beta1[],
                  double samples_beta2[],
                  double samples_beta3[],
                  double samples_alpha1[],
                  double samples_alpha2[],
                  double samples_alpha3[],
                  double samples_kappa1[],
                  double samples_kappa2[],
                  double samples_kappa3[],
                  double samples_theta[],
                  double samples_gamma[],
                    double samples_misc[],
                    double moveVec[])
{
    GetRNGstate();
    
    time_t now;
    
    int i, j, M;
    
    /* Survival Data */
    
    gsl_vector *survTime1    = gsl_vector_alloc(*n);
    gsl_vector *survTime2    = gsl_vector_alloc(*n);
    gsl_vector *survEvent1   = gsl_vector_alloc(*n);
    gsl_vector *survEvent2   = gsl_vector_alloc(*n);
    for(i = 0; i < *n; i++)
    {
        gsl_vector_set(survTime1, i, survData[(0 * *n) + i]);
        gsl_vector_set(survEvent1, i, survData[(1* *n) + i]);
        gsl_vector_set(survTime2, i, survData[(2 * *n) + i]);
        gsl_vector_set(survEvent2, i, survData[(3* *n) + i]);
    }

    int nP1, nP2, nP3;
    
    if(*p1 > 0) nP1 = *p1;
    if(*p1 == 0) nP1 = 1;
    if(*p2 > 0) nP2 = *p2;
    if(*p2 == 0) nP2 = 1;
    if(*p3 > 0) nP3 = *p3;
    if(*p3 == 0) nP3 = 1;
    
    gsl_matrix *survCov1     = gsl_matrix_calloc(*n, nP1);
    gsl_matrix *survCov2     = gsl_matrix_calloc(*n, nP2);
    gsl_matrix *survCov3     = gsl_matrix_calloc(*n, nP3);
    
    if(*p1 >0)
    {
        for(i = 0; i < *n; i++)
        {
            for(j = 0; j < *(p1); j++)
            {
                gsl_matrix_set(survCov1, i, j, survData[((4+j)* *n) + i]);
            }
        }
    }
    
    if(*p2 >0)
    {
        for(i = 0; i < *n; i++)
        {
            for(j = 0; j < *(p2); j++)
            {
                gsl_matrix_set(survCov2, i, j, survData[((4+(*p1)+j)* *n) + i]);
            }
        }
    }
    
    if(*p3 >0)
    {
        for(i = 0; i < *n; i++)
        {
            for(j = 0; j < *(p3); j++)
            {
                gsl_matrix_set(survCov3, i, j, survData[((4+(*p1)+(*p2)+j)* *n) + i]);
            }
        }
    }
    
    
    gsl_vector *case01   = gsl_vector_alloc(*n);
    gsl_vector *case11   = gsl_vector_alloc(*n);
    
    gsl_vector_memcpy(case01, survEvent1);
    gsl_vector_scale(case01, -1);
    gsl_vector_add_constant(case01, 1);
    gsl_vector_mul(case01, survEvent2);
    
    gsl_vector_memcpy(case11, survEvent1);
    gsl_vector_mul(case11, survEvent2);
    
    gsl_vector *yStar = gsl_vector_calloc(*n);
    gsl_vector_memcpy(yStar, survTime2);
    gsl_vector_sub(yStar, survTime1);
    
    
    
    /* Hyperparameters */
    
    double a1       = hyperParams[0];
    double b1       = hyperParams[1];
    double a2       = hyperParams[2];
    double b2       = hyperParams[3];
    double a3       = hyperParams[4];
    double b3       = hyperParams[5];
    double c1       = hyperParams[6];
    double d1       = hyperParams[7];
    double c2       = hyperParams[8];
    double d2       = hyperParams[9];
    double c3       = hyperParams[10];
    double d3       = hyperParams[11];
    double psi      = hyperParams[12];
    double omega    = hyperParams[13];
    

    /* varialbes for M-H step */
    
    double mhProp_alpha1_var = mcmcParams[0];
    double mhProp_alpha2_var = mcmcParams[1];
    double mhProp_alpha3_var = mcmcParams[2];
    double mhProp_theta_var  = mcmcParams[3];
    
    
    
    /* Starting values */
    
    gsl_vector *beta1 = gsl_vector_calloc(nP1);
    gsl_vector *beta2 = gsl_vector_calloc(nP2);
    gsl_vector *beta3 = gsl_vector_calloc(nP3);
    
    if(*p1 > 0)
    {
        for(j = 0; j < *p1; j++) gsl_vector_set(beta1, j, startValues[j]);
    }
    if(*p2 > 0)
    {
        for(j = 0; j < *p2; j++) gsl_vector_set(beta2, j, startValues[j + *p1]);
    }
    if(*p3 > 0)
    {
        for(j = 0; j < *p3; j++) gsl_vector_set(beta3, j, startValues[j + *p1 + *p2]);
    }
    
    double alpha1 = startValues[*p1 + *p2 + *p3];
    double alpha2 = startValues[*p1 + *p2 + *p3 + 1];
    double alpha3 = startValues[*p1 + *p2 + *p3 + 2];
    double kappa1 = startValues[*p1 + *p2 + *p3 + 3];
    double kappa2 = startValues[*p1 + *p2 + *p3 + 4];
    double kappa3 = startValues[*p1 + *p2 + *p3 + 5];
    double theta  = startValues[*p1 + *p2 + *p3 + 6];
    
    gsl_vector *gamma = gsl_vector_calloc(*n);
    
    for(i = 0; i < *n; i++)
    {
        gsl_vector_set(gamma, i, startValues[*p1 + *p2 + *p3 + 7 + i]);
    }
    

    /* Variables required for storage of samples */
    
    int StoreInx;
    
    gsl_vector *accept_beta1 = gsl_vector_calloc(nP1);
    gsl_vector *accept_beta2 = gsl_vector_calloc(nP2);
    gsl_vector *accept_beta3 = gsl_vector_calloc(nP3);
    
    int accept_alpha1 = 0;
    int accept_alpha2 = 0;
    int accept_alpha3 = 0;
    int accept_theta  = 0;
    

    
    /* Compute probabilities for various types of moves */
    
    double pRP1, pRP2, pRP3, pSH1, pSH2, pSH3, pSC1, pSC2, pSC3, pDP, pFP, choice;
    int move, numUpdate;
    
    numUpdate = 8;
    if(*p1 > 0) numUpdate += 1;
    if(*p2 > 0) numUpdate += 1;
    if(*p3 > 0) numUpdate += 1;
    
    pRP1 = (*p1 > 0) ? (double) 1/numUpdate : 0;
    pRP2 = (*p2 > 0) ? (double) 1/numUpdate : 0;
    pRP3 = (*p3 > 0) ? (double) 1/numUpdate : 0;
    pSH1 = (double) 1/numUpdate;
    pSH2 = (double) 1/numUpdate;
    pSH3 = (double) 1/numUpdate;
    pSC1 = (double) 1/numUpdate;
    pSC2 = (double) 1/numUpdate;
    pSC3 = (double) 1/numUpdate;
    pDP  = (double) 1/numUpdate;
    pFP  = 1-(pRP1 + pRP2 + pRP3 + pSH1 + pSH2 + pSH3 + pSC1 + pSC2 + pSC3 + pDP);
    
    
    
    
    
    
    for(M = 0; M < *numReps; M++)
    {
        /* selecting a move */
        /* move: 1=RP1, 2=RP2, 3=RP3, 4=SH1, 5=SH2, 6=SH3 */
        /* move: 7=SC1, 8=SC2, 9=SC3, 10=DP, 11=FP */
        
        choice  = runif(0, 1);
        move    = 1;
        if(choice > pRP1) move = 2;
        if(choice > pRP1 + pRP2) move = 3;
        if(choice > pRP1 + pRP2 + pRP3) move = 4;
        if(choice > pRP1 + pRP2 + pRP3 + pSH1) move = 5;
        if(choice > pRP1 + pRP2 + pRP3 + pSH1 + pSH2) move = 6;
        if(choice > pRP1 + pRP2 + pRP3 + pSH1 + pSH2 + pSH3) move = 7;
        if(choice > pRP1 + pRP2 + pRP3 + pSH1 + pSH2 + pSH3 + pSC1) move = 8;
        if(choice > pRP1 + pRP2 + pRP3 + pSH1 + pSH2 + pSH3 + pSC1 + pSC2) move = 9;
        if(choice > pRP1 + pRP2 + pRP3 + pSH1 + pSH2 + pSH3 + pSC1 + pSC2 + pSC3) move = 10;
        if(choice > pRP1 + pRP2 + pRP3 + pSH1 + pSH2 + pSH3 + pSC1 + pSC2 + pSC3 + pDP) move = 11;
        
        
        moveVec[M] = (double) move;        
        
        
        /* updating regression parameter: beta1 */
        
        if(move == 1)
        {
            BweibScrSM_updateRP1(beta1, &alpha1, &kappa1, gamma, survTime1, survEvent1, survCov1, accept_beta1);
        }
        
        
        /* updating regression parameter: beta2 */
        
        if(move == 2)
        {
            BweibScrSM_updateRP2(beta2, &alpha2, &kappa2, gamma, survTime1, case01, survCov2, accept_beta2);
        }
        
        
        /* updating regression parameter: beta3 */
        
        if(move == 3)
        {
            BweibScrSM_updateRP3(beta3, &alpha3, &kappa3, gamma, yStar, case11, survCov3, accept_beta3);
        }
        
        
        /* updating shape parameter: alpha1 */
        
        if(move == 4)
        {
            BweibScrSM_updateSC1(beta1, &alpha1, &kappa1, gamma, survTime1, survEvent1, survCov1, mhProp_alpha1_var, a1, b1, &accept_alpha1);
        }
        
        
        /* updating shape parameter: alpha2 */
        
        if(move == 5)
        {
            BweibScrSM_updateSC2(beta2, &alpha2, &kappa2, gamma, survTime1, survTime2, case01, survCov2, mhProp_alpha2_var, a2, b2, &accept_alpha2);
        }
        

        /* updating shape parameter: alpha3 */
        
        if(move == 6)
        {
            BweibScrSM_updateSC3(beta3, &alpha3, &kappa3, gamma, yStar, case11, survCov3, mhProp_alpha3_var, a3, b3, &accept_alpha3);
        }


        /* updating shape parameter: kappa1 */
        
        if(move == 7)
        {
            BweibScrSM_updateSH1(beta1, &alpha1, &kappa1, gamma, survTime1, survEvent1, survCov1, c1, d1);
        }
        
        
        /* updating shape parameter: kappa2 */
        
        if(move == 8)
        {
            BweibScrSM_updateSH2(beta2, &alpha2, &kappa2, gamma, survTime1, case01, survCov2, c2, d2);
        }
        
        
        /* updating shape parameter: kappa3 */
        
        if(move == 9)
        {
            BweibScrSM_updateSH3(beta3, &alpha3, &kappa3, gamma, yStar, case11, survCov3, c3, d3);
        }
        
        
        /* updating variance parameter: theta */
        
        if(move == 10)
        {
            BweibScrSM_updateDP(gamma, &theta, mhProp_theta_var, psi, omega, &accept_theta);
        }
        
        
        /* updating frailty parameter: gamma */
        
        if(move == 11)
        {
            BweibScrSM_updateFP(beta1, beta2, beta3, alpha1, alpha2, alpha3, kappa1, kappa2, kappa3, theta, gamma, survTime1, yStar, survEvent1, survEvent2, survCov1, survCov2, survCov3);
        }
        
        
        
        
        /* Storing posterior samples */
        
        
        if( ( (M+1) % *thin ) == 0 && (M+1) > (*numReps * *burninPerc))
        {
            StoreInx = (M+1)/(*thin)- (*numReps * *burninPerc)/(*thin);
 
            samples_alpha1[StoreInx - 1] = alpha1;
            samples_alpha2[StoreInx - 1] = alpha2;
            samples_alpha3[StoreInx - 1] = alpha3;
            
            samples_kappa1[StoreInx - 1] = kappa1;
            samples_kappa2[StoreInx - 1] = kappa2;
            samples_kappa3[StoreInx - 1] = kappa3;
            
            samples_theta[StoreInx - 1] = theta;
 
            if(*p1 >0)
            {
                for(j = 0; j < *p1; j++) samples_beta1[(StoreInx - 1) * (*p1) + j] = gsl_vector_get(beta1, j);
            }
            if(*p2 >0)
            {
                for(j = 0; j < *p2; j++) samples_beta2[(StoreInx - 1) * (*p2) + j] = gsl_vector_get(beta2, j);
            }
            if(*p3 >0)
            {
                for(j = 0; j < *p3; j++) samples_beta3[(StoreInx - 1) * (*p3) + j] = gsl_vector_get(beta3, j);
            }
            
            if(*nGam_save == *n)
            {
                for(i = 0; i < *n; i++)
                {
                    samples_gamma[(StoreInx - 1) * (*n) + i] = gsl_vector_get(gamma, i);
                }
            }
            if(*nGam_save < *n)
            {
                for(i = 0; i < *nGam_save; i++)
                {
                    samples_gamma[(StoreInx - 1) * (*nGam_save) + i] = gsl_vector_get(gamma, i);
                }
            }
            
            if(M == (*numReps - 1))
            {
                if(*p1 >0)
                {
                    for(j = 0; j < *p1; j++) samples_misc[j] = (int) gsl_vector_get(accept_beta1, j);
                }
                if(*p2 >0)
                {
                    for(j = 0; j < *p2; j++) samples_misc[*p1 + j] = (int) gsl_vector_get(accept_beta2, j);
                }
                if(*p3 >0)
                {
                for(j = 0; j < *p3; j++) samples_misc[*p1 + *p2 + j] = (int) gsl_vector_get(accept_beta3, j);
                }
                samples_misc[*p1 + *p2 + *p3]       = accept_alpha1;
                samples_misc[*p1 + *p2 + *p3 + 1]   = accept_alpha2;
                samples_misc[*p1 + *p2 + *p3 + 2]   = accept_alpha3;
                samples_misc[*p1 + *p2 + *p3 + 3]   = accept_theta;
            }
            
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





















