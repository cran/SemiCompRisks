/*
 TO COMPILE USE THE CODE:
 
 R CMD SHLIB BAFTscr.c BAFTscr_Updates.c BAFTscr_Utilities.c -lgsl -lgslcblas
 
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
#include "BAFT_LNscr.h"


/* */
void BAFTscrmcmc(double Ymat[],
                 double y1LInf[],
                 double y1UInf[],
                 double y2LInf[],
                 double y2UInf[],
                 double c0Inf[],
                 double X1mat[],
                 double X2mat[],
                 double X3mat[],
                 int *n,
                 int *p1,
                 int *p2,
                 int *p3,
                 double hyperP[],
                 double mcmcP[],
                 double startValues[],
                 int *numReps,
                 int *thin,
                 double *burninPerc,
                 int *nGam_save,
                 double samples_y1[],
                 double samples_y2[],
                 double samples_y1_NA[],
                 double samples_beta1[],
                 double samples_beta2[],
                 double samples_beta3[],
                 double samples_beta01[],
                 double samples_beta02[],
                 double samples_beta03[],
                 double samples_sigSq1[],
                 double samples_sigSq2[],
                 double samples_sigSq3[],
                 double samples_theta[],
                 double samples_gamma[],
                 double samples_misc[],
                 double samples_logLH[],
                 double LH_i_mean[],
                 double invLH_i_mean[])
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
    
    gsl_vector *y1L = gsl_vector_calloc(*n);
    gsl_vector *y1U = gsl_vector_calloc(*n);
    gsl_vector *y2L = gsl_vector_calloc(*n);
    gsl_vector *y2U = gsl_vector_calloc(*n);
    gsl_vector *y1L_neginf = gsl_vector_calloc(*n);
    gsl_vector *y1U_posinf = gsl_vector_calloc(*n);
    gsl_vector *y2L_neginf = gsl_vector_calloc(*n);
    gsl_vector *y2U_posinf = gsl_vector_calloc(*n);
    gsl_vector *c0_neginf = gsl_vector_calloc(*n);
    gsl_vector *y1_NA = gsl_vector_calloc(*n);
    gsl_vector *c0 = gsl_vector_calloc(*n);
    
    int nP1, nP2, nP3;
    
    if(*p1 > 0) nP1 = *p1;
    if(*p1 == 0) nP1 = 1;
    if(*p2 > 0) nP2 = *p2;
    if(*p2 == 0) nP2 = 1;
    if(*p3 > 0) nP3 = *p3;
    if(*p3 == 0) nP3 = 1;
    
    gsl_matrix *X1 = gsl_matrix_calloc(*n, (nP1));
    gsl_matrix *X2 = gsl_matrix_calloc(*n, (nP2));
    gsl_matrix *X3 = gsl_matrix_calloc(*n, (nP3));
    
    for(i = 0; i < *n; i++)
    {
        gsl_vector_set(y1L, i, Ymat[(0 * *n) + i]);
        gsl_vector_set(y1U, i, Ymat[(1 * *n) + i]);
        gsl_vector_set(y2L, i, Ymat[(2 * *n) + i]);
        gsl_vector_set(y2U, i, Ymat[(3 * *n) + i]);
        gsl_vector_set(c0, i, Ymat[(4 * *n) + i]);
        gsl_vector_set(y1L_neginf, i, y1LInf[i]);
        gsl_vector_set(y1U_posinf, i, y1UInf[i]);
        gsl_vector_set(y2L_neginf, i, y2LInf[i]);
        gsl_vector_set(y2U_posinf, i, y2UInf[i]);
        gsl_vector_set(c0_neginf, i, c0Inf[i]);
        if(*p1 >0)
        {
            for(j = 0; j < *p1; j++)
            {
                gsl_matrix_set(X1, i, j, X1mat[(j* *n) + i]);
            }
        }
        if(*p2 >0)
        {
            for(j = 0; j < *p2; j++)
            {
                gsl_matrix_set(X2, i, j, X2mat[(j* *n) + i]);
            }
        }
        if(*p3 >0)
        {
            for(j = 0; j < *p3; j++)
            {
                gsl_matrix_set(X3, i, j, X3mat[(j* *n) + i]);
            }
        }
    }
    
    /* Hyperparameters */

    double a_theta = hyperP[0];
    double b_theta = hyperP[1];
    
    double a_sigSq1 = hyperP[2];
    double b_sigSq1 = hyperP[3];
    double a_sigSq2 = hyperP[4];
    double b_sigSq2 = hyperP[5];
    double a_sigSq3 = hyperP[6];
    double b_sigSq3 = hyperP[7];
    
    double beta1_prop_var = mcmcP[0];
    double beta2_prop_var = mcmcP[1];
    double beta3_prop_var = mcmcP[2];
    double beta01_prop_var = mcmcP[3];
    double beta02_prop_var = mcmcP[4];
    double beta03_prop_var = mcmcP[5];
    double sigSq1_prop_var = mcmcP[6];
    double sigSq2_prop_var = mcmcP[7];
    double sigSq3_prop_var = mcmcP[8];
    double gamma_prop_var = mcmcP[9];
    
    /* Starting values */
    
    gsl_vector *y1 = gsl_vector_calloc(*n);
    gsl_vector *y2 = gsl_vector_calloc(*n);
    gsl_vector *gamma = gsl_vector_calloc(*n);
    gsl_vector *beta1 = gsl_vector_calloc(nP1);
    gsl_vector *beta2 = gsl_vector_calloc(nP2);
    gsl_vector *beta3 = gsl_vector_calloc(nP3);
    
    double beta01 = startValues[*n+*n+*p1+*p2+*p3];
    double beta02 = startValues[*n+*n+*p1+*p2+*p3+1];
    double beta03 = startValues[*n+*n+*p1+*p2+*p3+2];
    double sigSq1 = startValues[*n+*n+*p1+*p2+*p3+3];
    double sigSq2 = startValues[*n+*n+*p1+*p2+*p3+4];
    double sigSq3 = startValues[*n+*n+*p1+*p2+*p3+5];
    double theta = startValues[*n+*n+*p1+*p2+*p3+3+3+*n];
    
    for(i = 0; i < *n; i++)
    {
        gsl_vector_set(y1, i, startValues[i]);
        gsl_vector_set(y2, i, startValues[(1 * *n) + i]);
        gsl_vector_set(gamma, i, startValues[*n+*n+*p1+*p2+*p3+3+3+ i]);
    }
    
    if(*p1 > 0)
    {
        for(i = 0; i < *p1; i++) gsl_vector_set(beta1, i, startValues[*n+*n+i]);
    }
    if(*p2 > 0)
    {
        for(i = 0; i < *p2; i++) gsl_vector_set(beta2, i, startValues[*n+*n+*p1+i]);
    }
    if(*p3 > 0)
    {
        for(i = 0; i < *p3; i++) gsl_vector_set(beta3, i, startValues[*n+*n+*p1+*p2+i]);
    }
    
    
    /* Variables required for storage of samples */
    
    int StoreInx;
    gsl_vector *accept_beta1 = gsl_vector_calloc(nP1);
    gsl_vector *accept_beta2 = gsl_vector_calloc(nP2);
    gsl_vector *accept_beta3 = gsl_vector_calloc(nP3);
    gsl_vector *accept_gamma = gsl_vector_calloc(*n);
    int accept_beta01 = 0;
    int accept_beta02 = 0;
    int accept_beta03 = 0;
    int accept_sigSq1 = 0;
    int accept_sigSq2 = 0;
    int accept_sigSq3 = 0;
    
    /* For posterior predictive checks */
    
    gsl_vector *invf_i_mean = gsl_vector_calloc(*n);
    gsl_vector *invf_i = gsl_vector_calloc(*n);
    gsl_vector *f_i_mean = gsl_vector_calloc(*n);
    gsl_vector *f_i = gsl_vector_calloc(*n);
    
    double val, logLH;
    

    /* selecting a move */
    /* move: 1= y, 2= beta1, 3= beta2, 4= beta3, 5= beta01, 6= beta02, 7= beta03, 8= sigSq1, 9= sigSq2, 10= sigSq3, 11= gamma, 12= theta */
    
    double pY = 0.05;
    double pBeta1 = 0.1;
    double pBeta2 = 0.1;
    double pBeta3 = 0.1;
    double pBeta01 = 0.05;
    double pBeta02 = 0.05;
    double pBeta03 = 0.05;
    double pSigSq1 = 0.05;
    double pSigSq2 = 0.05;
    double pSigSq3 = 0.05;
    double pGamma = 0.1;
    double pTheta = 1 - (pY+pBeta1+pBeta2+pBeta3+pBeta01+pBeta02+pBeta03+pSigSq1+pSigSq2+pSigSq3+pGamma);
    int move;
    gsl_vector *moveProb = gsl_vector_calloc(12);
    gsl_vector_set(moveProb, 0, pY);
    gsl_vector_set(moveProb, 1, pBeta1);
    gsl_vector_set(moveProb, 2, pBeta2);
    gsl_vector_set(moveProb, 3, pBeta3);
    gsl_vector_set(moveProb, 4, pBeta01);
    gsl_vector_set(moveProb, 5, pBeta02);
    gsl_vector_set(moveProb, 6, pBeta03);
    gsl_vector_set(moveProb, 7, pSigSq1);
    gsl_vector_set(moveProb, 8, pSigSq2);
    gsl_vector_set(moveProb, 9, pSigSq3);
    gsl_vector_set(moveProb, 10, pGamma);
    gsl_vector_set(moveProb, 11, pTheta);
    
    
    
    for(M = 0; M < *numReps; M++)
    {
        
        move = c_multinom_sample(rr, moveProb, 12);
        
        if(move == 1)
        {
            /* Updating y1 and y2         */
            BAFT_LNscr_update_y(y1L, y1U, y2L, y2U, y1L_neginf, y2L_neginf, y1U_posinf, y2U_posinf, y1_NA, X1, X2, X3, y1, y2, beta1, beta2, beta3, gamma, beta01, beta02, beta03, sigSq1, sigSq2, sigSq3);
        }
        
        if(move == 2)
        {
            /* Updating beta1*/
            if(*p1 > 0)
            {
                BAFT_LNscr_update_beta1(y1_NA, c0, c0_neginf, X1, X2, X3, y1, y2, beta1, beta2, beta3, gamma, beta01, beta02, beta03, sigSq1, sigSq2, sigSq3, beta1_prop_var, accept_beta1);
            }
        }
        
        if(move == 3)
        {
            /* Updating beta2         */
            if(*p2 > 0)
            {
                BAFT_LNscr_update_beta2(y1_NA, c0, c0_neginf, X1, X2, X3, y1, y2, beta1, beta2, beta3, gamma, beta01, beta02, beta03, sigSq1, sigSq2, sigSq3, beta2_prop_var, accept_beta2);
            }
        }
        
        if(move == 4)
        {
            /* Updating beta3         */
            if(*p3 > 0)
            {
                BAFT_LNscr_update_beta3(y1_NA, c0, c0_neginf, X1, X2, X3, y1, y2, beta1, beta2, beta3, gamma, beta01, beta02, beta03, sigSq1, sigSq2, sigSq3, beta3_prop_var, accept_beta3);
            }
        }
        
        if(move == 5)
        {
            /* Updating beta01 */
            BAFT_LNscr_update_beta01(y1_NA, c0, c0_neginf, X1, X2, X3, y1, y2, beta1, beta2, beta3, gamma, &beta01, beta02, beta03, sigSq1, sigSq2, sigSq3, beta01_prop_var, &accept_beta01);
        }
        
        if(move == 6)
        {
            /* Updating beta02*/
            BAFT_LNscr_update_beta02(y1_NA, c0, c0_neginf, X1, X2, X3, y1, y2, beta1, beta2, beta3, gamma, beta01, &beta02, beta03, sigSq1, sigSq2, sigSq3, beta02_prop_var, &accept_beta02);
        }
        
        if(move == 7)
        {
            /* Updating beta03*/
            BAFT_LNscr_update_beta03(y1_NA, c0, c0_neginf, X1, X2, X3, y1, y2, beta1, beta2, beta3, gamma, beta01, beta02, &beta03, sigSq1, sigSq2, sigSq3, beta03_prop_var, &accept_beta03);
        }
        
        if(move == 8)
        {
            /* Updating sigSq1*/
            BAFT_LNscr_update_sigSq1(y1_NA, c0, c0_neginf, X1, X2, X3, y1, y2, beta1, beta2, beta3, gamma, beta01, beta02, beta03, &sigSq1, sigSq2, sigSq3, a_sigSq1, b_sigSq1, sigSq1_prop_var, &accept_sigSq1);
        }
        
        if(move == 9)
        {
            /* Updating sigSq2*/
            BAFT_LNscr_update_sigSq2(y1_NA, c0, c0_neginf, X1, X2, X3, y1, y2, beta1, beta2, beta3, gamma, beta01, beta02, beta03, sigSq1, &sigSq2, sigSq3, a_sigSq2, b_sigSq2, sigSq2_prop_var, &accept_sigSq2);
        }
        
        if(move == 10)
        {
            /* Updating sigSq3*/
            BAFT_LNscr_update_sigSq3(y1_NA, c0, c0_neginf, X1, X2, X3, y1, y2, beta1, beta2, beta3, gamma, beta01, beta02, beta03, sigSq1, sigSq2, &sigSq3, a_sigSq3, b_sigSq3, sigSq3_prop_var, &accept_sigSq3);
        }
        
        if(move == 11)
        {
            /* Updating gamma*/
            BAFT_LNscr_update_gamma(y1_NA, c0, c0_neginf, X1, X2, X3, y1, y2, beta1, beta2, beta3, gamma, beta01, beta02, beta03, sigSq1, sigSq2, sigSq3, theta, gamma_prop_var, accept_gamma);
        }

        if(move == 12)
        {
            /* Updating theta*/
            BAFT_LNscr_update_theta(gamma, &theta, a_theta, b_theta);
        }

        
        /* Storing posterior samples */
        
        if( ( (M+1) % *thin ) == 0 && (M+1) > (*numReps * *burninPerc))
        {
            StoreInx = (M+1)/(*thin)- (*numReps * *burninPerc)/(*thin);
            
            for(i = 0; i < *n; i++)
            {
                samples_y1[(StoreInx - 1) * (*n) + i] = gsl_vector_get(y1, i);
                samples_y2[(StoreInx - 1) * (*n) + i] = gsl_vector_get(y2, i);
                samples_y1_NA[(StoreInx - 1) * (*n) + i] = gsl_vector_get(y1_NA, i);
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
            samples_beta01[StoreInx - 1] = beta01;
            samples_beta02[StoreInx - 1] = beta02;
            samples_beta03[StoreInx - 1] = beta03;
            samples_sigSq1[StoreInx - 1] = sigSq1;
            samples_sigSq2[StoreInx - 1] = sigSq2;
            samples_sigSq3[StoreInx - 1] = sigSq3;
            samples_theta[StoreInx - 1] = theta;
            
            /* Deviance */
            
            logLH = 0;
            for(i = 0; i < *n; i++)
            {
                if(gsl_vector_get(y1_NA, i) == 0)
                {
                    log_Jpdf_Upper_BAFT_LN(i, gsl_vector_get(y1, i), gsl_vector_get(y2, i), gsl_vector_get(c0, i), c0_neginf, X1, X2, X3, beta1, beta2, beta3, gamma, beta01, beta02, beta03, sigSq1, sigSq2, sigSq3, &val);
                }else
                {
                    log_Jpdf_Lower_BAFT_LN(i, gsl_vector_get(y2, i), gsl_vector_get(c0, i), c0_neginf, X1, X2, beta1, beta2, gamma, beta01, beta02, sigSq1, sigSq2, &val);
                }
                logLH += val;
                gsl_vector_set(f_i, i, exp(val));
                gsl_vector_set(invf_i, i, 1/exp(val));
            }
            
            gsl_vector_scale(f_i_mean, (double) StoreInx - 1);
            gsl_vector_add(f_i_mean, f_i);
            gsl_vector_scale(f_i_mean, (double)1/StoreInx);
            
            gsl_vector_scale(invf_i_mean, (double) StoreInx - 1);
            gsl_vector_add(invf_i_mean, invf_i);
            gsl_vector_scale(invf_i_mean, (double)1/StoreInx);
            
            samples_logLH[StoreInx - 1] = logLH;
            
            
            
        }
        
        if(M == (*numReps - 1))
        {
            for(i = 0; i < *p1; i++)
            {
                samples_misc[i] = (int) gsl_vector_get(accept_beta1, i);
            }
            for(i = 0; i < *p2; i++)
            {
                samples_misc[*p1+i] = (int) gsl_vector_get(accept_beta2, i);
            }
            for(i = 0; i < *p3; i++)
            {
                samples_misc[*p1+*p2+i] = (int) gsl_vector_get(accept_beta3, i);
            }
            
            samples_misc[*p1+*p2+*p3] = (int) accept_beta01;
            samples_misc[*p1+*p2+*p3+1] = (int) accept_beta02;
            samples_misc[*p1+*p2+*p3+2] = (int) accept_beta03;
            samples_misc[*p1+*p2+*p3+3] = (int) accept_sigSq1;
            samples_misc[*p1+*p2+*p3+4] = (int) accept_sigSq2;
            samples_misc[*p1+*p2+*p3+5] = (int) accept_sigSq3;
            
            for(i = 0; i < *n; i++)
            {
                samples_misc[*p1+*p2+*p3+6+i] = (int) gsl_vector_get(accept_gamma, i);
            }            
            
            for(i = 0; i < *n; i++)
            {
                LH_i_mean[i] = gsl_vector_get(f_i, i);
                invLH_i_mean[i] = gsl_vector_get(invf_i, i);
            }
            
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





















