/*
 TO COMPILE USE THE CODE:
 
 R CMD SHLIB BAFT_DPscr.c BAFT_DPscr_Updates.c BAFT_DPscr_Utilities.c -lgsl -lgslcblas
 
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
#include "BAFT_DPscr.h"


#ifdef NAN
/* NAN is supported */
#endif
#ifdef INFINITY
/* INFINITY is supported */
#endif


/* */
void BAFT_DPscrmcmc(double Ymat[],
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
                    int *nY1_save,
                    int *nY2_save,
                    int *nY1_NA_save,
                    double samples_y1[],
                    double samples_y2[],
                    double samples_y1_NA[],
                    double samples_beta1[],
                    double samples_beta2[],
                    double samples_beta3[],
                    double samples_r1[],
                    double samples_r2[],
                    double samples_r3[],
                    double samples_beta01[],
                    double samples_beta02[],
                    double samples_beta03[],
                    double samples_sigSq1[],
                    double samples_sigSq2[],
                    double samples_sigSq3[],
                    double samples_tau1[],
                    double samples_tau2[],
                    double samples_tau3[],
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
    
    double a01 = hyperP[2];
    double b01 = hyperP[3];
    double a02 = hyperP[4];
    double b02 = hyperP[5];
    double a03 = hyperP[6];
    double b03 = hyperP[7];
    
    double aTau1 = hyperP[8];
    double bTau1 = hyperP[9];
    double aTau2 = hyperP[10];
    double bTau2 = hyperP[11];
    double aTau3 = hyperP[12];
    double bTau3 = hyperP[13];
    double mu01  = hyperP[14];
    double mu02  = hyperP[15];
    double mu03  = hyperP[16];
    double sigSq01  = hyperP[17];
    double sigSq02  = hyperP[18];
    double sigSq03  = hyperP[19];
    
    double beta1_prop_var = mcmcP[0];
    double beta2_prop_var = mcmcP[1];
    double beta3_prop_var = mcmcP[2];
    double beta01_prop_var = mcmcP[3];
    double beta02_prop_var = mcmcP[4];
    double beta03_prop_var = mcmcP[5];
    double zeta1_prop_var = mcmcP[6];
    double zeta2_prop_var = mcmcP[7];
    double zeta3_prop_var = mcmcP[8];
    double gamma_prop_var = mcmcP[9];
    
    /* Starting values */
    
    gsl_vector *y1 = gsl_vector_calloc(*n);
    gsl_vector *y2 = gsl_vector_calloc(*n);
    gsl_vector *beta1 = gsl_vector_calloc(nP1);
    gsl_vector *beta2 = gsl_vector_calloc(nP2);
    gsl_vector *beta3 = gsl_vector_calloc(nP3);
    gsl_vector *gamma = gsl_vector_calloc(*n);
    gsl_vector *r1 = gsl_vector_calloc(*n);
    gsl_vector *r2 = gsl_vector_calloc(*n);
    gsl_vector *r3 = gsl_vector_calloc(*n);
    
    
    for(i = 0; i < *n; i++)
    {
        gsl_vector_set(y1, i, startValues[i]);
        gsl_vector_set(y2, i, startValues[(1 * *n) + i]);
        gsl_vector_set(r1, i, startValues[*n+*n+*p1+*p2+*p3+i]);
        gsl_vector_set(r2, i, startValues[*n+*n+*p1+*p2+*p3+*n+i]);
        gsl_vector_set(r3, i, startValues[*n+*n+*p1+*p2+*p3+*n+*n+i]);
        gsl_vector_set(gamma, i, startValues[*n+*n+*p1+*p2+*p3+*n+*n+*n+3+i]);
    }
    
    for(i = 0; i < *n; i++)
    {
        if(gsl_vector_get(y1, i) < gsl_vector_get(y2, i))
        {
            gsl_vector_set(y1_NA, i, 0);
        }else
        {
            gsl_vector_set(y1_NA, i, 1);
        }
    }
    
    
    if(*p1 > 0)
    {
        for(j = 0; j < *p1; j++) gsl_vector_set(beta1, j, startValues[*n+*n+j]);
    }
    if(*p2 > 0)
    {
        for(j = 0; j < *p2; j++) gsl_vector_set(beta2, j, startValues[j + *n+*n+*p1]);
    }
    if(*p3 > 0)
    {
        for(j = 0; j < *p3; j++) gsl_vector_set(beta3, j, startValues[j + *n+*n+*p1+*p2]);
    }
    
    double tau1 = startValues[*n+*n+*p1+*p2+*p3+*n+*n+*n];
    double tau2 = startValues[*n+*n+*p1+*p2+*p3+*n+*n+*n+1];
    double tau3 = startValues[*n+*n+*p1+*p2+*p3+*n+*n+*n+2];
    double theta = startValues[*n+*n+*p1+*p2+*p3+*n+*n+*n+3+*n];
    
    
    
    /* Variables required for storage of samples */
    
    int StoreInx;
    gsl_vector *accept_beta1 = gsl_vector_calloc(nP1);
    gsl_vector *accept_beta2 = gsl_vector_calloc(nP2);
    gsl_vector *accept_beta3 = gsl_vector_calloc(nP3);
    gsl_vector *accept_gamma = gsl_vector_calloc(*n);
    
    int accept_beta01 = 0;
    int accept_beta02 = 0;
    int accept_beta03 = 0;
    int accept_zeta1 = 0;
    int accept_zeta2 = 0;
    int accept_zeta3 = 0;
    
    gsl_vector *mu1_all = gsl_vector_calloc(*n);
    gsl_vector *zeta1_all = gsl_vector_calloc(*n);
    gsl_vector *mu2_all = gsl_vector_calloc(*n);
    gsl_vector *zeta2_all = gsl_vector_calloc(*n);
    gsl_vector *mu3_all = gsl_vector_calloc(*n);
    gsl_vector *zeta3_all = gsl_vector_calloc(*n);
    gsl_vector *mu3_vec = gsl_vector_calloc(*n);
    gsl_vector *zeta3_vec = gsl_vector_calloc(*n);
    
    gsl_vector *r1Uniq = gsl_vector_calloc(*n);
    gsl_vector *r1Uniq_count = gsl_vector_calloc(*n);
    gsl_vector *r2Uniq = gsl_vector_calloc(*n);
    gsl_vector *r2Uniq_count = gsl_vector_calloc(*n);
    gsl_vector *r3Uniq = gsl_vector_calloc(*n);
    gsl_vector *r3Uniq_count = gsl_vector_calloc(*n);
    
    int nClass_DP1;
    int nClass_DP2;
    int nClass_DP3;
    
    
    /* For posterior predictive checks */
    
    gsl_vector *invf_i_mean = gsl_vector_calloc(*n);
    gsl_vector *invf_i = gsl_vector_calloc(*n);
    gsl_vector *f_i_mean = gsl_vector_calloc(*n);
    gsl_vector *f_i = gsl_vector_calloc(*n);
    
    double val, logLH;
    
    
    
    for(i = 0; i < *n; i++)
    {
        if(gsl_vector_get(y1, i) >= gsl_vector_get(y2, i))
        {
            gsl_vector_set(r3, i, 0);
        }
    }
    
    c_uniq1(r1, r1Uniq, r1Uniq_count, &nClass_DP1);
    c_uniq1(r2, r2Uniq, r2Uniq_count, &nClass_DP2);
    c_uniq1_h3(r3, r3Uniq, r3Uniq_count, y1_NA, &nClass_DP3);
    
    for(i = 0; i < *n; i++)
    {
        if(gsl_vector_get(y1_NA, i) == 0)
        {
            gsl_vector_set(mu3_vec, i, 0.1);
            gsl_vector_set(zeta3_vec, i, 2);
        }
    }
    
    for(i = 0; i < nClass_DP1; i++)
    {
        gsl_vector_set(mu1_all, i, startValues[*n+*n+*p1+*p2+*p3+*n+*n+*n+3+*n+1+i]);
        gsl_vector_set(zeta1_all, i, startValues[*n+*n+*p1+*p2+*p3+*n+*n+*n+3+*n+*n+*n+*n+1+i]);
    }
    
    for(i = 0; i < nClass_DP2; i++)
    {
        gsl_vector_set(mu2_all, i, startValues[*n+*n+*p1+*p2+*p3+*n+*n+*n+3+*n+*n+1+i]);
        gsl_vector_set(zeta2_all, i, startValues[*n+*n+*p1+*p2+*p3+*n+*n+*n+3+*n+*n+*n+*n+*n+1+i]);
    }
    
    c_uniq_h3(r3, r3Uniq, r3Uniq_count, mu3_all, zeta3_all, mu3_vec, zeta3_vec, y1_NA, &nClass_DP3);
    
    /*     */
    BAFT_DPscr_update_y(y1L, y1U, y2L, y2U, y1L_neginf, y2L_neginf, y1U_posinf, y2U_posinf, c0_neginf, y1_NA, X1, X2, X3, y1, y2, beta1, beta2, beta3, r1, r2, r3, mu1_all, mu2_all, mu3_all, mu3_vec, zeta1_all, zeta2_all, zeta3_all, zeta3_vec, r1Uniq, r2Uniq, r3Uniq, r1Uniq_count, r2Uniq_count, r3Uniq_count, gamma, &nClass_DP1, &nClass_DP2, &nClass_DP3, mu03, sigSq03, a03, b03, tau3, rr);
    
    c_uniq_h3(r3, r3Uniq, r3Uniq_count, mu3_all, zeta3_all, mu3_vec, zeta3_vec, y1_NA, &nClass_DP3);
    
    
    /* selecting a move */
    /* move: 1= DPM1, 2= DPM2, 3= DPM3, 4= y, 5= tau1, 6= tau2, 7= tau3, 8= beta1, 9= beta2, 10= beta3, 11= gamma, 12= theta */
    
    double pDPM1 = 0.1;
    double pDPM2 = 0.1;
    double pDPM3 = 0.05;
    double pY = 0.05;
    double pTau1 = 0.05;
    double pTau2 = 0.05;
    double pTau3 = 0.05;
    double pBeta1 = 0.08;
    double pBeta2 = 0.12;
    double pBeta3 = 0.08;
    double pGamma = 0.1;
    double pTheta = 1 - (pDPM1+pDPM2+pDPM3+pY+pTau1+pTau2+pTau3+pBeta1+pBeta2+pBeta3+pGamma);
    int move;
    gsl_vector *moveProb = gsl_vector_calloc(12);
    gsl_vector_set(moveProb, 0, pDPM1);
    gsl_vector_set(moveProb, 1, pDPM2);
    gsl_vector_set(moveProb, 2, pDPM3);
    gsl_vector_set(moveProb, 3, pY);
    gsl_vector_set(moveProb, 4, pTau1);
    gsl_vector_set(moveProb, 5, pTau2);
    gsl_vector_set(moveProb, 6, pTau3);
    gsl_vector_set(moveProb, 7, pBeta1);
    gsl_vector_set(moveProb, 8, pBeta2);
    gsl_vector_set(moveProb, 9, pBeta3);
    gsl_vector_set(moveProb, 10, pGamma);
    gsl_vector_set(moveProb, 11, pTheta);
    
    
    for(M = 0; M < *numReps; M++)
    {
        move = c_multinom_sample(rr, moveProb, 12);
        
        if(move == 1)
        {
            /* update mu1 and zeta1         */
            BAFT_DPscr_update_mu_zeta1(c0, c0_neginf, X1, y1, y2, y1_NA, beta1, gamma, r1, mu1_all, zeta1_all, r1Uniq, r1Uniq_count, tau1, a01, b01, mu01, sigSq01, beta01_prop_var, zeta1_prop_var, &accept_beta01, &accept_zeta1, &nClass_DP1, rr);
        }
        
        
        
        if(move == 2)
        {
            /* update mu2 and zeta2         */
            BAFT_DPscr_update_mu_zeta2(c0, c0_neginf, X2, y1, y2, y1_NA, beta2, gamma, r2, mu2_all, zeta2_all, r2Uniq, r2Uniq_count, tau2, a02, b02, mu02, sigSq02, beta02_prop_var, zeta2_prop_var, &accept_beta02, &accept_zeta2, &nClass_DP2, rr);
        }
        
        
        
        if(move == 3)
        {
            /* update mu3 and zeta3         */
            BAFT_DPscr_update_mu_zeta3(c0, c0_neginf, X3, y1, y2, y1_NA, beta3, gamma, r3, mu3_all, mu3_vec, zeta3_all, zeta3_vec, r3Uniq, r3Uniq_count, tau3, a03, b03, mu03,  sigSq03, beta03_prop_var, zeta3_prop_var, &accept_beta03, &accept_zeta3, &nClass_DP3, rr);
        }
        
        if(move == 4)
        {
            /* Updating y1 and y2          */
            BAFT_DPscr_update_y(y1L, y1U, y2L, y2U, y1L_neginf, y2L_neginf, y1U_posinf, y2U_posinf, c0_neginf, y1_NA, X1, X2, X3, y1, y2, beta1, beta2, beta3, r1, r2, r3, mu1_all, mu2_all, mu3_all, mu3_vec, zeta1_all, zeta2_all, zeta3_all, zeta3_vec, r1Uniq, r2Uniq, r3Uniq, r1Uniq_count, r2Uniq_count, r3Uniq_count, gamma, &nClass_DP1, &nClass_DP2, &nClass_DP3, mu03, sigSq03, a03, b03, tau3, rr);
        }
        
        
        
        if(move == 5)
        {
            /* update tau1         */
            BAFT_DPscr_update_tau(n, &tau1, aTau1, bTau1, &nClass_DP1);
        }
        
        if(move == 6)
        {
            /* update tau2             */
            BAFT_DPscr_update_tau(n, &tau2, aTau2, bTau2, &nClass_DP2);
        }
        
        if(move == 7)
        {
            /* update tau3         */
            BAFT_DPscr_update_tau(n, &tau3, aTau3, bTau3, &nClass_DP3);
        }
        
        if(move == 8)
        {
            /* Updating beta1         */
            if(*p1 > 0)
            {
                BAFT_DPscr_update_beta1(y1_NA, c0, c0_neginf, X1, y1, y2, beta1, gamma, r1, mu1_all, zeta1_all, r1Uniq, r1Uniq_count, &nClass_DP1, beta1_prop_var, accept_beta1);
            }
        }
        
        
        if(move == 9)
        {
            /* Updating beta2          */
            if(*p2 > 0)
            {
                BAFT_DPscr_update_beta2(y1_NA, c0, c0_neginf, X2, y1, y2, beta2, gamma, r2, mu2_all, zeta2_all, r2Uniq, r2Uniq_count, &nClass_DP2, beta2_prop_var, accept_beta2);
            }
        }
        
        if(move == 10)
        {
            /* Updating beta3         */
            if(*p3 > 0)
            {
                BAFT_DPscr_update_beta3(y1_NA, c0, c0_neginf, X3, y1, y2, beta3, gamma, r3, mu3_all, zeta3_all, r3Uniq, r3Uniq_count, &nClass_DP3, beta3_prop_var, accept_beta3);
            }
        }
        
        if(move == 11)
        {
            /* Updating gamma                  */
            BAFT_DPscr_update_gamma(y1_NA, c0, c0_neginf, X1, X2, X3, y1, y2, beta1, beta2, beta3, gamma, r1, r2, r3, mu1_all, mu2_all, mu3_all, zeta1_all, zeta2_all, zeta3_all, r1Uniq, r2Uniq, r3Uniq, r1Uniq_count, r2Uniq_count, r3Uniq_count, &nClass_DP1, &nClass_DP2, &nClass_DP3, theta, gamma_prop_var, accept_gamma);
        }
        
        if(move == 12)
        {
            /* Updating theta         */
            BAFT_DPscr_update_theta(gamma, &theta, a_theta, b_theta);
        }
        
        
        /* Storing posterior samples */
        
        if( ( (M+1) % *thin ) == 0 && (M+1) > (*numReps * *burninPerc))
        {
            StoreInx = (M+1)/(*thin)- (*numReps * *burninPerc)/(*thin);
            
            if(*nY1_save == *n)
            {
                for(i = 0; i < *n; i++)
                {
                    samples_y1[(StoreInx - 1) * (*n) + i] = gsl_vector_get(y1, i);
                }
            }
            if(*nY1_save < *n)
            {
                for(i = 0; i < *nY1_save; i++)
                {
                    samples_y1[(StoreInx - 1) * (*nY1_save) + i] = gsl_vector_get(y1, i);
                }
            }
            if(*nY2_save == *n)
            {
                for(i = 0; i < *n; i++)
                {
                    samples_y2[(StoreInx - 1) * (*n) + i] = gsl_vector_get(y2, i);
                }
            }
            if(*nY2_save < *n)
            {
                for(i = 0; i < *nY2_save; i++)
                {
                    samples_y2[(StoreInx - 1) * (*nY2_save) + i] = gsl_vector_get(y2, i);
                }
            }
            if(*nY1_NA_save == *n)
            {
                for(i = 0; i < *n; i++)
                {
                    samples_y1_NA[(StoreInx - 1) * (*n) + i] = gsl_vector_get(y1_NA, i);
                }
            }
            if(*nY1_NA_save < *n)
            {
                for(i = 0; i < *nY1_NA_save; i++)
                {
                    samples_y1_NA[(StoreInx - 1) * (*nY1_NA_save) + i] = gsl_vector_get(y1_NA, i);
                }
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
            
            for(j = 0; j < *n; j++)
            {
                samples_r1[(StoreInx - 1) * (*n) + j] = gsl_vector_get(r1, j);
                samples_r2[(StoreInx - 1) * (*n) + j] = gsl_vector_get(r2, j);
                samples_r3[(StoreInx - 1) * (*n) + j] = gsl_vector_get(r3, j);
                
                samples_beta01[(StoreInx - 1) * (*n) + j] = gsl_vector_get(mu1_all, j);
                samples_beta02[(StoreInx - 1) * (*n) + j] = gsl_vector_get(mu2_all, j);
                samples_beta03[(StoreInx - 1) * (*n) + j] = gsl_vector_get(mu3_all, j);
                
                samples_sigSq1[(StoreInx - 1) * (*n) + j] = 1/gsl_vector_get(zeta1_all, j);
                samples_sigSq2[(StoreInx - 1) * (*n) + j] = 1/gsl_vector_get(zeta2_all, j);
                samples_sigSq3[(StoreInx - 1) * (*n) + j] = 1/gsl_vector_get(zeta3_all, j);
            }
            
            samples_tau1[StoreInx - 1] = tau1;
            samples_tau2[StoreInx - 1] = tau2;
            samples_tau3[StoreInx - 1] = tau3;
            samples_theta[StoreInx - 1] = theta;
            
            
            /* Deviance */
            
            logLH = 0;
            for(i = 0; i < *n; i++)
            {
                if(gsl_vector_get(y1_NA, i) == 0)
                {
                    log_Jpdf_Upper_BAFT_DP(i, gsl_vector_get(y1, i), gsl_vector_get(y2, i), gsl_vector_get(c0, i), c0_neginf, X1, X2, X3, beta1, beta2, beta3, gamma, mu1_all, mu2_all, mu3_all, zeta1_all, zeta2_all, zeta3_all, r1Uniq, r2Uniq, r3Uniq, r1Uniq_count, r2Uniq_count, r3Uniq_count, nClass_DP1, nClass_DP2, nClass_DP3, &val);
                }else
                {
                    log_Jpdf_Lower_BAFT_DP(i, gsl_vector_get(y2, i), gsl_vector_get(c0, i), c0_neginf, X1, X2, beta1, beta2, gamma, mu1_all, mu2_all, zeta1_all, zeta2_all, r1Uniq, r2Uniq, r1Uniq_count, r2Uniq_count, nClass_DP1, nClass_DP2, &val);
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
            
            for(i = 0; i < *n; i++)
            {
                samples_misc[*p1+*p2+*p3+i] = (int) gsl_vector_get(accept_gamma, i);
            }
            samples_misc[*p1+*p2+*p3+*n] = (int) accept_beta01;
            samples_misc[*p1+*p2+*p3+*n+1] = (int) accept_beta02;
            samples_misc[*p1+*p2+*p3+*n+2] = (int) accept_beta03;
            samples_misc[*p1+*p2+*p3+*n+3] = (int) accept_zeta1;
            samples_misc[*p1+*p2+*p3+*n+4] = (int) accept_zeta2;
            samples_misc[*p1+*p2+*p3+*n+5] = (int) accept_zeta3;
            
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





















