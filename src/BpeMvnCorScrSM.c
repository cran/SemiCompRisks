/*
 TO COMPILE USE THE CODE:
 
 R CMD SHLIB BpeMvnCorScrSM.c BpeMvnCorScrSM_Updates.c BpeMvnCorScrSM_Utilities.c -lgsl -lgslcblas
 
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

#include "BpeMvnCorScrSM.h"



/* */
void BpeMvnCorScrSMmcmc(double survData[],
                      int *n,
                      int *p1,
                      int *p2,
                      int *p3,
                      int *J,
                      double nj[],
                      double hyperParams[],
                      double startValues[],
                      double mcmcParams[],
                      int *numReps,
                      int *thin,
                      double *burninPerc,
                      int *nGam_save,
                      double samples_beta1[],
                      double samples_beta2[],
                      double samples_beta3[],
                      double samples_mu_lam1[],
                      double samples_mu_lam2[],
                      double samples_mu_lam3[],
                      double samples_sigSq_lam1[],
                      double samples_sigSq_lam2[],
                      double samples_sigSq_lam3[],
                      double samples_K1[],
                      double samples_K2[],
                      double samples_K3[],
                      double samples_s1[],
                      double samples_s2[],
                      double samples_s3[],
                      double samples_nu2[],
                      double samples_nu3[],
                      double samples_theta[],
                      double samples_V1[],
                      double samples_V2[],
                      double samples_V3[],
                      double samples_Sigma_V[],
                      double samples_gamma[],
                      double samples_gamma_last[],
                      double samples_misc[],
                      double lambda1_fin[],
                      double lambda2_fin[],
                      double lambda3_fin[],
                        double gammaP[],
                        double dev[],
                        double invLH[],
                        double *logLH_fin,
                        double lpml[],
                        double moveVec[])
{
    GetRNGstate();
    
    time_t now;    
    
    int i, j, MM;

    
    /* Survival Data */
    
    gsl_vector *survTime1    = gsl_vector_alloc(*n);
    gsl_vector *survTime2    = gsl_vector_alloc(*n);
    gsl_vector *survEvent1   = gsl_vector_alloc(*n);
    gsl_vector *survEvent2   = gsl_vector_alloc(*n);
    gsl_vector *cluster      = gsl_vector_alloc(*n);
    for(i = 0; i < *n; i++)
    {
        gsl_vector_set(survTime1, i, survData[(0 * *n) + i]);
        gsl_vector_set(survEvent1, i, survData[(1* *n) + i]);
        gsl_vector_set(survTime2, i, survData[(2 * *n) + i]);
        gsl_vector_set(survEvent2, i, survData[(3* *n) + i]);
        gsl_vector_set(cluster, i, survData[(4* *n) + i]);
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
                gsl_matrix_set(survCov1, i, j, survData[((5+j)* *n) + i]);
            }
        }
    }
    
    if(*p2 >0)
    {
        for(i = 0; i < *n; i++)
        {
            for(j = 0; j < *(p2); j++)
            {
                gsl_matrix_set(survCov2, i, j, survData[((5+(*p1)+j)* *n) + i]);
            }
        }
    }
    
    if(*p3 >0)
    {
        for(i = 0; i < *n; i++)
        {
            for(j = 0; j < *(p3); j++)
            {
                gsl_matrix_set(survCov3, i, j, survData[((5+(*p1)+(*p2)+j)* *n) + i]);
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
    
    
    gsl_vector *n_j = gsl_vector_calloc(*J);
    
    for(j = 0; j < *J; j++)
    {
        gsl_vector_set(n_j, j, nj[j]);
    }
    
    /*     
    
    for(i = 0; i < 20; i++)
    {
        printf("cluster = %.1f\n", gsl_vector_get(cluster, i));
    }
    for(i = 0; i < 20; i++)
    {
        printf("n_j = %d\n", (int) nj[i]);
    }
     */

    
    
    /* Hyperparameters */
    
    double a1       = hyperParams[0];
    double b1       = hyperParams[1];
    double a2       = hyperParams[2];
    double b2       = hyperParams[3];
    double a3       = hyperParams[4];
    double b3       = hyperParams[5];
    double alpha1    = hyperParams[6];
    double alpha2    = hyperParams[7];
    double alpha3    = hyperParams[8];
    double c_lam1    = hyperParams[9];
    double c_lam2    = hyperParams[10];
    double c_lam3    = hyperParams[11];
    double psi      = hyperParams[12];
    double omega    = hyperParams[13];
    double rho_v    = hyperParams[14];
    
    gsl_matrix *Psi_v = gsl_matrix_calloc(3,3);
    for(i = 0; i < 3; i++)
    {
        for(j = 0; j < 3; j++)
        {
            gsl_matrix_set(Psi_v, j, i, hyperParams[15 + i*3 + j]);
        }
    }
    
    
    /* varialbes for birth and death moves */
    
    double C1               = mcmcParams[0];
    double C2               = mcmcParams[1];
    double C3               = mcmcParams[2];
    double delPert1         = mcmcParams[3];
    double delPert2         = mcmcParams[4];
    double delPert3         = mcmcParams[5];
    int num_s_propBI1       = mcmcParams[6];
    int num_s_propBI2       = mcmcParams[7];
    int num_s_propBI3       = mcmcParams[8];
    int K1_max              = mcmcParams[9];
    int K2_max              = mcmcParams[10];
    int K3_max              = mcmcParams[11];
    double s1_max           = mcmcParams[12];
    double s2_max           = mcmcParams[13];
    double s3_max           = mcmcParams[14];
    
    gsl_vector *s_propBI1    = gsl_vector_calloc(num_s_propBI1);
    gsl_vector *s_propBI2    = gsl_vector_calloc(num_s_propBI2);
    gsl_vector *s_propBI3    = gsl_vector_calloc(num_s_propBI3);
    for(j = 0; j < num_s_propBI1; j++) gsl_vector_set(s_propBI1, j, mcmcParams[18+j]);
    for(j = 0; j < num_s_propBI2; j++) gsl_vector_set(s_propBI2, j, mcmcParams[18+num_s_propBI1+j]);
    for(j = 0; j < num_s_propBI3; j++) gsl_vector_set(s_propBI3, j, mcmcParams[18+num_s_propBI1+num_s_propBI2+j]);
    
    
    
    /* time points where lambda values are stored  */
    
    int nTime_lambda1 = mcmcParams[15];
    int nTime_lambda2 = mcmcParams[16];
    int nTime_lambda3 = mcmcParams[17];
    
    gsl_vector *time_lambda1 = gsl_vector_calloc(nTime_lambda1);
    gsl_vector *time_lambda2 = gsl_vector_calloc(nTime_lambda2);
    gsl_vector *time_lambda3 = gsl_vector_calloc(nTime_lambda3);
    
    for(i = 0; i < nTime_lambda1; i++)
    {
        gsl_vector_set(time_lambda1, i, mcmcParams[18+num_s_propBI1+num_s_propBI2+num_s_propBI3+i]);
    }
    for(i = 0; i < nTime_lambda2; i++)
    {
        gsl_vector_set(time_lambda2, i, mcmcParams[18+num_s_propBI1+num_s_propBI2+num_s_propBI3+nTime_lambda1+i]);
    }
    for(i = 0; i < nTime_lambda3; i++)
    {
        gsl_vector_set(time_lambda3, i, mcmcParams[18+num_s_propBI1+num_s_propBI2+num_s_propBI3+nTime_lambda1+nTime_lambda2+i]);
    }
    
    
    /* varialbes for M-H step */
    
    
    double mhProp_theta_var  = mcmcParams[18+num_s_propBI1+num_s_propBI2+num_s_propBI3+nTime_lambda1+nTime_lambda2+nTime_lambda3];
    

    double mhProp_V1_var  = mcmcParams[18+num_s_propBI1+num_s_propBI2+num_s_propBI3+nTime_lambda1+nTime_lambda2+nTime_lambda3+2];
    double mhProp_V2_var  = mcmcParams[18+num_s_propBI1+num_s_propBI2+num_s_propBI3+nTime_lambda1+nTime_lambda2+nTime_lambda3+3];
    double mhProp_V3_var  = mcmcParams[18+num_s_propBI1+num_s_propBI2+num_s_propBI3+nTime_lambda1+nTime_lambda2+nTime_lambda3+4];

    

    
    
    
    
    
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
    
    int K1              = startValues[*p1 + *p2 + *p3];
    int K2              = startValues[*p1 + *p2 + *p3 + 1];
    int K3              = startValues[*p1 + *p2 + *p3 + 2];
    double mu_lam1      = startValues[*p1 + *p2 + *p3 + 3];
    double mu_lam2      = startValues[*p1 + *p2 + *p3 + 4];
    double mu_lam3      = startValues[*p1 + *p2 + *p3 + 5];
    double sigSq_lam1   = startValues[*p1 + *p2 + *p3 + 6];
    double sigSq_lam2   = startValues[*p1 + *p2 + *p3 + 7];
    double sigSq_lam3   = startValues[*p1 + *p2 + *p3 + 8];
    double theta        = startValues[*p1 + *p2 + *p3 + 9];
    
    gsl_vector *gamma = gsl_vector_calloc(*n);
    for(i = 0; i < *n; i++)
    {
        gsl_vector_set(gamma, i, startValues[*p1 + *p2 + *p3 + 10 + i]);
    }
    
    
    gsl_vector *lambda1  = gsl_vector_calloc(K1_max+1);
    gsl_vector *lambda2  = gsl_vector_calloc(K2_max+1);
    gsl_vector *lambda3  = gsl_vector_calloc(K3_max+1);
    for(j = 0; j < (K1+1); j++) gsl_vector_set(lambda1, j, startValues[*p1 + *p2 + *p3 + 10 + *n + j]);
    for(j = 0; j < (K2+1); j++) gsl_vector_set(lambda2, j, startValues[*p1 + *p2 + *p3 + 10 + *n + (K1+1) + j]);
    for(j = 0; j < (K3+1); j++) gsl_vector_set(lambda3, j, startValues[*p1 + *p2 + *p3 + 10 + *n + (K1+1) + (K2+1) + j]);
    
    gsl_vector *s1       = gsl_vector_calloc(K1_max+1);
    gsl_vector *s2       = gsl_vector_calloc(K2_max+1);
    gsl_vector *s3       = gsl_vector_calloc(K3_max+1);
    for(j = 0; j < (K1+1); j++) gsl_vector_set(s1, j, startValues[*p1 + *p2 + *p3 + 10 + *n + (K1+1) + (K2+1) + (K3+1) + j]);
    for(j = 0; j < (K2+1); j++) gsl_vector_set(s2, j, startValues[*p1 + *p2 + *p3 + 10 + *n + (K1+1)*2 + (K2+1) + (K3+1) + j]);
    for(j = 0; j < (K3+1); j++) gsl_vector_set(s3, j, startValues[*p1 + *p2 + *p3 + 10 + *n + (K1+1)*2 + (K2+1)*2 + (K3+1) + j]);
    
    double nu2 = startValues[*p1 + *p2 + *p3 + 10 + *n + (K1+1)*2 + (K2+1)*2 + (K3+1)*2];
    double nu3 = startValues[*p1 + *p2 + *p3 + 10 + *n + (K1+1)*2 + (K2+1)*2 + (K3+1)*2 + 1];
    
    gsl_vector *V1 = gsl_vector_calloc(*J);
    gsl_vector *V2 = gsl_vector_calloc(*J);
    gsl_vector *V3 = gsl_vector_calloc(*J);
    
    for(j = 0; j < *J; j++)
    {
        gsl_vector_set(V1, j, startValues[*p1 + *p2 + *p3 + 10 + *n + (K1+1)*2 + (K2+1)*2 + (K3+1)*2 + 2 + j]);
        gsl_vector_set(V2, j, startValues[*p1 + *p2 + *p3 + 10 + *n + (K1+1)*2 + (K2+1)*2 + (K3+1)*2 + 2 + *J + j]);
        gsl_vector_set(V3, j, startValues[*p1 + *p2 + *p3 + 10 + *n + (K1+1)*2 + (K2+1)*2 + (K3+1)*2 + 2 + *J*2 + j]);
    }
    
    gsl_matrix *Sigma_V = gsl_matrix_calloc(3,3);
    for(i = 0; i < 3; i++)
    {
        for(j = 0; j < 3; j++)
        {
            gsl_matrix_set(Sigma_V, j, i, startValues[*p1 + *p2 + *p3 + 10 + *n + (K1+1)*2 + (K2+1)*2 + (K3+1)*2 + 2 + *J*3 + i*3 + j]);
        }
    }
    
    gsl_vector *xbeta1 = gsl_vector_calloc(*n);
    gsl_vector *xbeta2 = gsl_vector_calloc(*n);
    gsl_vector *xbeta3 = gsl_vector_calloc(*n);
    gsl_blas_dgemv(CblasNoTrans, 1, survCov1, beta1, 0, xbeta1);
    gsl_blas_dgemv(CblasNoTrans, 1, survCov2, beta2, 0, xbeta2);
    gsl_blas_dgemv(CblasNoTrans, 1, survCov3, beta3, 0, xbeta3);
    
    
    
    
    /* Calculating Sigma_lam (from W and Q) */
    
    
    gsl_matrix *Sigma_lam1       = gsl_matrix_calloc(K1_max+1, K1_max+1);
    gsl_matrix *Sigma_lam2       = gsl_matrix_calloc(K2_max+1, K2_max+1);
    gsl_matrix *Sigma_lam3       = gsl_matrix_calloc(K3_max+1, K3_max+1);
    gsl_matrix *invSigma_lam1    = gsl_matrix_calloc(K1_max+1, K1_max+1);
    gsl_matrix *invSigma_lam2    = gsl_matrix_calloc(K2_max+1, K2_max+1);
    gsl_matrix *invSigma_lam3    = gsl_matrix_calloc(K3_max+1, K3_max+1);
    gsl_matrix *W1               = gsl_matrix_calloc(K1_max+1, K1_max+1);
    gsl_matrix *Q1               = gsl_matrix_calloc(K1_max+1, K1_max+1);
    gsl_matrix *W2               = gsl_matrix_calloc(K2_max+1, K2_max+1);
    gsl_matrix *Q2               = gsl_matrix_calloc(K2_max+1, K2_max+1);
    gsl_matrix *W3               = gsl_matrix_calloc(K3_max+1, K3_max+1);
    gsl_matrix *Q3               = gsl_matrix_calloc(K3_max+1, K3_max+1);
    
    cal_Sigma(Sigma_lam1, invSigma_lam1, W1, Q1, s1, c_lam1, K1);
    cal_Sigma(Sigma_lam2, invSigma_lam2, W2, Q2, s2, c_lam2, K2);
    cal_Sigma(Sigma_lam3, invSigma_lam3, W3, Q3, s3, c_lam3, K3);
    
    
    
    
    /* Variables required for storage of samples */
    
    int StoreInx;
    
    gsl_vector *accept_beta1 = gsl_vector_calloc(nP1);
    gsl_vector *accept_beta2 = gsl_vector_calloc(nP2);
    gsl_vector *accept_beta3 = gsl_vector_calloc(nP3);
    

    gsl_vector *accept_V1       = gsl_vector_calloc(*J);
    gsl_vector *accept_V2       = gsl_vector_calloc(*J);
    gsl_vector *accept_V3       = gsl_vector_calloc(*J);
    
    int accept_BI1      = 0;
    int accept_DI1      = 0;
    int accept_BI2      = 0;
    int accept_DI2      = 0;
    int accept_BI3      = 0;
    int accept_DI3      = 0;
    int accept_theta    = 0;


    

    
    
    
    /* For posterior predictive checks */
    
    gsl_vector *beta1_mean = gsl_vector_calloc(nP1);
    gsl_vector *beta2_mean = gsl_vector_calloc(nP2);
    gsl_vector *beta3_mean = gsl_vector_calloc(nP3);
    
    gsl_vector *xbeta1_mean = gsl_vector_calloc(*n);
    gsl_vector *xbeta2_mean = gsl_vector_calloc(*n);
    gsl_vector *xbeta3_mean = gsl_vector_calloc(*n);
    
    gsl_vector *gamma_mean = gsl_vector_calloc(*n);
    
    gsl_vector *lambda1_mean = gsl_vector_calloc(nTime_lambda1);
    gsl_vector *lambda2_mean = gsl_vector_calloc(nTime_lambda2);
    gsl_vector *lambda3_mean = gsl_vector_calloc(nTime_lambda3);
    
    gsl_vector *lambda1_vec = gsl_vector_calloc(nTime_lambda1);
    gsl_vector *lambda2_vec = gsl_vector_calloc(nTime_lambda2);
    gsl_vector *lambda3_vec = gsl_vector_calloc(nTime_lambda3);
    
    gsl_vector *V1_mean = gsl_vector_calloc(*J);
    gsl_vector *V2_mean = gsl_vector_calloc(*J);
    gsl_vector *V3_mean = gsl_vector_calloc(*J);
    
    double logLH;
    
    gsl_vector *invLH_vec = gsl_vector_calloc(*n);
    gsl_vector *invLH_mean = gsl_vector_calloc(*n);
    gsl_vector *cpo_vec = gsl_vector_calloc(*n);
    
    double invLHval, lpml_temp;
    
    double theta_mean;
    
    
    
    
    
 
    
    /* Compute probabilities for various types of moves (22 moves)*/
    
    double pRP1, pRP2, pRP3, pBH1, pBH2, pBH3, pSP1, pSP2, pSP3, pBI1, pBI2, pBI3, pDI1, pDI2, pDI3, pFP, pDP;
    int move, numUpdate;
    
    gsl_vector *pB1          = gsl_vector_calloc(K1_max);
    gsl_vector *pD1          = gsl_vector_calloc(K1_max);
    gsl_vector *rho_lam1_vec = gsl_vector_calloc(K1_max);
    gsl_vector *pB2          = gsl_vector_calloc(K2_max);
    gsl_vector *pD2          = gsl_vector_calloc(K2_max);
    gsl_vector *rho_lam2_vec = gsl_vector_calloc(K2_max);
    gsl_vector *pB3          = gsl_vector_calloc(K3_max);
    gsl_vector *pD3          = gsl_vector_calloc(K3_max);
    gsl_vector *rho_lam3_vec = gsl_vector_calloc(K3_max);
    double  rho_lam1, rho_lam2, rho_lam3;
    
    for(j = 0; j < K1_max; j++)
    {
        gsl_vector_set(pB1, j, c_min(1, alpha1/(j+1 + 1)));
        gsl_vector_set(pD1, j, c_min(1, (j+1)/alpha1));
    }
    for(j = 0; j < K2_max; j++)
    {
        gsl_vector_set(pB2, j, c_min(1, alpha2/(j+1 + 1)));
        gsl_vector_set(pD2, j, c_min(1, (j+1)/alpha2));
    }
    for(j = 0; j < K3_max; j++)
    {
        gsl_vector_set(pB3, j, c_min(1, alpha3/(j+1 + 1)));
        gsl_vector_set(pD3, j, c_min(1, (j+1)/alpha3));
    }
    
    for(j = 0; j < K1_max; j++) gsl_vector_set(rho_lam1_vec, j, C1/(gsl_vector_get(pB1, j) + gsl_vector_get(pD1, j)));
    for(j = 0; j < K2_max; j++) gsl_vector_set(rho_lam2_vec, j, C2/(gsl_vector_get(pB2, j) + gsl_vector_get(pD2, j)));
    for(j = 0; j < K3_max; j++) gsl_vector_set(rho_lam3_vec, j, C3/(gsl_vector_get(pB3, j) + gsl_vector_get(pD3, j)));
    
    rho_lam1 = gsl_vector_min(rho_lam1_vec);
    rho_lam2 = gsl_vector_min(rho_lam2_vec);
    rho_lam3 = gsl_vector_min(rho_lam3_vec);
    
    
    numUpdate = 19;
    if(*p1 > 0) numUpdate += 1;
    if(*p2 > 0) numUpdate += 1;
    if(*p3 > 0) numUpdate += 1;
    
    
    
    
    
    
    
    for(MM = 0; MM < *numReps; MM++)
    {
        

        
        if(K1 < K1_max)
        {
            pBI1 = rho_lam1 * c_min(1, alpha1/(K1+1));
            pDI1 = rho_lam1 * c_min(1, K1/alpha1);
        }
        if(K2 < K2_max)
        {
            pBI2 = rho_lam2 * c_min(1, alpha2/(K2+1));
            pDI2 = rho_lam2 * c_min(1, K2/alpha2);
        }
        if(K3 < K3_max)
        {
            pBI3 = rho_lam3 * c_min(1, alpha3/(K3+1));
            pDI3 = rho_lam3 * c_min(1, K3/alpha3);
        }
        if(K1 >= K1_max)
        {
            pBI1 = 0;
            pDI1 = rho_lam1 * 2;
        }
        if(K2 >= K2_max)
        {
            pBI2 = 0;
            pDI2 = rho_lam2 * 2;
        }
        if(K3 >= K3_max)
        {
            pBI3 = 0;
            pDI3 = rho_lam3 * 2;
        }
         
         
        
        
        
        
       

        
        double pMP, pCP1, pCP2, pCP3, pVP, choice;        
        
        
        pDP = 0.2;
        pMP = 0;
        
        double probSub = (1 - pMP - pDP - pBI1 - pBI2 - pBI3 - pDI1 - pDI2 - pDI3)/(numUpdate-8);
        
        pRP1 = (*p1 > 0) ? probSub : 0;
        pRP2 = (*p2 > 0) ? probSub : 0;
        pRP3 = (*p3 > 0) ? probSub : 0;
        pBH1 = probSub;
        pBH2 = probSub;
        pBH3 = probSub;
        pSP1 = probSub;
        pSP2 = probSub;
        pSP3 = probSub;
        pFP  = probSub;
        pCP1 = probSub;
        pCP2 = probSub;
        pCP3 = probSub;
        pVP  = probSub;
         
         

        
        /* selecting a move */
        /* move: 1=RP1, 2=RP2, 3=RP3, 4=BH1, 5=BH2, 6=BH3, 7=SP1, 8=SP2, 9=SP3, 10=FP,
         11=DP, 12=BI1, 13=BI2, 14=BI3, 15=DI1, 16=DI2, 17=DI3, 18=CP1, 19=CP2, 20=CP3, 21=VP, 22=MP */
        
        

        
        choice  = runif(0, 1);
        move    = 1;
        if(choice > pRP1) move = 2;
        if(choice > pRP1+pRP2) move = 3;
        if(choice > pRP1+pRP2+pRP3) move = 4;
        if(choice > pRP1+pRP2+pRP3+pBH1) move = 5;
        if(choice > pRP1+pRP2+pRP3+pBH1+pBH2) move = 6;
        if(choice > pRP1+pRP2+pRP3+pBH1+pBH2+pBH3) move = 7;
        if(choice > pRP1+pRP2+pRP3+pBH1+pBH2+pBH3+pSP1) move = 8;
        if(choice > pRP1+pRP2+pRP3+pBH1+pBH2+pBH3+pSP1+pSP2) move = 9;
        if(choice > pRP1+pRP2+pRP3+pBH1+pBH2+pBH3+pSP1+pSP2+pSP3) move = 10;
        if(choice > pRP1+pRP2+pRP3+pBH1+pBH2+pBH3+pSP1+pSP2+pSP3+pFP) move = 11;
        if(choice > pRP1+pRP2+pRP3+pBH1+pBH2+pBH3+pSP1+pSP2+pSP3+pFP+pDP) move = 12;
        if(choice > pRP1+pRP2+pRP3+pBH1+pBH2+pBH3+pSP1+pSP2+pSP3+pFP+pDP+pBI1) move = 13;
        if(choice > pRP1+pRP2+pRP3+pBH1+pBH2+pBH3+pSP1+pSP2+pSP3+pFP+pDP+pBI1+pBI2) move = 14;
        if(choice > pRP1+pRP2+pRP3+pBH1+pBH2+pBH3+pSP1+pSP2+pSP3+pFP+pDP+pBI1+pBI2+pBI3) move = 15;
        if(choice > pRP1+pRP2+pRP3+pBH1+pBH2+pBH3+pSP1+pSP2+pSP3+pFP+pDP+pBI1+pBI2+pBI3+pDI1) move = 16;
        if(choice > pRP1+pRP2+pRP3+pBH1+pBH2+pBH3+pSP1+pSP2+pSP3+pFP+pDP+pBI1+pBI2+pBI3+pDI1+pDI2) move = 17;
        if(choice > pRP1+pRP2+pRP3+pBH1+pBH2+pBH3+pSP1+pSP2+pSP3+pFP+pDP+pBI1+pBI2+pBI3+pDI1+pDI2+pDI3) move = 18;
        if(choice > pRP1+pRP2+pRP3+pBH1+pBH2+pBH3+pSP1+pSP2+pSP3+pFP+pDP+pBI1+pBI2+pBI3+pDI1+pDI2+pDI3+pCP1) move = 19;
        if(choice > pRP1+pRP2+pRP3+pBH1+pBH2+pBH3+pSP1+pSP2+pSP3+pFP+pDP+pBI1+pBI2+pBI3+pDI1+pDI2+pDI3+pCP1+pCP2) move = 20;
        if(choice > pRP1+pRP2+pRP3+pBH1+pBH2+pBH3+pSP1+pSP2+pSP3+pFP+pDP+pBI1+pBI2+pBI3+pDI1+pDI2+pDI3+pCP1+pCP2+pCP3) move = 21;
        if(choice > pRP1+pRP2+pRP3+pBH1+pBH2+pBH3+pSP1+pSP2+pSP3+pFP+pDP+pBI1+pBI2+pBI3+pDI1+pDI2+pDI3+pCP1+pCP2+pCP3+pVP) move = 22;
        
         
         
        moveVec[MM] = (double) move;
        
        
        /* updating regression parameter: beta1         
        
        move = 1;
         
          */
        

        
        if(move == 1)
        {
            BpeMvnCorScrSM_updateRP1(beta1, xbeta1, gamma, V1, lambda1, s1, survTime1, survEvent1, cluster, survCov1, K1, accept_beta1);
        }
        
        
        /* updating regression parameter: beta2
         
         move = 2;
          */
        
        

         
        if(move == 2)
        {
            BpeMvnCorScrSM_updateRP2(beta2, xbeta2, nu2, gamma, V2, lambda2, s2, survTime1, case01, cluster, survCov2, K2, accept_beta2);
        }
 
        
        
        
        /* updating regression parameter: beta3         
         
         move = 3;

        */
        
        
        
        
        
        if(move == 3)
        {
            BpeMvnCorScrSM_updateRP3(beta3, xbeta3, nu3, gamma, V3, lambda3, s3, yStar, case11, cluster, survCov3, K3, accept_beta3);
        }
        
        
        
        
        /* updating log-baseline hazard function parameter: lambda1 
        
        move = 4;*/
        
        


        if(move == 4)
        {
            BpeMvnCorScrSM_updateBH1(lambda1, s1, xbeta1, gamma, V1, survTime1, survEvent1, cluster, Sigma_lam1, invSigma_lam1, W1, Q1, mu_lam1, sigSq_lam1, K1);
        }
        
        
        
        

        /* updating log-baseline hazard function parameter: lambda2 
        
        move = 5;*/
        
        

        
        
        if(move == 5)
        {
            BpeMvnCorScrSM_updateBH2(lambda2, s2, xbeta2, gamma, V2, nu2, survTime1, survTime2, case01, cluster, Sigma_lam2, invSigma_lam2, W2, Q2, mu_lam2, sigSq_lam2, K2);
        }
        
        
        /* updating log-baseline hazard function parameter: lambda3 
        
        move = 6;*/
        
        

        
        if(move == 6)
        {
            BpeMvnCorScrSM_updateBH3(lambda3, s3, xbeta3, gamma, V3, nu3, yStar, case11, cluster, Sigma_lam3, invSigma_lam3, W3, Q3, mu_lam3, sigSq_lam3, K3);
        }
        
        
        
        
        /* updating second stage survival components: mu_lam1 and sigSq_lam1 
         
        move = 7;*/
        
        
        if(move == 7)
        {
            BpeMvnCorScrSM_updateSP1(&mu_lam1, &sigSq_lam1, lambda1, Sigma_lam1, invSigma_lam1, a1, b1, K1);
        }
        
        
        
        /* updating second stage survival components: mu_lam2 and sigSq_lam2         
         
         move = 8;*/


        
        if(move == 8)
        {
            BpeMvnCorScrSM_updateSP2(&mu_lam2, &sigSq_lam2, lambda2, Sigma_lam2, invSigma_lam2, a2, b2, K2);
        }
        
        
        /* updating second stage survival components: mu_lam3 and sigSq_lam3 
        
        move = 9;*/
        

        
        if(move == 9)
        {
            BpeMvnCorScrSM_updateSP3(&mu_lam3, &sigSq_lam3, lambda3, Sigma_lam3, invSigma_lam3, a3, b3, K3);
        }
        
        
        /* updating frailty parameter: gamma
         
         move = 10;*/
        

        
        if(move == 10)
        {
            
            BpeMvnCorScrSM_updateFP_Gibbs(gamma, theta, xbeta1, xbeta2, xbeta3, lambda1, lambda2, lambda3, s1, s2, s3, K1, K2, K3, V1, V2, V3, survTime1, yStar, survEvent1, survEvent2, cluster);
        }
        
        
        
        
        /* updating variance parameter: theta
         
         move = 11;        */
        
        
        
        if(move == 11)
        {
            BpeMvnCorScrSM_updateDP(gamma, &theta, mhProp_theta_var, psi, omega, &accept_theta);
        }
        
        
        /* Updating the number of splits and their positions: K1 and s1 (Birth move) 
        
        move = 12;*/
        
        if(move == 12)
        {
            BpeMvnCorScrSM_updateBI1(s1, &K1, &accept_BI1, survTime1, survEvent1, gamma, xbeta1, V1, cluster, Sigma_lam1, invSigma_lam1, W1, Q1, lambda1, s_propBI1, num_s_propBI1, delPert1, alpha1, c_lam1, mu_lam1, sigSq_lam1, s1_max);
        }
        
        
        /* Updating the number of splits and their positions: K2 and s2 (Birth move)
         
         move = 13;*/
        
        if(move == 13)
        {
            BpeMvnCorScrSM_updateBI2(s2, &K2, &accept_BI2, survTime1, survTime2, case01, gamma, xbeta2, V2, cluster, Sigma_lam2, invSigma_lam2, W2, Q2, lambda2, s_propBI2, num_s_propBI2, delPert2, alpha2, c_lam2, mu_lam2, sigSq_lam2, s2_max);
        }
        
        
        /* Updating the number of splits and their positions: K3 and s3 (Birth move) 
        
        move = 14;*/
        
        if(move == 14)
        {
            BpeMvnCorScrSM_updateBI3(s3, &K3, &accept_BI3, survTime1, yStar, case11, gamma, xbeta3, V3, cluster, Sigma_lam3, invSigma_lam3, W3, Q3, lambda3, s_propBI3, num_s_propBI3, delPert3, alpha3, c_lam3, mu_lam3, sigSq_lam3, s3_max);
        }
        
        
        
        /* Updating the number of splits and their positions: K1 and s1 (Death move) 
        
        move = 15;*/
        
        
        
        if(move == 15)
        {
            BpeMvnCorScrSM_updateDI1(s1, &K1, &accept_DI1, survTime1, survEvent1, gamma, xbeta1, V1, cluster, Sigma_lam1, invSigma_lam1, W1, Q1, lambda1, s_propBI1, num_s_propBI1, delPert1, alpha1, c_lam1, mu_lam1, sigSq_lam1, s1_max, K1_max);
        }
        
        
        
        /* Updating the number of splits and their positions: K2 and s2 (Death move) 
        
        move = 16;*/
        
        
        if(move == 16)
        {
            BpeMvnCorScrSM_updateDI2(s2, &K2, &accept_DI2, survTime1, survTime2, case01, gamma, xbeta2, V2, cluster, Sigma_lam2, invSigma_lam2, W2, Q2, lambda2, s_propBI2, num_s_propBI2, delPert2, alpha2, c_lam2, mu_lam2, sigSq_lam2, s2_max, K2_max);
        }
        
        
        /* Updating the number of splits and their positions: K3 and K3 (Death move) 
        
        move = 17;*/
        
        if(move == 17)
        {
            BpeMvnCorScrSM_updateDI3(s3, &K3, &accept_DI3, survTime1, yStar, case11, gamma, xbeta3, V3, cluster, Sigma_lam3, invSigma_lam3, W3, Q3, lambda3, s_propBI3, num_s_propBI3, delPert3, alpha3, c_lam3, mu_lam3, sigSq_lam3, s3_max, K3_max);
        }
        
  
        /* updating cluster-specific random effect: V1
         
         move = 18;*/
        
        if(move == 18)
        {
            BpeMvnCorScrSM_updateCP1_amcmc(beta1, gamma, V1, V2, V3, lambda1, s1, K1, Sigma_V, survTime1, survEvent1, cluster, survCov1, n_j, accept_V1, mhProp_V1_var);
            
        }
        
        
        /* updating cluster-specific random effect: V2
        
        move = 19;*/
        
        if(move == 19)
        {
            BpeMvnCorScrSM_updateCP2_amcmc(beta2, nu2, gamma, V1, V2, V3, lambda2, s2, K2, Sigma_V, survTime1, case01, cluster, survCov2, n_j, accept_V2, mhProp_V2_var);
            
        }
        
        
        
        /* updating cluster-specific random effect: V3
        
        move = 20;*/
        
        if(move == 20)
        {
            BpeMvnCorScrSM_updateCP3_amcmc(beta3, nu3, gamma, V1, V2, V3, lambda3, s3, K3, Sigma_V, yStar, case11, cluster, survCov3, n_j, accept_V3, mhProp_V3_var);
            
        }
        
        /* updating variance covariance matrix of cluster-specific random effect: Sigma_V 
         
         
         move = 21;*/
        
        if(move == 21)
        {
            BpeMvnCorScrSM_updateVP(V1, V2, V3, Sigma_V, rho_v, Psi_v);
        }
        
        
        
        
        
             
     
        /* Storing posterior samples */
        
        
        if( ( (MM+1) % *thin ) == 0 && (MM+1) > (*numReps * *burninPerc))
        {
            StoreInx = (MM+1)/(*thin)- (*numReps * *burninPerc)/(*thin);
            
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
            
            for(j = 0; j < *J; j++) samples_V1[(StoreInx - 1) * (*J) + j] = gsl_vector_get(V1, j);
            for(j = 0; j < *J; j++) samples_V2[(StoreInx - 1) * (*J) + j] = gsl_vector_get(V2, j);
            for(j = 0; j < *J; j++) samples_V3[(StoreInx - 1) * (*J) + j] = gsl_vector_get(V3, j);
            
            
            for(i = 0; i < 3; i++)
            {
                for(j = 0; j < 3; j++)
                {
                    samples_Sigma_V[(StoreInx - 1) * 9 + i*3 + j] = gsl_matrix_get(Sigma_V, i, j);
                }
            }
            
            
            samples_theta[StoreInx - 1] = theta;
            
            for(i = 0; i < nTime_lambda1; i++)
            {
                j = 0;
                while(gsl_vector_get(time_lambda1, i) > gsl_vector_get(s1, j))
                {
                    j += 1;
                }
                lambda1_fin[(StoreInx - 1) * (nTime_lambda1) + i] = gsl_vector_get(lambda1, j);
                gsl_vector_set(lambda1_vec, i, gsl_vector_get(lambda1, j));
            }
            
            
            for(i = 0; i < nTime_lambda2; i++)
            {
                j = 0;
                while(gsl_vector_get(time_lambda2, i) > gsl_vector_get(s2, j))
                {
                    j += 1;
                }
                lambda2_fin[(StoreInx - 1) * (nTime_lambda2) + i] = gsl_vector_get(lambda2, j);
                gsl_vector_set(lambda2_vec, i, gsl_vector_get(lambda2, j));
            }
            
            
            
            for(i = 0; i < nTime_lambda3; i++)
            {
                j = 0;
                while(gsl_vector_get(time_lambda3, i) > gsl_vector_get(s3, j))
                {
                    j += 1;
                }
                lambda3_fin[(StoreInx - 1) * (nTime_lambda3) + i] = gsl_vector_get(lambda3, j);
                gsl_vector_set(lambda3_vec, i, gsl_vector_get(lambda3, j));
            }
            
            
            
            
            
            
            /* */
            
            
            samples_mu_lam1[StoreInx - 1] = mu_lam1;
            samples_mu_lam2[StoreInx - 1] = mu_lam2;
            samples_mu_lam3[StoreInx - 1] = mu_lam3;
            samples_sigSq_lam1[StoreInx - 1] = sigSq_lam1;
            samples_sigSq_lam2[StoreInx - 1] = sigSq_lam2;
            samples_sigSq_lam3[StoreInx - 1] = sigSq_lam3;
            samples_K1[StoreInx - 1] = K1;
            samples_K2[StoreInx - 1] = K2;
            samples_K3[StoreInx - 1] = K3;
            
            for(j = 0; j < K1+1; j++)
            {
                samples_s1[(StoreInx - 1) * (K1_max+1) + j] = gsl_vector_get(s1, j);
            }
            for(j = 0; j < K2+1; j++)
            {
                samples_s2[(StoreInx - 1) * (K2_max+1) + j] = gsl_vector_get(s2, j);
            }
            for(j = 0; j < K3+1; j++)
            {
                samples_s3[(StoreInx - 1) * (K3_max+1) + j] = gsl_vector_get(s3, j);
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
            
            
            
            /* deviance */
            
            
            BpeMvnCorScrSM_logMLH(beta1, beta2, beta3, xbeta1, xbeta2, xbeta3, theta, lambda1, lambda2, lambda3, s1, s2, s3, V1, V2, V3, survTime1, survTime2, yStar, survEvent1, survEvent2, case01, case11, survCov1, survCov2, survCov3, cluster, K1, K2, K3, &logLH);
            
            dev[StoreInx - 1] = -2*logLH;
            
            gsl_vector_scale(gamma_mean, (double) StoreInx - 1);
            gsl_vector_add(gamma_mean, gamma);
            gsl_vector_scale(gamma_mean, (double)1/StoreInx);
            
            if(*p1 >0)
            {
                gsl_vector_scale(beta1_mean, (double) StoreInx - 1);
                gsl_vector_add(beta1_mean, beta1);
                gsl_vector_scale(beta1_mean, (double)1/StoreInx);
            }
            if(*p2 >0)
            {
                gsl_vector_scale(beta2_mean, (double) StoreInx - 1);
                gsl_vector_add(beta2_mean, beta2);
                gsl_vector_scale(beta2_mean, (double)1/StoreInx);
            }
            if(*p3 >0)
            {
                gsl_vector_scale(beta3_mean, (double) StoreInx - 1);
                gsl_vector_add(beta3_mean, beta3);
                gsl_vector_scale(beta3_mean, (double)1/StoreInx);
            }
            
            gsl_vector_scale(lambda1_mean, (double) StoreInx - 1);
            gsl_vector_add(lambda1_mean, lambda1_vec);
            gsl_vector_scale(lambda1_mean, (double)1/StoreInx);
            
            gsl_vector_scale(lambda2_mean, (double) StoreInx - 1);
            gsl_vector_add(lambda2_mean, lambda2_vec);
            gsl_vector_scale(lambda2_mean, (double)1/StoreInx);
            
            gsl_vector_scale(lambda3_mean, (double) StoreInx - 1);
            gsl_vector_add(lambda3_mean, lambda3_vec);
            gsl_vector_scale(lambda3_mean, (double)1/StoreInx);
            
            
            gsl_vector_scale(V1_mean, (double) StoreInx - 1);
            gsl_vector_add(V1_mean, V1);
            gsl_vector_scale(V1_mean, (double)1/StoreInx);
            
            gsl_vector_scale(V2_mean, (double) StoreInx - 1);
            gsl_vector_add(V2_mean, V2);
            gsl_vector_scale(V2_mean, (double)1/StoreInx);
            
            gsl_vector_scale(V3_mean, (double) StoreInx - 1);
            gsl_vector_add(V3_mean, V3);
            gsl_vector_scale(V3_mean, (double)1/StoreInx);
            
            
            theta_mean = ((double) StoreInx - 1) * theta_mean;
            theta_mean = theta_mean + theta;
            theta_mean = theta_mean/((double) StoreInx);
            
            
            
            
            /* CPO_i */
            
            for(i = 0; i < *n; i++)
            {
                BpeMvnCorScrSM_logMLH_i(i, beta1, beta2, beta3, xbeta1, xbeta2, xbeta3, theta, lambda1, lambda2, lambda3, s1, s2, s3, V1, V2, V3, survTime1, survTime2, yStar, survEvent1, survEvent2, case01, case11, survCov1, survCov2, survCov3, cluster, K1, K2, K3, &invLHval);
                
                invLHval = 1/exp(invLHval);
                
                gsl_vector_set(invLH_vec, i, invLHval);
                
            }
            
            
            
            gsl_vector_scale(invLH_mean, (double) StoreInx - 1);
            
            
            
            gsl_vector_add(invLH_mean, invLH_vec);
            
            
            
            gsl_vector_scale(invLH_mean, (double)1/StoreInx);
            
            
            
            
            for(i = 0; i < *n; i++)
            {
                gsl_vector_set(cpo_vec, i, 1/gsl_vector_get(invLH_mean, i));
            }
            
            lpml_temp = 0;
            
            for(i = 0; i < *n; i++)
            {
                lpml_temp += log(gsl_vector_get(cpo_vec, i));
            }
            
            lpml[StoreInx - 1] = lpml_temp;
            
            
        }
        
        
        
        
        
        if(MM == (*numReps - 1))
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
            
            samples_misc[*p1 + *p2 + *p3] = accept_BI1;
            samples_misc[*p1 + *p2 + *p3 + 1] = accept_DI1;
            samples_misc[*p1 + *p2 + *p3 + 2] = accept_BI2;
            samples_misc[*p1 + *p2 + *p3 + 3] = accept_DI2;
            samples_misc[*p1 + *p2 + *p3 + 4] = accept_BI3;
            samples_misc[*p1 + *p2 + *p3 + 5] = accept_DI3;
            samples_misc[*p1 + *p2 + *p3 + 6] = accept_theta;
            
            for(i = 0; i < *J; i++) samples_misc[*p1 + *p2 + *p3 + 10 + *n + *n + i] = (int) gsl_vector_get(accept_V1, i);
            for(i = 0; i < *J; i++) samples_misc[*p1 + *p2 + *p3 + 10 + *n + *n + *J + i] = (int) gsl_vector_get(accept_V2, i);
            for(i = 0; i < *J; i++) samples_misc[*p1 + *p2 + *p3 + 10 + *n + *n + *J + *J + i] = (int) gsl_vector_get(accept_V3, i);
            
            for(i = 0; i < *n; i++)
            {
                samples_gamma_last[i] = gsl_vector_get(gamma, i);
            }
            
            for(i = 0; i < *n; i++)
            {
                gammaP[i] = gsl_vector_get(gamma_mean, i);
            }
            
            for(i = 0; i < *n; i++)
            {
                invLH[i] = gsl_vector_get(invLH_mean, i);
            }
            
            
            gsl_blas_dgemv(CblasNoTrans, 1, survCov1, beta1_mean, 0, xbeta1_mean);
            gsl_blas_dgemv(CblasNoTrans, 1, survCov2, beta2_mean, 0, xbeta2_mean);
            gsl_blas_dgemv(CblasNoTrans, 1, survCov3, beta3_mean, 0, xbeta3_mean);
            
            BpeMvnCorScrSM_logMLH(beta1_mean, beta2_mean, beta3_mean, xbeta1_mean, xbeta2_mean, xbeta3_mean, theta_mean, lambda1_mean, lambda2_mean, lambda3_mean, time_lambda1, time_lambda2, time_lambda3, V1_mean, V2_mean, V3_mean, survTime1, survTime2, yStar, survEvent1, survEvent2, case01, case11, survCov1, survCov2, survCov3, cluster, nTime_lambda1-1, nTime_lambda2-1, nTime_lambda3-1, logLH_fin);
            
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





















