/*
 TO COMPILE USE THE CODE:
 
 R CMD SHLIB BweibMvnCorScrSM.c BweibMvnCorScrSM_Updates.c BweibMvnCorScrSM_Utilities.c -lgsl -lgslcblas
 
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

#include "BweibMvnCorScrSM.h"





/* */
void BweibMvnCorScrSMmcmc(double survData[],
                  int *n,
                  int *p1,
                  int *p2,
                  int *p3,
                  int *J,
                  double nj[],
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
                  double samples_nu2[],
                  double samples_nu3[],
                  double samples_theta[],
                  double samples_V1[],
                  double samples_V2[],
                  double samples_V3[],
                  double samples_Sigma_V[],
                  double samples_gamma[],
                  double samples_misc[],                          
                          double gammaP[],
                          double dev[],
                          double invLH[],
                          double *logLH_fin,
                          double lpml[],
                          double lpml2[],
                          double moveVec[])
{
    GetRNGstate();
    
    time_t now;    
    
    int i, j, MM;
    int lastChgProp = 0;

    
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
    double rho_v    = hyperParams[14];

    

    gsl_matrix *Psi_v = gsl_matrix_calloc(3,3);
    for(i = 0; i < 3; i++)
    {
        for(j = 0; j < 3; j++)
        {
            gsl_matrix_set(Psi_v, j, i, hyperParams[15 + i*3 + j]);
        }
    }
    

    

    /* varialbes for M-H step */
    
    double mhProp_alpha1_var = mcmcParams[0];
    double mhProp_alpha2_var = mcmcParams[1];
    double mhProp_alpha3_var = mcmcParams[2];
    double mhProp_theta_var  = mcmcParams[3];

    double mhProp_V1_var  = mcmcParams[5];
    double mhProp_V2_var  = mcmcParams[6];
    double mhProp_V3_var  = mcmcParams[7];

    
    gsl_vector *mhGam_chk = gsl_vector_calloc(*n); /* 0=RW, 1=IP */
    
    
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
    
    double nu2 = startValues[*p1 + *p2 + *p3 + *n + 7];
    double nu3 = startValues[*p1 + *p2 + *p3 + *n + 8];
    
    gsl_vector *V1 = gsl_vector_calloc(*J);
    gsl_vector *V2 = gsl_vector_calloc(*J);
    gsl_vector *V3 = gsl_vector_calloc(*J);
    
    for(j = 0; j < *J; j++)
    {
        gsl_vector_set(V1, j, startValues[*p1 + *p2 + *p3 + *n + 9 + j]);
        gsl_vector_set(V2, j, startValues[*p1 + *p2 + *p3 + *n + 9 + *J + j]);
        gsl_vector_set(V3, j, startValues[*p1 + *p2 + *p3 + *n + 9 + *J*2 + j]);
    }
    
    gsl_matrix *Sigma_V = gsl_matrix_calloc(3,3);
    for(i = 0; i < 3; i++)
    {
        for(j = 0; j < 3; j++)
        {
            gsl_matrix_set(Sigma_V, j, i, startValues[*p1 + *p2 + *p3 + *n + 9 + *J*3 + i*3 + j]);
        }
    }
    
   
    
    

    /* Variables required for storage of samples */
    
    int StoreInx;
    
    gsl_vector *accept_beta1 = gsl_vector_calloc(nP1);
    gsl_vector *accept_beta2 = gsl_vector_calloc(nP2);
    gsl_vector *accept_beta3 = gsl_vector_calloc(nP3);
    
    gsl_vector *accept_gamma    = gsl_vector_calloc(*n);
    gsl_vector *accept_V1       = gsl_vector_calloc(*J);
    gsl_vector *accept_V2       = gsl_vector_calloc(*J);
    gsl_vector *accept_V3       = gsl_vector_calloc(*J);    
    
    int accept_alpha1 = 0;
    int accept_alpha2 = 0;
    int accept_alpha3 = 0;
    int accept_theta  = 0;
    int accept_nu2  = 0;
    int accept_nu3  = 0;
    


    /* For posterior predictive checks */
    
    gsl_vector *beta1_mean = gsl_vector_calloc(nP1);
    gsl_vector *beta2_mean = gsl_vector_calloc(nP2);
    gsl_vector *beta3_mean = gsl_vector_calloc(nP3);
    



    
    gsl_vector *gamma_mean = gsl_vector_calloc(*n);
    
    gsl_vector *V1_mean = gsl_vector_calloc(*J);
    gsl_vector *V2_mean = gsl_vector_calloc(*J);
    gsl_vector *V3_mean = gsl_vector_calloc(*J);
    
    double logLH;
    
    gsl_vector *invLH_vec = gsl_vector_calloc(*n);
    gsl_vector *invLH_mean = gsl_vector_calloc(*n);
    gsl_vector *cpo_vec = gsl_vector_calloc(*n);
    



    
    double invLHval, lpml_temp;
    
    double alpha1_mean=0;
    double alpha2_mean=0;
    double alpha3_mean=0;
    double kappa1_mean=0;
    double kappa2_mean=0;
    double kappa3_mean=0;
    double theta_mean=0;
    
    
    /* Compute probabilities for various types of moves */
    
    double pRP1, pRP2, pRP3, pSH1, pSH2, pSH3, pSC1, pSC2, pSC3, pDP, pFP, pMP, pCP1, pCP2, pCP3, pVP, choice;
    int move, numUpdate;
    
    
    
    
    numUpdate = 13;
    if(*p1 > 0) numUpdate += 1;
    if(*p2 > 0) numUpdate += 1;
    if(*p3 > 0) numUpdate += 1;
    
    pDP  = (double) 0.15;
    pSC1 = (double) 0.15;
    pSC2 = (double) 0.15;
    pSC3 = (double) 0.15;
    pMP = (double) 0;      
    
    
    double probSub = (1 - pDP - pSC1 - pSC2 - pSC3 -pMP)/(numUpdate-5);
    
    
    pRP1 = (*p1 > 0) ? probSub : 0;
    pRP2 = (*p2 > 0) ? probSub : 0;
    pRP3 = (*p3 > 0) ? probSub : 0;
    pSH1 = probSub;
    pSH2 = probSub;
    pSH3 = probSub;
    pFP  = probSub;
    /*pMP  = probSub;*/
    pCP1 = probSub;
    pCP2 = probSub;
    pCP3 = probSub;
    pVP  = 1-(pRP1 + pRP2 + pRP3 + pSH1 + pSH2 + pSH3 + pSC1 + pSC2 + pSC3 + pDP + pFP + pMP + pCP1 + pCP2 + pCP3);
    
    
   
    
        
    for(MM = 0; MM < *numReps; MM++)
    {
        /* selecting a move */
        /* move: 1=RP1, 2=RP2, 3=RP3, 4=SH1, 5=SH2, 6=SH3 */
        /* move: 7=SC1, 8=SC2, 9=SC3, 10=DP, 11=FP 12=MP*/
        /* move: 13=CP1, 14=CP2, 15=CP3, 16=VP*/
        
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
        if(choice > pRP1 + pRP2 + pRP3 + pSH1 + pSH2 + pSH3 + pSC1 + pSC2 + pSC3 + pDP + pFP) move = 12;
        if(choice > pRP1 + pRP2 + pRP3 + pSH1 + pSH2 + pSH3 + pSC1 + pSC2 + pSC3 + pDP + pFP + pMP) move = 13;
        if(choice > pRP1 + pRP2 + pRP3 + pSH1 + pSH2 + pSH3 + pSC1 + pSC2 + pSC3 + pDP + pFP + pMP + pCP1) move = 14;
        if(choice > pRP1 + pRP2 + pRP3 + pSH1 + pSH2 + pSH3 + pSC1 + pSC2 + pSC3 + pDP + pFP + pMP + pCP1 + pCP2) move = 15;
        if(choice > pRP1 + pRP2 + pRP3 + pSH1 + pSH2 + pSH3 + pSC1 + pSC2 + pSC3 + pDP + pFP + pMP + pCP1 + pCP2 + pCP3) move = 16;
        
        
        moveVec[MM] = (double) move;
        
        /* updating regression parameter: beta1 
        
        move = 100;*/  
        
        if(move == 1)
        {
            BweibMvnCorScrSM_updateRP1(beta1, &alpha1, &kappa1, gamma, V1, survTime1, survEvent1, cluster, survCov1, accept_beta1);
        }
        
        
        
        /* updating regression parameter: beta2  
        
        move = 200; */  
        
        if(move == 2)
        {
            BweibMvnCorScrSM_updateRP2(beta2, &alpha2, &kappa2, &nu2, gamma, V2, survTime1, case01, cluster, survCov2, accept_beta2);
        }

         
        
        /* updating regression parameter: beta3 
        
        move = 30;  */  
        
        if(move == 3)
        {
            BweibMvnCorScrSM_updateRP3(beta3, &alpha3, &kappa3, &nu3, gamma, V3, yStar, case11, cluster, survCov3, accept_beta3);
        }

        
        /* updating shape parameter: alpha1
        
        move = 40;       */  
        
        if(move == 4)
        {
            BweibMvnCorScrSM_updateSC1_rw2(beta1, &alpha1, &kappa1, gamma, V1, survTime1, survEvent1, cluster, survCov1, mhProp_alpha1_var, a1, b1, &accept_alpha1);
        }

         
         
        
        /* updating shape parameter: alpha2
        
        move = 50;       */  
        
        if(move == 5)
        {
            BweibMvnCorScrSM_updateSC2_rw2(beta2, &alpha2, &kappa2, &nu2, gamma, V2, survTime1, survTime2, case01, cluster, survCov2, mhProp_alpha2_var, a2, b2, &accept_alpha2);
        }
       
  


        /* updating shape parameter: alpha3  
        
        move = 60;        */       
    
        
        if(move == 6)
        {
            BweibMvnCorScrSM_updateSC3_rw2(beta3, &alpha3, &kappa3, &nu3, gamma, V3, yStar, case11, cluster, survCov3, mhProp_alpha3_var, a3, b3, &accept_alpha3);
        }
  

         

        /* updating shape parameter: kappa1 
  
        move = 70;     */
        
        if(move == 7)
        {
            BweibMvnCorScrSM_updateSH1(beta1, &alpha1, &kappa1, gamma, V1, survTime1, survEvent1, cluster, survCov1, c1, d1);
        }

  
        
         
        /* updating shape parameter: kappa2 
        
        move = 80;         */        
        
        if(move == 8)
        {
            BweibMvnCorScrSM_updateSH2(beta2, &alpha2, &kappa2, &nu2, gamma, V2, survTime1, case01, cluster, survCov2, c2, d2);
        }
         


        
        /* updating shape parameter: kappa3 
        
        move = 90;        */ 
        
        if(move == 9)
        {
            BweibMvnCorScrSM_updateSH3(beta3, &alpha3, &kappa3, &nu3, gamma, V3, yStar, case11, cluster, survCov3, c3, d3);
        }


    
        
        
        
        
        /* updating variance parameter: theta
        
        move = 100;        */  
        
        if(move == 10)
        {
            BweibMvnCorScrSM_updateDP(gamma, &theta, mhProp_theta_var, psi, omega, &accept_theta);
        }
         


        
        
        /* updating frailty parameter: gamma      
        
        move = 110;*/           
         
        if(move == 11)
        {

             BweibMvnCorScrSM_updateFP_Gibbs(beta1, beta2, beta3, alpha1, alpha2, alpha3, kappa1, kappa2, kappa3, theta, gamma, V1, V2, V3, survTime1, yStar, survEvent1, survEvent2, cluster, survCov1, survCov2, survCov3);

        }

        
        



        
        
        /* updating nu2 and nu3 
        
        move = 120;   */  
        
        if(move == 12000)
        {
            BweibMvnCorScrSM_updateMP(beta2, beta3, &alpha2, &alpha3, &kappa2, &kappa3, &nu2, &nu3, gamma, V2, V3, survTime1, survTime2, case01, case11, cluster, survCov2, survCov3, &accept_nu2, &accept_nu3);
        }
        

        
        
        
        /* updating cluster-specific random effect: V1   
         
         move = 130;*/
    
        if(move == 13)
        {
            BweibMvnCorScrSM_updateCP1_amcmc(beta1, alpha1, kappa1, gamma, V1, V2, V3, Sigma_V, survTime1, survEvent1, cluster, survCov1, n_j, accept_V1, mhProp_V1_var);
 
        }
        
        
        
        /* updating cluster-specific random effect: V2 
        
         move = 140; */
        
        if(move == 14)
        {
            BweibMvnCorScrSM_updateCP2_amcmc(beta2, alpha2, kappa2, nu2, gamma, V1, V2, V3, Sigma_V, survTime1, case01, cluster, survCov2, n_j, accept_V2, mhProp_V2_var);
        }
        

        /* updating cluster-specific random effect: V3 
        
        
        move = 150;*/
        
        if(move == 15)
        {
            BweibMvnCorScrSM_updateCP3_amcmc(beta3, alpha3, kappa3, nu3, gamma, V1, V2, V3, Sigma_V, yStar, case11, cluster, survCov3, n_j, accept_V3, mhProp_V3_var);
        }

        
        
        /* updating variance covariance matrix of cluster-specific random effect: Sigma_V x
        
        
        move = 160;*/
        
        if(move == 16)
        {
             BweibMvnCorScrSM_updateVP(V1, V2, V3, Sigma_V, rho_v, Psi_v);   
        }
        
        
        
        
        
        /*        */
       


        
        
        /* Storing posterior samples */
        
        
        if( ( (MM+1) % *thin ) == 0 && (MM+1) > (*numReps * *burninPerc))
        {
            StoreInx = (MM+1)/(*thin)- (*numReps * *burninPerc)/(*thin);
 
            samples_alpha1[StoreInx - 1] = alpha1;
            samples_alpha2[StoreInx - 1] = alpha2;
            samples_alpha3[StoreInx - 1] = alpha3;
            
            samples_kappa1[StoreInx - 1] = kappa1;
            samples_kappa2[StoreInx - 1] = kappa2;
            samples_kappa3[StoreInx - 1] = kappa3;
            
            samples_nu2[StoreInx - 1] = nu2;
            samples_nu3[StoreInx - 1] = nu3;
            
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
            
            
            BweibMvnCorScrSM_logMLH(beta1, beta2, beta3, alpha1, alpha2, alpha3, kappa1, kappa2, kappa3, theta, V1, V2, V3, survTime1, survTime2, yStar, survEvent1, survEvent2, case01, case11, survCov1, survCov2, survCov3, cluster, &logLH);
            
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
            
            gsl_vector_scale(V1_mean, (double) StoreInx - 1);
            gsl_vector_add(V1_mean, V1);
            gsl_vector_scale(V1_mean, (double)1/StoreInx);
            
            gsl_vector_scale(V2_mean, (double) StoreInx - 1);
            gsl_vector_add(V2_mean, V2);
            gsl_vector_scale(V2_mean, (double)1/StoreInx);
            
            gsl_vector_scale(V3_mean, (double) StoreInx - 1);
            gsl_vector_add(V3_mean, V3);
            gsl_vector_scale(V3_mean, (double)1/StoreInx);
            
            
            alpha1_mean = ((double) StoreInx - 1) * alpha1_mean;
            alpha1_mean = alpha1_mean + alpha1;
            alpha1_mean = alpha1_mean/((double) StoreInx);
            
            alpha2_mean = ((double) StoreInx - 1) * alpha2_mean;
            alpha2_mean = alpha2_mean + alpha2;
            alpha2_mean = alpha2_mean/((double) StoreInx);
            
            alpha3_mean = ((double) StoreInx - 1) * alpha3_mean;
            alpha3_mean = alpha3_mean + alpha3;
            alpha3_mean = alpha3_mean/((double) StoreInx);
            
            kappa1_mean = ((double) StoreInx - 1) * kappa1_mean;
            kappa1_mean = kappa1_mean + kappa1;
            kappa1_mean = kappa1_mean/((double) StoreInx);
            
            kappa2_mean = ((double) StoreInx - 1) * kappa2_mean;
            kappa2_mean = kappa2_mean + kappa2;
            kappa2_mean = kappa2_mean/((double) StoreInx);
            
            kappa3_mean = ((double) StoreInx - 1) * kappa3_mean;
            kappa3_mean = kappa3_mean + kappa3;
            kappa3_mean = kappa3_mean/((double) StoreInx);
            
            theta_mean = ((double) StoreInx - 1) * theta_mean;
            theta_mean = theta_mean + theta;
            theta_mean = theta_mean/((double) StoreInx);
            
            
            /* CPO_i */
            
            for(i = 0; i < *n; i++)
            {
                BweibMvnCorScrSM_logMLH_i(i, beta1, beta2, beta3, alpha1, alpha2, alpha3, kappa1, kappa2, kappa3, theta, V1, V2, V3, survTime1, survTime2, yStar, survEvent1, survEvent2, case01, case11, survCov1, survCov2, survCov3, cluster, &invLHval);
                
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
                samples_misc[*p1 + *p2 + *p3]       = accept_alpha1;
                samples_misc[*p1 + *p2 + *p3 + 1]   = accept_alpha2;
                samples_misc[*p1 + *p2 + *p3 + 2]   = accept_alpha3;
                samples_misc[*p1 + *p2 + *p3 + 3]   = accept_theta;
                samples_misc[*p1 + *p2 + *p3 + 4]   = accept_nu2;
                samples_misc[*p1 + *p2 + *p3 + 5]   = accept_nu3;
                
                for(j = 0; j < *n; j++) samples_misc[*p1 + *p2 + *p3 + 6 + j] = (int) gsl_vector_get(accept_gamma, j);
                
                samples_misc[*p1 + *p2 + *p3 + 5 + *n + 1]   = lastChgProp;
                
                for(i = 0; i < *n; i++) samples_misc[*p1 + *p2 + *p3 + 7 + *n + i] = (int) gsl_vector_get(mhGam_chk, i);
                
                for(i = 0; i < *J; i++) samples_misc[*p1 + *p2 + *p3 + 7 + *n + *n + i] = (int) gsl_vector_get(accept_V1, i);
                for(i = 0; i < *J; i++) samples_misc[*p1 + *p2 + *p3 + 7 + *n + *n + *J + i] = (int) gsl_vector_get(accept_V2, i);
                for(i = 0; i < *J; i++) samples_misc[*p1 + *p2 + *p3 + 7 + *n + *n + *J + *J + i] = (int) gsl_vector_get(accept_V3, i);
            
            for(i = 0; i < *n; i++)
            {
                gammaP[i] = gsl_vector_get(gamma_mean, i);
            }
            
            for(i = 0; i < *n; i++)
            {
                invLH[i] = gsl_vector_get(invLH_mean, i);
            }
            
            
            BweibMvnCorScrSM_logMLH(beta1_mean, beta2_mean, beta3_mean, alpha1_mean, alpha2_mean, alpha3_mean, kappa1_mean, kappa2_mean, kappa3_mean, theta_mean, V1_mean, V2_mean, V3_mean, survTime1, survTime2, yStar, survEvent1, survEvent2, case01, case11, survCov1, survCov2, survCov3, cluster, logLH_fin);
            
            
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





















