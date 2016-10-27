


#include <stdio.h>
#include <math.h>

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_sort_vector.h"
#include "gsl/gsl_heapsort.h"
#include "gsl/gsl_sf.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

#include "R.h"
#include "Rmath.h"

#include "BpeScr.h"
#include "BpeScrSM.h"
#include "BpeSurv.h"
#include "BweibScr.h"
#include "BweibScrSM.h"
#include "BweibSurv.h"
#include "BweibMvnCorScr.h"
#include "BweibMvnCorScrSM.h"
#include "BweibDpCorScr.h"
#include "BweibDpCorScrSM.h"
#include "BpeMvnCorScr.h"
#include "BpeMvnCorScrSM.h"
#include "BpeDpCorScr.h"
#include "BpeDpCorScrSM.h"
#include "BweibCorSurv.h"
#include "BweibDpCorSurv.h"
#include "BpeMvnCorSurv.h"
#include "BpeDpCorSurv.h"

#include "BAFT_DPscr.h"
#include "BAFT_LNscr.h"
#include "BAFT_DPsurv.h"
#include "BAFT_LNsurv.h"




#define Pi 3.141592653589793238462643383280



/********* For PEM-DPM (univariate) model *************/


/* evaluating log-likelihood function */

/**/

void BpeDpCorSurv_logLH(gsl_vector *beta,
                        gsl_vector *xbeta,
                        gsl_vector *lambda,
                        gsl_vector *s,
                        gsl_vector *V,
                        gsl_vector *survTime,
                        gsl_vector *survEvent,
                        gsl_matrix *survCov,
                        gsl_vector *cluster,
                        int K,
                        double *val)
{
    double logLH = 0;
    

    int n = survTime -> size;
    
    int i, j, jj;
    
    double Del, cumHaz;
    
    for(i = 0; i < n; i++)
    {
        cumHaz = 0;
        
        jj = (int) gsl_vector_get(cluster, i) - 1;
        
        if(gsl_vector_get(survEvent, i) == 1)
        {
            for(j = 0; j < K+1; j++)
            {
                if(j == 0 && gsl_vector_get(survTime, i) <= gsl_vector_get(s, 0))
                {
                    logLH += gsl_vector_get(lambda, j);
                }
                if(j != 0 && gsl_vector_get(survTime, i) > gsl_vector_get(s, j-1) && gsl_vector_get(survTime, i) <= gsl_vector_get(s, j))
                {
                    logLH += gsl_vector_get(lambda, j);
                }
            }
            logLH += gsl_vector_get(xbeta, i);
            logLH += gsl_vector_get(V, jj);
        }
        
        for(j = 0; j < K+1; j++)
        {
            if(j > 0)
            {
                Del = c_max(0, (c_min(gsl_vector_get(s, j), gsl_vector_get(survTime, i)) - gsl_vector_get(s, j-1)));
            }
            if(j == 0)
            {
                Del = c_max(0, c_min(gsl_vector_get(s, j), gsl_vector_get(survTime, i)) - 0);
            }
            cumHaz += Del* exp(gsl_vector_get(lambda, j));
        }
        
        cumHaz *= exp(gsl_vector_get(xbeta, i)+gsl_vector_get(V, jj));
        
        
        logLH += -cumHaz;
    }
    
    *val = logLH;
    
    return;
    
}





/********* For PEM-Normal (univariate) model *************/


/* evaluating log-likelihood function */

/**/

void BpeMvnCorSurv_logLH(gsl_vector *beta,
                         gsl_vector *xbeta,
                         gsl_vector *lambda,
                         gsl_vector *s,
                         gsl_vector *V,
                         gsl_vector *survTime,
                         gsl_vector *survEvent,
                         gsl_matrix *survCov,
                         gsl_vector *cluster,
                         int K,
                         double *val)
{
    double logLH = 0;
    

    int n = survTime -> size;
    
    int i, j, jj;
    
    double Del, cumHaz;
    
    for(i = 0; i < n; i++)
    {
        cumHaz = 0;
        
        jj = (int) gsl_vector_get(cluster, i) - 1;
        
        if(gsl_vector_get(survEvent, i) == 1)
        {
            for(j = 0; j < K+1; j++)
            {
                if(j == 0 && gsl_vector_get(survTime, i) <= gsl_vector_get(s, 0))
                {
                    logLH += gsl_vector_get(lambda, j);
                }
                if(j != 0 && gsl_vector_get(survTime, i) > gsl_vector_get(s, j-1) && gsl_vector_get(survTime, i) <= gsl_vector_get(s, j))
                {
                    logLH += gsl_vector_get(lambda, j);
                }
            }
            logLH += gsl_vector_get(xbeta, i);
            logLH += gsl_vector_get(V, jj);
        }
        
        for(j = 0; j < K+1; j++)
        {
            if(j > 0)
            {
                Del = c_max(0, (c_min(gsl_vector_get(s, j), gsl_vector_get(survTime, i)) - gsl_vector_get(s, j-1)));
            }
            if(j == 0)
            {
                Del = c_max(0, c_min(gsl_vector_get(s, j), gsl_vector_get(survTime, i)) - 0);
            }
            cumHaz += Del* exp(gsl_vector_get(lambda, j));
        }
        
        cumHaz *= exp(gsl_vector_get(xbeta, i)+gsl_vector_get(V, jj));
        
        
        logLH += -cumHaz;
    }
    
    *val = logLH;
    
    return;
    
}




/********* For Weibull-DPM (univariate) model *************/

double Qfunc_univ(double V,
             double mu0,
             double zeta0,
             double a0,
             double b0)
{
    double lterm1, term2, term3, val;
    
    lterm1 = gsl_sf_lngamma(a0+0.5) - gsl_sf_lngamma(a0);
    term2 = sqrt(zeta0/(2*Pi*b0*(zeta0 + 1)));
    term3 = pow(zeta0*pow(V - mu0, 2)/(2*b0*(zeta0 + 1)) + 1, -a0-0.5);
    
    val = exp(lterm1)*term2*term3;
    
    return val;
}





/********* For Weibull-Normal (univariate) model *************/

/********* For PEM-DPM-SM model *************/

/* evaluating log-likelihood function */

/**/

void BpeDpCorScrSM_logMLH(gsl_vector *beta1,
                          gsl_vector *beta2,
                          gsl_vector *beta3,
                          gsl_vector *xbeta1,
                          gsl_vector *xbeta2,
                          gsl_vector *xbeta3,
                          double theta,
                          gsl_vector *lambda1,
                          gsl_vector *lambda2,
                          gsl_vector *lambda3,
                          gsl_vector *s1,
                          gsl_vector *s2,
                          gsl_vector *s3,
                          gsl_vector *V1,
                          gsl_vector *V2,
                          gsl_vector *V3,
                          gsl_vector *survTime1,
                          gsl_vector *survTime2,
                          gsl_vector *yStar,
                          gsl_vector *survEvent1,
                          gsl_vector *survEvent2,
                          gsl_vector *case01,
                          gsl_vector *case11,
                          gsl_matrix *survCov1,
                          gsl_matrix *survCov2,
                          gsl_matrix *survCov3,
                          gsl_vector *cluster,
                          int K1,
                          int K2,
                          int K3,
                          double *val)
{
    double gfunc;
    double logLH = 0;
    

    int n = survTime1 -> size;
    
    int i, j, jj;
    
    for(i = 0; i < n; i++)
    {
        
        
        jj = (int) gsl_vector_get(cluster, i) - 1;
        
        if(gsl_vector_get(survEvent1, i) == 1)
        {
            for(j = 0; j < K1+1; j++)
            {
                if(j == 0 && gsl_vector_get(survTime1, i) <= gsl_vector_get(s1, 0))
                {
                    logLH += gsl_vector_get(lambda1, j);
                }
                if(j != 0 && gsl_vector_get(survTime1, i) > gsl_vector_get(s1, j-1) && gsl_vector_get(survTime1, i) <= gsl_vector_get(s1, j))
                {
                    logLH += gsl_vector_get(lambda1, j);
                }
            }
            
            logLH += gsl_vector_get(xbeta1, i);
            logLH += gsl_vector_get(V1, jj);
        }
        
        if(gsl_vector_get(case01, i) ==  1)
        {
            for(j = 0; j < K2+1; j++)
            {
                if(j == 0 && gsl_vector_get(survTime2, i) <= gsl_vector_get(s2, 0))
                {
                    logLH += gsl_vector_get(lambda2, j);
                }
                if(j != 0 && gsl_vector_get(survTime2, i) > gsl_vector_get(s2, j-1) && gsl_vector_get(survTime2, i) <= gsl_vector_get(s2, j))
                {
                    logLH += gsl_vector_get(lambda2, j);
                }
            }
            
            logLH += gsl_vector_get(xbeta2, i);
            logLH += gsl_vector_get(V2, jj);
        }
        
        if(gsl_vector_get(case11, i) == 1)
        {
            for(j = 0; j < K3+1; j++)
            {
                if(j == 0 && gsl_vector_get(survTime2, i) <= gsl_vector_get(s3, 0))
                {
                    logLH += gsl_vector_get(lambda3, j);
                }
                if(j != 0 && gsl_vector_get(survTime2, i) > gsl_vector_get(s3, j-1) && gsl_vector_get(survTime2, i) <= gsl_vector_get(s3, j))
                {
                    logLH += gsl_vector_get(lambda3, j);
                }
            }
            
            logLH += gsl_vector_get(xbeta3, i);
            logLH += gsl_vector_get(V3, jj);
        }
        
        gfunc = BpeDpCorScrSM_wFunc(i, xbeta1, xbeta2, xbeta3, lambda1, lambda2, lambda3, jj, V1, V2, V3, s1, s2, s3, K1, K2, K3, survTime1, survTime2);
        
        logLH += (- 1/theta - gsl_vector_get(survEvent1, i) - gsl_vector_get(survEvent2, i)) * log(1 + theta * gfunc);
    }
    
    *val = logLH;
    
    return;
    
}






/* evaluating log-likelihood function for subject i */

/**/

void BpeDpCorScrSM_logMLH_i(int i,
                            gsl_vector *beta1,
                            gsl_vector *beta2,
                            gsl_vector *beta3,
                            gsl_vector *xbeta1,
                            gsl_vector *xbeta2,
                            gsl_vector *xbeta3,
                            double theta,
                            gsl_vector *lambda1,
                            gsl_vector *lambda2,
                            gsl_vector *lambda3,
                            gsl_vector *s1,
                            gsl_vector *s2,
                            gsl_vector *s3,
                            gsl_vector *V1,
                            gsl_vector *V2,
                            gsl_vector *V3,
                            gsl_vector *survTime1,
                            gsl_vector *survTime2,
                            gsl_vector *yStar,
                            gsl_vector *survEvent1,
                            gsl_vector *survEvent2,
                            gsl_vector *case01,
                            gsl_vector *case11,
                            gsl_matrix *survCov1,
                            gsl_matrix *survCov2,
                            gsl_matrix *survCov3,
                            gsl_vector *cluster,
                            int K1,
                            int K2,
                            int K3,
                            double *val)
{
    double gfunc;
    double logLH = 0;
    


    
    int j, jj;
    
    
    
    jj = (int) gsl_vector_get(cluster, i) - 1;
    
    if(gsl_vector_get(survEvent1, i) == 1)
    {
        for(j = 0; j < K1+1; j++)
        {
            if(j == 0 && gsl_vector_get(survTime1, i) <= gsl_vector_get(s1, 0))
            {
                logLH += gsl_vector_get(lambda1, j);
            }
            if(j != 0 && gsl_vector_get(survTime1, i) > gsl_vector_get(s1, j-1) && gsl_vector_get(survTime1, i) <= gsl_vector_get(s1, j))
            {
                logLH += gsl_vector_get(lambda1, j);
            }
        }
        
        logLH += gsl_vector_get(xbeta1, i);
        logLH += gsl_vector_get(V1, jj);
    }
    
    if(gsl_vector_get(case01, i) ==  1)
    {
        for(j = 0; j < K2+1; j++)
        {
            if(j == 0 && gsl_vector_get(survTime2, i) <= gsl_vector_get(s2, 0))
            {
                logLH += gsl_vector_get(lambda2, j);
            }
            if(j != 0 && gsl_vector_get(survTime2, i) > gsl_vector_get(s2, j-1) && gsl_vector_get(survTime2, i) <= gsl_vector_get(s2, j))
            {
                logLH += gsl_vector_get(lambda2, j);
            }
        }
        
        logLH += gsl_vector_get(xbeta2, i);
        logLH += gsl_vector_get(V2, jj);
    }
    
    if(gsl_vector_get(case11, i) == 1)
    {
        for(j = 0; j < K3+1; j++)
        {
            if(j == 0 && gsl_vector_get(yStar, i) <= gsl_vector_get(s3, 0))
            {
                logLH += gsl_vector_get(lambda3, j);
            }
            if(j != 0 && gsl_vector_get(yStar, i) > gsl_vector_get(s3, j-1) && gsl_vector_get(yStar, i) <= gsl_vector_get(s3, j))
            {
                logLH += gsl_vector_get(lambda3, j);
            }
        }
        
        logLH += gsl_vector_get(xbeta3, i);
        logLH += gsl_vector_get(V3, jj);
    }
    
    gfunc = BpeDpCorScrSM_wFunc(i, xbeta1, xbeta2, xbeta3, lambda1, lambda2, lambda3, jj, V1, V2, V3, s1, s2, s3, K1, K2, K3, survTime1, yStar);
    
    logLH += (- 1/theta - gsl_vector_get(survEvent1, i) - gsl_vector_get(survEvent2, i)) * log(1 + theta * gfunc);
    
    *val = logLH;
    
    return;
    
}





/* evaluating log-likelihood function */

/**/

void BpeDpCorScrSM_logLH(gsl_vector *beta1,
                         gsl_vector *beta2,
                         gsl_vector *beta3,
                         gsl_vector *xbeta1,
                         gsl_vector *xbeta2,
                         gsl_vector *xbeta3,
                         gsl_vector *gamma,
                         gsl_vector *lambda1,
                         gsl_vector *lambda2,
                         gsl_vector *lambda3,
                         gsl_vector *s1,
                         gsl_vector *s2,
                         gsl_vector *s3,
                         gsl_vector *V1,
                         gsl_vector *V2,
                         gsl_vector *V3,
                         gsl_vector *survTime1,
                         gsl_vector *survTime2,
                         gsl_vector *yStar,
                         gsl_vector *survEvent1,
                         gsl_vector *case01,
                         gsl_vector *case11,
                         gsl_matrix *survCov1,
                         gsl_matrix *survCov2,
                         gsl_matrix *survCov3,
                         gsl_vector *cluster,
                         int K1,
                         int K2,
                         int K3,
                         double *val)
{
    double gam;
    double logLH = 0;
    

    int n = survTime1 -> size;
    
    int i, j, jj;
    
    for(i = 0; i < n; i++)
    {
        gam = gsl_vector_get(gamma, i);
        
        jj = (int) gsl_vector_get(cluster, i) - 1;
        
        if(gsl_vector_get(survEvent1, i) == 1)
        {
            for(j = 0; j < K1+1; j++)
            {
                if(j == 0 && gsl_vector_get(survTime1, i) <= gsl_vector_get(s1, 0))
                {
                    logLH += gsl_vector_get(lambda1, j);
                }
                if(j != 0 && gsl_vector_get(survTime1, i) > gsl_vector_get(s1, j-1) && gsl_vector_get(survTime1, i) <= gsl_vector_get(s1, j))
                {
                    logLH += gsl_vector_get(lambda1, j);
                }
            }
            logLH += log(gam);
            logLH += gsl_vector_get(xbeta1, i);
            logLH += gsl_vector_get(V1, jj);
        }
        
        if(gsl_vector_get(case01, i) ==  1)
        {
            for(j = 0; j < K2+1; j++)
            {
                if(j == 0 && gsl_vector_get(survTime2, i) <= gsl_vector_get(s2, 0))
                {
                    logLH += gsl_vector_get(lambda2, j);
                }
                if(j != 0 && gsl_vector_get(survTime2, i) > gsl_vector_get(s2, j-1) && gsl_vector_get(survTime2, i) <= gsl_vector_get(s2, j))
                {
                    logLH += gsl_vector_get(lambda2, j);
                }
            }
            logLH += log(gam);
            logLH += gsl_vector_get(xbeta2, i);
            logLH += gsl_vector_get(V2, jj);
        }
        
        if(gsl_vector_get(case11, i) == 1)
        {
            for(j = 0; j < K3+1; j++)
            {
                if(j == 0 && gsl_vector_get(survTime2, i) <= gsl_vector_get(s3, 0))
                {
                    logLH += gsl_vector_get(lambda3, j);
                }
                if(j != 0 && gsl_vector_get(survTime2, i) > gsl_vector_get(s3, j-1) && gsl_vector_get(survTime2, i) <= gsl_vector_get(s3, j))
                {
                    logLH += gsl_vector_get(lambda3, j);
                }
            }
            logLH += log(gam);
            logLH += gsl_vector_get(xbeta3, i);
            logLH += gsl_vector_get(V3, jj);
        }
        logLH += -gam * BpeDpCorScrSM_wFunc(i, xbeta1, xbeta2, xbeta3, lambda1, lambda2, lambda3, jj, V1, V2, V3, s1, s2, s3, K1, K2, K3, survTime1, survTime2);
    }
    
    *val = logLH;
    
    return;
    
}






/* evaluating log-likelihood function for subject i */

/**/

void BpeDpCorScrSM_logLH_i(int i,
                           gsl_vector *beta1,
                           gsl_vector *beta2,
                           gsl_vector *beta3,
                           gsl_vector *xbeta1,
                           gsl_vector *xbeta2,
                           gsl_vector *xbeta3,
                           gsl_vector *gamma,
                           gsl_vector *lambda1,
                           gsl_vector *lambda2,
                           gsl_vector *lambda3,
                           gsl_vector *s1,
                           gsl_vector *s2,
                           gsl_vector *s3,
                           gsl_vector *V1,
                           gsl_vector *V2,
                           gsl_vector *V3,
                           gsl_vector *survTime1,
                           gsl_vector *survTime2,
                           gsl_vector *yStar,
                           gsl_vector *survEvent1,
                           gsl_vector *case01,
                           gsl_vector *case11,
                           gsl_matrix *survCov1,
                           gsl_matrix *survCov2,
                           gsl_matrix *survCov3,
                           gsl_vector *cluster,
                           int K1,
                           int K2,
                           int K3,
                           double *val)
{
    double gam;
    double logLH = 0;
    


    
    int j, jj;
    
    gam = gsl_vector_get(gamma, i);
    
    jj = (int) gsl_vector_get(cluster, i) - 1;
    
    if(gsl_vector_get(survEvent1, i) == 1)
    {
        for(j = 0; j < K1+1; j++)
        {
            if(j == 0 && gsl_vector_get(survTime1, i) <= gsl_vector_get(s1, 0))
            {
                logLH += gsl_vector_get(lambda1, j);
            }
            if(j != 0 && gsl_vector_get(survTime1, i) > gsl_vector_get(s1, j-1) && gsl_vector_get(survTime1, i) <= gsl_vector_get(s1, j))
            {
                logLH += gsl_vector_get(lambda1, j);
            }
        }
        logLH += log(gam);
        logLH += gsl_vector_get(xbeta1, i);
        logLH += gsl_vector_get(V1, jj);
    }
    
    if(gsl_vector_get(case01, i) ==  1)
    {
        for(j = 0; j < K2+1; j++)
        {
            if(j == 0 && gsl_vector_get(survTime2, i) <= gsl_vector_get(s2, 0))
            {
                logLH += gsl_vector_get(lambda2, j);
            }
            if(j != 0 && gsl_vector_get(survTime2, i) > gsl_vector_get(s2, j-1) && gsl_vector_get(survTime2, i) <= gsl_vector_get(s2, j))
            {
                logLH += gsl_vector_get(lambda2, j);
            }
        }
        logLH += log(gam);
        logLH += gsl_vector_get(xbeta2, i);
        logLH += gsl_vector_get(V2, jj);
    }
    
    if(gsl_vector_get(case11, i) == 1)
    {
        for(j = 0; j < K3+1; j++)
        {
            if(j == 0 && gsl_vector_get(yStar, i) <= gsl_vector_get(s3, 0))
            {
                logLH += gsl_vector_get(lambda3, j);
            }
            if(j != 0 && gsl_vector_get(yStar, i) > gsl_vector_get(s3, j-1) && gsl_vector_get(yStar, i) <= gsl_vector_get(s3, j))
            {
                logLH += gsl_vector_get(lambda3, j);
            }
        }
        logLH += log(gam);
        logLH += gsl_vector_get(xbeta3, i);
        logLH += gsl_vector_get(V3, jj);
    }
    
    logLH += -gam * BpeDpCorScrSM_wFunc(i, xbeta1, xbeta2, xbeta3, lambda1, lambda2, lambda3, jj, V1, V2, V3, s1, s2, s3, K1, K2, K3, survTime1, yStar);
    
    *val = logLH;
    
    return;
    
}





/*
 Evaluate w(y1, y2) function
 */
double BpeDpCorScrSM_wFunc(int subjInx,
                           gsl_vector *xbeta1,
                           gsl_vector *xbeta2,
                           gsl_vector *xbeta3,
                           gsl_vector *lambda1,
                           gsl_vector *lambda2,
                           gsl_vector *lambda3,
                           int jj,
                           gsl_vector *V1,
                           gsl_vector *V2,
                           gsl_vector *V3,
                           gsl_vector *s1,
                           gsl_vector *s2,
                           gsl_vector *s3,
                           int J1,
                           int J2,
                           int J3,
                           gsl_vector *survTime1,
                           gsl_vector *yStar)
{
    int i = subjInx;
    double cumHaz1, cumHaz2, cumHaz3diff;
    double Del, wVal;
    int j;
    
    cumHaz1 = 0; cumHaz2 = 0; cumHaz3diff = 0;
    
    for(j = 0; j < J1+1; j++)
    {
        if(j > 0)
        {
            Del = c_max(0, (c_min(gsl_vector_get(s1, j), gsl_vector_get(survTime1, i)) - gsl_vector_get(s1, j-1)));
        }
        if(j == 0)
        {
            Del = c_max(0, c_min(gsl_vector_get(s1, j), gsl_vector_get(survTime1, i)) - 0);
        }
        cumHaz1 += Del* exp(gsl_vector_get(lambda1, j)) * exp(gsl_vector_get(xbeta1, i)+gsl_vector_get(V1, jj));
    }
    
    for(j = 0; j < J2+1; j++)
    {
        if(j > 0)
        {
            Del = c_max(0, (c_min(gsl_vector_get(s2, j), gsl_vector_get(survTime1, i)) - gsl_vector_get(s2, j-1)));
        }
        if(j == 0)
        {
            Del = c_max(0, c_min(gsl_vector_get(s2, j), gsl_vector_get(survTime1, i)) - 0);
        }
        cumHaz2 += Del* exp(gsl_vector_get(lambda2, j)) * exp(gsl_vector_get(xbeta2, i)+gsl_vector_get(V2, jj));
    }
    
    for(j = 0; j < J3+1; j++)
    {
        if(j > 0)
        {
            Del = c_max(0, (c_min(gsl_vector_get(s3, j), gsl_vector_get(yStar, i)) - gsl_vector_get(s3, j-1)));
        }
        if(j == 0)
        {
            Del = c_max(0, c_min(gsl_vector_get(s3, j), gsl_vector_get(yStar, i)) - 0);
        }
        cumHaz3diff += Del* exp(gsl_vector_get(lambda3, j)) * exp(gsl_vector_get(xbeta3, i)+gsl_vector_get(V3, jj));
    }
    
    

    wVal = cumHaz1 + cumHaz2 + cumHaz3diff;
    
    return wVal;
}





















































/********* For PEM-DPM-M model *************/




/* evaluating log-likelihood function */

/**/

void BpeDpCorScr_logMLH(gsl_vector *beta1,
                        gsl_vector *beta2,
                        gsl_vector *beta3,
                        gsl_vector *xbeta1,
                        gsl_vector *xbeta2,
                        gsl_vector *xbeta3,
                        double theta,
                        gsl_vector *lambda1,
                        gsl_vector *lambda2,
                        gsl_vector *lambda3,
                        gsl_vector *s1,
                        gsl_vector *s2,
                        gsl_vector *s3,
                        gsl_vector *V1,
                        gsl_vector *V2,
                        gsl_vector *V3,
                        gsl_vector *survTime1,
                        gsl_vector *survTime2,
                        gsl_vector *survEvent1,
                        gsl_vector *survEvent2,
                        gsl_vector *case01,
                        gsl_vector *case11,
                        gsl_matrix *survCov1,
                        gsl_matrix *survCov2,
                        gsl_matrix *survCov3,
                        gsl_vector *cluster,
                        int K1,
                        int K2,
                        int K3,
                        double *val)
{
    double gfunc;
    double logLH = 0;
    

    int n = survTime1 -> size;
    
    int i, j, jj;
    
    for(i = 0; i < n; i++)
    {
        
        jj = (int) gsl_vector_get(cluster, i) - 1;
        
        if(gsl_vector_get(survEvent1, i) == 1)
        {
            for(j = 0; j < K1+1; j++)
            {
                if(j == 0 && gsl_vector_get(survTime1, i) <= gsl_vector_get(s1, 0))
                {
                    logLH += gsl_vector_get(lambda1, j);
                }
                if(j != 0 && gsl_vector_get(survTime1, i) > gsl_vector_get(s1, j-1) && gsl_vector_get(survTime1, i) <= gsl_vector_get(s1, j))
                {
                    logLH += gsl_vector_get(lambda1, j);
                }
            }
            logLH += gsl_vector_get(xbeta1, i);
            logLH += gsl_vector_get(V1, jj);
        }
        
        if(gsl_vector_get(case01, i) ==  1)
        {
            for(j = 0; j < K2+1; j++)
            {
                if(j == 0 && gsl_vector_get(survTime2, i) <= gsl_vector_get(s2, 0))
                {
                    logLH += gsl_vector_get(lambda2, j);
                }
                if(j != 0 && gsl_vector_get(survTime2, i) > gsl_vector_get(s2, j-1) && gsl_vector_get(survTime2, i) <= gsl_vector_get(s2, j))
                {
                    logLH += gsl_vector_get(lambda2, j);
                }
            }
            logLH += gsl_vector_get(xbeta2, i);
            logLH += gsl_vector_get(V2, jj);
        }
        
        if(gsl_vector_get(case11, i) == 1)
        {
            for(j = 0; j < K3+1; j++)
            {
                if(j == 0 && gsl_vector_get(survTime2, i) <= gsl_vector_get(s3, 0))
                {
                    logLH += gsl_vector_get(lambda3, j);
                }
                if(j != 0 && gsl_vector_get(survTime2, i) > gsl_vector_get(s3, j-1) && gsl_vector_get(survTime2, i) <= gsl_vector_get(s3, j))
                {
                    logLH += gsl_vector_get(lambda3, j);
                }
            }
            logLH += gsl_vector_get(xbeta3, i);
            logLH += gsl_vector_get(V3, jj);
        }
        gfunc = BpeDpCorScr_wFunc(i, xbeta1, xbeta2, xbeta3, lambda1, lambda2, lambda3, jj, V1, V2, V3, s1, s2, s3, K1, K2, K3, survTime1, survTime2);
        
        logLH += (- 1/theta - gsl_vector_get(survEvent1, i) - gsl_vector_get(survEvent2, i)) * log(1 + theta * gfunc);
    }
    
    *val = logLH;
    
    return;
    
}





/* evaluating log-likelihood function for subject i */

/**/

void BpeDpCorScr_logMLH_i(int i,
                          gsl_vector *beta1,
                          gsl_vector *beta2,
                          gsl_vector *beta3,
                          gsl_vector *xbeta1,
                          gsl_vector *xbeta2,
                          gsl_vector *xbeta3,
                          double theta,
                          gsl_vector *lambda1,
                          gsl_vector *lambda2,
                          gsl_vector *lambda3,
                          gsl_vector *s1,
                          gsl_vector *s2,
                          gsl_vector *s3,
                          gsl_vector *V1,
                          gsl_vector *V2,
                          gsl_vector *V3,
                          gsl_vector *survTime1,
                          gsl_vector *survTime2,
                          gsl_vector *survEvent1,
                          gsl_vector *survEvent2,
                          gsl_vector *case01,
                          gsl_vector *case11,
                          gsl_matrix *survCov1,
                          gsl_matrix *survCov2,
                          gsl_matrix *survCov3,
                          gsl_vector *cluster,
                          int K1,
                          int K2,
                          int K3,
                          double *val)
{
    double gfunc;
    double logLH = 0;
    


    
    int j, jj;
    
    jj = (int) gsl_vector_get(cluster, i) - 1;
    
    if(gsl_vector_get(survEvent1, i) == 1)
    {
        for(j = 0; j < K1+1; j++)
        {
            if(j == 0 && gsl_vector_get(survTime1, i) <= gsl_vector_get(s1, 0))
            {
                logLH += gsl_vector_get(lambda1, j);
            }
            if(j != 0 && gsl_vector_get(survTime1, i) > gsl_vector_get(s1, j-1) && gsl_vector_get(survTime1, i) <= gsl_vector_get(s1, j))
            {
                logLH += gsl_vector_get(lambda1, j);
            }
        }
        logLH += gsl_vector_get(xbeta1, i);
        logLH += gsl_vector_get(V1, jj);
    }
    
    if(gsl_vector_get(case01, i) ==  1)
    {
        for(j = 0; j < K2+1; j++)
        {
            if(j == 0 && gsl_vector_get(survTime2, i) <= gsl_vector_get(s2, 0))
            {
                logLH += gsl_vector_get(lambda2, j);
            }
            if(j != 0 && gsl_vector_get(survTime2, i) > gsl_vector_get(s2, j-1) && gsl_vector_get(survTime2, i) <= gsl_vector_get(s2, j))
            {
                logLH += gsl_vector_get(lambda2, j);
            }
        }
        logLH += gsl_vector_get(xbeta2, i);
        logLH += gsl_vector_get(V2, jj);
    }
    
    if(gsl_vector_get(case11, i) == 1)
    {
        for(j = 0; j < K3+1; j++)
        {
            if(j == 0 && gsl_vector_get(survTime2, i) <= gsl_vector_get(s3, 0))
            {
                logLH += gsl_vector_get(lambda3, j);
            }
            if(j != 0 && gsl_vector_get(survTime2, i) > gsl_vector_get(s3, j-1) && gsl_vector_get(survTime2, i) <= gsl_vector_get(s3, j))
            {
                logLH += gsl_vector_get(lambda3, j);
            }
        }
        logLH += gsl_vector_get(xbeta3, i);
        logLH += gsl_vector_get(V3, jj);
    }
    
    
    gfunc = BpeDpCorScr_wFunc(i, xbeta1, xbeta2, xbeta3, lambda1, lambda2, lambda3, jj, V1, V2, V3, s1, s2, s3, K1, K2, K3, survTime1, survTime2);
    
    logLH += (- 1/theta - gsl_vector_get(survEvent1, i) - gsl_vector_get(survEvent2, i)) * log(1 + theta * gfunc);
    
    *val = logLH;
    
    return;
    
}





/* evaluating log-likelihood function */

/**/

void BpeDpCorScr_logLH(gsl_vector *beta1,
                       gsl_vector *beta2,
                       gsl_vector *beta3,
                       gsl_vector *xbeta1,
                       gsl_vector *xbeta2,
                       gsl_vector *xbeta3,
                       gsl_vector *gamma,
                       gsl_vector *lambda1,
                       gsl_vector *lambda2,
                       gsl_vector *lambda3,
                       gsl_vector *s1,
                       gsl_vector *s2,
                       gsl_vector *s3,
                       gsl_vector *V1,
                       gsl_vector *V2,
                       gsl_vector *V3,
                       gsl_vector *survTime1,
                       gsl_vector *survTime2,
                       gsl_vector *survEvent1,
                       gsl_vector *case01,
                       gsl_vector *case11,
                       gsl_matrix *survCov1,
                       gsl_matrix *survCov2,
                       gsl_matrix *survCov3,
                       gsl_vector *cluster,
                       int K1,
                       int K2,
                       int K3,
                       double *val)
{
    double gam;
    double logLH = 0;
    

    int n = survTime1 -> size;
    
    int i, j, jj;
    
    for(i = 0; i < n; i++)
    {
        gam = gsl_vector_get(gamma, i);
        
        jj = (int) gsl_vector_get(cluster, i) - 1;
        
        if(gsl_vector_get(survEvent1, i) == 1)
        {
            for(j = 0; j < K1+1; j++)
            {
                if(j == 0 && gsl_vector_get(survTime1, i) <= gsl_vector_get(s1, 0))
                {
                    logLH += gsl_vector_get(lambda1, j);
                }
                if(j != 0 && gsl_vector_get(survTime1, i) > gsl_vector_get(s1, j-1) && gsl_vector_get(survTime1, i) <= gsl_vector_get(s1, j))
                {
                    logLH += gsl_vector_get(lambda1, j);
                }
            }
            logLH += log(gam);
            logLH += gsl_vector_get(xbeta1, i);
            logLH += gsl_vector_get(V1, jj);
        }
        
        if(gsl_vector_get(case01, i) ==  1)
        {
            for(j = 0; j < K2+1; j++)
            {
                if(j == 0 && gsl_vector_get(survTime2, i) <= gsl_vector_get(s2, 0))
                {
                    logLH += gsl_vector_get(lambda2, j);
                }
                if(j != 0 && gsl_vector_get(survTime2, i) > gsl_vector_get(s2, j-1) && gsl_vector_get(survTime2, i) <= gsl_vector_get(s2, j))
                {
                    logLH += gsl_vector_get(lambda2, j);
                }
            }
            logLH += log(gam);
            logLH += gsl_vector_get(xbeta2, i);
            logLH += gsl_vector_get(V2, jj);
        }
        
        if(gsl_vector_get(case11, i) == 1)
        {
            for(j = 0; j < K3+1; j++)
            {
                if(j == 0 && gsl_vector_get(survTime2, i) <= gsl_vector_get(s3, 0))
                {
                    logLH += gsl_vector_get(lambda3, j);
                }
                if(j != 0 && gsl_vector_get(survTime2, i) > gsl_vector_get(s3, j-1) && gsl_vector_get(survTime2, i) <= gsl_vector_get(s3, j))
                {
                    logLH += gsl_vector_get(lambda3, j);
                }
            }
            logLH += log(gam);
            logLH += gsl_vector_get(xbeta3, i);
            logLH += gsl_vector_get(V3, jj);
        }
        logLH += -gam * BpeDpCorScr_wFunc(i, xbeta1, xbeta2, xbeta3, lambda1, lambda2, lambda3, jj, V1, V2, V3, s1, s2, s3, K1, K2, K3, survTime1, survTime2);
    }
    
    *val = logLH;
    
    return;
    
}





/* evaluating log-likelihood function for subject i */

/**/

void BpeDpCorScr_logLH_i(int i,
                         gsl_vector *beta1,
                         gsl_vector *beta2,
                         gsl_vector *beta3,
                         gsl_vector *xbeta1,
                         gsl_vector *xbeta2,
                         gsl_vector *xbeta3,
                         gsl_vector *gamma,
                         gsl_vector *lambda1,
                         gsl_vector *lambda2,
                         gsl_vector *lambda3,
                         gsl_vector *s1,
                         gsl_vector *s2,
                         gsl_vector *s3,
                         gsl_vector *V1,
                         gsl_vector *V2,
                         gsl_vector *V3,
                         gsl_vector *survTime1,
                         gsl_vector *survTime2,
                         gsl_vector *survEvent1,
                         gsl_vector *case01,
                         gsl_vector *case11,
                         gsl_matrix *survCov1,
                         gsl_matrix *survCov2,
                         gsl_matrix *survCov3,
                         gsl_vector *cluster,
                         int K1,
                         int K2,
                         int K3,
                         double *val)
{
    double gam;
    double logLH = 0;
    


    
    int j, jj;
    
    gam = gsl_vector_get(gamma, i);
    
    jj = (int) gsl_vector_get(cluster, i) - 1;
    
    if(gsl_vector_get(survEvent1, i) == 1)
    {
        for(j = 0; j < K1+1; j++)
        {
            if(j == 0 && gsl_vector_get(survTime1, i) <= gsl_vector_get(s1, 0))
            {
                logLH += gsl_vector_get(lambda1, j);
            }
            if(j != 0 && gsl_vector_get(survTime1, i) > gsl_vector_get(s1, j-1) && gsl_vector_get(survTime1, i) <= gsl_vector_get(s1, j))
            {
                logLH += gsl_vector_get(lambda1, j);
            }
        }
        logLH += log(gam);
        logLH += gsl_vector_get(xbeta1, i);
        logLH += gsl_vector_get(V1, jj);
    }
    
    if(gsl_vector_get(case01, i) ==  1)
    {
        for(j = 0; j < K2+1; j++)
        {
            if(j == 0 && gsl_vector_get(survTime2, i) <= gsl_vector_get(s2, 0))
            {
                logLH += gsl_vector_get(lambda2, j);
            }
            if(j != 0 && gsl_vector_get(survTime2, i) > gsl_vector_get(s2, j-1) && gsl_vector_get(survTime2, i) <= gsl_vector_get(s2, j))
            {
                logLH += gsl_vector_get(lambda2, j);
            }
        }
        logLH += log(gam);
        logLH += gsl_vector_get(xbeta2, i);
        logLH += gsl_vector_get(V2, jj);
    }
    
    if(gsl_vector_get(case11, i) == 1)
    {
        for(j = 0; j < K3+1; j++)
        {
            if(j == 0 && gsl_vector_get(survTime2, i) <= gsl_vector_get(s3, 0))
            {
                logLH += gsl_vector_get(lambda3, j);
            }
            if(j != 0 && gsl_vector_get(survTime2, i) > gsl_vector_get(s3, j-1) && gsl_vector_get(survTime2, i) <= gsl_vector_get(s3, j))
            {
                logLH += gsl_vector_get(lambda3, j);
            }
        }
        logLH += log(gam);
        logLH += gsl_vector_get(xbeta3, i);
        logLH += gsl_vector_get(V3, jj);
    }
    
    logLH += -gam * BpeDpCorScr_wFunc(i, xbeta1, xbeta2, xbeta3, lambda1, lambda2, lambda3, jj, V1, V2, V3, s1, s2, s3, K1, K2, K3, survTime1, survTime2);
    
    *val = logLH;
    
    return;
    
}





/*
 Evaluate w(y1, y2) function
 */
double BpeDpCorScr_wFunc(int subjInx,
                         gsl_vector *xbeta1,
                         gsl_vector *xbeta2,
                         gsl_vector *xbeta3,
                         gsl_vector *lambda1,
                         gsl_vector *lambda2,
                         gsl_vector *lambda3,
                         int jj,
                         gsl_vector *V1,
                         gsl_vector *V2,
                         gsl_vector *V3,
                         gsl_vector *s1,
                         gsl_vector *s2,
                         gsl_vector *s3,
                         int J1,
                         int J2,
                         int J3,
                         gsl_vector *survTime1,
                         gsl_vector *survTime2)
{
    int i = subjInx;
    double cumHaz1, cumHaz2, cumHaz3diff;
    double Del, wVal;
    int j;
    
    cumHaz1 = 0; cumHaz2 = 0; cumHaz3diff = 0;
    
    for(j = 0; j < J1+1; j++)
    {
        if(j > 0)
        {
            Del = c_max(0, (c_min(gsl_vector_get(s1, j), gsl_vector_get(survTime1, i)) - gsl_vector_get(s1, j-1)));
        }
        if(j == 0)
        {
            Del = c_max(0, c_min(gsl_vector_get(s1, j), gsl_vector_get(survTime1, i)) - 0);
        }
        cumHaz1 += Del* exp(gsl_vector_get(lambda1, j)) * exp(gsl_vector_get(xbeta1, i)+gsl_vector_get(V1, jj));
    }
    
    for(j = 0; j < J2+1; j++)
    {
        if(j > 0)
        {
            Del = c_max(0, (c_min(gsl_vector_get(s2, j), gsl_vector_get(survTime1, i)) - gsl_vector_get(s2, j-1)));
        }
        if(j == 0)
        {
            Del = c_max(0, c_min(gsl_vector_get(s2, j), gsl_vector_get(survTime1, i)) - 0);
        }
        cumHaz2 += Del* exp(gsl_vector_get(lambda2, j)) * exp(gsl_vector_get(xbeta2, i)+gsl_vector_get(V2, jj));
    }
    
    for(j = 0; j < J3+1; j++)
    {
        if(j > 0)
        {
            Del = c_max(0, (c_min(gsl_vector_get(s3, j), gsl_vector_get(survTime2, i)) - c_max(gsl_vector_get(s3, j-1), gsl_vector_get(survTime1, i))));
        }
        if(j == 0)
        {
            Del = c_max(0, (c_min(gsl_vector_get(s3, j), gsl_vector_get(survTime2, i)) - c_max(0, gsl_vector_get(survTime1, i))));
        }
        cumHaz3diff += Del* exp(gsl_vector_get(lambda3, j)) * exp(gsl_vector_get(xbeta3, i)+gsl_vector_get(V3, jj));
    }
    

    
    wVal = cumHaz1 + cumHaz2 + cumHaz3diff;
    
    return wVal;
}













































/********* For PEM-MVN-SM model *************/


/* evaluating log-likelihood function for subject i */

/**/

void BpeMvnCorScrSM_logMLH_i(int i,
                             gsl_vector *beta1,
                             gsl_vector *beta2,
                             gsl_vector *beta3,
                             gsl_vector *xbeta1,
                             gsl_vector *xbeta2,
                             gsl_vector *xbeta3,
                             double theta,
                             gsl_vector *lambda1,
                             gsl_vector *lambda2,
                             gsl_vector *lambda3,
                             gsl_vector *s1,
                             gsl_vector *s2,
                             gsl_vector *s3,
                             gsl_vector *V1,
                             gsl_vector *V2,
                             gsl_vector *V3,
                             gsl_vector *survTime1,
                             gsl_vector *survTime2,
                             gsl_vector *yStar,
                             gsl_vector *survEvent1,
                             gsl_vector *survEvent2,
                             gsl_vector *case01,
                             gsl_vector *case11,
                             gsl_matrix *survCov1,
                             gsl_matrix *survCov2,
                             gsl_matrix *survCov3,
                             gsl_vector *cluster,
                             int K1,
                             int K2,
                             int K3,
                             double *val)
{
    double gfunc;
    double logLH = 0;
    


    
    int j, jj;
    
    
    
    jj = (int) gsl_vector_get(cluster, i) - 1;
    
    if(gsl_vector_get(survEvent1, i) == 1)
    {
        for(j = 0; j < K1+1; j++)
        {
            if(j == 0 && gsl_vector_get(survTime1, i) <= gsl_vector_get(s1, 0))
            {
                logLH += gsl_vector_get(lambda1, j);
            }
            if(j != 0 && gsl_vector_get(survTime1, i) > gsl_vector_get(s1, j-1) && gsl_vector_get(survTime1, i) <= gsl_vector_get(s1, j))
            {
                logLH += gsl_vector_get(lambda1, j);
            }
        }
        
        logLH += gsl_vector_get(xbeta1, i);
        logLH += gsl_vector_get(V1, jj);
    }
    
    if(gsl_vector_get(case01, i) ==  1)
    {
        for(j = 0; j < K2+1; j++)
        {
            if(j == 0 && gsl_vector_get(survTime2, i) <= gsl_vector_get(s2, 0))
            {
                logLH += gsl_vector_get(lambda2, j);
            }
            if(j != 0 && gsl_vector_get(survTime2, i) > gsl_vector_get(s2, j-1) && gsl_vector_get(survTime2, i) <= gsl_vector_get(s2, j))
            {
                logLH += gsl_vector_get(lambda2, j);
            }
        }
        
        logLH += gsl_vector_get(xbeta2, i);
        logLH += gsl_vector_get(V2, jj);
    }
    
    if(gsl_vector_get(case11, i) == 1)
    {
        for(j = 0; j < K3+1; j++)
        {
            if(j == 0 && gsl_vector_get(yStar, i) <= gsl_vector_get(s3, 0))
            {
                logLH += gsl_vector_get(lambda3, j);
            }
            if(j != 0 && gsl_vector_get(yStar, i) > gsl_vector_get(s3, j-1) && gsl_vector_get(yStar, i) <= gsl_vector_get(s3, j))
            {
                logLH += gsl_vector_get(lambda3, j);
            }
        }
        
        logLH += gsl_vector_get(xbeta3, i);
        logLH += gsl_vector_get(V3, jj);
    }
    
    gfunc = BpeMvnCorScrSM_wFunc(i, xbeta1, xbeta2, xbeta3, lambda1, lambda2, lambda3, jj, V1, V2, V3, s1, s2, s3, K1, K2, K3, survTime1, yStar);
    
    
    logLH += (- 1/theta - gsl_vector_get(survEvent1, i) - gsl_vector_get(survEvent2, i)) * log(1 + theta * gfunc);
    
    *val = logLH;
    
    return;
    
}





/* evaluating log-likelihood function */

/**/

void BpeMvnCorScrSM_logMLH(gsl_vector *beta1,
                           gsl_vector *beta2,
                           gsl_vector *beta3,
                           gsl_vector *xbeta1,
                           gsl_vector *xbeta2,
                           gsl_vector *xbeta3,
                           double theta,
                           gsl_vector *lambda1,
                           gsl_vector *lambda2,
                           gsl_vector *lambda3,
                           gsl_vector *s1,
                           gsl_vector *s2,
                           gsl_vector *s3,
                           gsl_vector *V1,
                           gsl_vector *V2,
                           gsl_vector *V3,
                           gsl_vector *survTime1,
                           gsl_vector *survTime2,
                           gsl_vector *yStar,
                           gsl_vector *survEvent1,
                           gsl_vector *survEvent2,
                           gsl_vector *case01,
                           gsl_vector *case11,
                           gsl_matrix *survCov1,
                           gsl_matrix *survCov2,
                           gsl_matrix *survCov3,
                           gsl_vector *cluster,
                           int K1,
                           int K2,
                           int K3,
                           double *val)
{
    double gfunc;
    double logLH = 0;
    

    int n = survTime1 -> size;
    
    int i, j, jj;
    
    for(i = 0; i < n; i++)
    {
        
        
        jj = (int) gsl_vector_get(cluster, i) - 1;
        
        if(gsl_vector_get(survEvent1, i) == 1)
        {
            for(j = 0; j < K1+1; j++)
            {
                if(j == 0 && gsl_vector_get(survTime1, i) <= gsl_vector_get(s1, 0))
                {
                    logLH += gsl_vector_get(lambda1, j);
                }
                if(j != 0 && gsl_vector_get(survTime1, i) > gsl_vector_get(s1, j-1) && gsl_vector_get(survTime1, i) <= gsl_vector_get(s1, j))
                {
                    logLH += gsl_vector_get(lambda1, j);
                }
            }
            
            logLH += gsl_vector_get(xbeta1, i);
            logLH += gsl_vector_get(V1, jj);
        }
        
        if(gsl_vector_get(case01, i) ==  1)
        {
            for(j = 0; j < K2+1; j++)
            {
                if(j == 0 && gsl_vector_get(survTime2, i) <= gsl_vector_get(s2, 0))
                {
                    logLH += gsl_vector_get(lambda2, j);
                }
                if(j != 0 && gsl_vector_get(survTime2, i) > gsl_vector_get(s2, j-1) && gsl_vector_get(survTime2, i) <= gsl_vector_get(s2, j))
                {
                    logLH += gsl_vector_get(lambda2, j);
                }
            }
            
            logLH += gsl_vector_get(xbeta2, i);
            logLH += gsl_vector_get(V2, jj);
        }
        
        if(gsl_vector_get(case11, i) == 1)
        {
            for(j = 0; j < K3+1; j++)
            {
                if(j == 0 && gsl_vector_get(yStar, i) <= gsl_vector_get(s3, 0))
                {
                    logLH += gsl_vector_get(lambda3, j);
                }
                if(j != 0 && gsl_vector_get(yStar, i) > gsl_vector_get(s3, j-1) && gsl_vector_get(yStar, i) <= gsl_vector_get(s3, j))
                {
                    logLH += gsl_vector_get(lambda3, j);
                }
            }
            
            logLH += gsl_vector_get(xbeta3, i);
            logLH += gsl_vector_get(V3, jj);
        }
        gfunc = BpeMvnCorScrSM_wFunc(i, xbeta1, xbeta2, xbeta3, lambda1, lambda2, lambda3, jj, V1, V2, V3, s1, s2, s3, K1, K2, K3, survTime1, yStar);
        
        logLH +=  (- 1/theta - gsl_vector_get(survEvent1, i) - gsl_vector_get(survEvent2, i)) * log(1 + theta * gfunc);
    }
    
    *val = logLH;
    
    return;
    
}







/* evaluating log-likelihood function for subject i */

/**/

void BpeMvnCorScrSM_logLH_i(int i,
                            gsl_vector *beta1,
                            gsl_vector *beta2,
                            gsl_vector *beta3,
                            gsl_vector *xbeta1,
                            gsl_vector *xbeta2,
                            gsl_vector *xbeta3,
                            gsl_vector *gamma,
                            gsl_vector *lambda1,
                            gsl_vector *lambda2,
                            gsl_vector *lambda3,
                            gsl_vector *s1,
                            gsl_vector *s2,
                            gsl_vector *s3,
                            gsl_vector *V1,
                            gsl_vector *V2,
                            gsl_vector *V3,
                            gsl_vector *survTime1,
                            gsl_vector *survTime2,
                            gsl_vector *yStar,
                            gsl_vector *survEvent1,
                            gsl_vector *case01,
                            gsl_vector *case11,
                            gsl_matrix *survCov1,
                            gsl_matrix *survCov2,
                            gsl_matrix *survCov3,
                            gsl_vector *cluster,
                            int K1,
                            int K2,
                            int K3,
                            double *val)
{
    double gam;
    double logLH = 0;
    


    
    int j, jj;
    
    gam = gsl_vector_get(gamma, i);
    
    jj = (int) gsl_vector_get(cluster, i) - 1;
    
    if(gsl_vector_get(survEvent1, i) == 1)
    {
        for(j = 0; j < K1+1; j++)
        {
            if(j == 0 && gsl_vector_get(survTime1, i) <= gsl_vector_get(s1, 0))
            {
                logLH += gsl_vector_get(lambda1, j);
            }
            if(j != 0 && gsl_vector_get(survTime1, i) > gsl_vector_get(s1, j-1) && gsl_vector_get(survTime1, i) <= gsl_vector_get(s1, j))
            {
                logLH += gsl_vector_get(lambda1, j);
            }
        }
        logLH += log(gam);
        logLH += gsl_vector_get(xbeta1, i);
        logLH += gsl_vector_get(V1, jj);
    }
    
    if(gsl_vector_get(case01, i) ==  1)
    {
        for(j = 0; j < K2+1; j++)
        {
            if(j == 0 && gsl_vector_get(survTime2, i) <= gsl_vector_get(s2, 0))
            {
                logLH += gsl_vector_get(lambda2, j);
            }
            if(j != 0 && gsl_vector_get(survTime2, i) > gsl_vector_get(s2, j-1) && gsl_vector_get(survTime2, i) <= gsl_vector_get(s2, j))
            {
                logLH += gsl_vector_get(lambda2, j);
            }
        }
        logLH += log(gam);
        logLH += gsl_vector_get(xbeta2, i);
        logLH += gsl_vector_get(V2, jj);
    }
    
    if(gsl_vector_get(case11, i) == 1)
    {
        for(j = 0; j < K3+1; j++)
        {
            if(j == 0 && gsl_vector_get(yStar, i) <= gsl_vector_get(s3, 0))
            {
                logLH += gsl_vector_get(lambda3, j);
            }
            if(j != 0 && gsl_vector_get(yStar, i) > gsl_vector_get(s3, j-1) && gsl_vector_get(yStar, i) <= gsl_vector_get(s3, j))
            {
                logLH += gsl_vector_get(lambda3, j);
            }
        }
        logLH += log(gam);
        logLH += gsl_vector_get(xbeta3, i);
        logLH += gsl_vector_get(V3, jj);
    }
    
    logLH += -gam * BpeMvnCorScrSM_wFunc(i, xbeta1, xbeta2, xbeta3, lambda1, lambda2, lambda3, jj, V1, V2, V3, s1, s2, s3, K1, K2, K3, survTime1, yStar);
    
    *val = logLH;
    
    return;
    
}





/* evaluating log-likelihood function */

/**/

void BpeMvnCorScrSM_logLH(gsl_vector *beta1,
                          gsl_vector *beta2,
                          gsl_vector *beta3,
                          gsl_vector *xbeta1,
                          gsl_vector *xbeta2,
                          gsl_vector *xbeta3,
                          gsl_vector *gamma,
                          gsl_vector *lambda1,
                          gsl_vector *lambda2,
                          gsl_vector *lambda3,
                          gsl_vector *s1,
                          gsl_vector *s2,
                          gsl_vector *s3,
                          gsl_vector *V1,
                          gsl_vector *V2,
                          gsl_vector *V3,
                          gsl_vector *survTime1,
                          gsl_vector *survTime2,
                          gsl_vector *yStar,
                          gsl_vector *survEvent1,
                          gsl_vector *case01,
                          gsl_vector *case11,
                          gsl_matrix *survCov1,
                          gsl_matrix *survCov2,
                          gsl_matrix *survCov3,
                          gsl_vector *cluster,
                          int K1,
                          int K2,
                          int K3,
                          double *val)
{
    double gam;
    double logLH = 0;
    

    int n = survTime1 -> size;
    
    int i, j, jj;
    
    for(i = 0; i < n; i++)
    {
        gam = gsl_vector_get(gamma, i);
        
        jj = (int) gsl_vector_get(cluster, i) - 1;
        
        if(gsl_vector_get(survEvent1, i) == 1)
        {
            for(j = 0; j < K1+1; j++)
            {
                if(j == 0 && gsl_vector_get(survTime1, i) <= gsl_vector_get(s1, 0))
                {
                    logLH += gsl_vector_get(lambda1, j);
                }
                if(j != 0 && gsl_vector_get(survTime1, i) > gsl_vector_get(s1, j-1) && gsl_vector_get(survTime1, i) <= gsl_vector_get(s1, j))
                {
                    logLH += gsl_vector_get(lambda1, j);
                }
            }
            logLH += log(gam);
            logLH += gsl_vector_get(xbeta1, i);
            logLH += gsl_vector_get(V1, jj);
        }
        
        if(gsl_vector_get(case01, i) ==  1)
        {
            for(j = 0; j < K2+1; j++)
            {
                if(j == 0 && gsl_vector_get(survTime2, i) <= gsl_vector_get(s2, 0))
                {
                    logLH += gsl_vector_get(lambda2, j);
                }
                if(j != 0 && gsl_vector_get(survTime2, i) > gsl_vector_get(s2, j-1) && gsl_vector_get(survTime2, i) <= gsl_vector_get(s2, j))
                {
                    logLH += gsl_vector_get(lambda2, j);
                }
            }
            logLH += log(gam);
            logLH += gsl_vector_get(xbeta2, i);
            logLH += gsl_vector_get(V2, jj);
        }
        
        if(gsl_vector_get(case11, i) == 1)
        {
            for(j = 0; j < K3+1; j++)
            {
                if(j == 0 && gsl_vector_get(yStar, i) <= gsl_vector_get(s3, 0))
                {
                    logLH += gsl_vector_get(lambda3, j);
                }
                if(j != 0 && gsl_vector_get(yStar, i) > gsl_vector_get(s3, j-1) && gsl_vector_get(yStar, i) <= gsl_vector_get(s3, j))
                {
                    logLH += gsl_vector_get(lambda3, j);
                }
            }
            logLH += log(gam);
            logLH += gsl_vector_get(xbeta3, i);
            logLH += gsl_vector_get(V3, jj);
        }
        logLH += -gam * BpeMvnCorScrSM_wFunc(i, xbeta1, xbeta2, xbeta3, lambda1, lambda2, lambda3, jj, V1, V2, V3, s1, s2, s3, K1, K2, K3, survTime1, yStar);
    }
    
    *val = logLH;
    
    return;
    
}







/*
 Evaluate w(y1, y2) function
 */
double BpeMvnCorScrSM_wFunc(int subjInx,
                            gsl_vector *xbeta1,
                            gsl_vector *xbeta2,
                            gsl_vector *xbeta3,
                            gsl_vector *lambda1,
                            gsl_vector *lambda2,
                            gsl_vector *lambda3,
                            int jj,
                            gsl_vector *V1,
                            gsl_vector *V2,
                            gsl_vector *V3,
                            gsl_vector *s1,
                            gsl_vector *s2,
                            gsl_vector *s3,
                            int J1,
                            int J2,
                            int J3,
                            gsl_vector *survTime1,
                            gsl_vector *yStar)
{
    int i = subjInx;
    double cumHaz1, cumHaz2, cumHaz3diff;
    double Del, wVal;
    int j;
    
    cumHaz1 = 0; cumHaz2 = 0; cumHaz3diff = 0;
    
    for(j = 0; j < J1+1; j++)
    {
        if(j > 0)
        {
            Del = c_max(0, (c_min(gsl_vector_get(s1, j), gsl_vector_get(survTime1, i)) - gsl_vector_get(s1, j-1)));
        }
        if(j == 0)
        {
            Del = c_max(0, c_min(gsl_vector_get(s1, j), gsl_vector_get(survTime1, i)) - 0);
        }
        cumHaz1 += Del* exp(gsl_vector_get(lambda1, j)) * exp(gsl_vector_get(xbeta1, i)+gsl_vector_get(V1, jj));
    }
    
    for(j = 0; j < J2+1; j++)
    {
        if(j > 0)
        {
            Del = c_max(0, (c_min(gsl_vector_get(s2, j), gsl_vector_get(survTime1, i)) - gsl_vector_get(s2, j-1)));
        }
        if(j == 0)
        {
            Del = c_max(0, c_min(gsl_vector_get(s2, j), gsl_vector_get(survTime1, i)) - 0);
        }
        cumHaz2 += Del* exp(gsl_vector_get(lambda2, j)) * exp(gsl_vector_get(xbeta2, i)+gsl_vector_get(V2, jj));
    }
    
    for(j = 0; j < J3+1; j++)
    {
        if(j > 0)
        {
            Del = c_max(0, (c_min(gsl_vector_get(s3, j), gsl_vector_get(yStar, i)) - gsl_vector_get(s3, j-1)));
        }
        if(j == 0)
        {
            Del = c_max(0, c_min(gsl_vector_get(s3, j), gsl_vector_get(yStar, i)) - 0);
        }
        cumHaz3diff += Del* exp(gsl_vector_get(lambda3, j)) * exp(gsl_vector_get(xbeta3, i)+gsl_vector_get(V3, jj));
    }
    

    
    wVal = cumHaz1 + cumHaz2 + cumHaz3diff;
    
    return wVal;
}







































/********* For PEM-MVN-M model *************/



/* evaluating log-likelihood function */

/**/

void BpeMvnCorScr_logLH(gsl_vector *beta1,
                        gsl_vector *beta2,
                        gsl_vector *beta3,
                        gsl_vector *xbeta1,
                        gsl_vector *xbeta2,
                        gsl_vector *xbeta3,
                        gsl_vector *gamma,
                        gsl_vector *lambda1,
                        gsl_vector *lambda2,
                        gsl_vector *lambda3,
                        gsl_vector *s1,
                        gsl_vector *s2,
                        gsl_vector *s3,
                        gsl_vector *V1,
                        gsl_vector *V2,
                        gsl_vector *V3,
                        gsl_vector *survTime1,
                        gsl_vector *survTime2,
                        gsl_vector *survEvent1,
                        gsl_vector *case01,
                        gsl_vector *case11,
                        gsl_matrix *survCov1,
                        gsl_matrix *survCov2,
                        gsl_matrix *survCov3,
                        gsl_vector *cluster,
                        int K1,
                        int K2,
                        int K3,
                        double *val)
{
    double gam;
    double logLH = 0;
    

    int n = survTime1 -> size;
    
    int i, j, jj;
    
    for(i = 0; i < n; i++)
    {
        gam = gsl_vector_get(gamma, i);
        
        jj = (int) gsl_vector_get(cluster, i) - 1;
        
        if(gsl_vector_get(survEvent1, i) == 1)
        {
            for(j = 0; j < K1+1; j++)
            {
                if(j == 0 && gsl_vector_get(survTime1, i) <= gsl_vector_get(s1, 0))
                {
                    logLH += gsl_vector_get(lambda1, j);
                }
                if(j != 0 && gsl_vector_get(survTime1, i) > gsl_vector_get(s1, j-1) && gsl_vector_get(survTime1, i) <= gsl_vector_get(s1, j))
                {
                    logLH += gsl_vector_get(lambda1, j);
                }
            }
            logLH += log(gam);
            logLH += gsl_vector_get(xbeta1, i);
            logLH += gsl_vector_get(V1, jj);
        }
        
        if(gsl_vector_get(case01, i) ==  1)
        {
            for(j = 0; j < K2+1; j++)
            {
                if(j == 0 && gsl_vector_get(survTime2, i) <= gsl_vector_get(s2, 0))
                {
                    logLH += gsl_vector_get(lambda2, j);
                }
                if(j != 0 && gsl_vector_get(survTime2, i) > gsl_vector_get(s2, j-1) && gsl_vector_get(survTime2, i) <= gsl_vector_get(s2, j))
                {
                    logLH += gsl_vector_get(lambda2, j);
                }
            }
            logLH += log(gam);
            logLH += gsl_vector_get(xbeta2, i);
            logLH += gsl_vector_get(V2, jj);
        }
        
        if(gsl_vector_get(case11, i) == 1)
        {
            for(j = 0; j < K3+1; j++)
            {
                if(j == 0 && gsl_vector_get(survTime2, i) <= gsl_vector_get(s3, 0))
                {
                    logLH += gsl_vector_get(lambda3, j);
                }
                if(j != 0 && gsl_vector_get(survTime2, i) > gsl_vector_get(s3, j-1) && gsl_vector_get(survTime2, i) <= gsl_vector_get(s3, j))
                {
                    logLH += gsl_vector_get(lambda3, j);
                }
            }
            logLH += log(gam);
            logLH += gsl_vector_get(xbeta3, i);
            logLH += gsl_vector_get(V3, jj);
        }
        logLH += -gam * BpeMvnCorScr_wFunc(i, xbeta1, xbeta2, xbeta3, lambda1, lambda2, lambda3, jj, V1, V2, V3, s1, s2, s3, K1, K2, K3, survTime1, survTime2);
    }
    
    *val = logLH;
    
    return;
    
}





/* evaluating log-likelihood function for subject i */

/**/

void BpeMvnCorScr_logLH_i(int i,
                          gsl_vector *beta1,
                          gsl_vector *beta2,
                          gsl_vector *beta3,
                          gsl_vector *xbeta1,
                          gsl_vector *xbeta2,
                          gsl_vector *xbeta3,
                          gsl_vector *gamma,
                          gsl_vector *lambda1,
                          gsl_vector *lambda2,
                          gsl_vector *lambda3,
                          gsl_vector *s1,
                          gsl_vector *s2,
                          gsl_vector *s3,
                          gsl_vector *V1,
                          gsl_vector *V2,
                          gsl_vector *V3,
                          gsl_vector *survTime1,
                          gsl_vector *survTime2,
                          gsl_vector *survEvent1,
                          gsl_vector *case01,
                          gsl_vector *case11,
                          gsl_matrix *survCov1,
                          gsl_matrix *survCov2,
                          gsl_matrix *survCov3,
                          gsl_vector *cluster,
                          int K1,
                          int K2,
                          int K3,
                          double *val)
{
    double gam;
    double logLH = 0;
    


    
    int j, jj;
    
    gam = gsl_vector_get(gamma, i);
    
    jj = (int) gsl_vector_get(cluster, i) - 1;
    
    if(gsl_vector_get(survEvent1, i) == 1)
    {
        for(j = 0; j < K1+1; j++)
        {
            if(j == 0 && gsl_vector_get(survTime1, i) <= gsl_vector_get(s1, 0))
            {
                logLH += gsl_vector_get(lambda1, j);
            }
            if(j != 0 && gsl_vector_get(survTime1, i) > gsl_vector_get(s1, j-1) && gsl_vector_get(survTime1, i) <= gsl_vector_get(s1, j))
            {
                logLH += gsl_vector_get(lambda1, j);
            }
        }
        logLH += log(gam);
        logLH += gsl_vector_get(xbeta1, i);
        logLH += gsl_vector_get(V1, jj);
    }
    
    if(gsl_vector_get(case01, i) ==  1)
    {
        for(j = 0; j < K2+1; j++)
        {
            if(j == 0 && gsl_vector_get(survTime2, i) <= gsl_vector_get(s2, 0))
            {
                logLH += gsl_vector_get(lambda2, j);
            }
            if(j != 0 && gsl_vector_get(survTime2, i) > gsl_vector_get(s2, j-1) && gsl_vector_get(survTime2, i) <= gsl_vector_get(s2, j))
            {
                logLH += gsl_vector_get(lambda2, j);
            }
        }
        logLH += log(gam);
        logLH += gsl_vector_get(xbeta2, i);
        logLH += gsl_vector_get(V2, jj);
    }
    
    if(gsl_vector_get(case11, i) == 1)
    {
        for(j = 0; j < K3+1; j++)
        {
            if(j == 0 && gsl_vector_get(survTime2, i) <= gsl_vector_get(s3, 0))
            {
                logLH += gsl_vector_get(lambda3, j);
            }
            if(j != 0 && gsl_vector_get(survTime2, i) > gsl_vector_get(s3, j-1) && gsl_vector_get(survTime2, i) <= gsl_vector_get(s3, j))
            {
                logLH += gsl_vector_get(lambda3, j);
            }
        }
        logLH += log(gam);
        logLH += gsl_vector_get(xbeta3, i);
        logLH += gsl_vector_get(V3, jj);
    }
    
    logLH += -gam * BpeMvnCorScr_wFunc(i, xbeta1, xbeta2, xbeta3, lambda1, lambda2, lambda3, jj, V1, V2, V3, s1, s2, s3, K1, K2, K3, survTime1, survTime2);
    
    *val = logLH;
    
    return;
    
}





/* evaluating density function for subject i */

/**/

void BpeMvnCorScr_logf_i(int i,
                         gsl_vector *beta1,
                         gsl_vector *beta2,
                         gsl_vector *beta3,
                         gsl_vector *xbeta1,
                         gsl_vector *xbeta2,
                         gsl_vector *xbeta3,
                         gsl_vector *gamma,
                         gsl_vector *lambda1,
                         gsl_vector *lambda2,
                         gsl_vector *lambda3,
                         gsl_vector *s1,
                         gsl_vector *s2,
                         gsl_vector *s3,
                         gsl_vector *V1,
                         gsl_vector *V2,
                         gsl_vector *V3,
                         gsl_vector *survTime1,
                         gsl_vector *survTime2,
                         gsl_vector *survEvent1,
                         gsl_vector *case01,
                         gsl_vector *case11,
                         gsl_matrix *survCov1,
                         gsl_matrix *survCov2,
                         gsl_matrix *survCov3,
                         gsl_vector *cluster,
                         int K1,
                         int K2,
                         int K3,
                         double *val)
{
    double gam;
    double logf = 0;
    


    
    int j, jj;
    
    gam = gsl_vector_get(gamma, i);
    
    jj = (int) gsl_vector_get(cluster, i) - 1;
    
    
    if(gsl_vector_get(survTime1, i) != gsl_vector_get(survTime2, i))
    {
        for(j = 0; j < K1+1; j++)
        {
            if(j == 0 && gsl_vector_get(survTime1, i) <= gsl_vector_get(s1, 0))
            {
                logf += gsl_vector_get(lambda1, j);
            }
            if(j != 0 && gsl_vector_get(survTime1, i) > gsl_vector_get(s1, j-1) && gsl_vector_get(survTime1, i) <= gsl_vector_get(s1, j))
            {
                logf += gsl_vector_get(lambda1, j);
            }
        }
        logf += log(gam);
        logf += gsl_vector_get(xbeta1, i);
        logf += gsl_vector_get(V1, jj);
        
        for(j = 0; j < K3+1; j++)
        {
            if(j == 0 && gsl_vector_get(survTime2, i) <= gsl_vector_get(s3, 0))
            {
                logf += gsl_vector_get(lambda3, j);
            }
            if(j != 0 && gsl_vector_get(survTime2, i) > gsl_vector_get(s3, j-1) && gsl_vector_get(survTime2, i) <= gsl_vector_get(s3, j))
            {
                logf += gsl_vector_get(lambda3, j);
            }
        }
        logf += log(gam);
        logf += gsl_vector_get(xbeta3, i);
        logf += gsl_vector_get(V3, jj);
        
    }
    
    if(gsl_vector_get(survTime1, i) == gsl_vector_get(survTime2, i))
    {
        for(j = 0; j < K2+1; j++)
        {
            if(j == 0 && gsl_vector_get(survTime2, i) <= gsl_vector_get(s2, 0))
            {
                logf += gsl_vector_get(lambda2, j);
            }
            if(j != 0 && gsl_vector_get(survTime2, i) > gsl_vector_get(s2, j-1) && gsl_vector_get(survTime2, i) <= gsl_vector_get(s2, j))
            {
                logf += gsl_vector_get(lambda2, j);
            }
        }
        logf += log(gam);
        logf += gsl_vector_get(xbeta2, i);
        logf += gsl_vector_get(V2, jj);
    }
    
    logf += -gam * BpeMvnCorScr_wFunc(i, xbeta1, xbeta2, xbeta3, lambda1, lambda2, lambda3, jj, V1, V2, V3, s1, s2, s3, K1, K2, K3, survTime1, survTime2);
    
    *val = logf;
    
    return;
    
}








/*
 Density calculation for a multivariate normal distribution
 */
void c_dmvnormSH(gsl_vector *x,
                 double     mu,
                 double     sigma,
                 gsl_matrix *AInv,
                 double     *value)
{
    int signum, K = x->size;
	double sigmaSqInv = pow(sigma, -2);
    
	gsl_vector *muVec      = gsl_vector_alloc(K);
    gsl_vector *diff       = gsl_vector_alloc(K);
	gsl_matrix *SigmaInv   = gsl_matrix_alloc(K, K);
    gsl_matrix *SigmaInvLU = gsl_matrix_alloc(K, K);
    gsl_permutation *p     = gsl_permutation_alloc(K);
    
	gsl_vector_set_all(muVec, mu);
	gsl_vector_memcpy(diff, x);
	gsl_vector_sub(diff, muVec);
    
	gsl_matrix_memcpy(SigmaInv, AInv);
	gsl_matrix_scale(SigmaInv, sigmaSqInv);
    gsl_matrix_memcpy(SigmaInvLU, SigmaInv);
    gsl_linalg_LU_decomp(SigmaInvLU, p, &signum);
    
	c_quadform_vMv(diff, SigmaInv, value);
	*value = (log(gsl_linalg_LU_det(SigmaInvLU, signum)) - log(pow(2*Pi, K)) - *value) / 2;
    
	gsl_vector_free(muVec);
	gsl_vector_free(diff);
	gsl_matrix_free(SigmaInv);
	gsl_matrix_free(SigmaInvLU);
	gsl_permutation_free(p);
    return;
}




/*
 Evaluate w(y1, y2) function
 */
double BpeMvnCorScr_wFunc(int subjInx,
                          gsl_vector *xbeta1,
                          gsl_vector *xbeta2,
                          gsl_vector *xbeta3,
                          gsl_vector *lambda1,
                          gsl_vector *lambda2,
                          gsl_vector *lambda3,
                          int jj,
                          gsl_vector *V1,
                          gsl_vector *V2,
                          gsl_vector *V3,
                          gsl_vector *s1,
                          gsl_vector *s2,
                          gsl_vector *s3,
                          int J1,
                          int J2,
                          int J3,
                          gsl_vector *survTime1,
                          gsl_vector *survTime2)
{
    int i = subjInx;
    double cumHaz1, cumHaz2, cumHaz3diff;
    double Del, wVal;
    int j;
    
    cumHaz1 = 0; cumHaz2 = 0; cumHaz3diff = 0;
    
    for(j = 0; j < J1+1; j++)
    {
        if(j > 0)
        {
            Del = c_max(0, (c_min(gsl_vector_get(s1, j), gsl_vector_get(survTime1, i)) - gsl_vector_get(s1, j-1)));
        }
        if(j == 0)
        {
            Del = c_max(0, c_min(gsl_vector_get(s1, j), gsl_vector_get(survTime1, i)) - 0);
        }
        cumHaz1 += Del* exp(gsl_vector_get(lambda1, j)) * exp(gsl_vector_get(xbeta1, i)+gsl_vector_get(V1, jj));
    }
    
    for(j = 0; j < J2+1; j++)
    {
        if(j > 0)
        {
            Del = c_max(0, (c_min(gsl_vector_get(s2, j), gsl_vector_get(survTime1, i)) - gsl_vector_get(s2, j-1)));
        }
        if(j == 0)
        {
            Del = c_max(0, c_min(gsl_vector_get(s2, j), gsl_vector_get(survTime1, i)) - 0);
        }
        cumHaz2 += Del* exp(gsl_vector_get(lambda2, j)) * exp(gsl_vector_get(xbeta2, i)+gsl_vector_get(V2, jj));
    }
    
    for(j = 0; j < J3+1; j++)
    {
        if(j > 0)
        {
            Del = c_max(0, (c_min(gsl_vector_get(s3, j), gsl_vector_get(survTime2, i)) - c_max(gsl_vector_get(s3, j-1), gsl_vector_get(survTime1, i))));
        }
        if(j == 0)
        {
            Del = c_max(0, (c_min(gsl_vector_get(s3, j), gsl_vector_get(survTime2, i)) - c_max(0, gsl_vector_get(survTime1, i))));
        }
        cumHaz3diff += Del* exp(gsl_vector_get(lambda3, j)) * exp(gsl_vector_get(xbeta3, i)+gsl_vector_get(V3, jj));
    }
    
    
    
    /*
     printf("ch1 = %.3f\n", cumHaz1);
     printf("ch2 = %.3f\n", cumHaz2);
     printf("ch3 = %.3f\n", cumHaz3diff);
     
     */
    
    
    wVal = cumHaz1 + cumHaz2 + cumHaz3diff;
    
    return wVal;
}
















/********* For Weibull-DPM-SM model *************/



/* evaluating log-likelihood function */

/**/

void BweibDpCorScrSM_logMLH(gsl_vector *beta1,
                            gsl_vector *beta2,
                            gsl_vector *beta3,
                            double alpha1,
                            double alpha2,
                            double alpha3,
                            double kappa1,
                            double kappa2,
                            double kappa3,
                            double theta,
                            gsl_vector *V1,
                            gsl_vector *V2,
                            gsl_vector *V3,
                            gsl_vector *survTime1,
                            gsl_vector *survTime2,
                            gsl_vector *yStar,
                            gsl_vector *survEvent1,
                            gsl_vector *survEvent2,
                            gsl_vector *case01,
                            gsl_vector *case11,
                            gsl_matrix *survCov1,
                            gsl_matrix *survCov2,
                            gsl_matrix *survCov3,
                            gsl_vector *cluster,
                            double *val)
{
    double gfunc, LP1, LP2, LP3;
    double logLH = 0;
    

    int n = survTime1 -> size;
    
    int i, jj;
    
    for(i = 0; i < n; i++)
    {
        
        
        gsl_vector_view Xi1 = gsl_matrix_row(survCov1, i);
        gsl_blas_ddot(&Xi1.vector, beta1, &LP1);
        gsl_vector_view Xi2 = gsl_matrix_row(survCov2, i);
        gsl_blas_ddot(&Xi2.vector, beta2, &LP2);
        gsl_vector_view Xi3 = gsl_matrix_row(survCov3, i);
        gsl_blas_ddot(&Xi3.vector, beta3, &LP3);
        
        jj = (int) gsl_vector_get(cluster, i) - 1;
        
        if(gsl_vector_get(survEvent1, i) == 1)
        {
            logLH += log(alpha1) + log(kappa1) + (alpha1-1)*log(gsl_vector_get(survTime1, i));
            
            logLH += LP1;
            logLH += gsl_vector_get(V1, jj);
        }
        
        if(gsl_vector_get(case01, i) ==  1)
        {
            logLH += log(alpha2) + log(kappa2) + (alpha2-1)*log(gsl_vector_get(survTime2, i));
            
            logLH += LP2;
            logLH += gsl_vector_get(V2, jj);
        }
        
        if(gsl_vector_get(case11, i) == 1 && gsl_vector_get(yStar, i) != 0)
        {
            logLH += log(alpha3) + log(kappa3) + (alpha3-1)*log(gsl_vector_get(yStar, i));
            
            logLH += LP3;
            logLH += gsl_vector_get(V3, jj);
        }
        
        gfunc = BweibDpCorScrSM_wFunc_old(i, beta1, beta2, beta3, alpha1, alpha2, alpha3, kappa1, kappa2, kappa3, V1, V2, V3, survTime1, yStar, cluster, survCov1, survCov2, survCov3);
        
        logLH += (- 1/theta - gsl_vector_get(survEvent1, i) - gsl_vector_get(survEvent2, i)) * log(1 + theta * gfunc);
    }
    
    *val = logLH;
    
    return;
    
    
}





/* evaluating log-likelihood function for subject i */

/**/

void BweibDpCorScrSM_logMLH_i(int i,
                              gsl_vector *beta1,
                              gsl_vector *beta2,
                              gsl_vector *beta3,
                              double alpha1,
                              double alpha2,
                              double alpha3,
                              double kappa1,
                              double kappa2,
                              double kappa3,
                              double theta,
                              gsl_vector *V1,
                              gsl_vector *V2,
                              gsl_vector *V3,
                              gsl_vector *survTime1,
                              gsl_vector *survTime2,
                              gsl_vector *yStar,
                              gsl_vector *survEvent1,
                              gsl_vector *survEvent2,
                              gsl_vector *case01,
                              gsl_vector *case11,
                              gsl_matrix *survCov1,
                              gsl_matrix *survCov2,
                              gsl_matrix *survCov3,
                              gsl_vector *cluster,
                              double *val)
{
    double gfunc, LP1, LP2, LP3;
    double logLH = 0;
    


    
    int jj;
    
    
    
    gsl_vector_view Xi1 = gsl_matrix_row(survCov1, i);
    gsl_blas_ddot(&Xi1.vector, beta1, &LP1);
    gsl_vector_view Xi2 = gsl_matrix_row(survCov2, i);
    gsl_blas_ddot(&Xi2.vector, beta2, &LP2);
    gsl_vector_view Xi3 = gsl_matrix_row(survCov3, i);
    gsl_blas_ddot(&Xi3.vector, beta3, &LP3);
    
    jj = (int) gsl_vector_get(cluster, i) - 1;
    
    if(gsl_vector_get(survEvent1, i) == 1)
    {
        logLH += log(alpha1) + log(kappa1) + (alpha1-1)*log(gsl_vector_get(survTime1, i));
        
        logLH += LP1;
        logLH += gsl_vector_get(V1, jj);
    }
    
    if(gsl_vector_get(case01, i) ==  1)
    {
        logLH += log(alpha2) + log(kappa2) + (alpha2-1)*log(gsl_vector_get(survTime2, i));
        
        logLH += LP2;
        logLH += gsl_vector_get(V2, jj);
    }
    
    if(gsl_vector_get(case11, i) == 1 && gsl_vector_get(yStar, i) != 0)
    {
        logLH += log(alpha3) + log(kappa3) + (alpha3-1)*log(gsl_vector_get(yStar, i));
        
        logLH += LP3;
        logLH += gsl_vector_get(V3, jj);
    }
    
    gfunc = BweibDpCorScrSM_wFunc_old(i, beta1, beta2, beta3, alpha1, alpha2, alpha3, kappa1, kappa2, kappa3, V1, V2, V3, survTime1, yStar, cluster, survCov1, survCov2, survCov3);
    
    logLH += (- 1/theta - gsl_vector_get(survEvent1, i) - gsl_vector_get(survEvent2, i)) * log(1 + theta * gfunc);
    
    *val = logLH;
    
    return;
    
    
}






/* evaluating log-likelihood function */

/**/

void BweibDpCorScrSM_logLH(gsl_vector *beta1,
                           gsl_vector *beta2,
                           gsl_vector *beta3,
                           double alpha1,
                           double alpha2,
                           double alpha3,
                           double kappa1,
                           double kappa2,
                           double kappa3,
                           gsl_vector *gamma,
                           gsl_vector *V1,
                           gsl_vector *V2,
                           gsl_vector *V3,
                           gsl_vector *survTime1,
                           gsl_vector *survTime2,
                           gsl_vector *yStar,
                           gsl_vector *survEvent1,
                           gsl_vector *case01,
                           gsl_vector *case11,
                           gsl_matrix *survCov1,
                           gsl_matrix *survCov2,
                           gsl_matrix *survCov3,
                           gsl_vector *cluster,
                           double *val)
{
    double gam, LP1, LP2, LP3;
    double logLH = 0;
    

    int n = survTime1 -> size;
    
    int i, jj;
    
    for(i = 0; i < n; i++)
    {
        gam = gsl_vector_get(gamma, i);
        
        gsl_vector_view Xi1 = gsl_matrix_row(survCov1, i);
        gsl_blas_ddot(&Xi1.vector, beta1, &LP1);
        gsl_vector_view Xi2 = gsl_matrix_row(survCov2, i);
        gsl_blas_ddot(&Xi2.vector, beta2, &LP2);
        gsl_vector_view Xi3 = gsl_matrix_row(survCov3, i);
        gsl_blas_ddot(&Xi3.vector, beta3, &LP3);
        
        jj = (int) gsl_vector_get(cluster, i) - 1;
        
        if(gsl_vector_get(survEvent1, i) == 1)
        {
            logLH += log(alpha1) + log(kappa1) + (alpha1-1)*log(gsl_vector_get(survTime1, i));
            logLH += log(gam);
            logLH += LP1;
            logLH += gsl_vector_get(V1, jj);
        }
        
        if(gsl_vector_get(case01, i) ==  1)
        {
            logLH += log(alpha2) + log(kappa2) + (alpha2-1)*log(gsl_vector_get(survTime2, i));
            logLH += log(gam);
            logLH += LP2;
            logLH += gsl_vector_get(V2, jj);
        }
        
        if(gsl_vector_get(case11, i) == 1 && gsl_vector_get(yStar, i) != 0)
        {
            logLH += log(alpha3) + log(kappa3) + (alpha3-1)*log(gsl_vector_get(yStar, i));
            logLH += log(gam);
            logLH += LP3;
            logLH += gsl_vector_get(V3, jj);
        }
        logLH += -gam * BweibDpCorScrSM_wFunc_old(i, beta1, beta2, beta3, alpha1, alpha2, alpha3, kappa1, kappa2, kappa3, V1, V2, V3, survTime1, yStar, cluster, survCov1, survCov2, survCov3);
    }
    
    *val = logLH;
    
    return;
    
    
}





/* evaluating log-likelihood function for subject i */

/**/

void BweibDpCorScrSM_logLH_i(int i,
                             gsl_vector *beta1,
                             gsl_vector *beta2,
                             gsl_vector *beta3,
                             double alpha1,
                             double alpha2,
                             double alpha3,
                             double kappa1,
                             double kappa2,
                             double kappa3,
                             gsl_vector *gamma,
                             gsl_vector *V1,
                             gsl_vector *V2,
                             gsl_vector *V3,
                             gsl_vector *survTime1,
                             gsl_vector *survTime2,
                             gsl_vector *yStar,
                             gsl_vector *survEvent1,
                             gsl_vector *case01,
                             gsl_vector *case11,
                             gsl_matrix *survCov1,
                             gsl_matrix *survCov2,
                             gsl_matrix *survCov3,
                             gsl_vector *cluster,
                             double *val)
{
    double gam, LP1, LP2, LP3;
    double logLH = 0;
    


    
    int jj;
    
    gam = gsl_vector_get(gamma, i);
    
    gsl_vector_view Xi1 = gsl_matrix_row(survCov1, i);
    gsl_blas_ddot(&Xi1.vector, beta1, &LP1);
    gsl_vector_view Xi2 = gsl_matrix_row(survCov2, i);
    gsl_blas_ddot(&Xi2.vector, beta2, &LP2);
    gsl_vector_view Xi3 = gsl_matrix_row(survCov3, i);
    gsl_blas_ddot(&Xi3.vector, beta3, &LP3);
    
    jj = (int) gsl_vector_get(cluster, i) - 1;
    
    if(gsl_vector_get(survEvent1, i) == 1)
    {
        logLH += log(alpha1) + log(kappa1) + (alpha1-1)*log(gsl_vector_get(survTime1, i));
        logLH += log(gam);
        logLH += LP1;
        logLH += gsl_vector_get(V1, jj);
    }
    
    if(gsl_vector_get(case01, i) ==  1)
    {
        logLH += log(alpha2) + log(kappa2) + (alpha2-1)*log(gsl_vector_get(survTime2, i));
        logLH += log(gam);
        logLH += LP2;
        logLH += gsl_vector_get(V2, jj);
    }
    
    if(gsl_vector_get(case11, i) == 1 && gsl_vector_get(yStar, i) != 0)
    {
        logLH += log(alpha3) + log(kappa3) + (alpha3-1)*log(gsl_vector_get(yStar, i));
        logLH += log(gam);
        logLH += LP3;
        logLH += gsl_vector_get(V3, jj);
    }
    logLH += -gam * BweibDpCorScrSM_wFunc_old(i, beta1, beta2, beta3, alpha1, alpha2, alpha3, kappa1, kappa2, kappa3, V1, V2, V3, survTime1, yStar, cluster, survCov1, survCov2, survCov3);
    
    
    *val = logLH;
    
    return;
    
    
}










/*
 Evaluate w(y1, y2) function with nu2 and nu3
 */
double BweibDpCorScrSM_wFunc(int subjInx,
                             gsl_vector *beta1,
                             gsl_vector *beta2,
                             gsl_vector *beta3,
                             double alpha1,
                             double alpha2,
                             double alpha3,
                             double kappa1,
                             double kappa2,
                             double kappa3,
                             double nu2,
                             double nu3,
                             double gam,
                             gsl_vector *V1,
                             gsl_vector *V2,
                             gsl_vector *V3,
                             gsl_vector *survTime1,
                             gsl_vector *survTime2,
                             gsl_vector *cluster,
                             gsl_matrix *survCov1,
                             gsl_matrix *survCov2,
                             gsl_matrix *survCov3)
{
    int i = subjInx;
    double LP1, LP2, LP3, cumHaz1, cumHaz2, cumHaz3diff;
    double wVal;
    int jj = (int) gsl_vector_get(cluster, i) - 1;
    
    gsl_vector_view Xi1 = gsl_matrix_row(survCov1, i);
    gsl_blas_ddot(&Xi1.vector, beta1 ,&LP1);
    
    gsl_vector_view Xi2 = gsl_matrix_row(survCov2, i);
    gsl_blas_ddot(&Xi2.vector, beta2 ,&LP2);
    
    gsl_vector_view Xi3 = gsl_matrix_row(survCov3, i);
    gsl_blas_ddot(&Xi3.vector, beta3 ,&LP3);
    
    cumHaz1     = kappa1 * pow(gsl_vector_get(survTime1, i), alpha1) * exp(LP1 + gsl_vector_get(V1, jj));
    cumHaz2     = kappa2 * pow(gsl_vector_get(survTime1, i), alpha2) * exp(LP2 + gsl_vector_get(V2, jj));
    cumHaz3diff = kappa3 * (pow(gsl_vector_get(survTime2, i), alpha3)-pow(gsl_vector_get(survTime1, i), alpha3)) * exp(LP3 + gsl_vector_get(V3, jj));
    
    wVal = gam * cumHaz1 + pow(gam, nu2) * cumHaz2 + pow(gam, nu3) * cumHaz3diff;
    
    return wVal;
}



/*
 Evaluate w(y1, y2) function
 */
double BweibDpCorScrSM_wFunc_old(int subjInx,
                                 gsl_vector *beta1,
                                 gsl_vector *beta2,
                                 gsl_vector *beta3,
                                 double alpha1,
                                 double alpha2,
                                 double alpha3,
                                 double kappa1,
                                 double kappa2,
                                 double kappa3,
                                 gsl_vector *V1,
                                 gsl_vector *V2,
                                 gsl_vector *V3,
                                 gsl_vector *survTime1,
                                 gsl_vector *yStar,
                                 gsl_vector *cluster,
                                 gsl_matrix *survCov1,
                                 gsl_matrix *survCov2,
                                 gsl_matrix *survCov3)
{
    int i = subjInx;
    int j = (int) gsl_vector_get(cluster, i) - 1;
    double LP1, LP2, LP3, cumHaz1, cumHaz2, cumHaz3diff;
    double wVal;
    
    gsl_vector_view Xi1 = gsl_matrix_row(survCov1, i);
    gsl_blas_ddot(&Xi1.vector, beta1 ,&LP1);
    
    gsl_vector_view Xi2 = gsl_matrix_row(survCov2, i);
    gsl_blas_ddot(&Xi2.vector, beta2 ,&LP2);
    
    gsl_vector_view Xi3 = gsl_matrix_row(survCov3, i);
    gsl_blas_ddot(&Xi3.vector, beta3 ,&LP3);
    
    cumHaz1     = kappa1 * pow(gsl_vector_get(survTime1, i), alpha1) * exp(LP1 + gsl_vector_get(V1, j));
    cumHaz2     = kappa2 * pow(gsl_vector_get(survTime1, i), alpha2) * exp(LP2 + gsl_vector_get(V2, j));
    cumHaz3diff = kappa3 * pow(gsl_vector_get(yStar, i), alpha3) * exp(LP3 + gsl_vector_get(V3, j));
    
    wVal = cumHaz1 + cumHaz2 + cumHaz3diff;
    
    return wVal;
}









/********* For Weibull-DPM-M model *************/



/* evaluating log-likelihood function */

/**/

void BweibDpCorScr_logLH(gsl_vector *beta1,
                         gsl_vector *beta2,
                         gsl_vector *beta3,
                         double alpha1,
                         double alpha2,
                         double alpha3,
                         double kappa1,
                         double kappa2,
                         double kappa3,
                         gsl_vector *gamma,
                         gsl_vector *V1,
                         gsl_vector *V2,
                         gsl_vector *V3,
                         gsl_vector *survTime1,
                         gsl_vector *survTime2,
                         gsl_vector *survEvent1,
                         gsl_vector *case01,
                         gsl_vector *case11,
                         gsl_matrix *survCov1,
                         gsl_matrix *survCov2,
                         gsl_matrix *survCov3,
                         gsl_vector *cluster,
                         double *val)
{
    double gam, LP1, LP2, LP3;
    double logLH = 0;
    

    int n = survTime1 -> size;
    
    int i, jj;
    
    for(i = 0; i < n; i++)
    {
        gam = gsl_vector_get(gamma, i);
        
        gsl_vector_view Xi1 = gsl_matrix_row(survCov1, i);
        gsl_blas_ddot(&Xi1.vector, beta1, &LP1);
        gsl_vector_view Xi2 = gsl_matrix_row(survCov2, i);
        gsl_blas_ddot(&Xi2.vector, beta2, &LP2);
        gsl_vector_view Xi3 = gsl_matrix_row(survCov3, i);
        gsl_blas_ddot(&Xi3.vector, beta3, &LP3);
        
        jj = (int) gsl_vector_get(cluster, i) - 1;
        
        if(gsl_vector_get(survEvent1, i) == 1)
        {
            logLH += log(alpha1) + log(kappa1) + (alpha1-1)*log(gsl_vector_get(survTime1, i));
            logLH += log(gam);
            logLH += LP1;
            logLH += gsl_vector_get(V1, jj);
        }
        
        if(gsl_vector_get(case01, i) ==  1)
        {
            logLH += log(alpha2) + log(kappa2) + (alpha2-1)*log(gsl_vector_get(survTime2, i));
            logLH += log(gam);
            logLH += LP2;
            logLH += gsl_vector_get(V2, jj);
        }
        
        if(gsl_vector_get(case11, i) == 1)
        {
            logLH += log(alpha3) + log(kappa3) + (alpha3-1)*log(gsl_vector_get(survTime2, i));
            logLH += log(gam);
            logLH += LP3;
            logLH += gsl_vector_get(V3, jj);
        }
        logLH += -gam * BweibDpCorScr_wFunc_old(i, beta1, beta2, beta3, alpha1, alpha2, alpha3, kappa1, kappa2, kappa3, V1, V2, V3, survTime1, survTime2, cluster, survCov1, survCov2, survCov3);
    }
    
    *val = logLH;
    
    return;
    
    
}





/* evaluating log-likelihood function for subject i */

/**/

void BweibDpCorScr_logLH_i(int i,
                           gsl_vector *beta1,
                           gsl_vector *beta2,
                           gsl_vector *beta3,
                           double alpha1,
                           double alpha2,
                           double alpha3,
                           double kappa1,
                           double kappa2,
                           double kappa3,
                           gsl_vector *gamma,
                           gsl_vector *V1,
                           gsl_vector *V2,
                           gsl_vector *V3,
                           gsl_vector *survTime1,
                           gsl_vector *survTime2,
                           gsl_vector *survEvent1,
                           gsl_vector *case01,
                           gsl_vector *case11,
                           gsl_matrix *survCov1,
                           gsl_matrix *survCov2,
                           gsl_matrix *survCov3,
                           gsl_vector *cluster,
                           double *val)
{
    double gam, LP1, LP2, LP3;
    double logLH = 0;
    


    
    int jj;
    
    gam = gsl_vector_get(gamma, i);
    
    gsl_vector_view Xi1 = gsl_matrix_row(survCov1, i);
    gsl_blas_ddot(&Xi1.vector, beta1, &LP1);
    gsl_vector_view Xi2 = gsl_matrix_row(survCov2, i);
    gsl_blas_ddot(&Xi2.vector, beta2, &LP2);
    gsl_vector_view Xi3 = gsl_matrix_row(survCov3, i);
    gsl_blas_ddot(&Xi3.vector, beta3, &LP3);
    
    jj = (int) gsl_vector_get(cluster, i) - 1;
    
    if(gsl_vector_get(survEvent1, i) == 1)
    {
        logLH += log(alpha1) + log(kappa1) + (alpha1-1)*log(gsl_vector_get(survTime1, i));
        logLH += log(gam);
        logLH += LP1;
        logLH += gsl_vector_get(V1, jj);
    }
    
    if(gsl_vector_get(case01, i) ==  1)
    {
        logLH += log(alpha2) + log(kappa2) + (alpha2-1)*log(gsl_vector_get(survTime2, i));
        logLH += log(gam);
        logLH += LP2;
        logLH += gsl_vector_get(V2, jj);
    }
    
    if(gsl_vector_get(case11, i) == 1)
    {
        logLH += log(alpha3) + log(kappa3) + (alpha3-1)*log(gsl_vector_get(survTime2, i));
        logLH += log(gam);
        logLH += LP3;
        logLH += gsl_vector_get(V3, jj);
    }
    logLH += -gam * BweibDpCorScr_wFunc_old(i, beta1, beta2, beta3, alpha1, alpha2, alpha3, kappa1, kappa2, kappa3, V1, V2, V3, survTime1, survTime2, cluster, survCov1, survCov2, survCov3);
    
    
    *val = logLH;
    
    return;
    
    
}





/* evaluating log-likelihood function */

/**/

void BweibDpCorScr_logMLH(gsl_vector *beta1,
                          gsl_vector *beta2,
                          gsl_vector *beta3,
                          double alpha1,
                          double alpha2,
                          double alpha3,
                          double kappa1,
                          double kappa2,
                          double kappa3,
                          double theta,
                          gsl_vector *V1,
                          gsl_vector *V2,
                          gsl_vector *V3,
                          gsl_vector *survTime1,
                          gsl_vector *survTime2,
                          gsl_vector *survEvent1,
                          gsl_vector *survEvent2,
                          gsl_vector *case01,
                          gsl_vector *case11,
                          gsl_matrix *survCov1,
                          gsl_matrix *survCov2,
                          gsl_matrix *survCov3,
                          gsl_vector *cluster,
                          double *val)
{
    double gfunc, LP1, LP2, LP3;
    double logLH = 0;
    

    int n = survTime1 -> size;
    
    int i, jj;
    
    for(i = 0; i < n; i++)
    {
        gsl_vector_view Xi1 = gsl_matrix_row(survCov1, i);
        gsl_blas_ddot(&Xi1.vector, beta1, &LP1);
        gsl_vector_view Xi2 = gsl_matrix_row(survCov2, i);
        gsl_blas_ddot(&Xi2.vector, beta2, &LP2);
        gsl_vector_view Xi3 = gsl_matrix_row(survCov3, i);
        gsl_blas_ddot(&Xi3.vector, beta3, &LP3);
        
        jj = (int) gsl_vector_get(cluster, i) - 1;
        
        if(gsl_vector_get(survEvent1, i) == 1)
        {
            logLH += log(alpha1) + log(kappa1) + (alpha1-1)*log(gsl_vector_get(survTime1, i));
            logLH += LP1;
            logLH += gsl_vector_get(V1, jj);
        }
        
        if(gsl_vector_get(case01, i) ==  1)
        {
            logLH += log(alpha2) + log(kappa2) + (alpha2-1)*log(gsl_vector_get(survTime2, i));
            logLH += LP2;
            logLH += gsl_vector_get(V2, jj);
        }
        
        if(gsl_vector_get(case11, i) == 1)
        {
            logLH += log(alpha3) + log(kappa3) + (alpha3-1)*log(gsl_vector_get(survTime2, i));
            logLH += LP3;
            logLH += gsl_vector_get(V3, jj);
        }
        gfunc = BweibDpCorScr_wFunc_old(i, beta1, beta2, beta3, alpha1, alpha2, alpha3, kappa1, kappa2, kappa3, V1, V2, V3, survTime1, survTime2, cluster, survCov1, survCov2, survCov3);
        
        logLH += (- 1/theta - gsl_vector_get(survEvent1, i) - gsl_vector_get(survEvent2, i)) * log(1 + theta * gfunc);
    }
    
    *val = logLH;
    
    return;
    
    
}





/* evaluating log-likelihood function for subject i */

/**/

void BweibDpCorScr_logMLH_i(int i,
                            gsl_vector *beta1,
                            gsl_vector *beta2,
                            gsl_vector *beta3,
                            double alpha1,
                            double alpha2,
                            double alpha3,
                            double kappa1,
                            double kappa2,
                            double kappa3,
                            double theta,
                            gsl_vector *V1,
                            gsl_vector *V2,
                            gsl_vector *V3,
                            gsl_vector *survTime1,
                            gsl_vector *survTime2,
                            gsl_vector *survEvent1,
                            gsl_vector *survEvent2,
                            gsl_vector *case01,
                            gsl_vector *case11,
                            gsl_matrix *survCov1,
                            gsl_matrix *survCov2,
                            gsl_matrix *survCov3,
                            gsl_vector *cluster,
                            double *val)
{
    double gfunc, LP1, LP2, LP3;
    double logLH = 0;
    


    
    int jj;
    
    gsl_vector_view Xi1 = gsl_matrix_row(survCov1, i);
    gsl_blas_ddot(&Xi1.vector, beta1, &LP1);
    gsl_vector_view Xi2 = gsl_matrix_row(survCov2, i);
    gsl_blas_ddot(&Xi2.vector, beta2, &LP2);
    gsl_vector_view Xi3 = gsl_matrix_row(survCov3, i);
    gsl_blas_ddot(&Xi3.vector, beta3, &LP3);
    
    jj = (int) gsl_vector_get(cluster, i) - 1;
    
    if(gsl_vector_get(survEvent1, i) == 1)
    {
        logLH += log(alpha1) + log(kappa1) + (alpha1-1)*log(gsl_vector_get(survTime1, i));
        logLH += LP1;
        logLH += gsl_vector_get(V1, jj);
    }
    
    if(gsl_vector_get(case01, i) ==  1)
    {
        logLH += log(alpha2) + log(kappa2) + (alpha2-1)*log(gsl_vector_get(survTime2, i));
        logLH += LP2;
        logLH += gsl_vector_get(V2, jj);
    }
    
    if(gsl_vector_get(case11, i) == 1)
    {
        logLH += log(alpha3) + log(kappa3) + (alpha3-1)*log(gsl_vector_get(survTime2, i));
        logLH += LP3;
        logLH += gsl_vector_get(V3, jj);
    }
    
    gfunc = BweibDpCorScr_wFunc_old(i, beta1, beta2, beta3, alpha1, alpha2, alpha3, kappa1, kappa2, kappa3, V1, V2, V3, survTime1, survTime2, cluster, survCov1, survCov2, survCov3);
    
    logLH += (- 1/theta - gsl_vector_get(survEvent1, i) - gsl_vector_get(survEvent2, i)) * log(1 + theta * gfunc);
    
    *val = logLH;
    
    return;
    
    
}




int c_multinom_sample(gsl_rng *rr,
                      gsl_vector *prob,
                      int length_prob)
{
    int ii, val;
    int KK = length_prob;
    double probK[KK];
    
    for(ii = 0; ii < KK; ii++)
    {
        probK[ii] = gsl_vector_get(prob, ii);
    }
    
    unsigned int samples[KK];
    
    gsl_ran_multinomial(rr, KK, 1, probK, samples);
    
    for(ii = 0; ii < KK; ii++)
    {
        if(samples[ii] == 1) val = ii + 1;
    }
    
    return val;
}


double Qfunc(gsl_vector *V,
             gsl_vector *mu0,
             double zeta0,
             gsl_matrix *Psi0,
             double rho0)
{
    double expr1, expr2num, expr2den, expr3num, expr3den, val2;
    
    gsl_matrix *expr = gsl_matrix_calloc(3,3);
    gsl_matrix *temp_expr = gsl_matrix_calloc(3,3);
    gsl_vector *temp_vec = gsl_vector_calloc(3);
    
    expr1 = pow(Pi*sqrt(2*(1+zeta0)), -3);
    c_mgamma3(rho0/2+0.5, &expr2num);
    c_mgamma3(rho0/2, &expr2den);
    c_det(Psi0, &expr3num);
    expr3num = rho0/2 * log(expr3num);
    
    gsl_matrix_memcpy(expr, Psi0);
    gsl_blas_dger(1, V, V, temp_expr);
    gsl_matrix_add(expr, temp_expr);
    gsl_matrix_set_zero(temp_expr);
    gsl_blas_dger(1, mu0, mu0, temp_expr);
    gsl_matrix_scale(temp_expr, 1/zeta0);
    gsl_matrix_add(expr, temp_expr);
    gsl_matrix_set_zero(temp_expr);
    
    gsl_vector_memcpy(temp_vec, mu0);
    gsl_vector_scale(temp_vec, 1/zeta0);
    gsl_vector_add(temp_vec, V);
    gsl_blas_dger(1, temp_vec, temp_vec, temp_expr);
    gsl_matrix_scale(temp_expr, (double) 1/(1+1/zeta0));
    gsl_matrix_sub(expr, temp_expr);
    
    c_det(expr, &expr3den);
    expr3den = (rho0/2 + 0.5) * log(expr3den);
    
    val2 = expr1 * exp(expr2num - expr2den) * exp(expr3num - expr3den);
    
    gsl_matrix_free(expr);
    gsl_matrix_free(temp_expr);
    gsl_vector_free(temp_vec);
    
    return val2;
}




void c_det(gsl_matrix *A,
           double *result)
{
    int signum, K = A->size1;
    gsl_matrix *matLU = gsl_matrix_alloc(K, K);
    gsl_permutation *permu    = gsl_permutation_alloc(K);
    gsl_matrix_memcpy(matLU, A);
    
    gsl_linalg_LU_decomp(matLU, permu, &signum);
    
    *result = gsl_linalg_LU_det(matLU, signum);
    
    gsl_matrix_free(matLU);
    gsl_permutation_free(permu);
}


/* multivarite log - gamma function */

void c_mgamma3(double a, double *val)
{
    *val = 1.5 * log(Pi) + lgamma(a) + lgamma(a-0.5) + lgamma(a-1);
    
    return;
}




/*
 Random number generation for multivariate normal distribution
 mean (n)
 Var (n x n)
 sample (numSpl x n)
 */


void c_rmvnorm(gsl_matrix *sample,
               gsl_vector *mean,
               gsl_matrix *Var)
{
    int n = sample->size2;
    int numSpl = sample->size1;
    int i, j;
    double spl;
    
    gsl_matrix *temp = gsl_matrix_alloc(n, n);
    
    gsl_matrix_memcpy(temp, Var);
    gsl_linalg_cholesky_decomp(temp);
    
    for(i = 0; i < n; i ++){
        for(j = 0; j < n; j++){
            if(i > j){
                gsl_matrix_set(temp, i, j, 0);
            }
        }
    }
    
    for(i = 0; i < numSpl; i ++){
        for(j = 0; j < n; j ++){
            spl = rnorm(0, 1);
            gsl_matrix_set(sample, i, j, spl);
        }
    }
    
    gsl_blas_dtrmm(CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, 1, temp, sample);
    
    for(i = 0; i < numSpl; i++){
        gsl_vector_view sampleRow = gsl_matrix_row(sample, i);
        gsl_vector_add(&sampleRow.vector, mean);
    }
    
    gsl_matrix_free(temp);
    
    return;
}





/*
 Evaluate w(y1, y2) function with nu2 and nu3
 */
double BweibDpCorScr_wFunc(int subjInx,
                           gsl_vector *beta1,
                           gsl_vector *beta2,
                           gsl_vector *beta3,
                           double alpha1,
                           double alpha2,
                           double alpha3,
                           double kappa1,
                           double kappa2,
                           double kappa3,
                           double nu2,
                           double nu3,
                           double gam,
                           gsl_vector *V1,
                           gsl_vector *V2,
                           gsl_vector *V3,
                           gsl_vector *survTime1,
                           gsl_vector *survTime2,
                           gsl_vector *cluster,
                           gsl_matrix *survCov1,
                           gsl_matrix *survCov2,
                           gsl_matrix *survCov3)
{
    int i = subjInx;
    double LP1, LP2, LP3, cumHaz1, cumHaz2, cumHaz3diff;
    double wVal;
    int jj = (int) gsl_vector_get(cluster, i) - 1;
    
    gsl_vector_view Xi1 = gsl_matrix_row(survCov1, i);
    gsl_blas_ddot(&Xi1.vector, beta1 ,&LP1);
    
    gsl_vector_view Xi2 = gsl_matrix_row(survCov2, i);
    gsl_blas_ddot(&Xi2.vector, beta2 ,&LP2);
    
    gsl_vector_view Xi3 = gsl_matrix_row(survCov3, i);
    gsl_blas_ddot(&Xi3.vector, beta3 ,&LP3);
    
    cumHaz1     = kappa1 * pow(gsl_vector_get(survTime1, i), alpha1) * exp(LP1 + gsl_vector_get(V1, jj));
    cumHaz2     = kappa2 * pow(gsl_vector_get(survTime1, i), alpha2) * exp(LP2 + gsl_vector_get(V2, jj));
    cumHaz3diff = kappa3 * (pow(gsl_vector_get(survTime2, i), alpha3)-pow(gsl_vector_get(survTime1, i), alpha3)) * exp(LP3 + gsl_vector_get(V3, jj));
    
    wVal = gam * cumHaz1 + pow(gam, nu2) * cumHaz2 + pow(gam, nu3) * cumHaz3diff;
    
    return wVal;
}



/*
 Evaluate w(y1, y2) function
 */
double BweibDpCorScr_wFunc_old(int subjInx,
                               gsl_vector *beta1,
                               gsl_vector *beta2,
                               gsl_vector *beta3,
                               double alpha1,
                               double alpha2,
                               double alpha3,
                               double kappa1,
                               double kappa2,
                               double kappa3,
                               gsl_vector *V1,
                               gsl_vector *V2,
                               gsl_vector *V3,
                               gsl_vector *survTime1,
                               gsl_vector *survTime2,
                               gsl_vector *cluster,
                               gsl_matrix *survCov1,
                               gsl_matrix *survCov2,
                               gsl_matrix *survCov3)
{
    int i = subjInx;
    int j = (int) gsl_vector_get(cluster, i) - 1;
    double LP1, LP2, LP3, cumHaz1, cumHaz2, cumHaz3diff;
    double wVal;
    
    gsl_vector_view Xi1 = gsl_matrix_row(survCov1, i);
    gsl_blas_ddot(&Xi1.vector, beta1 ,&LP1);
    
    gsl_vector_view Xi2 = gsl_matrix_row(survCov2, i);
    gsl_blas_ddot(&Xi2.vector, beta2 ,&LP2);
    
    gsl_vector_view Xi3 = gsl_matrix_row(survCov3, i);
    gsl_blas_ddot(&Xi3.vector, beta3 ,&LP3);
    
    cumHaz1     = kappa1 * pow(gsl_vector_get(survTime1, i), alpha1) * exp(LP1 + gsl_vector_get(V1, j));
    cumHaz2     = kappa2 * pow(gsl_vector_get(survTime1, i), alpha2) * exp(LP2 + gsl_vector_get(V2, j));
    cumHaz3diff = kappa3 * (pow(gsl_vector_get(survTime2, i), alpha3)-pow(gsl_vector_get(survTime1, i), alpha3)) * exp(LP3 + gsl_vector_get(V3, j));
    
    wVal = cumHaz1 + cumHaz2 + cumHaz3diff;
    
    return wVal;
}
















/********* For Weibull-MVN-SM model *************/


/* evaluating log-likelihood function */

/**/

void BweibMvnCorScrSM_logMLH(gsl_vector *beta1,
                             gsl_vector *beta2,
                             gsl_vector *beta3,
                             double alpha1,
                             double alpha2,
                             double alpha3,
                             double kappa1,
                             double kappa2,
                             double kappa3,
                             double theta,
                             gsl_vector *V1,
                             gsl_vector *V2,
                             gsl_vector *V3,
                             gsl_vector *survTime1,
                             gsl_vector *survTime2,
                             gsl_vector *yStar,
                             gsl_vector *survEvent1,
                             gsl_vector *survEvent2,
                             gsl_vector *case01,
                             gsl_vector *case11,
                             gsl_matrix *survCov1,
                             gsl_matrix *survCov2,
                             gsl_matrix *survCov3,
                             gsl_vector *cluster,
                             double *val)
{
    double gfunc, LP1, LP2, LP3;
    double logLH = 0;
    

    int n = survTime1 -> size;
    
    int i, jj;
    
    for(i = 0; i < n; i++)
    {
        
        gsl_vector_view Xi1 = gsl_matrix_row(survCov1, i);
        gsl_blas_ddot(&Xi1.vector, beta1, &LP1);
        gsl_vector_view Xi2 = gsl_matrix_row(survCov2, i);
        gsl_blas_ddot(&Xi2.vector, beta2, &LP2);
        gsl_vector_view Xi3 = gsl_matrix_row(survCov3, i);
        gsl_blas_ddot(&Xi3.vector, beta3, &LP3);
        
        jj = (int) gsl_vector_get(cluster, i) - 1;
        
        if(gsl_vector_get(survEvent1, i) == 1)
        {
            logLH += log(alpha1) + log(kappa1) + (alpha1-1)*log(gsl_vector_get(survTime1, i));
            logLH += LP1;
            logLH += gsl_vector_get(V1, jj);
        }
        
        if(gsl_vector_get(case01, i) ==  1)
        {
            logLH += log(alpha2) + log(kappa2) + (alpha2-1)*log(gsl_vector_get(survTime2, i));
            logLH += LP2;
            logLH += gsl_vector_get(V2, jj);
        }
        
        if(gsl_vector_get(case11, i) == 1 && gsl_vector_get(yStar, i) != 0)
        {
            logLH += log(alpha3) + log(kappa3) + (alpha3-1)*log(gsl_vector_get(yStar, i));
            logLH += LP3;
            logLH += gsl_vector_get(V3, jj);
        }
        gfunc = BweibMvnCorScrSM_wFunc_old(i, beta1, beta2, beta3, alpha1, alpha2, alpha3, kappa1, kappa2, kappa3, V1, V2, V3, survTime1, yStar, cluster, survCov1, survCov2, survCov3);
        
        logLH += (- 1/theta - gsl_vector_get(survEvent1, i) - gsl_vector_get(survEvent2, i)) * log(1 + theta * gfunc);
    }
    
    *val = logLH;
    
    return;
    
    
}





/* evaluating log-likelihood function for subject i */

/**/

void BweibMvnCorScrSM_logMLH_i(int i,
                               gsl_vector *beta1,
                               gsl_vector *beta2,
                               gsl_vector *beta3,
                               double alpha1,
                               double alpha2,
                               double alpha3,
                               double kappa1,
                               double kappa2,
                               double kappa3,
                               double theta,
                               gsl_vector *V1,
                               gsl_vector *V2,
                               gsl_vector *V3,
                               gsl_vector *survTime1,
                               gsl_vector *survTime2,
                               gsl_vector *yStar,
                               gsl_vector *survEvent1,
                               gsl_vector *survEvent2,
                               gsl_vector *case01,
                               gsl_vector *case11,
                               gsl_matrix *survCov1,
                               gsl_matrix *survCov2,
                               gsl_matrix *survCov3,
                               gsl_vector *cluster,
                               double *val)
{
    double gfunc, LP1, LP2, LP3;
    double logLH = 0;
    


    
    int jj;
    
    
    
    gsl_vector_view Xi1 = gsl_matrix_row(survCov1, i);
    gsl_blas_ddot(&Xi1.vector, beta1, &LP1);
    gsl_vector_view Xi2 = gsl_matrix_row(survCov2, i);
    gsl_blas_ddot(&Xi2.vector, beta2, &LP2);
    gsl_vector_view Xi3 = gsl_matrix_row(survCov3, i);
    gsl_blas_ddot(&Xi3.vector, beta3, &LP3);
    
    jj = (int) gsl_vector_get(cluster, i) - 1;
    
    if(gsl_vector_get(survEvent1, i) == 1)
    {
        logLH += log(alpha1) + log(kappa1) + (alpha1-1)*log(gsl_vector_get(survTime1, i));
        logLH += LP1;
        logLH += gsl_vector_get(V1, jj);
    }
    
    if(gsl_vector_get(case01, i) ==  1)
    {
        logLH += log(alpha2) + log(kappa2) + (alpha2-1)*log(gsl_vector_get(survTime2, i));
        logLH += LP2;
        logLH += gsl_vector_get(V2, jj);
    }
    
    if(gsl_vector_get(case11, i) == 1 && gsl_vector_get(yStar, i) != 0)
    {
        logLH += log(alpha3) + log(kappa3) + (alpha3-1)*log(gsl_vector_get(yStar, i));
        logLH += LP3;
        logLH += gsl_vector_get(V3, jj);
    }
    
    gfunc = BweibMvnCorScrSM_wFunc_old(i, beta1, beta2, beta3, alpha1, alpha2, alpha3, kappa1, kappa2, kappa3, V1, V2, V3, survTime1, yStar, cluster, survCov1, survCov2, survCov3);
    
    logLH += (- 1/theta - gsl_vector_get(survEvent1, i) - gsl_vector_get(survEvent2, i)) * log(1 + theta * gfunc);
    
    
    *val = logLH;
    
    return;
    
    
}



/* evaluating log-likelihood function */

/**/

void BweibMvnCorScrSM_logLH(gsl_vector *beta1,
                            gsl_vector *beta2,
                            gsl_vector *beta3,
                            double alpha1,
                            double alpha2,
                            double alpha3,
                            double kappa1,
                            double kappa2,
                            double kappa3,
                            gsl_vector *gamma,
                            gsl_vector *V1,
                            gsl_vector *V2,
                            gsl_vector *V3,
                            gsl_vector *survTime1,
                            gsl_vector *survTime2,
                            gsl_vector *yStar,
                            gsl_vector *survEvent1,
                            gsl_vector *case01,
                            gsl_vector *case11,
                            gsl_matrix *survCov1,
                            gsl_matrix *survCov2,
                            gsl_matrix *survCov3,
                            gsl_vector *cluster,
                            double *val)
{
    double gam, LP1, LP2, LP3;
    double logLH = 0;
    

    int n = survTime1 -> size;
    
    int i, jj;
    
    for(i = 0; i < n; i++)
    {
        gam = gsl_vector_get(gamma, i);
        
        gsl_vector_view Xi1 = gsl_matrix_row(survCov1, i);
        gsl_blas_ddot(&Xi1.vector, beta1, &LP1);
        gsl_vector_view Xi2 = gsl_matrix_row(survCov2, i);
        gsl_blas_ddot(&Xi2.vector, beta2, &LP2);
        gsl_vector_view Xi3 = gsl_matrix_row(survCov3, i);
        gsl_blas_ddot(&Xi3.vector, beta3, &LP3);
        
        jj = (int) gsl_vector_get(cluster, i) - 1;
        
        if(gsl_vector_get(survEvent1, i) == 1)
        {
            logLH += log(alpha1) + log(kappa1) + (alpha1-1)*log(gsl_vector_get(survTime1, i));
            logLH += log(gam);
            logLH += LP1;
            logLH += gsl_vector_get(V1, jj);
        }
        
        if(gsl_vector_get(case01, i) ==  1)
        {
            logLH += log(alpha2) + log(kappa2) + (alpha2-1)*log(gsl_vector_get(survTime2, i));
            logLH += log(gam);
            logLH += LP2;
            logLH += gsl_vector_get(V2, jj);
        }
        
        if(gsl_vector_get(case11, i) == 1 && gsl_vector_get(yStar, i) != 0)
        {
            logLH += log(alpha3) + log(kappa3) + (alpha3-1)*log(gsl_vector_get(yStar, i));
            logLH += log(gam);
            logLH += LP3;
            logLH += gsl_vector_get(V3, jj);
        }
        logLH += -gam * BweibMvnCorScrSM_wFunc_old(i, beta1, beta2, beta3, alpha1, alpha2, alpha3, kappa1, kappa2, kappa3, V1, V2, V3, survTime1, yStar, cluster, survCov1, survCov2, survCov3);
    }
    
    *val = logLH;
    
    return;
    
    
}





/* evaluating log-likelihood function for subject i */

/**/

void BweibMvnCorScrSM_logLH_i(int i,
                              gsl_vector *beta1,
                              gsl_vector *beta2,
                              gsl_vector *beta3,
                              double alpha1,
                              double alpha2,
                              double alpha3,
                              double kappa1,
                              double kappa2,
                              double kappa3,
                              gsl_vector *gamma,
                              gsl_vector *V1,
                              gsl_vector *V2,
                              gsl_vector *V3,
                              gsl_vector *survTime1,
                              gsl_vector *survTime2,
                              gsl_vector *yStar,
                              gsl_vector *survEvent1,
                              gsl_vector *case01,
                              gsl_vector *case11,
                              gsl_matrix *survCov1,
                              gsl_matrix *survCov2,
                              gsl_matrix *survCov3,
                              gsl_vector *cluster,
                              double *val)
{
    double gam, LP1, LP2, LP3;
    double logLH = 0;
    


    
    int jj;
    
    gam = gsl_vector_get(gamma, i);
    
    gsl_vector_view Xi1 = gsl_matrix_row(survCov1, i);
    gsl_blas_ddot(&Xi1.vector, beta1, &LP1);
    gsl_vector_view Xi2 = gsl_matrix_row(survCov2, i);
    gsl_blas_ddot(&Xi2.vector, beta2, &LP2);
    gsl_vector_view Xi3 = gsl_matrix_row(survCov3, i);
    gsl_blas_ddot(&Xi3.vector, beta3, &LP3);
    
    jj = (int) gsl_vector_get(cluster, i) - 1;
    
    if(gsl_vector_get(survEvent1, i) == 1)
    {
        logLH += log(alpha1) + log(kappa1) + (alpha1-1)*log(gsl_vector_get(survTime1, i));
        logLH += log(gam);
        logLH += LP1;
        logLH += gsl_vector_get(V1, jj);
    }
    
    if(gsl_vector_get(case01, i) ==  1)
    {
        logLH += log(alpha2) + log(kappa2) + (alpha2-1)*log(gsl_vector_get(survTime2, i));
        logLH += log(gam);
        logLH += LP2;
        logLH += gsl_vector_get(V2, jj);
    }
    
    if(gsl_vector_get(case11, i) == 1 && gsl_vector_get(yStar, i) != 0)
    {
        logLH += log(alpha3) + log(kappa3) + (alpha3-1)*log(gsl_vector_get(yStar, i));
        logLH += log(gam);
        logLH += LP3;
        logLH += gsl_vector_get(V3, jj);
    }
    logLH += -gam * BweibMvnCorScrSM_wFunc_old(i, beta1, beta2, beta3, alpha1, alpha2, alpha3, kappa1, kappa2, kappa3, V1, V2, V3, survTime1, yStar, cluster, survCov1, survCov2, survCov3);
    
    
    *val = logLH;
    
    return;
    
    
}





/*
 Evaluate w(y1, y2) function with nu2 and nu3
 */
double BweibMvnCorScrSM_wFunc(int subjInx,
                              gsl_vector *beta1,
                              gsl_vector *beta2,
                              gsl_vector *beta3,
                              double alpha1,
                              double alpha2,
                              double alpha3,
                              double kappa1,
                              double kappa2,
                              double kappa3,
                              double nu2,
                              double nu3,
                              double gam,
                              gsl_vector *V1,
                              gsl_vector *V2,
                              gsl_vector *V3,
                              gsl_vector *survTime1,
                              gsl_vector *survTime2,
                              gsl_vector *cluster,
                              gsl_matrix *survCov1,
                              gsl_matrix *survCov2,
                              gsl_matrix *survCov3)
{
    int i = subjInx;
    double LP1, LP2, LP3, cumHaz1, cumHaz2, cumHaz3diff;
    double wVal;
    int jj = (int) gsl_vector_get(cluster, i) - 1;
    
    gsl_vector_view Xi1 = gsl_matrix_row(survCov1, i);
    gsl_blas_ddot(&Xi1.vector, beta1 ,&LP1);
    
    gsl_vector_view Xi2 = gsl_matrix_row(survCov2, i);
    gsl_blas_ddot(&Xi2.vector, beta2 ,&LP2);
    
    gsl_vector_view Xi3 = gsl_matrix_row(survCov3, i);
    gsl_blas_ddot(&Xi3.vector, beta3 ,&LP3);
    
    cumHaz1     = kappa1 * pow(gsl_vector_get(survTime1, i), alpha1) * exp(LP1 + gsl_vector_get(V1, jj));
    cumHaz2     = kappa2 * pow(gsl_vector_get(survTime1, i), alpha2) * exp(LP2 + gsl_vector_get(V2, jj));
    cumHaz3diff = kappa3 * (pow(gsl_vector_get(survTime2, i), alpha3)-pow(gsl_vector_get(survTime1, i), alpha3)) * exp(LP3 + gsl_vector_get(V3, jj));
    
    wVal = gam * cumHaz1 + pow(gam, nu2) * cumHaz2 + pow(gam, nu3) * cumHaz3diff;
    
    return wVal;
}



/*
 Evaluate w(y1, y2) function
 */
double BweibMvnCorScrSM_wFunc_old(int subjInx,
                                  gsl_vector *beta1,
                                  gsl_vector *beta2,
                                  gsl_vector *beta3,
                                  double alpha1,
                                  double alpha2,
                                  double alpha3,
                                  double kappa1,
                                  double kappa2,
                                  double kappa3,
                                  gsl_vector *V1,
                                  gsl_vector *V2,
                                  gsl_vector *V3,
                                  gsl_vector *survTime1,
                                  gsl_vector *yStar,
                                  gsl_vector *cluster,
                                  gsl_matrix *survCov1,
                                  gsl_matrix *survCov2,
                                  gsl_matrix *survCov3)
{
    int i = subjInx;
    int j = (int) gsl_vector_get(cluster, i) - 1;
    double LP1, LP2, LP3, cumHaz1, cumHaz2, cumHaz3diff;
    double wVal;
    
    gsl_vector_view Xi1 = gsl_matrix_row(survCov1, i);
    gsl_blas_ddot(&Xi1.vector, beta1 ,&LP1);
    
    gsl_vector_view Xi2 = gsl_matrix_row(survCov2, i);
    gsl_blas_ddot(&Xi2.vector, beta2 ,&LP2);
    
    gsl_vector_view Xi3 = gsl_matrix_row(survCov3, i);
    gsl_blas_ddot(&Xi3.vector, beta3 ,&LP3);
    
    cumHaz1     = kappa1 * pow(gsl_vector_get(survTime1, i), alpha1) * exp(LP1 + gsl_vector_get(V1, j));
    cumHaz2     = kappa2 * pow(gsl_vector_get(survTime1, i), alpha2) * exp(LP2 + gsl_vector_get(V2, j));
    cumHaz3diff = kappa3 * pow(gsl_vector_get(yStar, i), alpha3) * exp(LP3 + gsl_vector_get(V3, j));
    
    wVal = cumHaz1 + cumHaz2 + cumHaz3diff;
    
    return wVal;
}







/********* For Weibull-MVN-M model *************/

/* evaluating log-likelihood function */

/**/

void BweibMvnCorScr_logLH(gsl_vector *beta1,
                          gsl_vector *beta2,
                          gsl_vector *beta3,
                          double alpha1,
                          double alpha2,
                          double alpha3,
                          double kappa1,
                          double kappa2,
                          double kappa3,
                          gsl_vector *gamma,
                          gsl_vector *V1,
                          gsl_vector *V2,
                          gsl_vector *V3,
                          gsl_vector *survTime1,
                          gsl_vector *survTime2,
                          gsl_vector *survEvent1,
                          gsl_vector *case01,
                          gsl_vector *case11,
                          gsl_matrix *survCov1,
                          gsl_matrix *survCov2,
                          gsl_matrix *survCov3,
                          gsl_vector *cluster,
                          double *val)
{
    double gam, LP1, LP2, LP3;
    double logLH = 0;
    

    int n = survTime1 -> size;
    
    int i, jj;
    
    for(i = 0; i < n; i++)
    {
        gam = gsl_vector_get(gamma, i);
        
        gsl_vector_view Xi1 = gsl_matrix_row(survCov1, i);
        gsl_blas_ddot(&Xi1.vector, beta1, &LP1);
        gsl_vector_view Xi2 = gsl_matrix_row(survCov2, i);
        gsl_blas_ddot(&Xi2.vector, beta2, &LP2);
        gsl_vector_view Xi3 = gsl_matrix_row(survCov3, i);
        gsl_blas_ddot(&Xi3.vector, beta3, &LP3);
        
        jj = (int) gsl_vector_get(cluster, i) - 1;
        
        if(gsl_vector_get(survEvent1, i) == 1)
        {
            logLH += log(alpha1) + log(kappa1) + (alpha1-1)*log(gsl_vector_get(survTime1, i));
            logLH += log(gam);
            logLH += LP1;
            logLH += gsl_vector_get(V1, jj);
        }
        
        if(gsl_vector_get(case01, i) ==  1)
        {
            logLH += log(alpha2) + log(kappa2) + (alpha2-1)*log(gsl_vector_get(survTime2, i));
            logLH += log(gam);
            logLH += LP2;
            logLH += gsl_vector_get(V2, jj);
        }
        
        if(gsl_vector_get(case11, i) == 1)
        {
            logLH += log(alpha3) + log(kappa3) + (alpha3-1)*log(gsl_vector_get(survTime2, i));
            logLH += log(gam);
            logLH += LP3;
            logLH += gsl_vector_get(V3, jj);
        }
        logLH += -gam * BweibMvnCorScr_wFunc_old(i, beta1, beta2, beta3, alpha1, alpha2, alpha3, kappa1, kappa2, kappa3, V1, V2, V3, survTime1, survTime2, cluster, survCov1, survCov2, survCov3);
    }
    
    *val = logLH;
    
    return;
    
    
}





/* evaluating log-likelihood function for subject i */

/**/

void BweibMvnCorScr_logLH_i(int i,
                            gsl_vector *beta1,
                            gsl_vector *beta2,
                            gsl_vector *beta3,
                            double alpha1,
                            double alpha2,
                            double alpha3,
                            double kappa1,
                            double kappa2,
                            double kappa3,
                            gsl_vector *gamma,
                            gsl_vector *V1,
                            gsl_vector *V2,
                            gsl_vector *V3,
                            gsl_vector *survTime1,
                            gsl_vector *survTime2,
                            gsl_vector *survEvent1,
                            gsl_vector *case01,
                            gsl_vector *case11,
                            gsl_matrix *survCov1,
                            gsl_matrix *survCov2,
                            gsl_matrix *survCov3,
                            gsl_vector *cluster,
                            double *val)
{
    double gam, LP1, LP2, LP3;
    double logLH = 0;
    


    
    int jj;
    
    gam = gsl_vector_get(gamma, i);
    
    gsl_vector_view Xi1 = gsl_matrix_row(survCov1, i);
    gsl_blas_ddot(&Xi1.vector, beta1, &LP1);
    gsl_vector_view Xi2 = gsl_matrix_row(survCov2, i);
    gsl_blas_ddot(&Xi2.vector, beta2, &LP2);
    gsl_vector_view Xi3 = gsl_matrix_row(survCov3, i);
    gsl_blas_ddot(&Xi3.vector, beta3, &LP3);
    
    jj = (int) gsl_vector_get(cluster, i) - 1;
    
    if(gsl_vector_get(survEvent1, i) == 1)
    {
        logLH += log(alpha1) + log(kappa1) + (alpha1-1)*log(gsl_vector_get(survTime1, i));
        logLH += log(gam);
        logLH += LP1;
        logLH += gsl_vector_get(V1, jj);
    }
    
    if(gsl_vector_get(case01, i) ==  1)
    {
        logLH += log(alpha2) + log(kappa2) + (alpha2-1)*log(gsl_vector_get(survTime2, i));
        logLH += log(gam);
        logLH += LP2;
        logLH += gsl_vector_get(V2, jj);
    }
    
    if(gsl_vector_get(case11, i) == 1)
    {
        logLH += log(alpha3) + log(kappa3) + (alpha3-1)*log(gsl_vector_get(survTime2, i));
        logLH += log(gam);
        logLH += LP3;
        logLH += gsl_vector_get(V3, jj);
    }
    logLH += -gam * BweibMvnCorScr_wFunc_old(i, beta1, beta2, beta3, alpha1, alpha2, alpha3, kappa1, kappa2, kappa3, V1, V2, V3, survTime1, survTime2, cluster, survCov1, survCov2, survCov3);
    
    
    *val = logLH;
    
    return;
    
    
}








/* evaluating log-likelihood function */

/**/

void BweibMvnCorScr_logMLH(gsl_vector *beta1,
                           gsl_vector *beta2,
                           gsl_vector *beta3,
                           double alpha1,
                           double alpha2,
                           double alpha3,
                           double kappa1,
                           double kappa2,
                           double kappa3,
                           double theta,
                           gsl_vector *V1,
                           gsl_vector *V2,
                           gsl_vector *V3,
                           gsl_vector *survTime1,
                           gsl_vector *survTime2,
                           gsl_vector *survEvent1,
                           gsl_vector *survEvent2,
                           gsl_vector *case01,
                           gsl_vector *case11,
                           gsl_matrix *survCov1,
                           gsl_matrix *survCov2,
                           gsl_matrix *survCov3,
                           gsl_vector *cluster,
                           double *val)
{
    double gfunc, LP1, LP2, LP3;
    double logLH = 0;
    

    int n = survTime1 -> size;
    
    int i, jj;
    
    for(i = 0; i < n; i++)
    {
        gsl_vector_view Xi1 = gsl_matrix_row(survCov1, i);
        gsl_blas_ddot(&Xi1.vector, beta1, &LP1);
        gsl_vector_view Xi2 = gsl_matrix_row(survCov2, i);
        gsl_blas_ddot(&Xi2.vector, beta2, &LP2);
        gsl_vector_view Xi3 = gsl_matrix_row(survCov3, i);
        gsl_blas_ddot(&Xi3.vector, beta3, &LP3);
        
        jj = (int) gsl_vector_get(cluster, i) - 1;
        
        if(gsl_vector_get(survEvent1, i) == 1)
        {
            logLH += log(alpha1) + log(kappa1) + (alpha1-1)*log(gsl_vector_get(survTime1, i));
            logLH += LP1;
            logLH += gsl_vector_get(V1, jj);
        }
        
        if(gsl_vector_get(case01, i) ==  1)
        {
            logLH += log(alpha2) + log(kappa2) + (alpha2-1)*log(gsl_vector_get(survTime2, i));
            logLH += LP2;
            logLH += gsl_vector_get(V2, jj);
        }
        
        if(gsl_vector_get(case11, i) == 1)
        {
            logLH += log(alpha3) + log(kappa3) + (alpha3-1)*log(gsl_vector_get(survTime2, i));
            logLH += LP3;
            logLH += gsl_vector_get(V3, jj);
        }
        gfunc = BweibMvnCorScr_wFunc_old(i, beta1, beta2, beta3, alpha1, alpha2, alpha3, kappa1, kappa2, kappa3, V1, V2, V3, survTime1, survTime2, cluster, survCov1, survCov2, survCov3);
        
        logLH += (- 1/theta - gsl_vector_get(survEvent1, i) - gsl_vector_get(survEvent2, i)) * log(1 + theta * gfunc);
    }
    
    *val = logLH;
    
    return;
    
    
}





/* evaluating log-likelihood function for subject i */

/**/

void BweibMvnCorScr_logMLH_i(int i,
                             gsl_vector *beta1,
                             gsl_vector *beta2,
                             gsl_vector *beta3,
                             double alpha1,
                             double alpha2,
                             double alpha3,
                             double kappa1,
                             double kappa2,
                             double kappa3,
                             double theta,
                             gsl_vector *V1,
                             gsl_vector *V2,
                             gsl_vector *V3,
                             gsl_vector *survTime1,
                             gsl_vector *survTime2,
                             gsl_vector *survEvent1,
                             gsl_vector *survEvent2,
                             gsl_vector *case01,
                             gsl_vector *case11,
                             gsl_matrix *survCov1,
                             gsl_matrix *survCov2,
                             gsl_matrix *survCov3,
                             gsl_vector *cluster,
                             double *val)
{
    double gfunc, LP1, LP2, LP3;
    double logLH = 0;
    


    
    int jj;
    
    gsl_vector_view Xi1 = gsl_matrix_row(survCov1, i);
    gsl_blas_ddot(&Xi1.vector, beta1, &LP1);
    gsl_vector_view Xi2 = gsl_matrix_row(survCov2, i);
    gsl_blas_ddot(&Xi2.vector, beta2, &LP2);
    gsl_vector_view Xi3 = gsl_matrix_row(survCov3, i);
    gsl_blas_ddot(&Xi3.vector, beta3, &LP3);
    
    jj = (int) gsl_vector_get(cluster, i) - 1;
    
    if(gsl_vector_get(survEvent1, i) == 1)
    {
        logLH += log(alpha1) + log(kappa1) + (alpha1-1)*log(gsl_vector_get(survTime1, i));
        logLH += LP1;
        logLH += gsl_vector_get(V1, jj);
    }
    
    if(gsl_vector_get(case01, i) ==  1)
    {
        logLH += log(alpha2) + log(kappa2) + (alpha2-1)*log(gsl_vector_get(survTime2, i));
        logLH += LP2;
        logLH += gsl_vector_get(V2, jj);
    }
    
    if(gsl_vector_get(case11, i) == 1)
    {
        logLH += log(alpha3) + log(kappa3) + (alpha3-1)*log(gsl_vector_get(survTime2, i));
        logLH += LP3;
        logLH += gsl_vector_get(V3, jj);
    }
    
    gfunc = BweibMvnCorScr_wFunc_old(i, beta1, beta2, beta3, alpha1, alpha2, alpha3, kappa1, kappa2, kappa3, V1, V2, V3, survTime1, survTime2, cluster, survCov1, survCov2, survCov3);
    
    logLH += (- 1/theta - gsl_vector_get(survEvent1, i) - gsl_vector_get(survEvent2, i)) * log(1 + theta * gfunc);
    
    *val = logLH;
    
    return;
    
    
}


/*
 Random generation from the Inverse Wishart distribution
 */

void c_riwishart(int v,
                 gsl_matrix *X_ori,
                 gsl_matrix *sample)
{
    int i, df;
    double normVal;
    
    gsl_matrix *X = gsl_matrix_calloc(3, 3);
    matrixInv(X_ori, X);
    
    gsl_matrix *cholX = gsl_matrix_calloc(3, 3);
    gsl_matrix *ZZ = gsl_matrix_calloc(3, 3);
    gsl_matrix *XX = gsl_matrix_calloc(3, 3);
    gsl_matrix *KK = gsl_matrix_calloc(3, 3);
    
    gsl_matrix_memcpy(cholX, X);
    gsl_linalg_cholesky_decomp(cholX);
    
    gsl_matrix_set(cholX, 1, 0, 0);
    gsl_matrix_set(cholX, 2, 0, 0);
    gsl_matrix_set(cholX, 2, 1, 0);
    
    for(i = 0; i < 3; i++)
    {
        df = v - i;
        gsl_matrix_set(ZZ, i, i, sqrt(rchisq(df)));
    }
    
    normVal = rnorm(0, 1);
    gsl_matrix_set(ZZ, 0, 1, normVal);
    
    normVal = rnorm(0, 1);
    gsl_matrix_set(ZZ, 0, 2, normVal);
    
    normVal = rnorm(0, 1);
    gsl_matrix_set(ZZ, 1, 2, normVal);
    
    
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, ZZ, cholX, 0, XX);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, XX, XX, 0, KK);
    matrixInv(KK, sample);
    
    gsl_matrix_free(X);
    gsl_matrix_free(cholX);
    gsl_matrix_free(XX);
    gsl_matrix_free(ZZ);
    gsl_matrix_free(KK);
    
}




/*
 Evaluate w(y1, y2) function with nu2 and nu3
 */
double BweibMvnCorScr_wFunc(int subjInx,
                            gsl_vector *beta1,
                            gsl_vector *beta2,
                            gsl_vector *beta3,
                            double alpha1,
                            double alpha2,
                            double alpha3,
                            double kappa1,
                            double kappa2,
                            double kappa3,
                            double nu2,
                            double nu3,
                            double gam,
                            gsl_vector *V1,
                            gsl_vector *V2,
                            gsl_vector *V3,
                            gsl_vector *survTime1,
                            gsl_vector *survTime2,
                            gsl_vector *cluster,
                            gsl_matrix *survCov1,
                            gsl_matrix *survCov2,
                            gsl_matrix *survCov3)
{
    int i = subjInx;
    double LP1, LP2, LP3, cumHaz1, cumHaz2, cumHaz3diff;
    double wVal;
    int jj = (int) gsl_vector_get(cluster, i) - 1;
    
    gsl_vector_view Xi1 = gsl_matrix_row(survCov1, i);
    gsl_blas_ddot(&Xi1.vector, beta1 ,&LP1);
    
    gsl_vector_view Xi2 = gsl_matrix_row(survCov2, i);
    gsl_blas_ddot(&Xi2.vector, beta2 ,&LP2);
    
    gsl_vector_view Xi3 = gsl_matrix_row(survCov3, i);
    gsl_blas_ddot(&Xi3.vector, beta3 ,&LP3);
    
    cumHaz1     = kappa1 * pow(gsl_vector_get(survTime1, i), alpha1) * exp(LP1 + gsl_vector_get(V1, jj));
    cumHaz2     = kappa2 * pow(gsl_vector_get(survTime1, i), alpha2) * exp(LP2 + gsl_vector_get(V2, jj));
    cumHaz3diff = kappa3 * (pow(gsl_vector_get(survTime2, i), alpha3)-pow(gsl_vector_get(survTime1, i), alpha3)) * exp(LP3 + gsl_vector_get(V3, jj));
    
    wVal = gam * cumHaz1 + pow(gam, nu2) * cumHaz2 + pow(gam, nu3) * cumHaz3diff;
    
    return wVal;
}



/*
 Evaluate w(y1, y2) function
 */
double BweibMvnCorScr_wFunc_old(int subjInx,
                                gsl_vector *beta1,
                                gsl_vector *beta2,
                                gsl_vector *beta3,
                                double alpha1,
                                double alpha2,
                                double alpha3,
                                double kappa1,
                                double kappa2,
                                double kappa3,
                                gsl_vector *V1,
                                gsl_vector *V2,
                                gsl_vector *V3,
                                gsl_vector *survTime1,
                                gsl_vector *survTime2,
                                gsl_vector *cluster,
                                gsl_matrix *survCov1,
                                gsl_matrix *survCov2,
                                gsl_matrix *survCov3)
{
    int i = subjInx;
    int j = (int) gsl_vector_get(cluster, i) - 1;
    double LP1, LP2, LP3, cumHaz1, cumHaz2, cumHaz3diff;
    double wVal;
    
    gsl_vector_view Xi1 = gsl_matrix_row(survCov1, i);
    gsl_blas_ddot(&Xi1.vector, beta1 ,&LP1);
    
    gsl_vector_view Xi2 = gsl_matrix_row(survCov2, i);
    gsl_blas_ddot(&Xi2.vector, beta2 ,&LP2);
    
    gsl_vector_view Xi3 = gsl_matrix_row(survCov3, i);
    gsl_blas_ddot(&Xi3.vector, beta3 ,&LP3);
    
    cumHaz1     = kappa1 * pow(gsl_vector_get(survTime1, i), alpha1) * exp(LP1 + gsl_vector_get(V1, j));
    cumHaz2     = kappa2 * pow(gsl_vector_get(survTime1, i), alpha2) * exp(LP2 + gsl_vector_get(V2, j));
    cumHaz3diff = kappa3 * (pow(gsl_vector_get(survTime2, i), alpha3)-pow(gsl_vector_get(survTime1, i), alpha3)) * exp(LP3 + gsl_vector_get(V3, j));
    
    wVal = cumHaz1 + cumHaz2 + cumHaz3diff;
    
    return wVal;
}







/*
 Density calculation for a multivariate normal distribution
 */
void c_dmvnorm2(gsl_vector *x,
               gsl_vector *mu,
               double     sigma,
               gsl_matrix *AInv,
               double     *value)
{
    int signum, K = x->size;
	double sigmaSqInv = pow(sigma, -2);
    
    gsl_vector *diff       = gsl_vector_alloc(K);
	gsl_matrix *SigmaInv   = gsl_matrix_alloc(K, K);
    gsl_matrix *SigmaInvLU = gsl_matrix_alloc(K, K);
    gsl_permutation *p     = gsl_permutation_alloc(K);
    
	gsl_vector_memcpy(diff, x);
	gsl_vector_sub(diff, mu);
    
	gsl_matrix_memcpy(SigmaInv, AInv);
	gsl_matrix_scale(SigmaInv, sigmaSqInv);
    gsl_matrix_memcpy(SigmaInvLU, SigmaInv);
    gsl_linalg_LU_decomp(SigmaInvLU, p, &signum);
    
	c_quadform_vMv(diff, SigmaInv, value);
    
    *value = (log(gsl_linalg_LU_det(SigmaInvLU, signum)) - log(pow(2*Pi, K)) - *value) / 2;
    
	gsl_vector_free(diff);
	gsl_matrix_free(SigmaInv);
	gsl_matrix_free(SigmaInvLU);
	gsl_permutation_free(p);
    return;
}









/*********************************************************************************
 *********************************************************************************
 *********************************************************************************
 *********************************************************************************
 *********************************************************************************
 *********************************************************************************
 *********************************************************************************
 *********************************************************************************
 *********************************************************************************
 *********************************************************************************
 *********************************************************************************
 Functions for models in BayesID  */


/*
 Evaluate w(y1, y2) function
 */
double Bscr_wFunc(int subjInx,
             gsl_vector *xbeta1,
             gsl_vector *xbeta2,
             gsl_vector *xbeta3,
             gsl_vector *lambda1,
             gsl_vector *lambda2,
             gsl_vector *lambda3,
             gsl_vector *s1,
             gsl_vector *s2,
             gsl_vector *s3,
             int J1,
             int J2,
             int J3,
             gsl_vector *survTime1,
             gsl_vector *survTime2)
{
    int i = subjInx;
    double cumHaz1, cumHaz2, cumHaz3diff;
    double Del, wVal;
    int j;
    
    cumHaz1 = 0; cumHaz2 = 0; cumHaz3diff = 0;
    
    for(j = 0; j < J1+1; j++)
    {
        if(j > 0)
        {
            Del = c_max(0, (c_min(gsl_vector_get(s1, j), gsl_vector_get(survTime1, i)) - gsl_vector_get(s1, j-1)));
        }
        if(j == 0)
        {
            Del = c_max(0, c_min(gsl_vector_get(s1, j), gsl_vector_get(survTime1, i)) - 0);
        }
        cumHaz1 += Del* exp(gsl_vector_get(lambda1, j)) * exp(gsl_vector_get(xbeta1, i));
    }
    
    for(j = 0; j < J2+1; j++)
    {
        if(j > 0)
        {
            Del = c_max(0, (c_min(gsl_vector_get(s2, j), gsl_vector_get(survTime1, i)) - gsl_vector_get(s2, j-1)));
        }
        if(j == 0)
        {
            Del = c_max(0, c_min(gsl_vector_get(s2, j), gsl_vector_get(survTime1, i)) - 0);
        }
        cumHaz2 += Del* exp(gsl_vector_get(lambda2, j)) * exp(gsl_vector_get(xbeta2, i));
    }
    
    for(j = 0; j < J3+1; j++)
    {
        if(j > 0)
        {
            Del = c_max(0, (c_min(gsl_vector_get(s3, j), gsl_vector_get(survTime2, i)) - c_max(gsl_vector_get(s3, j-1), gsl_vector_get(survTime1, i))));
        }
        if(j == 0)
        {
            Del = c_max(0, (c_min(gsl_vector_get(s3, j), gsl_vector_get(survTime2, i)) - c_max(0, gsl_vector_get(survTime1, i))));
        }
        cumHaz3diff += Del* exp(gsl_vector_get(lambda3, j)) * exp(gsl_vector_get(xbeta3, i));
    }
        
    wVal = cumHaz1 + cumHaz2 + cumHaz3diff;
    
    return wVal;
}






/*
 Evaluate w(y1, y2) function
 */
double BscrSM_wFunc(int subjInx,
                    gsl_vector *xbeta1,
                    gsl_vector *xbeta2,
                    gsl_vector *xbeta3,
                    gsl_vector *lambda1,
                    gsl_vector *lambda2,
                    gsl_vector *lambda3,
                    gsl_vector *s1,
                    gsl_vector *s2,
                    gsl_vector *s3,
                    int J1,
                    int J2,
                    int J3,
                    gsl_vector *survTime1,
                    gsl_vector *yStar)
{
    int i = subjInx;
    double cumHaz1, cumHaz2, cumHaz3diff;
    double Del, wVal;
    int j;
    
    cumHaz1 = 0; cumHaz2 = 0; cumHaz3diff = 0;
    
    for(j = 0; j < J1+1; j++)
    {
        if(j > 0)
        {
            Del = c_max(0, (c_min(gsl_vector_get(s1, j), gsl_vector_get(survTime1, i)) - gsl_vector_get(s1, j-1)));
        }
        if(j == 0)
        {
            Del = c_max(0, c_min(gsl_vector_get(s1, j), gsl_vector_get(survTime1, i)) - 0);
        }
        cumHaz1 += Del* exp(gsl_vector_get(lambda1, j)) * exp(gsl_vector_get(xbeta1, i));
    }
    
    for(j = 0; j < J2+1; j++)
    {
        if(j > 0)
        {
            Del = c_max(0, (c_min(gsl_vector_get(s2, j), gsl_vector_get(survTime1, i)) - gsl_vector_get(s2, j-1)));
        }
        if(j == 0)
        {
            Del = c_max(0, c_min(gsl_vector_get(s2, j), gsl_vector_get(survTime1, i)) - 0);
        }
        cumHaz2 += Del* exp(gsl_vector_get(lambda2, j)) * exp(gsl_vector_get(xbeta2, i));
    }
    
    for(j = 0; j < J3+1; j++)
    {
        if(j > 0)
        {
            Del = c_max(0, (c_min(gsl_vector_get(s3, j), gsl_vector_get(yStar, i)) - gsl_vector_get(s3, j-1)));
        }
        if(j == 0)
        {
            Del = c_max(0, c_min(gsl_vector_get(s3, j), gsl_vector_get(yStar, i)) - 0);
        }
        cumHaz3diff += Del* exp(gsl_vector_get(lambda3, j)) * exp(gsl_vector_get(xbeta3, i));
    }
        
    wVal = cumHaz1 + cumHaz2 + cumHaz3diff;
    
    return wVal;
}







/*
 Evaluate w(y1, y2) function
 */
double BweibScr_wFunc(int subjInx,
                      gsl_vector *beta1,
                      gsl_vector *beta2,
                      gsl_vector *beta3,
                      double alpha1,
                      double alpha2,
                      double alpha3,
                      double kappa1,
                      double kappa2,
                      double kappa3,
                      gsl_vector *survTime1,
                      gsl_vector *survTime2,
                      gsl_matrix *survCov1,
                      gsl_matrix *survCov2,
                      gsl_matrix *survCov3)
{
    int i = subjInx;
    double LP1, LP2, LP3, cumHaz1, cumHaz2, cumHaz3diff;
    double wVal;
    
    gsl_vector_view Xi1 = gsl_matrix_row(survCov1, i);
    gsl_blas_ddot(&Xi1.vector, beta1 ,&LP1);
    
    gsl_vector_view Xi2 = gsl_matrix_row(survCov2, i);
    gsl_blas_ddot(&Xi2.vector, beta2 ,&LP2);
    
    gsl_vector_view Xi3 = gsl_matrix_row(survCov3, i);
    gsl_blas_ddot(&Xi3.vector, beta3 ,&LP3);
    
    cumHaz1     = kappa1 * pow(gsl_vector_get(survTime1, i), alpha1) * exp(LP1);
    cumHaz2     = kappa2 * pow(gsl_vector_get(survTime1, i), alpha2) * exp(LP2);
    cumHaz3diff = kappa3 * (pow(gsl_vector_get(survTime2, i), alpha3)-pow(gsl_vector_get(survTime1, i), alpha3)) * exp(LP3);
    
    wVal = cumHaz1 + cumHaz2 + cumHaz3diff;
    
    return wVal;
}


/*
 Evaluate w(y1, y2) function
 */
double BweibScrSM_wFunc(int subjInx,
                      gsl_vector *beta1,
                      gsl_vector *beta2,
                      gsl_vector *beta3,
                      double alpha1,
                      double alpha2,
                      double alpha3,
                      double kappa1,
                      double kappa2,
                      double kappa3,
                      gsl_vector *survTime1,
                      gsl_vector *yStar,
                      gsl_matrix *survCov1,
                      gsl_matrix *survCov2,
                      gsl_matrix *survCov3)
{
    int i = subjInx;
    double LP1, LP2, LP3, cumHaz1, cumHaz2, cumHaz3diff;
    double wVal;
    
    gsl_vector_view Xi1 = gsl_matrix_row(survCov1, i);
    gsl_blas_ddot(&Xi1.vector, beta1 ,&LP1);
    
    gsl_vector_view Xi2 = gsl_matrix_row(survCov2, i);
    gsl_blas_ddot(&Xi2.vector, beta2 ,&LP2);
    
    gsl_vector_view Xi3 = gsl_matrix_row(survCov3, i);
    gsl_blas_ddot(&Xi3.vector, beta3 ,&LP3);
    
    cumHaz1     = kappa1 * pow(gsl_vector_get(survTime1, i), alpha1) * exp(LP1);
    cumHaz2     = kappa2 * pow(gsl_vector_get(survTime1, i), alpha2) * exp(LP2);
    cumHaz3diff = kappa3 * pow(gsl_vector_get(yStar, i), alpha3) * exp(LP3);
    
    wVal = cumHaz1 + cumHaz2 + cumHaz3diff;
    
    return wVal;
}







void cal_Sigma(gsl_matrix *Sigma_lam,
               gsl_matrix *invSigma_lam,
               gsl_matrix *W,
               gsl_matrix *Q,
               gsl_vector *s,
               double c_lam,
               int J)
{
    int j;
    gsl_matrix_view W_sub               = gsl_matrix_submatrix(W, 0, 0, J+1, J+1);
    gsl_matrix_view Q_sub               = gsl_matrix_submatrix(Q, 0, 0, J+1, J+1);
    gsl_matrix_view Sigma_lam_sub       = gsl_matrix_submatrix(Sigma_lam, 0, 0, J+1, J+1);
    gsl_matrix_view invSigma_lam_sub    = gsl_matrix_submatrix(invSigma_lam, 0, 0, J+1, J+1);
    gsl_matrix *ImW                     = gsl_matrix_calloc(J+1, J+1);
    gsl_matrix *invImW                  = gsl_matrix_calloc(J+1, J+1);
    
    gsl_vector *diff_s = gsl_vector_calloc(J+1);
    
    if(J+1 == 1)
    {
        gsl_vector_set(diff_s, 0, gsl_vector_get(s, 0));
        
        gsl_matrix_set(&W_sub.matrix, 0, 0, 0);
        
        gsl_matrix_set(&Q_sub.matrix, 0, 0, 2/(2*gsl_vector_get(diff_s, 0)));
        
        gsl_matrix_set(&Sigma_lam_sub.matrix, 0, 0, gsl_matrix_get(&Q_sub.matrix, 0, 0));
        gsl_matrix_set(&invSigma_lam_sub.matrix, 0, 0, 1/gsl_matrix_get(&Sigma_lam_sub.matrix, 0, 0));

    }
    
    if(J+1 ==2)
    {
        gsl_vector_set(diff_s, 0, gsl_vector_get(s, 0));
        gsl_vector_set(diff_s, 1, gsl_vector_get(s, 1)-gsl_vector_get(s, 0));
        
        gsl_matrix_set(&W_sub.matrix, 0, 1, c_lam * (gsl_vector_get(diff_s, 0) + gsl_vector_get(diff_s, 1)) / (2*gsl_vector_get(diff_s, 0) + gsl_vector_get(diff_s, 1)));
        
        gsl_matrix_set(&W_sub.matrix, J, J-1, c_lam * (gsl_vector_get(diff_s, J-1) + gsl_vector_get(diff_s, J)) / (gsl_vector_get(diff_s, J-1) + 2*gsl_vector_get(diff_s, J)));
        
        gsl_matrix_set(&Q_sub.matrix, 0, 0, 2/(2*gsl_vector_get(diff_s, 0) + gsl_vector_get(diff_s, 1)));
        
        gsl_matrix_set(&Q_sub.matrix, J, J, 2/(gsl_vector_get(diff_s, J-1) + 2*gsl_vector_get(diff_s, J)));
        
        gsl_matrix_set_identity(ImW);
        gsl_matrix_sub(ImW, &W_sub.matrix);
        
        c_solve(ImW, invImW);
        
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, invImW, &Q_sub.matrix, 0, &Sigma_lam_sub.matrix);
        
        c_solve(&Sigma_lam_sub.matrix, &invSigma_lam_sub.matrix);
        
        gsl_matrix_free(ImW);
        gsl_matrix_free(invImW);
        
    }
    
    if(J+1 >= 3){
        gsl_vector_set(diff_s, 0, gsl_vector_get(s, 0));
        for(j = 1; j < J+1; j++)
        {
            gsl_vector_set(diff_s, j, gsl_vector_get(s, j)-gsl_vector_get(s, j-1));
        }
        
        
        gsl_matrix_set(&W_sub.matrix, 0, 1, c_lam * (gsl_vector_get(diff_s, 0) + gsl_vector_get(diff_s, 1)) / (2*gsl_vector_get(diff_s, 0) + gsl_vector_get(diff_s, 1)));
        
        gsl_matrix_set(&W_sub.matrix, J, J-1, c_lam * (gsl_vector_get(diff_s, J-1) + gsl_vector_get(diff_s, J)) / (gsl_vector_get(diff_s, J-1) + 2*gsl_vector_get(diff_s, J)));
        
        gsl_matrix_set(&Q_sub.matrix, 0, 0, 2/(2*gsl_vector_get(diff_s, 0) + gsl_vector_get(diff_s, 1)));
        
        gsl_matrix_set(&Q_sub.matrix, J, J, 2/(gsl_vector_get(diff_s, J-1) + 2*gsl_vector_get(diff_s, J)));
        
        for(j = 1; j < J; j++)
        {
            gsl_matrix_set(&Q_sub.matrix, j, j, 2/(gsl_vector_get(diff_s, j-1) + 2*gsl_vector_get(diff_s, j) + gsl_vector_get(diff_s, j+1)));
            gsl_matrix_set(&W_sub.matrix, j, j-1, c_lam * (gsl_vector_get(diff_s, j-1) + gsl_vector_get(diff_s, j)) / (gsl_vector_get(diff_s, j-1) + 2*gsl_vector_get(diff_s, j) + gsl_vector_get(diff_s, j+1)));
            gsl_matrix_set(&W_sub.matrix, j, j+1, c_lam * (gsl_vector_get(diff_s, j) + gsl_vector_get(diff_s, j+1)) / (gsl_vector_get(diff_s, j-1) + 2*gsl_vector_get(diff_s, j) + gsl_vector_get(diff_s, j+1)));
        }
        
        gsl_matrix_set_identity(ImW);
        gsl_matrix_sub(ImW, &W_sub.matrix);
        
        c_solve(ImW, invImW);
        
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, invImW, &Q_sub.matrix, 0, &Sigma_lam_sub.matrix);
        
        c_solve(&Sigma_lam_sub.matrix, &invSigma_lam_sub.matrix);
        
        gsl_matrix_free(ImW);
        gsl_matrix_free(invImW);
    }
    
    return;
}





void set_Ind(gsl_matrix *ind_d,
             gsl_matrix *ind_r,
             gsl_vector *nEvent,
             gsl_vector *s,
             gsl_vector *survTime,
             gsl_vector *survEvent,
             gsl_vector *case0yleq,
             gsl_vector *case0ygeq,
             gsl_vector *case1yleq,
             gsl_vector *case1ygeq,
             double s_max,
             int J)
{
    int i, j;
    int n = survTime -> size;
    
    for(i = 0; i < n; i++)
    {
        if(gsl_vector_get(survEvent, i)==0 & gsl_vector_get(survTime, i)<= s_max) gsl_vector_set(case0yleq, i, 1);
        if(gsl_vector_get(survEvent, i)==0 & gsl_vector_get(survTime, i)> s_max) gsl_vector_set(case0ygeq, i, 1);
        if(gsl_vector_get(survEvent, i)==1 & gsl_vector_get(survTime, i)<= s_max) gsl_vector_set(case1yleq, i, 1);
        if(gsl_vector_get(survEvent, i)==1 & gsl_vector_get(survTime, i)> s_max) gsl_vector_set(case1ygeq, i, 1);
    }
    
    
    for(i = 0; i < n; i++)
    {
        if(gsl_vector_get(case1yleq, i) == 1)
        {
            for(j = 0; j < J; j++)
            {
                if(gsl_vector_get(survTime, i) > gsl_vector_get(s, j) & gsl_vector_get(survTime, i) <= gsl_vector_get(s, j+1)) gsl_matrix_set(ind_d, i, j+1, 1);
                if(gsl_vector_get(survTime, i) > gsl_vector_get(s, j)) gsl_matrix_set(ind_r, i, j+1, 1);
            }
            if(gsl_vector_get(survTime, i) > 0 & gsl_vector_get(survTime, i) <= gsl_vector_get(s, 0)) gsl_matrix_set(ind_d, i, 0, 1);
        }
        
        if(gsl_vector_get(case0yleq, i) == 1)
        {
            for(j = 0; j < J; j++)
            {
                if(gsl_vector_get(survTime, i) > gsl_vector_get(s, j)) gsl_matrix_set(ind_r, i, j+1, 1);
            }
        }
        
        if(gsl_vector_get(case0ygeq, i) == 1 | gsl_vector_get(case1ygeq, i) == 1)
        {
            for(j = 0; j < J+1; j++)
            {
                gsl_matrix_set(ind_r, i, j, 1);
            }
        }
        gsl_matrix_set(ind_r, i, 0, 1);
    }
    
    for(j = 0; j < J+1; j++)
    {
        for(i = 0; i < n; i++)
        {
            gsl_vector_set(nEvent, j, gsl_vector_get(nEvent, j)+gsl_matrix_get(ind_d, i, j));
        }
    }
    
    return;
}






void cal_Delta(gsl_matrix *Delta,
               gsl_vector *survTime,
               gsl_vector *s,
               int J)
{
    int i, j;
    int n = survTime -> size;
    
    for(i = 0; i < n; i++)
    {
        for(j = 1; j < J+1; j++)
        {
            gsl_matrix_set(Delta, i, j, c_max( (c_min(gsl_vector_get(survTime, i), gsl_vector_get(s, j)) - gsl_vector_get(s, j-1)), 0));
        }
        gsl_matrix_set(Delta, i, 0, c_max( (c_min(gsl_vector_get(survTime, i), gsl_vector_get(s, 0)) - 0), 0));
    }
    
    return;
}





/*
 Evaluate the quadratic form: v^T M^{-1} v
 */
void c_quadform_vMv(gsl_vector *v,
                    gsl_matrix *Minv,
                    double     *value)
{
    int    d = v->size;
    gsl_vector *tempVec = gsl_vector_calloc(d);
    
    gsl_blas_dsymv(CblasUpper, 1, Minv, v, 0, tempVec);
    gsl_blas_ddot(v, tempVec, value);
    
    gsl_vector_free(tempVec);
    return;
}



/*
 Evaluate the quadratic form: v^T M^{-1} u
 - note v and u are assumed to be of the same length
 */
void c_quadform_vMu(gsl_vector *v,
                    gsl_matrix *Minv,
                    gsl_vector *u,
                    double     *value)
{
    int    d = v->size;
    gsl_vector *tempVec = gsl_vector_calloc(d);
    
    gsl_blas_dsymv(CblasUpper, 1, Minv, u, 0, tempVec);
    gsl_blas_ddot(v, tempVec, value);
    
    gsl_vector_free(tempVec);
    return;
}





/*
 Density calculation for a multivariate normal distribution
 */
void c_dmvnorm(gsl_vector *x,
               double     mu,
               double     sigma,
               gsl_matrix *AInv,
               double     *value)
{
    int signum, K = x->size;
	double sigmaSqInv = pow(sigma, -2);
    
	gsl_vector *muVec      = gsl_vector_alloc(K);
    gsl_vector *diff       = gsl_vector_alloc(K);
	gsl_matrix *SigmaInv   = gsl_matrix_alloc(K, K);
    gsl_matrix *SigmaInvLU = gsl_matrix_alloc(K, K);
    gsl_permutation *p     = gsl_permutation_alloc(K);
    
	gsl_vector_set_all(muVec, mu);
	gsl_vector_memcpy(diff, x);
	gsl_vector_sub(diff, muVec);
    
	gsl_matrix_memcpy(SigmaInv, AInv);
	gsl_matrix_scale(SigmaInv, sigmaSqInv);
    gsl_matrix_memcpy(SigmaInvLU, SigmaInv);
    gsl_linalg_LU_decomp(SigmaInvLU, p, &signum);
    
	c_quadform_vMv(diff, SigmaInv, value);
	*value = (log(gsl_linalg_LU_det(SigmaInvLU, signum)) - log(pow(2*Pi, K)) - *value) / 2;
    
	gsl_vector_free(muVec);
	gsl_vector_free(diff);
	gsl_matrix_free(SigmaInv);
	gsl_matrix_free(SigmaInvLU);
	gsl_permutation_free(p);
    return;
}



/*
 Evaluate the inverse of the matrix X
 */
void matrixInv(gsl_matrix *X, gsl_matrix *Xinv)
{
    int signum;
	int d = X->size1;
    gsl_matrix      *XLU = gsl_matrix_calloc(d, d);
    gsl_permutation *p   = gsl_permutation_alloc(d);
    
    gsl_matrix_memcpy(XLU, X);
    gsl_linalg_LU_decomp(XLU, p, &signum);
    gsl_linalg_LU_invert(XLU, p, Xinv);
    
    gsl_matrix_free(XLU);
    gsl_permutation_free(p);
    return;
}


/*
 Calculating column sums of matrix X
 */
void c_colSums(gsl_matrix *X, gsl_vector *v)
{
    int numCol = X->size2;
    int numRow = X->size1;    
    int i, j;
    double sum = 0;
    for(j = 0; j < numCol; j++)
    {
        i = 0;
        while(i < numRow)
        {
            sum = sum + gsl_matrix_get(X, i, j);
            i++;
        }
        gsl_vector_set(v, j, sum);
        sum = 0;
    }
    return;
}


/*
 Calculating row sums of matrix X
 */
void c_rowSums(gsl_matrix *X, gsl_vector *v)
{
    int numCol = X->size2;
    int numRow = X->size1;    
    int i, j;
    double sum = 0;
    for(i = 0; i < numRow; i++)
    {
        j = 0;
        while(j < numCol)
        {
            sum = sum + gsl_matrix_get(X, i, j);
            j++;
        }
        gsl_vector_set(v, i, sum);
        sum = 0;
    }
    return;
}


/*
 Replicate a vector v into rows of a matrix X
 */
void c_repVec_Rowmat(gsl_vector *v, gsl_matrix *X)
{
    int length = v->size;
    int numRep = X->size1;
    int i, j;
    for(i = 0; i < numRep; i++)
    {
        for(j = 0; j < length; j++)
        {
            gsl_matrix_set(X, i, j, gsl_vector_get(v, j));
        }
    }
    return;
}



/*
 Replicate a vector v into columns of a matrix X
 */
void c_repVec_Colmat(gsl_vector *v, gsl_matrix *X)
{
    int length = v->size;
    int numRep = X->size2;
    int i, j;
    for(j = 0; j < numRep; j++)
    {
        for(i = 0; i < length; i++)
        {
            gsl_matrix_set(X, i, j, gsl_vector_get(v, i));
        }
    }
    return;
}




/*
 Minimum of two numbers
 */
double c_min(double value1,
             double value2)
{
    double min = (value1 <= value2) ? value1 : value2;
    return min;
}



/*
 Maximum of two numbers
 */
double c_max(double value1,
             double value2)
{
    double max = (value1 >= value2) ? value1 : value2;
    return max;
}




/*
 Evaluate the inverse of the matrix M
 */
void c_solve(gsl_matrix *M,
             gsl_matrix *Minv)
{
    int signum;
	int d = M->size1;
    gsl_matrix      *MLU = gsl_matrix_calloc(d, d);
    gsl_permutation *p   = gsl_permutation_alloc(d);
    
    gsl_matrix_memcpy(MLU, M);
    gsl_linalg_LU_decomp(MLU, p, &signum);
    gsl_linalg_LU_invert(MLU, p, Minv);
    
    gsl_matrix_free(MLU);
    gsl_permutation_free(p);
    return;
}











/*********************************************************************************
 *********************************************************************************
 *********************************************************************************
 *********************************************************************************
 *********************************************************************************
 *********************************************************************************
 *********************************************************************************
 *********************************************************************************
 Functions for models in BayesID_AFT  */




double logistic(double x)
{
    double value = exp(x)/(1+exp(x));
    return value;
}





void set_r3_mu3_zeta3(gsl_vector *r3,
                      gsl_vector *mu3_vec,
                      gsl_vector *zeta3_vec,
                      gsl_vector *mu3_all,
                      gsl_vector *zeta3_all,
                      double y1,
                      double y2,
                      gsl_vector *c0_neginf,
                      gsl_vector *xbeta3,
                      gsl_vector *gamma,
                      gsl_vector *r3Uniq,
                      gsl_vector *r3Uniq_count,
                      int i,
                      int u3,
                      double mu03,
                      double sigSq03,
                      double a03,
                      double b03,
                      double tau3,
                      gsl_rng *rr)
{
    int j, k, r_ind, n_ir;
    double mean, sd, yStar, eta, val, sum_prob, b_mc;
    double zeta, mu;
    
    int n = r3 -> size;
    int u = u3;
    
    gsl_vector *prob = gsl_vector_calloc(u+1);
    
    zeta = INFINITY;
    while(zeta == INFINITY || isnan(zeta))
    {
        zeta = rgamma(a03, 1/b03);
    }
    
    mu = NAN;
    while(mu == INFINITY || isnan(mu))
    {
        mu = rnorm(mu03, sqrt(sigSq03));
    }
    
    eta = gsl_vector_get(xbeta3, i) + gsl_vector_get(gamma, i);
    yStar = y2 + log(1-exp(y1-y2));
    
    for(j = 0; j < u; j++)
    {
        n_ir = gsl_vector_get(r3Uniq_count, j);
        mean = gsl_vector_get(mu3_all, j) + eta;
        sd = sqrt(1/gsl_vector_get(zeta3_all, j));
        
        val = dnorm(yStar, mean, sd, 0);
        val *= (double) n_ir / (double)(n-1+tau3);
        gsl_vector_set(prob, j, val);
    }
    
    mean = mu + eta;
    sd = sqrt(1/zeta);
    
    val = dnorm(yStar, mean, sd, 0);
    val *= (double) tau3 / (double)(n-1+tau3);
    
    gsl_vector_set(prob, u, val);
    
    sum_prob = 0;
    for(k = 0; k < u+1; k++)
    {
        sum_prob += gsl_vector_get(prob, k);
    }
    
    if(sum_prob > pow(10, -300))
    {
        b_mc = 1/sum_prob;
        gsl_vector_scale(prob, b_mc);
    }else
    {
        for(k=0; k<u+1; k++)
        {
            gsl_vector_set(prob, k, (double) 1/(u+1));
        }
    }
    
    r_ind = c_multinom_sample(rr, prob, u+1);
    
    if(r_ind <= u)
    {
        gsl_vector_set(r3, i, gsl_vector_get(r3Uniq, r_ind-1));
        gsl_vector_set(mu3_vec, i, gsl_vector_get(mu3_all, r_ind-1));
        gsl_vector_set(zeta3_vec, i, gsl_vector_get(zeta3_all, r_ind-1));
    }else if(r_ind == u+1)
    {
        gsl_vector_set(r3, i, gsl_vector_max(r3Uniq)+1);
        gsl_vector_set(mu3_vec, i, mu);
        gsl_vector_set(zeta3_vec, i, zeta);
    }
    
    
    gsl_vector_free(prob);
    
    return;
}





/*
 Random generation from the Inverse Wishart distribution
 */

void c_riwishart_general(int v,
                         gsl_matrix *X_ori,
                         gsl_matrix *sample)
{
    int i, j, df;
    double normVal;
    int p = X_ori->size1;
    
    gsl_matrix *X = gsl_matrix_calloc(p, p);
    matrixInv(X_ori, X);
    
    gsl_matrix *cholX = gsl_matrix_calloc(p, p);
    gsl_matrix *ZZ = gsl_matrix_calloc(p, p);
    gsl_matrix *XX = gsl_matrix_calloc(p, p);
    gsl_matrix *KK = gsl_matrix_calloc(p, p);
    
    gsl_matrix_memcpy(cholX, X);
    gsl_linalg_cholesky_decomp(cholX);
    
    for(i = 0; i < p; i ++)
    {
        for(j = 0; j < i; j ++)
        {
            gsl_matrix_set(cholX, i, j, 0);
        }
    }
    
    for(i = 0; i < p; i++)
    {
        df = v - i;
        gsl_matrix_set(ZZ, i, i, sqrt(rchisq(df)));
    }
    
    for(i = 0; i < p; i++)
    {
        for(j = 0; j < i; j ++)
        {
            normVal = rnorm(0, 1);
            gsl_matrix_set(ZZ, i, j, normVal);
        }
    }
    
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, ZZ, cholX, 0, XX);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, XX, XX, 0, KK);
    matrixInv(KK, sample);
    
    gsl_matrix_free(X);
    gsl_matrix_free(cholX);
    gsl_matrix_free(XX);
    gsl_matrix_free(ZZ);
    gsl_matrix_free(KK);
    
}




/*
 Random generation from truncated normal distribution
 */
void c_rtnorm(double mean,
              double sd,
              double LL,
              double UL,
              int LL_neginf,
              int UL_posinf,
              double *value)
{
    int caseNO;
    int stop=0;
    double y, a, z, u, val, rho;
    
    LL = (LL-mean)/sd;
    UL = (UL-mean)/sd;
    
    if((LL_neginf == 1 && UL_posinf == 1) || (LL_neginf == 0 && LL < 0 && UL_posinf == 1) || (LL_neginf == 1 && UL >0 && UL_posinf == 0) || (UL_posinf == 0 && LL_neginf == 0 && LL <0 && UL > 0  && (UL-LL) > sqrt(2* Pi)))
    {
        caseNO = 1;
    }else if((LL >= 0 && UL > LL + 2*sqrt(exp(1))/(LL+sqrt(pow(LL,2)+4))*exp((2*LL-LL*sqrt(pow(LL, 2)+4))/4) && UL_posinf == 0 && LL_neginf == 0) || (LL >=0 && UL_posinf == 1 && LL_neginf == 0))
    {
        caseNO = 2;
    }else if((UL <= 0 && -LL > -UL + 2*sqrt(exp(1))/(-UL+sqrt(pow(UL, 2)+4))*exp((2*UL+UL*sqrt(pow(UL, 2)+4))/4) && LL_neginf == 0 && UL_posinf == 0) || (LL_neginf == 1 && UL <= 0 && UL_posinf == 0))
    {
        caseNO = 3;
    }else
    {
        caseNO = 4;
    }
    
    
    if(caseNO == 1)
    {
        while(stop == 0)
        {
            y = rnorm(0, 1);
            if(LL_neginf == 1 && UL_posinf == 1)
            {
                stop = 1;
                val = y;
            }else if(LL_neginf == 0 && UL_posinf == 1)
            {
                if(y > LL)
                {
                    stop = 1;
                    val = y;
                }
            }else if(LL_neginf == 1 && UL_posinf == 0)
            {
                if(y < UL)
                {
                    stop = 1;
                    val = y;
                }
            }else if(LL_neginf == 0 && UL_posinf == 0)
            {
                if(y > LL && y < UL)
                {
                    stop = 1;
                    val = y;
                }
            }
        }
    }
    if(caseNO == 2)
    {
        while(stop == 0)
        {
            a = (LL + sqrt(pow(LL, 2)+4))/2;
            z = rexp((double) 1/a) + LL;
            u = runif(0, 1);
            
            if(LL_neginf == 0 && UL_posinf == 0)
            {
                if(u <= exp(-pow(z-a, 2)/2) && z <= UL)
                {
                    stop = 1;
                    val = z;
                }
            }else if(LL_neginf == 0 && UL_posinf == 1)
            {
                if(u <= exp(-pow(z-a, 2)/2))
                {
                    stop = 1;
                    val = z;
                }
            }
        }
    }
    if(caseNO == 3)
    {
        while(stop == 0)
        {
            a = (-UL + sqrt(pow(UL, 2)+4))/2;
            z = rexp((double) 1/a) - UL;
            u = runif(0, 1);
            
            if(LL_neginf == 0 && UL_posinf == 0)
            {
                if(u <= exp(-pow(z-a, 2)/2) && z <= -LL)
                {
                    stop = 1;
                    val = -z;
                }
            }else if(LL_neginf == 1 && UL_posinf == 0)
            {
                if(u <= exp(-pow(z-a, 2)/2))
                {
                    stop = 1;
                    val = -z;
                }
            }
        }
    }
    if(caseNO == 4)
    {
        while(stop == 0)
        {
            if(LL_neginf == 0 && UL_posinf == 0)
            {
                z = runif(LL, UL);
                if(LL >0)
                {
                    rho = exp((pow(LL, 2) - pow(z, 2))/2);
                }else if(UL <0)
                {
                    rho = exp((pow(UL, 2) - pow(z, 2))/2);
                }else
                {
                    rho = exp(- pow(z, 2)/2);
                }
                u = runif(0, 1);
                
                if(u <= rho)
                {
                    stop = 1;
                    val = z;
                }
            }
        }
    }
    *value = mean + val * sd;
    return;
}





/* Evaluating the joint probability density function of (y1, y2) - upper wedge */

void log_Jpdf_Upper_BAFT_DP(int i,
                            double y1,
                            double y2,
                            double LT,
                            gsl_vector *c0_neginf,
                            gsl_matrix *X1,
                            gsl_matrix *X2,
                            gsl_matrix *X3,
                            gsl_vector *beta1,
                            gsl_vector *beta2,
                            gsl_vector *beta3,
                            gsl_vector *gamma,
                            gsl_vector *mu1_all,
                            gsl_vector *mu2_all,
                            gsl_vector *mu3_all,
                            gsl_vector *zeta1_all,
                            gsl_vector *zeta2_all,
                            gsl_vector *zeta3_all,
                            gsl_vector *r1Uniq,
                            gsl_vector *r2Uniq,
                            gsl_vector *r3Uniq,
                            gsl_vector *r1Uniq_count,
                            gsl_vector *r2Uniq_count,
                            gsl_vector *r3Uniq_count,
                            int u1,
                            int u2,
                            int u3,
                            double *value)
{
    double gam;
    double logf1_t1=0;
    double logS2_t1=0;
    double logS1_LT=0;
    double logS2_LT=0;
    double logf3_gap = 0;
    double eta1, eta2, eta3;
    double val = 0;
    double dummy = 0;
    
    double yStar=0;
    
    
    gam = gsl_vector_get(gamma, i);
    
    gsl_vector_view X1row = gsl_matrix_row(X1, i);
    gsl_vector_view X2row = gsl_matrix_row(X2, i);
    gsl_vector_view X3row = gsl_matrix_row(X3, i);
    
    gsl_blas_ddot(&X1row.vector, beta1, &eta1);
    gsl_blas_ddot(&X2row.vector, beta2, &eta2);
    gsl_blas_ddot(&X3row.vector, beta3, &eta3);
    
    log_fg_BAFT_DP(y1, u1, eta1, gam, mu1_all, zeta1_all, r1Uniq, r1Uniq_count, 1, 0, &logf1_t1, &dummy);
    
    if(y1 == y2)
    {
        logf3_gap = 0;
    }else if(y2 > y1)
    {
        yStar =y2+log(1-exp(y1-y2));
        log_fg_BAFT_DP(yStar, u3, eta3, gam, mu3_all, zeta3_all, r3Uniq, r3Uniq_count, 1, 0, &logf3_gap, &dummy);
    }
    
    log_fg_BAFT_DP(y1, u2, eta2, gam, mu2_all, zeta2_all, r2Uniq, r2Uniq_count, 0, 1, &dummy, &logS2_t1);
    
    if(gsl_vector_get(c0_neginf, i) == 0)
    {
        log_fg_BAFT_DP(LT, u1, eta1, gam, mu1_all, zeta1_all, r1Uniq, r1Uniq_count, 0, 1, &dummy, &logS1_LT);
        
        log_fg_BAFT_DP(LT, u2, eta2, gam, mu2_all, zeta2_all, r2Uniq, r2Uniq_count, 0, 1, &dummy, &logS2_LT);
    }
    
    val = (double) logf1_t1 + (double) logf3_gap + (double) logS2_t1 - (double) logS1_LT - (double) logS2_LT + y2 - yStar;
    
    *value = val;
    return;
}





/* Evaluating the joint probability density function of (y1, y2) - lower wedge */

void log_Jpdf_Lower_BAFT_DP(int i,
                            double y2,
                            double LT,
                            gsl_vector *c0_neginf,
                            gsl_matrix *X1,
                            gsl_matrix *X2,
                            gsl_vector *beta1,
                            gsl_vector *beta2,
                            gsl_vector *gamma,
                            gsl_vector *mu1_all,
                            gsl_vector *mu2_all,
                            gsl_vector *zeta1_all,
                            gsl_vector *zeta2_all,
                            gsl_vector *r1Uniq,
                            gsl_vector *r2Uniq,
                            gsl_vector *r1Uniq_count,
                            gsl_vector *r2Uniq_count,
                            int u1,
                            int u2,
                            double *value)
{
    double gam;
    double logf2_t2=0;
    double logS1_t2=0;
    double logS1_LT=0;
    double logS2_LT=0;
    double eta1, eta2;
    double val = 0;
    double dummy = 0;
    
    gam = gsl_vector_get(gamma, i);
    
    gsl_vector_view X1row = gsl_matrix_row(X1, i);
    gsl_vector_view X2row = gsl_matrix_row(X2, i);
    
    gsl_blas_ddot(&X1row.vector, beta1, &eta1);
    gsl_blas_ddot(&X2row.vector, beta2, &eta2);
    
    
    log_fg_BAFT_DP(y2, u2, eta2, gam, mu2_all, zeta2_all, r2Uniq, r2Uniq_count, 1, 0, &logf2_t2, &dummy);
    log_fg_BAFT_DP(y2, u1, eta1, gam, mu1_all, zeta1_all, r1Uniq, r1Uniq_count, 0, 1, &dummy, &logS1_t2);
    
    if(gsl_vector_get(c0_neginf, i) == 0)
    {
        log_fg_BAFT_DP(LT, u1, eta1, gam, mu1_all, zeta1_all, r1Uniq, r1Uniq_count, 0, 1, &dummy, &logS1_LT);
        log_fg_BAFT_DP(LT, u2, eta2, gam, mu2_all, zeta2_all, r2Uniq, r2Uniq_count, 0, 1, &dummy, &logS2_LT);
    }
    
    
    val = (double) logf2_t2 + (double) logS1_t2 - (double) logS1_LT - (double) logS2_LT;
    
    *value = val;
    return;
    
}






/*
 Random number generation for inverse gamma distribution
 alpha = shape, beta = rate
 
 */
void c_rigamma(double *temp,
               double alpha,
               double beta)
{
    double shape = alpha;
    double scale = (double) 1/beta;
    double gam=0;
    
    if(alpha > 0 && beta > 0){
        gam = rgamma(shape, scale);
    }
    *temp = (double) 1/gam;
    return;
}






double Qfunc_BAFT_DP(double V,
                     double mu0,
                     double zeta0,
                     double a0,
                     double b0)
{
    double lterm1, term2, term3, val;
    
    lterm1 = gsl_sf_lngamma(a0+0.5) - gsl_sf_lngamma(a0);
    term2 = sqrt(zeta0/(2*Pi*b0*(zeta0 + 1)));
    term3 = pow(zeta0*pow(V - mu0, 2)/(2*b0*(zeta0 + 1)) + 1, -a0-0.5);
    
    val = exp(lterm1)*term2*term3;
    
    return val;
}



int MFunc_BAFT_DP(double yL,
                  double yU,
                  int yL_neginf,
                  int yU_posinf,
                  gsl_vector *mu_all,
                  gsl_vector *zeta_all,
                  gsl_vector *rUniq_count,
                  int u,
                  double eta,
                  gsl_rng *rr)
{
    int k, selInd;
    double mixProb_den, mu_k, sigSq_k;
    
    gsl_vector *mixProb = gsl_vector_calloc(u);
    
    mixProb_den = 0;
    for(k=0; k<u; k++)
    {
        mu_k = gsl_vector_get(mu_all, k);
        sigSq_k = 1/gsl_vector_get(zeta_all, k);
        
        if(yU_posinf == 1)
        {
            if(yL_neginf == 1)
            {
                gsl_vector_set(mixProb, k, gsl_vector_get(rUniq_count, k));
            }else
            {
                gsl_vector_set(mixProb, k, gsl_vector_get(rUniq_count, k) * (1 - pnorm(yL, eta + mu_k, sqrt(sigSq_k), 1, 0)));
            }
        }else
        {
            if(yL_neginf == 1)
            {
                gsl_vector_set(mixProb, k, gsl_vector_get(rUniq_count, k) * (pnorm(yU, eta + mu_k, sqrt(sigSq_k), 1, 0)));
            }else
            {
                gsl_vector_set(mixProb, k, gsl_vector_get(rUniq_count, k) * (pnorm(yU, eta + mu_k, sqrt(sigSq_k), 1, 0) - pnorm(yL, eta + mu_k, sqrt(sigSq_k), 1, 0)));
            }
            
        }
        mixProb_den += gsl_vector_get(mixProb, k);
    }
    
    if(mixProb_den == 0)
    {
        for(k=0; k<u; k++)
        {
            gsl_vector_set(mixProb, k, (double) 1/u);
        }
    }else
    {
        gsl_vector_scale(mixProb, 1/mixProb_den);
    }
    
    selInd = c_multinom_sample(rr, mixProb, u);
    
    gsl_vector_free(mixProb);
    
    return selInd;
}




void log_fg_BAFT_DP(double y,
                    int u,
                    double xbeta,
                    double gam,
                    gsl_vector *mu_all,
                    gsl_vector *zeta_all,
                    gsl_vector *rUniq,
                    gsl_vector *rUniq_count,
                    int calf,
                    int calS,
                    double *logfg,
                    double *logSg)
{
    int k;
    double mu_k, sigSq_k, mixProb_den, val1, val2, log1, logk, temp;
    
    
    gsl_vector *mixProb = gsl_vector_calloc(u);
    
    mixProb_den=0;
    for(k=0;k<u;k++)
    {
        gsl_vector_set(mixProb, k, gsl_vector_get(rUniq_count, k));
        mixProb_den += gsl_vector_get(mixProb, k);
    }
    gsl_vector_scale(mixProb, (double) 1/mixProb_den);
    
    val1=0;
    val2=0;
    log1=0;
    logk=0;
    temp=0;
    for(k=0;k<u;k++)
    {
        mu_k = (double) gsl_vector_get(mu_all, k);
        sigSq_k = (double) pow(gsl_vector_get(zeta_all, k), -1);
        
        if(calf == 1)
        {
            val1 += gsl_vector_get(mixProb, k) * dnorm(y, xbeta+mu_k+gam, sqrt(sigSq_k), 0);
            
        }
        if(calS == 1)
        {
            val2 += gsl_vector_get(mixProb, k) * pnorm(y, xbeta+mu_k+gam, sqrt(sigSq_k), 0, 0);
        }
        
    }
    
    if(calf == 1)
    {
        if(val1 == 0)
        {
            *logfg = -600;
        }else
        {
            *logfg = log(val1);
        }
    }
    if(calS == 1)
    {
        if(val2 == 0)
        {
            *logSg = -600;
        }else
        {
            *logSg = log(val2);
        }
    }
    
    gsl_vector_free(mixProb);
    return;
}









/*
 Identify unique lables for cluster memberships:
 - r, mu_all, zeta_all musted be the updated one
 - rUniq, rUniq_count, u are calculated
 */
void c_uniq(gsl_vector *r,
            gsl_vector *rUniq,
            gsl_vector *rUniq_count,
            gsl_vector *mu_all,
            gsl_vector *zeta_all,
            int *u)
{
    int i, j, k;
    
    int u_old = *u;
    int n = r -> size;
    
    gsl_vector *rTemp = gsl_vector_calloc(n);
    gsl_vector *muTemp = gsl_vector_calloc(n);
    gsl_vector *zetaTemp = gsl_vector_calloc(n);
    gsl_vector *rUniq_old = gsl_vector_calloc(n);
    
    gsl_vector_memcpy(rUniq_old, rUniq);
    
    gsl_vector_set_zero(rUniq);
    gsl_vector_set_zero(rUniq_count);
    
    gsl_vector_memcpy(rTemp, r);
    *u = 0;
    for(i = 0; i < n; i++)
    {
        if(gsl_vector_get(rTemp, i) != 0)
        {
            *u += 1;
            gsl_vector_set(rUniq, *u-1, gsl_vector_get(rTemp, i));
            for(k=0; k < u_old; k++)
            {
                if(gsl_vector_get(rUniq_old, k) == gsl_vector_get(rUniq, *u-1))
                {
                    gsl_vector_set(muTemp, *u-1, gsl_vector_get(mu_all, k));
                    gsl_vector_set(zetaTemp, *u-1, gsl_vector_get(zeta_all, k));
                }
            }
            
            
            for(j = i; j < n; j++)
            {
                if(gsl_vector_get(rTemp, j) == gsl_vector_get(rUniq, *u-1))
                {
                    gsl_vector_set(rUniq_count, *u-1, gsl_vector_get(rUniq_count, *u-1)+1);
                    gsl_vector_set(rTemp, j, 0);
                }
            }
        }
    }
    
    gsl_vector_memcpy(mu_all, muTemp);
    gsl_vector_memcpy(zeta_all, zetaTemp);
    
    gsl_vector_free(rTemp);
    gsl_vector_free(muTemp);
    gsl_vector_free(zetaTemp);
    gsl_vector_free(rUniq_old);
    
    return;
}


/*
 Identify unique lables for cluster memberships:
 - r, mu_vec, zeta_vec musted be the updated one
 - rUniq, rUniq_count, u, mu_all, zeta_all are calculated
 */
void c_uniq_h3(gsl_vector *r,
               gsl_vector *rUniq,
               gsl_vector *rUniq_count,
               gsl_vector *mu_all,
               gsl_vector *zeta_all,
               gsl_vector *mu3_vec,
               gsl_vector *zeta3_vec,
               gsl_vector *y1_NA,
               int *u)
{
    int i, j;
    int n = r -> size;
    
    gsl_vector *rTemp = gsl_vector_calloc(n);
    gsl_vector *muTemp = gsl_vector_calloc(n);
    gsl_vector *zetaTemp = gsl_vector_calloc(n);
    
    gsl_vector_set_zero(rUniq);
    gsl_vector_set_zero(rUniq_count);
    
    gsl_vector_memcpy(rTemp, r);
    *u = 0;
    for(i = 0; i < n; i++)
    {
        if(gsl_vector_get(y1_NA, i) != 1)
        {
            if(gsl_vector_get(rTemp, i) != 0)
            {
                *u += 1;
                gsl_vector_set(rUniq, *u-1, gsl_vector_get(rTemp, i));
                gsl_vector_set(muTemp, *u-1, gsl_vector_get(mu3_vec, i));
                gsl_vector_set(zetaTemp, *u-1, gsl_vector_get(zeta3_vec, i));
                
                for(j = i; j < n; j++)
                {
                    if(gsl_vector_get(rTemp, j) == gsl_vector_get(rUniq, *u-1))
                    {
                        gsl_vector_set(rUniq_count, *u-1, gsl_vector_get(rUniq_count, *u-1)+1);
                        gsl_vector_set(rTemp, j, 0);
                    }
                }
            }
        }else
        {
            gsl_vector_set(r, i, 0);
            gsl_vector_set(mu3_vec, i, -exp(100000));
            gsl_vector_set(zeta3_vec, i, -exp(100000));
        }
    }
    
    gsl_vector_memcpy(mu_all, muTemp);
    gsl_vector_memcpy(zeta_all, zetaTemp);
    
    gsl_vector_free(rTemp);
    gsl_vector_free(muTemp);
    gsl_vector_free(zetaTemp);
    
    
    return;
}




/*
 Identify unique lables for cluster memberships for the first iteration
 : r must be the updated one --> set rUniq, rUniq_count, nClass
 */
void c_uniq1(gsl_vector *r,
             gsl_vector *rUniq,
             gsl_vector *rUniq_count,
             int *u)
{
    int i, j;
    int n = r -> size;
    
    gsl_vector *rTemp = gsl_vector_calloc(n);
    
    gsl_vector_set_zero(rUniq);
    gsl_vector_set_zero(rUniq_count);
    
    gsl_vector_memcpy(rTemp, r);
    *u = 1;
    for(i = 0; i < n; i++)
    {
        if(i == 0 && gsl_vector_get(rTemp, i) != 0)
        {
            gsl_vector_set(rUniq, 0, gsl_vector_get(rTemp, 0));
            for(j = i; j < n; j++)
            {
                if(gsl_vector_get(rTemp, j) == gsl_vector_get(rUniq, 0))
                {
                    gsl_vector_set(rUniq_count, 0, gsl_vector_get(rUniq_count, 0)+1);
                    gsl_vector_set(rTemp, j, 0);
                }
            }
        }
        
        if(i != 0 && gsl_vector_get(rTemp, i) != 0)
        {
            *u += 1;
            gsl_vector_set(rUniq, *u-1, gsl_vector_get(rTemp, i));
            
            for(j = i; j < n; j++)
            {
                if(gsl_vector_get(rTemp, j) == gsl_vector_get(rUniq, *u-1))
                {
                    gsl_vector_set(rUniq_count, *u-1, gsl_vector_get(rUniq_count, *u-1)+1);
                    gsl_vector_set(rTemp, j, 0);
                }
            }
        }
    }
    
    gsl_vector_free(rTemp);
    
    return;
}



/*
 Identify unique lables for cluster memberships for the first iteration
 : r and y1_NA must be the updated one --> set rUniq, rUniq_count, nClass
 */
void c_uniq1_h3(gsl_vector *r,
                gsl_vector *rUniq,
                gsl_vector *rUniq_count,
                gsl_vector *y1_NA,
                int *u)
{
    int i, j;
    int n = r -> size;
    
    gsl_vector *rTemp = gsl_vector_calloc(n);
    
    gsl_vector_set_zero(rUniq);
    gsl_vector_set_zero(rUniq_count);
    
    gsl_vector_memcpy(rTemp, r);
    *u = 0;
    for(i = 0; i < n; i++)
    {
        if(gsl_vector_get(y1_NA, i) == 0)
        {
            if(gsl_vector_get(rTemp, i) != 0)
            {
                *u += 1;
                gsl_vector_set(rUniq, *u-1, gsl_vector_get(rTemp, i));
                
                for(j = i; j < n; j++)
                {
                    if(gsl_vector_get(rTemp, j) == gsl_vector_get(rUniq, *u-1))
                    {
                        gsl_vector_set(rUniq_count, *u-1, gsl_vector_get(rUniq_count, *u-1)+1);
                        gsl_vector_set(rTemp, j, 0);
                    }
                }
            }
        }else
        {
            gsl_vector_set(r, i, 0);
        }
    }
    
    gsl_vector_free(rTemp);
    
    return;
}





















/* Evaluating the joint probability density function of (y1, y2) - upper wedge */

void log_Jpdf_Upper_BAFT_LN(int i,
                            double y1,
                            double y2,
                            double LT,
                            gsl_vector *c0_neginf,
                            gsl_matrix *X1,
                            gsl_matrix *X2,
                            gsl_matrix *X3,
                            gsl_vector *beta1,
                            gsl_vector *beta2,
                            gsl_vector *beta3,
                            gsl_vector *gamma,
                            double beta01,
                            double beta02,
                            double beta03,
                            double sigSq1,
                            double sigSq2,
                            double sigSq3,
                            double *value)
{
    double eta1, eta2, eta3, gam, y2_temp;
    double val = 0;
    y2_temp = 0;
    
    gam = gsl_vector_get(gamma, i);
    
    gsl_vector_view X1row = gsl_matrix_row(X1, i);
    gsl_vector_view X2row = gsl_matrix_row(X2, i);
    gsl_vector_view X3row = gsl_matrix_row(X3, i);
    
    gsl_blas_ddot(&X1row.vector, beta1, &eta1);
    gsl_blas_ddot(&X2row.vector, beta2, &eta2);
    gsl_blas_ddot(&X3row.vector, beta3, &eta3);
    
    eta1 += beta01 + gam;
    eta2 += beta02 + gam;
    eta3 += beta03 + gam;
    
    /* without left-truncation */
    val += dnorm(y1, eta1, sqrt(sigSq1), 1);
    val += pnorm(y1, eta2, sqrt(sigSq2), 0, 1);
    
    /*     */
    val += y2;
    
    if(y1 < y2)
    {
        val += dnorm(log(exp(y2)-exp(y1)), eta3, sqrt(sigSq3), 1);
        
        /*         */
        val -= y2 + log(1-exp(y1-y2));
        
    }
    
    if(gsl_vector_get(c0_neginf, i) == 0)
    {
        /* for left-truncation */
        val -= pnorm(LT, eta1, sqrt(sigSq1), 0, 1);
        val -= pnorm(LT, eta2, sqrt(sigSq2), 0, 1);
    }
    
    *value = val;
    return;
}


/* Evaluating the joint probability density function of (y1, y2) - lower wedge */

void log_Jpdf_Lower_BAFT_LN(int i,
                            double y2,
                            double LT,
                            gsl_vector *c0_neginf,
                            gsl_matrix *X1,
                            gsl_matrix *X2,
                            gsl_vector *beta1,
                            gsl_vector *beta2,
                            gsl_vector *gamma,
                            double beta01,
                            double beta02,
                            double sigSq1,
                            double sigSq2,
                            double *value)
{
    double eta1, eta2, gam;
    double val = 0;
    
    gam = gsl_vector_get(gamma, i);
    
    gsl_vector_view X1row = gsl_matrix_row(X1, i);
    gsl_vector_view X2row = gsl_matrix_row(X2, i);
    
    gsl_blas_ddot(&X1row.vector, beta1, &eta1);
    gsl_blas_ddot(&X2row.vector, beta2, &eta2);
    
    eta1 += beta01 + gam;
    eta2 += beta02 + gam;
    
    /* without left-truncation */
    val += dnorm(y2, eta2, sqrt(sigSq2), 1);
    val += pnorm(y2, eta1, sqrt(sigSq1), 0, 1);
    
    if(gsl_vector_get(c0_neginf, i) == 0)
    {
        /* for left-truncation */
        val -= pnorm(LT, eta1, sqrt(sigSq1), 0, 1);
        val -= pnorm(LT, eta2, sqrt(sigSq2), 0, 1);
    }
    
    *value = val;
    return;
}



/*
 Evaluate density function for mixture of truncated normals
 */
double fmixTN(double y,
              double yL,
              double yU,
              int yL_neginf,
              int yU_posinf,
              gsl_vector *mu_all,
              gsl_vector *zeta_all,
              gsl_vector *rUniq_count,
              int u,
              double eta)
{
    int k;
    double mixProb_den, mu_k, sigSq_k, temp;
    double TN_den, TN_num;
    
    gsl_vector *mixProb = gsl_vector_calloc(u);
    
    mixProb_den = 0;
    for(k=0;k<u;k++)
    {
        mu_k = (double) gsl_vector_get(mu_all, k);
        sigSq_k = (double) pow(gsl_vector_get(zeta_all, k), -1);
        
        if(yL_neginf == 0 && yU_posinf == 0)
        {
            gsl_vector_set(mixProb, k, gsl_vector_get(rUniq_count, k) * (pnorm(yU, eta + mu_k, sqrt(sigSq_k), 1, 0) - pnorm(yL, eta + mu_k, sqrt(sigSq_k), 1, 0)));
        }else if(yL_neginf == 0 && yU_posinf == 1)
        {
            gsl_vector_set(mixProb, k, gsl_vector_get(rUniq_count, k) * pnorm(yL, eta + mu_k, sqrt(sigSq_k), 0, 0));
        }else if(yL_neginf == 1 && yU_posinf == 0)
        {
            gsl_vector_set(mixProb, k, gsl_vector_get(rUniq_count, k) * pnorm(yU, eta + mu_k, sqrt(sigSq_k), 1, 0));
        }else if(yL_neginf == 1 && yU_posinf == 1)
        {
            gsl_vector_set(mixProb, k, gsl_vector_get(rUniq_count, k));
        }
        
        mixProb_den += gsl_vector_get(mixProb, k);
        
    }
    gsl_vector_scale(mixProb, (double) 1/mixProb_den);
    
    
    temp = 0;
    for(k=0;k<u;k++)
    {
        mu_k = (double) gsl_vector_get(mu_all, k);
        sigSq_k = (double) pow(gsl_vector_get(zeta_all, k), -1);
        
        
        if(yL_neginf == 0 && yU_posinf == 0)
        {
            TN_num = dnorm(y, eta + mu_k, sqrt(sigSq_k), 0);
            if(TN_num != 0)
            {
                TN_den = pnorm(yU, eta + mu_k, sqrt(sigSq_k), 1, 0)-pnorm(yL, eta + mu_k, sqrt(sigSq_k), 1, 0);
                if(TN_den !=0)
                {
                    temp += gsl_vector_get(mixProb, k) * TN_num/TN_den;
                }
            }
        }else if(yL_neginf == 0 && yU_posinf == 1)
        {
            TN_num = dnorm(y, eta + mu_k, sqrt(sigSq_k), 0);
            if(TN_num != 0)
            {
                TN_den = pnorm(yL, eta + mu_k, sqrt(sigSq_k), 0, 0);
                if(TN_den !=0)
                {
                    temp += gsl_vector_get(mixProb, k) * TN_num/TN_den;
                }
            }
        }else if(yL_neginf == 1 && yU_posinf == 0)
        {
            TN_num = dnorm(y, eta + mu_k, sqrt(sigSq_k), 0);
            if(TN_num != 0)
            {
                TN_den = pnorm(yU, eta + mu_k, sqrt(sigSq_k), 1, 0);
                if(TN_den !=0)
                {
                    temp += gsl_vector_get(mixProb, k) * TN_num/TN_den;
                }
            }
        }else if(yL_neginf == 1 && yU_posinf == 1)
        {
            temp += gsl_vector_get(mixProb, k) * dnorm(y, eta + mu_k, sqrt(sigSq_k), 0);
        }
        
    }
    
    gsl_vector_free(mixProb);
    return temp;
    
}



















