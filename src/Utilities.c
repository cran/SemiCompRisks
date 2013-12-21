


#include <stdio.h>
#include <math.h>

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_sort_vector.h"
#include "gsl/gsl_heapsort.h"

#include "R.h"
#include "Rmath.h"

#include "BpeScr.h"
#include "BpeScrSM.h"
#include "BpeSurv.h"
#include "BweibScr.h"
#include "BweibSurv.h"


#define Pi 3.141592653589793238462643383280


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
    
    /*
    printf("ch1 = %.3f\n", cumHaz1);
    printf("ch2 = %.3f\n", cumHaz2);
    printf("ch3 = %.3f\n", cumHaz3diff);
    */
    
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
    
    /*
     printf("ch1 = %.3f\n", cumHaz1);
     printf("ch2 = %.3f\n", cumHaz2);
     printf("ch3 = %.3f\n", cumHaz3diff);
     */
    
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








void cal_Sigma(gsl_matrix *Sigma_lam,
               gsl_matrix *invSigma_lam,
               gsl_matrix *W,
               gsl_matrix *Q,
               gsl_vector *s,
               double c_lam,
               int J)
{
    int i, j;
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








