//
//  BweibSurv.h
//  
//
//  Created by Kyu Ha Lee on 11/7/13.
//
//


extern void matrixInv(gsl_matrix *X, gsl_matrix *Xinv);

extern void c_colSums(gsl_matrix *X, gsl_vector *v);

extern void c_rowSums(gsl_matrix *X, gsl_vector *v);

extern void c_repVec_Rowmat(gsl_vector *v, gsl_matrix *X);

extern void c_repVec_Colmat(gsl_vector *v, gsl_matrix *X);

extern void BweibSurv_updateSH(gsl_vector *beta,
                     double *alpha,
                     double *kappa,
                     gsl_vector *survTime,
                     gsl_vector *survEvent,
                     gsl_matrix *survCov,
                     double c,
                     double d);

extern void BweibSurv_updateRP(gsl_vector *beta,
                     double *alpha,
                     double *kappa,
                     gsl_vector *survTime,
                     gsl_vector *survEvent,
                     gsl_matrix *survCov,
                     gsl_vector *accept_beta);

extern void BweibSurv_updateSC1(gsl_vector *beta,
                      double *alpha,
                      double *kappa,
                      gsl_vector *survTime,
                      gsl_vector *survEvent,
                      gsl_matrix *survCov,
                      double mhProp_alpha_var,
                      double a,
                      double b,
                      int *accept_alpha);

extern void BweibSurv_updateSC2(gsl_vector *beta,
                      double *alpha,
                      double *kappa,
                      gsl_vector *survTime,
                      gsl_vector *survEvent,
                      gsl_matrix *survCov,
                      double mhProp_alpha_var,
                      double a,
                      double b,
                      int *accept_alpha);
