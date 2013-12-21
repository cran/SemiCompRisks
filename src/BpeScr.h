//
//  BpeScr.h
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

extern double c_min(double value1,
                    double value2);


extern double c_max(double value1,
                    double value2);

extern void c_solve(gsl_matrix *M,
                    gsl_matrix *Minv);

extern double Bscr_wFunc(int subjInx,
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
                    gsl_vector *survTime2);


extern void cal_Sigma(gsl_matrix *Sigma_lam,
                      gsl_matrix *invSigma_lam,
                      gsl_matrix *W,
                      gsl_matrix *Q,
                      gsl_vector *s,
                      double c_lam,
                      int J);


extern void c_quadform_vMv(gsl_vector *v,
                           gsl_matrix *Minv,
                           double     *value);


extern void c_quadform_vMu(gsl_vector *v,
                           gsl_matrix *Minv,
                           gsl_vector *u,
                           double     *value);


extern void c_dmvnorm(gsl_vector *x,
                      double     mu,
                      double     sigma,
                      gsl_matrix *AInv,
                      double     *value);



extern void Bscr_updateRP1(gsl_vector *beta1,
                      gsl_vector *xbeta1,
                      gsl_vector *accept_beta1,
                      gsl_vector *gamma,
                      gsl_vector *lambda1,
                      gsl_vector *s1,
                      gsl_vector *survTime1,
                      gsl_vector *survEvent1,
                      gsl_matrix *survCov1,
                      int J1);


extern void Bscr_updateRP2(gsl_vector *beta2,
                      gsl_vector *xbeta2,
                      gsl_vector *accept_beta2,
                      gsl_vector *gamma,
                      gsl_vector *lambda2,
                      gsl_vector *s2,
                      gsl_vector *survTime1,
                      gsl_vector *case01,
                      gsl_matrix *survCov2,
                      int J2);


extern void Bscr_updateRP3(gsl_vector *beta3,
                      gsl_vector *xbeta3,
                      gsl_vector *accept_beta3,
                      gsl_vector *gamma,
                      gsl_vector *lambda3,
                      gsl_vector *s3,
                      gsl_vector *survTime1,
                      gsl_vector *survTime2,
                      gsl_vector *case11,
                      gsl_matrix *survCov3,
                      int J3);




extern void Bscr_updateBH1(gsl_vector *lambda1,
                      gsl_vector *s1,
                      gsl_vector *xbeta1,
                      gsl_vector *gamma,                      
                      gsl_vector *survTime1,
                      gsl_vector *survEvent1,
                      gsl_matrix *Sigma_lam1,
                      gsl_matrix *invSigma_lam1,
                      gsl_matrix *W1,
                      gsl_matrix *Q1,
                      double mu_lam1,
                      double sigSq_lam1,
                      int J1);



extern void Bscr_updateBH2(gsl_vector *lambda2,
                      gsl_vector *s2,
                      gsl_vector *xbeta2,
                      gsl_vector *gamma,
                      gsl_vector *survTime1,
                      gsl_vector *survTime2,
                      gsl_vector *case01,
                      gsl_matrix *Sigma_lam2,
                      gsl_matrix *invSigma_lam2,
                      gsl_matrix *W2,
                      gsl_matrix *Q2,
                      double mu_lam2,
                      double sigSq_lam2,
                      int J2);




extern void Bscr_updateBH3(gsl_vector *lambda3,
                      gsl_vector *s3,
                      gsl_vector *xbeta3,
                      gsl_vector *gamma,
                      gsl_vector *survTime1,
                      gsl_vector *survTime2,
                      gsl_vector *case11,
                      gsl_matrix *Sigma_lam3,
                      gsl_matrix *invSigma_lam3,
                      gsl_matrix *W3,
                      gsl_matrix *Q3,
                      double mu_lam3,
                      double sigSq_lam3,
                      int J3);






extern void Bscr_updateSP(double *mu_lam,
                     double *sigSq_lam,
                     gsl_vector *lambda,
                     gsl_matrix *Sigma_lam,
                     gsl_matrix *invSigma_lam,
                     double a,
                     double b,
                     int J);



extern void Bscr_updateSP1(double *mu_lam1,
                      double *sigSq_lam1,
                      gsl_vector *lambda1,
                      gsl_matrix *Sigma_lam1,
                      gsl_matrix *invSigma_lam1,
                      double a1,
                      double b1,
                      int J1);



extern void Bscr_updateSP2(double *mu_lam2,
                      double *sigSq_lam2,
                      gsl_vector *lambda2,
                      gsl_matrix *Sigma_lam2,
                      gsl_matrix *invSigma_lam2,
                      double a2,
                      double b2,
                      int J2);


extern void Bscr_updateSP3(double *mu_lam3,
                      double *sigSq_lam3,
                      gsl_vector *lambda3,
                      gsl_matrix *Sigma_lam3,
                      gsl_matrix *invSigma_lam3,
                      double a3,
                      double b3,
                      int J3);


extern void Bscr_updateFP(gsl_vector *gamma,
                     double theta,
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
                     gsl_vector *survTime2,
                     gsl_vector *survEvent1,
                     gsl_vector *survEvent2);






extern void Bscr_updateDP(gsl_vector *gamma,
                     double *theta,
                     double mhProp_theta_var,
                     double psi,
                     double omega,
                     int *accept_theta);






extern void Bscr_updateBI1(gsl_vector *s1,
                      int *J1,
                      int *accept_BI1,
                      gsl_vector *survTime1,
                      gsl_vector *survEvent1,
                      gsl_vector *gamma,
                      gsl_vector *xbeta1,
                      gsl_matrix *Sigma_lam1,
                      gsl_matrix *invSigma_lam1,
                      gsl_matrix *W1,
                      gsl_matrix *Q1,
                      gsl_vector *lambda1,
                      gsl_vector *s_propBI1,
                      int num_s_propBI1,
                      double delPert1,
                      int alpha1,
                      double c_lam1,
                      double mu_lam1,
                      double sigSq_lam1,
                      double s1_max);






extern void Bscr_updateDI1(gsl_vector *s1,
                      int *J1,
                      int *accept_DI1,
                      gsl_vector *survTime1,
                      gsl_vector *survEvent1,
                      gsl_vector *gamma,
                      gsl_vector *xbeta1,
                      gsl_matrix *Sigma_lam1,
                      gsl_matrix *invSigma_lam1,
                      gsl_matrix *W1,
                      gsl_matrix *Q1,
                      gsl_vector *lambda1,
                      gsl_vector *s_propBI1,
                      int num_s_propBI1,
                      double delPert1,
                      int alpha1,
                      double c_lam1,
                      double mu_lam1,
                      double sigSq_lam1,
                      double s1_max,
                      int J1_max);




extern void Bscr_updateBI2(gsl_vector *s2,
                      int *J2,
                      int *accept_BI2,
                      gsl_vector *survTime1,
                      gsl_vector *survTime2,
                      gsl_vector *case01,
                      gsl_vector *gamma,
                      gsl_vector *xbeta2,
                      gsl_matrix *Sigma_lam2,
                      gsl_matrix *invSigma_lam2,
                      gsl_matrix *W2,
                      gsl_matrix *Q2,
                      gsl_vector *lambda2,
                      gsl_vector *s_propBI2,
                      int num_s_propBI2,
                      double delPert2,
                      int alpha2,
                      double c_lam2,
                      double mu_lam2,
                      double sigSq_lam2,
                      double s2_max);




extern void Bscr_updateDI2(gsl_vector *s2,
                      int *J2,
                      int *accept_DI2,
                      gsl_vector *survTime1,
                      gsl_vector *survTime2,
                      gsl_vector *case01,
                      gsl_vector *gamma,
                      gsl_vector *xbeta2,
                      gsl_matrix *Sigma_lam2,
                      gsl_matrix *invSigma_lam2,
                      gsl_matrix *W2,
                      gsl_matrix *Q2,
                      gsl_vector *lambda2,
                      gsl_vector *s_propBI2,
                      int num_s_propBI2,
                      double delPert2,
                      int alpha2,
                      double c_lam2,
                      double mu_lam2,
                      double sigSq_lam2,
                      double s2_max,
                      int J2_max);




extern void Bscr_updateBI3(gsl_vector *s3,
                      int *J3,
                      int *accept_BI3,
                      gsl_vector *survTime1,
                      gsl_vector *survTime2,
                      gsl_vector *case11,
                      gsl_vector *gamma,
                      gsl_vector *xbeta3,
                      gsl_matrix *Sigma_lam3,
                      gsl_matrix *invSigma_lam3,
                      gsl_matrix *W3,
                      gsl_matrix *Q3,
                      gsl_vector *lambda3,
                      gsl_vector *s_propBI3,
                      int num_s_propBI3,
                      double delPert3,
                      int alpha3,
                      double c_lam3,
                      double mu_lam3,
                      double sigSq_lam3,
                      double s3_max);





extern void Bscr_updateDI3(gsl_vector *s3,
                      int *J3,
                      int *accept_DI3,
                      gsl_vector *survTime1,
                      gsl_vector *survTime2,
                      gsl_vector *case11,
                      gsl_vector *gamma,
                      gsl_vector *xbeta3,
                      gsl_matrix *Sigma_lam3,
                      gsl_matrix *invSigma_lam3,
                      gsl_matrix *W3,
                      gsl_matrix *Q3,
                      gsl_vector *lambda3,
                      gsl_vector *s_propBI3,
                      int num_s_propBI3,
                      double delPert3,
                      int alpha3,
                      double c_lam3,
                      double mu_lam3,
                      double sigSq_lam3,
                      double s3_max,
                      int J3_max);




























