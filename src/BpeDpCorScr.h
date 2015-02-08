extern void BpeDpCorScr_logMLH(gsl_vector *beta1,
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
                                double *val);


extern void BpeDpCorScr_logMLH_i(int i,
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
                                  double *val);


extern int c_multinom_sample(gsl_rng *rr,
                             gsl_vector *prob,
                             int length_prob);



extern double Qfunc(gsl_vector *V,
                    gsl_vector *mu0,
                    double zeta0,
                    gsl_matrix *Psi0,
                    double rho0);


extern void c_det(gsl_matrix *A, double *result);


extern void c_mgamma3(double a, double *val);



extern void BpeDpCorScr_logLH(gsl_vector *beta1,
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
                               double *val);


extern void BpeDpCorScr_logLH_i(int i,
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
                                double *val);

extern double c_min(double value1,
                    double value2);


extern double c_max(double value1,
                    double value2);

extern void matrixInv(gsl_matrix *X, gsl_matrix *Xinv);

extern void c_solve(gsl_matrix *M,
                    gsl_matrix *Minv);

extern void c_colSums(gsl_matrix *X, gsl_vector *v);

extern void c_rowSums(gsl_matrix *X, gsl_vector *v);

extern void c_repVec_Rowmat(gsl_vector *v, gsl_matrix *X);

extern void c_repVec_Colmat(gsl_vector *v, gsl_matrix *X);



extern void c_quadform_vMv(gsl_vector *v,
                           gsl_matrix *Minv,
                           double     *value);


extern void c_quadform_vMu(gsl_vector *v,
                           gsl_matrix *Minv,
                           gsl_vector *u,
                           double     *value);


extern void c_dmvnorm2(gsl_vector *x,
                      gsl_vector *mu,
                      double     sigma,
                      gsl_matrix *AInv,
                      double     *value);


extern void c_dmvnormSH(gsl_vector *x,
                      double     mu,
                      double     sigma,
                      gsl_matrix *AInv,
                      double     *value);


extern void c_rmvnorm(gsl_matrix *sample,
                      gsl_vector *mean,
                      gsl_matrix *Var);



extern void c_riwishart(int v,
                        gsl_matrix *X_ori,
                        gsl_matrix *sample);


extern void cal_Sigma(gsl_matrix *Sigma_lam,
                      gsl_matrix *invSigma_lam,
                      gsl_matrix *W,
                      gsl_matrix *Q,
                      gsl_vector *s,
                      double c_lam,
                      int J);



extern double BpeDpCorScr_wFunc(int subjInx,
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
                                 gsl_vector *survTime2);





extern void BpeDpCorScr_updateRP1(gsl_vector *beta1,
                                   gsl_vector *xbeta1,
                                   gsl_vector *gamma,
                                   gsl_vector *V1,
                                   gsl_vector *lambda1,
                                   gsl_vector *s1,
                                   gsl_vector *survTime1,
                                   gsl_vector *survEvent1,
                                   gsl_vector *cluster,
                                   gsl_matrix *survCov1,
                                   int K1,
                                   gsl_vector *accept_beta1);






extern void BpeDpCorScr_updateRP2(gsl_vector *beta2,
                                   gsl_vector *xbeta2,
                                   double nu2,
                                   gsl_vector *gamma,
                                   gsl_vector *V2,
                                   gsl_vector *lambda2,
                                   gsl_vector *s2,
                                   gsl_vector *survTime1,
                                   gsl_vector *case01,
                                   gsl_vector *cluster,
                                   gsl_matrix *survCov2,
                                   int K2,
                                   gsl_vector *accept_beta2);


extern void BpeDpCorScr_updateRP3(gsl_vector *beta3,
                                   gsl_vector *xbeta3,
                                   double nu3,
                                   gsl_vector *gamma,
                                   gsl_vector *V3,
                                   gsl_vector *lambda3,
                                   gsl_vector *s3,
                                   gsl_vector *survTime1,
                                   gsl_vector *survTime2,
                                   gsl_vector *case11,
                                   gsl_vector *cluster,
                                   gsl_matrix *survCov3,
                                   int K3,
                                   gsl_vector *accept_beta3);


extern void BpeDpCorScr_updateBH1(gsl_vector *lambda1,
                                   gsl_vector *s1,
                                   gsl_vector *xbeta1,
                                   gsl_vector *gamma,
                                   gsl_vector *V1,
                                   gsl_vector *survTime1,
                                   gsl_vector *survEvent1,
                                   gsl_vector *cluster,
                                   gsl_matrix *Sigma_lam1,
                                   gsl_matrix *invSigma_lam1,
                                   gsl_matrix *W1,
                                   gsl_matrix *Q1,
                                   double mu_lam1,
                                   double sigSq_lam1,
                                   int J1);


extern void BpeDpCorScr_updateBH2(gsl_vector *lambda2,
                                   gsl_vector *s2,
                                   gsl_vector *xbeta2,
                                   gsl_vector *gamma,
                                   gsl_vector *V2,
                                   double nu2,
                                   gsl_vector *survTime1,
                                   gsl_vector *survTime2,
                                   gsl_vector *case01,
                                   gsl_vector *cluster,
                                   gsl_matrix *Sigma_lam2,
                                   gsl_matrix *invSigma_lam2,
                                   gsl_matrix *W2,
                                   gsl_matrix *Q2,
                                   double mu_lam2,
                                   double sigSq_lam2,
                                   int J2);


extern void BpeDpCorScr_updateBH3(gsl_vector *lambda3,
                                   gsl_vector *s3,
                                   gsl_vector *xbeta3,
                                   gsl_vector *gamma,
                                   gsl_vector *V3,
                                   double nu3,
                                   gsl_vector *survTime1,
                                   gsl_vector *survTime2,
                                   gsl_vector *case11,
                                   gsl_vector *cluster,
                                   gsl_matrix *Sigma_lam3,
                                   gsl_matrix *invSigma_lam3,
                                   gsl_matrix *W3,
                                   gsl_matrix *Q3,
                                   double mu_lam3,
                                   double sigSq_lam3,
                                   int J3);



extern void BpeDpCorScr_updateSP1(double *mu_lam1,
                                   double *sigSq_lam1,
                                   gsl_vector *lambda1,
                                   gsl_matrix *Sigma_lam1,
                                   gsl_matrix *invSigma_lam1,
                                   double a1,
                                   double b1,
                                   int J1);




extern void BpeDpCorScr_updateSP2(double *mu_lam2,
                                   double *sigSq_lam2,
                                   gsl_vector *lambda2,
                                   gsl_matrix *Sigma_lam2,
                                   gsl_matrix *invSigma_lam2,
                                   double a2,
                                   double b2,
                                   int J2);


extern void BpeDpCorScr_updateSP3(double *mu_lam3,
                                   double *sigSq_lam3,
                                   gsl_vector *lambda3,
                                   gsl_matrix *Sigma_lam3,
                                   gsl_matrix *invSigma_lam3,
                                   double a3,
                                   double b3,
                                   int J3);




extern void BpeDpCorScr_updateFP_Gibbs(gsl_vector *gamma,
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
                                        gsl_vector *V1,
                                        gsl_vector *V2,
                                        gsl_vector *V3,
                                        gsl_vector *survTime1,
                                        gsl_vector *survTime2,
                                        gsl_vector *survEvent1,
                                        gsl_vector *survEvent2,
                                        gsl_vector *cluster);








extern void BpeDpCorScr_updateDP(gsl_vector *gamma,
                     double *theta,
                     double mhProp_theta_var,
                     double psi,
                     double omega,
                     int *accept_theta);







extern void BpeDpCorScr_updateBI1(gsl_vector *s1,
                                   int *J1,
                                   int *accept_BI1,
                                   gsl_vector *survTime1,
                                   gsl_vector *survEvent1,
                                   gsl_vector *gamma,
                                   gsl_vector *xbeta1,
                                   gsl_vector *V1,
                                   gsl_vector *cluster,
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





extern void BpeDpCorScr_updateDI1(gsl_vector *s1,
                            int *J1,
                            int *accept_DI1,
                            gsl_vector *survTime1,
                            gsl_vector *survEvent1,
                            gsl_vector *gamma,
                            gsl_vector *xbeta1,
                            gsl_vector *V1,
                            gsl_vector *cluster,
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


extern void BpeDpCorScr_updateBI2(gsl_vector *s2,
                                   int *J2,
                                   int *accept_BI2,
                                   gsl_vector *survTime1,
                                   gsl_vector *survTime2,
                                   gsl_vector *case01,
                                   gsl_vector *gamma,
                                   gsl_vector *xbeta2,
                                   gsl_vector *V2,
                                   gsl_vector *cluster,
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


extern void BpeDpCorScr_updateDI2(gsl_vector *s2,
                                   int *J2,
                                   int *accept_DI2,
                                   gsl_vector *survTime1,
                                   gsl_vector *survTime2,
                                   gsl_vector *case01,
                                   gsl_vector *gamma,
                                   gsl_vector *xbeta2,
                                   gsl_vector *V2,
                                   gsl_vector *cluster,
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






extern void BpeDpCorScr_updateBI3(gsl_vector *s3,
                                   int *J3,
                                   int *accept_BI3,
                                   gsl_vector *survTime1,
                                   gsl_vector *survTime2,
                                   gsl_vector *case11,
                                   gsl_vector *gamma,
                                   gsl_vector *xbeta3,
                                   gsl_vector *V3,
                                   gsl_vector *cluster,
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



void BpeDpCorScr_updateDI3(gsl_vector *s3,
                            int *J3,
                            int *accept_DI3,
                            gsl_vector *survTime1,
                            gsl_vector *survTime2,
                            gsl_vector *case11,
                            gsl_vector *gamma,
                            gsl_vector *xbeta3,
                            gsl_vector *V3,
                            gsl_vector *cluster,
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





extern void BpeDpCorScr_updateMP(gsl_vector *beta2,
                              gsl_vector *beta3,
                              double *alpha2,
                              double *alpha3,
                              double *kappa2,
                              double *kappa3,
                              double *nu2,
                              double *nu3,
                              gsl_vector *gamma,
                              gsl_vector *V2,
                              gsl_vector *V3,
                              gsl_vector *survTime1,
                              gsl_vector *survTime2,
                              gsl_vector *case01,
                              gsl_vector *case11,
                              gsl_vector *cluster,
                              gsl_matrix *survCov2,
                              gsl_matrix *survCov3,
                              int *accept_nu2,
                              int *accept_nu3);









extern void BpeDpCorScr_updateCP(gsl_vector *beta1,
                                 gsl_vector *beta2,
                                 gsl_vector *beta3,
                                 gsl_vector *lambda1,
                                 gsl_vector *lambda2,
                                 gsl_vector *lambda3,
                                 gsl_vector *s1,
                                 gsl_vector *s2,
                                 gsl_vector *s3,
                                 int K1,
                                 int K2,
                                 int K3,
                                 gsl_vector *gamma,
                                 gsl_vector *V1,
                                 gsl_vector *V2,
                                 gsl_vector *V3,
                                 double nu2,
                                 double nu3,
                                 gsl_vector *survTime1,
                                 gsl_vector *survTime2,
                                 gsl_vector *survEvent1,
                                 gsl_vector *case01,
                                 gsl_vector *case11,
                                 gsl_vector *cluster,
                                 gsl_matrix *survCov1,
                                 gsl_matrix *survCov2,
                                 gsl_matrix *survCov3,
                                 gsl_vector *n_j,
                                 gsl_vector *mu_all,
                                 gsl_matrix *Sigma_all,
                                 gsl_vector *c,
                                 gsl_vector *accept_V,
                                 gsl_vector *mu0,
                                 gsl_matrix *Psi0,
                                 double zeta0,
                                 double rho0,
                                 double tau,
                                 int *nClass_DP,
                                 gsl_rng *rr);















extern void BpeDpCorScr_updatePP(int *n,
                                   double *tau,
                                   double aTau,
                                   double bTau,
                                   int *nClass_DP);




