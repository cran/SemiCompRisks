


extern int c_multinom_sample(gsl_rng *rr,
                             gsl_vector *prob,
                             int length_prob);


extern double Qfunc_univ(double V,
                    double mu0,
                    double zeta0,
                    double a0,
                    double b0);



extern void matrixInv(gsl_matrix *X, gsl_matrix *Xinv);

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




extern void BweibDpCorSurv_updateRP(gsl_vector *beta,
                                  double *alpha,
                                  double *kappa,
                                  gsl_vector *V,
                                  gsl_vector *survTime,
                                  gsl_vector *survEvent,
                                  gsl_vector *cluster,
                                  gsl_matrix *survCov,
                                  gsl_vector *accept_beta);

extern void BweibDpCorSurv_updateSH_rw2(gsl_vector *beta,
                                   double *alpha,
                                   double *kappa,
                                   gsl_vector *V,
                                   gsl_vector *survTime,
                                   gsl_vector *survEvent,
                                   gsl_vector *cluster,
                                   gsl_matrix *survCov,
                                   double mhProp_alpha_var,
                                   double a,
                                   double b,
                                   int *accept_alpha);





extern void BweibDpCorSurv_updateSC(gsl_vector *beta,
                               double *alpha,
                               double *kappa,
                               gsl_vector *V,
                               gsl_vector *survTime,
                               gsl_vector *survEvent,
                               gsl_vector *cluster,
                               gsl_matrix *survCov,
                               double c,
                               double d);





extern void BweibDpCorSurv_updateCP(gsl_vector *beta,
                                  double alpha,
                                    double kappa,
                                  gsl_vector *V,
                                  gsl_vector *survTime,
                                  gsl_vector *survEvent,
                                  gsl_vector *cluster,
                                  gsl_matrix *survCov,
                                  gsl_vector *n_j,
                                  gsl_vector *mu_all,
                                  gsl_vector *zeta_all,
                                  gsl_vector *c,
                                  gsl_vector *accept_V,
                                    double mhProp_V_var,
                                  double mu0,
                                  double zeta0,
                                  double a0,
                                  double b0,
                                  double tau,
                                  int *nClass_DP,
                                  gsl_rng *rr);


extern void BweibDpCorSurv_updatePP(int *n,
                                  double *tau,
                                  double aTau,
                                  double bTau,
                                  int *nClass_DP);






