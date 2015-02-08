


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




extern void BweibCorSurv_updateRP(gsl_vector *beta,
                                  double *alpha,
                                  double *kappa,
                                  gsl_vector *V,
                                  gsl_vector *survTime,
                                  gsl_vector *survEvent,
                                  gsl_vector *cluster,
                                  gsl_matrix *survCov,
                                  gsl_vector *accept_beta);

extern void BweibSurv_updateSH_rw2(gsl_vector *beta,
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





extern void BweibSurv_updateSC(gsl_vector *beta,
                               double *alpha,
                               double *kappa,
                               gsl_vector *V,
                               gsl_vector *survTime,
                               gsl_vector *survEvent,
                               gsl_vector *cluster,
                               gsl_matrix *survCov,
                               double c,
                               double d);



extern void BweibSurv_updateCP(gsl_vector *beta,
                        double alpha,
                        double kappa,
                        gsl_vector *V,
                        double zeta,
                        gsl_vector *survTime,
                        gsl_vector *survEvent,
                        gsl_vector *cluster,
                        gsl_matrix *survCov,
                        gsl_vector *n_j,
                        gsl_vector *accept_V,
                        double mhProp_V_var);



extern void BweibSurv_updateVP(gsl_vector *V,
                               double *zeta,
                               double rho1,
                               double rho2);







