extern void BweibMvnCorScrSM_logMLH(gsl_vector *beta1,
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
                                    double *val);


extern void BweibMvnCorScrSM_logMLH_i(int i,
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
                                      double *val);



extern void BweibMvnCorScrSM_logLH(gsl_vector *beta1,
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
                                   double *val);




extern void BweibMvnCorScrSM_logLH_i(int i,
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
                                     double *val);

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




extern void c_riwishart(int v,
                        gsl_matrix *X_ori,
                        gsl_matrix *sample);



extern double BweibMvnCorScrSM_wFunc(int subjInx,
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
                             gsl_matrix *survCov3);



extern double BweibMvnCorScrSM_wFunc_old(int subjInx,
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
                                 gsl_matrix *survCov3);







extern void BweibMvnCorScrSM_updateRP1(gsl_vector *beta1,
                                  double *alpha1,
                                  double *kappa1,
                                  gsl_vector *gamma,
                                  gsl_vector *V1,
                                  gsl_vector *survTime1,
                                  gsl_vector *survEvent1,
                                  gsl_vector *cluster,
                                  gsl_matrix *survCov1,
                                  gsl_vector *accept_beta1);






extern void BweibMvnCorScrSM_updateRP2(gsl_vector *beta2,
                               double *alpha2,
                               double *kappa2,
                               double *nu2,
                               gsl_vector *gamma,
                               gsl_vector *V2,
                               gsl_vector *survTime1,
                               gsl_vector *case01,
                               gsl_vector *cluster,
                               gsl_matrix *survCov2,
                               gsl_vector *accept_beta2);


extern void BweibMvnCorScrSM_updateRP3(gsl_vector *beta3,
                               double *alpha3,
                               double *kappa3,
                               double *nu3,
                               gsl_vector *gamma,
                               gsl_vector *V3,
                               gsl_vector *yStar,
                               gsl_vector *case11,
                               gsl_vector *cluster,
                               gsl_matrix *survCov3,
                               gsl_vector *accept_beta3);



extern void BweibMvnCorScrSM_updateSH1(gsl_vector *beta1,
                               double *alpha1,
                               double *kappa1,
                               gsl_vector *gamma,
                               gsl_vector *V1,
                               gsl_vector *survTime1,
                               gsl_vector *survEvent1,
                               gsl_vector *cluster,
                               gsl_matrix *survCov1,
                               double c1,
                               double d1);


extern void BweibMvnCorScrSM_updateSH2(gsl_vector *beta2,
                               double *alpha2,
                               double *kappa2,
                               double *nu2,
                               gsl_vector *gamma,
                               gsl_vector *V2,
                               gsl_vector *survTime1,
                               gsl_vector *case01,
                               gsl_vector *cluster,
                               gsl_matrix *survCov2,
                               double c2,
                               double d2);


extern void BweibMvnCorScrSM_updateSH3(gsl_vector *beta3,
                               double *alpha3,
                               double *kappa3,
                               double *nu3,
                               gsl_vector *gamma,
                               gsl_vector *V3,
                               gsl_vector *yStar,
                               gsl_vector *case11,
                               gsl_vector *cluster,
                               gsl_matrix *survCov3,
                               double c3,
                               double d3);


extern void BweibMvnCorScrSM_updateFP(gsl_vector *beta1,
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
                              double theta,
                              gsl_vector *gamma,
                              gsl_vector *V1,
                              gsl_vector *V2,
                              gsl_vector *V3,
                              gsl_vector *survTime1,
                              gsl_vector *survTime2,
                              gsl_vector *survEvent1,
                              gsl_vector *survEvent2,
                              gsl_vector *cluster,                              
                              gsl_matrix *survCov1,
                              gsl_matrix *survCov2,
                              gsl_matrix *survCov3,
                              gsl_vector *accept_gamma,
                              double mhProp_gamma_var,
                              int numGamUpdate,
                              gsl_vector *mhGam_chk,
                              int *ChgProp);


extern void BweibMvnCorScrSM_updateFP_Gibbs(gsl_vector *beta1,
                                    gsl_vector *beta2,
                                    gsl_vector *beta3,
                                    double alpha1,
                                    double alpha2,
                                    double alpha3,
                                    double kappa1,
                                    double kappa2,
                                    double kappa3,
                                    double theta,
                                    gsl_vector *gamma,
                                    gsl_vector *V1,
                                    gsl_vector *V2,
                                    gsl_vector *V3,
                                    gsl_vector *survTime1,
                                    gsl_vector *yStar,
                                    gsl_vector *survEvent1,
                                    gsl_vector *survEvent2,
                                    gsl_vector *cluster,
                                    gsl_matrix *survCov1,
                                    gsl_matrix *survCov2,
                                    gsl_matrix *survCov3);







extern void BweibMvnCorScrSM_updateSC1(gsl_vector *beta1,
                               double *alpha1,
                               double *kappa1,
                               gsl_vector *gamma,
                               gsl_vector *V1,
                               gsl_vector *survTime1,
                               gsl_vector *survEvent1,
                               gsl_vector *cluster,
                               gsl_matrix *survCov1,
                               double mhProp_alpha1_var,
                               double a1,
                               double b1,
                               int *accept_alpha1);




extern void BweibMvnCorScrSM_updateSC2(gsl_vector *beta2,
                               double *alpha2,
                               double *kappa2,
                               double *nu2,
                               gsl_vector *gamma,
                               gsl_vector *V2,
                               gsl_vector *survTime1,
                               gsl_vector *survTime2,
                               gsl_vector *case01,
                               gsl_vector *cluster,
                               gsl_matrix *survCov2,
                               double mhProp_alpha2_var,
                               double a2,
                               double b2,
                               int *accept_alpha2);





extern void BweibMvnCorScrSM_updateSC3(gsl_vector *beta3,
                               double *alpha3,
                               double *kappa3,
                               double *nu3,
                               gsl_vector *gamma,
                               gsl_vector *V3,
                               gsl_vector *survTime1,
                               gsl_vector *survTime2,
                               gsl_vector *case11,
                               gsl_vector *cluster,
                               gsl_matrix *survCov3,
                               double mhProp_alpha3_var,
                               double a3,
                               double b3,
                               int *accept_alpha3);





extern void BweibMvnCorScrSM_updateDP(gsl_vector *gamma,
                     double *theta,
                     double mhProp_theta_var,
                     double psi,
                     double omega,
                     int *accept_theta);




extern void BweibMvnCorScrSM_updateMP(gsl_vector *beta2,
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




extern void BweibMvnCorScrSM_updateCP1(gsl_vector *beta1,
                               double alpha1,
                               double kappa1,
                               gsl_vector *gamma,
                               gsl_vector *V1,
                               gsl_vector *V2,
                               gsl_vector *V3,
                               gsl_matrix *Sigma_V,
                               gsl_vector *survTime1,
                               gsl_vector *survEvent1,
                               gsl_vector *cluster,
                               gsl_matrix *survCov1,
                               gsl_vector *n_j,
                               gsl_vector *accept_V1,
                               double mhProp_V1_var);




extern void BweibMvnCorScrSM_updateCP2(gsl_vector *beta2,
                               double alpha2,
                               double kappa2,
                               double nu2,
                               gsl_vector *gamma,
                               gsl_vector *V1,
                               gsl_vector *V2,
                               gsl_vector *V3,
                               gsl_matrix *Sigma_V,
                               gsl_vector *survTime1,
                               gsl_vector *case01,
                               gsl_vector *cluster,
                               gsl_matrix *survCov2,
                               gsl_vector *n_j,
                               gsl_vector *accept_V2,
                               double mhProp_V2_var);



extern void BweibMvnCorScrSM_updateCP3(gsl_vector *beta3,
                               double alpha3,
                               double kappa3,
                               double nu3,
                               gsl_vector *gamma,
                               gsl_vector *V1,
                               gsl_vector *V2,
                               gsl_vector *V3,
                               gsl_matrix *Sigma_V,
                               gsl_vector *survTime1,
                               gsl_vector *survTime2,
                               gsl_vector *case11,
                               gsl_vector *cluster,
                               gsl_matrix *survCov3,
                               gsl_vector *n_j,
                               gsl_vector *accept_V3,
                               double mhProp_V3_var);






extern void BweibMvnCorScrSM_updateSC1_amc(gsl_vector *beta1,
                            double *alpha1,
                            double *kappa1,
                            gsl_vector *gamma,
                            gsl_vector *V1,
                            gsl_vector *survTime1,
                            gsl_vector *survEvent1,
                            gsl_vector *cluster,
                            gsl_matrix *survCov1,
                            double mhProp_alpha1_var,
                            double a1,
                            double b1,
                            int *accept_alpha1);







extern void BweibMvnCorScrSM_updateSC1_amc(gsl_vector *beta1,
                            double *alpha1,
                            double *kappa1,
                            gsl_vector *gamma,
                            gsl_vector *V1,
                            gsl_vector *survTime1,
                            gsl_vector *survEvent1,
                            gsl_vector *cluster,
                            gsl_matrix *survCov1,
                            double mhProp_alpha1_var,
                            double a1,
                            double b1,
                            int *accept_alpha1);


extern void BweibMvnCorScrSM_updateSC1_rw2(gsl_vector *beta1,
                            double *alpha1,
                            double *kappa1,
                            gsl_vector *gamma,
                            gsl_vector *V1,
                            gsl_vector *survTime1,
                            gsl_vector *survEvent1,
                            gsl_vector *cluster,
                            gsl_matrix *survCov1,
                            double mhProp_alpha1_var,
                            double a1,
                            double b1,
                            int *accept_alpha1);


extern void BweibMvnCorScrSM_updateSC2_rw2(gsl_vector *beta2,
                            double *alpha2,
                            double *kappa2,
                            double *nu2,
                            gsl_vector *gamma,
                            gsl_vector *V2,
                            gsl_vector *survTime1,
                            gsl_vector *survTime2,
                            gsl_vector *case01,
                            gsl_vector *cluster,
                            gsl_matrix *survCov2,
                            double mhProp_alpha2_var,
                            double a2,
                            double b2,
                            int *accept_alpha2);


extern void BweibMvnCorScrSM_updateSC3_rw2(gsl_vector *beta3,
                            double *alpha3,
                            double *kappa3,
                            double *nu3,
                            gsl_vector *gamma,
                            gsl_vector *V3,
                            gsl_vector *yStar,
                            gsl_vector *case11,
                            gsl_vector *cluster,
                            gsl_matrix *survCov3,
                            double mhProp_alpha3_var,
                            double a3,
                            double b3,
                            int *accept_alpha3);




extern void BweibMvnCorScrSM_updateCP1_amcmc(gsl_vector *beta1,
                                           double alpha1,
                                           double kappa1,
                                           gsl_vector *gamma,
                                           gsl_vector *V1,
                                           gsl_vector *V2,
                                           gsl_vector *V3,
                                           gsl_matrix *Sigma_V,
                                           gsl_vector *survTime1,
                                           gsl_vector *survEvent1,
                                           gsl_vector *cluster,
                                           gsl_matrix *survCov1,
                                           gsl_vector *n_j,
                                           gsl_vector *accept_V1,
                                           double mhProp_V1_var);


extern void BweibMvnCorScrSM_updateCP2_amcmc(gsl_vector *beta2,
                              double alpha2,
                              double kappa2,
                              double nu2,
                              gsl_vector *gamma,
                              gsl_vector *V1,
                              gsl_vector *V2,
                              gsl_vector *V3,
                              gsl_matrix *Sigma_V,
                              gsl_vector *survTime1,
                              gsl_vector *case01,
                              gsl_vector *cluster,
                              gsl_matrix *survCov2,
                              gsl_vector *n_j,
                              gsl_vector *accept_V2,
                                     double mhProp_V2_var);


extern void BweibMvnCorScrSM_updateCP3_amcmc(gsl_vector *beta3,
                              double alpha3,
                              double kappa3,
                              double nu3,
                              gsl_vector *gamma,
                              gsl_vector *V1,
                              gsl_vector *V2,
                              gsl_vector *V3,
                              gsl_matrix *Sigma_V,
                              gsl_vector *yStar,
                              gsl_vector *case11,
                              gsl_vector *cluster,
                              gsl_matrix *survCov3,
                              gsl_vector *n_j,
                              gsl_vector *accept_V3,
                                     double mhProp_V3_var);


extern void BweibMvnCorScrSM_updateVP(gsl_vector *V1,
                                    gsl_vector *V2,
                                    gsl_vector *V3,
                                    gsl_matrix *Sigma_V,
                                    double rho_v,
                                    gsl_matrix *Psi_v);



