//
//  BweibSurv.h
//  
//
//  Created by Kyu Ha Lee on 11/7/13.
//
//


extern void BweibScrSM_logMLH(gsl_vector *beta1,
                              gsl_vector *beta2,
                              gsl_vector *beta3,
                              double alpha1,
                              double alpha2,
                              double alpha3,
                              double kappa1,
                              double kappa2,
                              double kappa3,
                              double theta,
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
                              double *val);



extern void BweibScrSM_logMLH_i(int i,
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
                                double *val);


extern double BweibScrSM_wFunc_old(int subjInx,
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
                                   gsl_matrix *survCov3);




extern void matrixInv(gsl_matrix *X, gsl_matrix *Xinv);

extern void c_colSums(gsl_matrix *X, gsl_vector *v);

extern void c_rowSums(gsl_matrix *X, gsl_vector *v);

extern void c_repVec_Rowmat(gsl_vector *v, gsl_matrix *X);

extern void c_repVec_Colmat(gsl_vector *v, gsl_matrix *X);

extern double BweibScrSM_wFunc(int subjInx,
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
                    gsl_matrix *survCov3);




extern void BweibScrSM_updateRP1(gsl_vector *beta1,
                      double *alpha1,
                      double *kappa1,
                      gsl_vector *gamma,
                      gsl_vector *survTime1,
                      gsl_vector *survEvent1,
                      gsl_matrix *survCov1,
                      gsl_vector *accept_beta1);



extern void BweibScrSM_updateRP2(gsl_vector *beta2,
                      double *alpha2,
                      double *kappa2,
                      gsl_vector *gamma,
                      gsl_vector *survTime1,
                      gsl_vector *case01,
                      gsl_matrix *survCov2,
                      gsl_vector *accept_beta2);


extern void BweibScrSM_updateRP3(gsl_vector *beta3,
                      double *alpha3,
                      double *kappa3,
                      gsl_vector *gamma,
                      gsl_vector *yStar,
                      gsl_vector *case11,
                      gsl_matrix *survCov3,
                      gsl_vector *accept_beta3);



extern void BweibScrSM_updateSH1(gsl_vector *beta1,
                      double *alpha1,
                      double *kappa1,
                      gsl_vector *gamma,
                      gsl_vector *survTime1,
                      gsl_vector *survEvent1,
                      gsl_matrix *survCov1,
                      double c1,
                      double d1);


extern void BweibScrSM_updateSH2(gsl_vector *beta2,
                      double *alpha2,
                      double *kappa2,
                      gsl_vector *gamma,
                      gsl_vector *survTime1,
                      gsl_vector *case01,
                      gsl_matrix *survCov2,
                      double c2,
                      double d2);


extern void BweibScrSM_updateSH3(gsl_vector *beta3,
                      double *alpha3,
                      double *kappa3,
                      gsl_vector *gamma,
                      gsl_vector *yStar,
                      gsl_vector *case11,
                      gsl_matrix *survCov3,
                      double c3,
                      double d3);




extern void BweibScrSM_updateFP(gsl_vector *beta1,
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
                     gsl_vector *survTime1,
                     gsl_vector *yStar,
                     gsl_vector *survEvent1,
                     gsl_vector *survEvent2,
                     gsl_matrix *survCov1,
                     gsl_matrix *survCov2,
                     gsl_matrix *survCov3);







extern void BweibScrSM_updateSC1(gsl_vector *beta1,
                      double *alpha1,
                      double *kappa1,
                      gsl_vector *gamma,
                      gsl_vector *survTime1,
                      gsl_vector *survEvent1,
                      gsl_matrix *survCov1,
                      double mhProp_alpha1_var,
                      double a1,
                      double b1,
                      int *accept_alpha1);




extern void BweibScrSM_updateSC2(gsl_vector *beta2,
                      double *alpha2,
                      double *kappa2,
                      gsl_vector *gamma,
                      gsl_vector *survTime1,
                      gsl_vector *survTime2,
                      gsl_vector *case01,
                      gsl_matrix *survCov2,
                      double mhProp_alpha2_var,
                      double a2,
                      double b2,
                      int *accept_alpha2);





extern void BweibScrSM_updateSC3(gsl_vector *beta3,
                      double *alpha3,
                      double *kappa3,
                      gsl_vector *gamma,
                      gsl_vector *yStar,
                      gsl_vector *case11,
                      gsl_matrix *survCov3,
                      double mhProp_alpha3_var,
                      double a3,
                      double b3,
                      int *accept_alpha3);





extern void BweibScrSM_updateDP(gsl_vector *gamma,
                     double *theta,
                     double mhProp_theta_var,
                     double psi,
                     double omega,
                     int *accept_theta);









