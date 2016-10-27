

extern void set_r3_mu3_zeta3(gsl_vector *r3,
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
                             gsl_rng *rr);

extern double c_min(double value1,
                    double value2);

extern double c_max(double value1,
                    double value2);

extern double logistic(double x);

extern void c_uniq(gsl_vector *r,
                   gsl_vector *rUniq,
                   gsl_vector *rUniq_count,
                   gsl_vector *mu_all,
                   gsl_vector *zeta_all,
                   int *u);

extern void c_uniq1(gsl_vector *r,
                    gsl_vector *rUniq,
                    gsl_vector *rUniq_count,
                    int *u);


extern void c_uniq_h3(gsl_vector *r,
                      gsl_vector *rUniq,
                      gsl_vector *rUniq_count,
                      gsl_vector *mu_all,
                      gsl_vector *zeta_all,
                      gsl_vector *mu3_vec,
                      gsl_vector *zeta3_vec,
                      gsl_vector *y1_NA,
                      int *u);

extern void c_uniq1_h3(gsl_vector *r,
                       gsl_vector *rUniq,
                       gsl_vector *rUniq_count,
                       gsl_vector *y1_NA,
                       int *u);


extern void BAFT_DPscr_update_mu_zeta3(gsl_vector *c0,
                                       gsl_vector *c0_neginf,
                                       gsl_matrix *X3,
                                       gsl_vector *y1,
                                       gsl_vector *y2,
                                       gsl_vector *y1_NA,
                                       gsl_vector *beta3,
                                       gsl_vector *gamma,
                                       gsl_vector *r3,
                                       gsl_vector *mu3_all,
                                       gsl_vector *mu3_vec,
                                       gsl_vector *zeta3_all,
                                       gsl_vector *zeta3_vec,
                                       gsl_vector *r3Uniq,
                                       gsl_vector *r3Uniq_count,
                                       double tau3,
                                       double a03,
                                       double b03,
                                       double mu03,
                                       double sigSq03,
                                       double beta03_prop_var,
                                       double zeta3_prop_var,
                                       int *accept_beta03,
                                       int *accept_zeta3,
                                       int *nClass_DP3,
                                       gsl_rng *rr);



extern void BAFT_DPscr_update_mu_zeta1(gsl_vector *c0,
                                       gsl_vector *c0_neginf,
                                       gsl_matrix *X1,
                                       gsl_vector *y1,
                                       gsl_vector *y2,
                                       gsl_vector *y1_NA,
                                       gsl_vector *beta1,
                                       gsl_vector *gamma,
                                       gsl_vector *r1,
                                       gsl_vector *mu1_all,
                                       gsl_vector *zeta1_all,
                                       gsl_vector *r1Uniq,
                                       gsl_vector *r1Uniq_count,
                                       double tau1,
                                       double a01,
                                       double b01,
                                       double mu01,
                                       double sigSq01,
                                       double beta01_prop_var,
                                       double zeta1_prop_var,
                                       int *accept_beta01,
                                       int *accept_zeta1,
                                       int *nClass_DP1,
                                       gsl_rng *rr);

extern void BAFT_DPscr_update_mu_zeta2(gsl_vector *c0,
                                       gsl_vector *c0_neginf,
                                       gsl_matrix *X2,
                                       gsl_vector *y1,
                                       gsl_vector *y2,
                                       gsl_vector *y1_NA,
                                       gsl_vector *beta2,
                                       gsl_vector *gamma,
                                       gsl_vector *r2,
                                       gsl_vector *mu2_all,
                                       gsl_vector *zeta2_all,
                                       gsl_vector *r2Uniq,
                                       gsl_vector *r2Uniq_count,
                                       double tau2,
                                       double a02,
                                       double b02,
                                       double mu02,
                                       double sigSq02,
                                       double beta02_prop_var,
                                       double zeta2_prop_var,
                                       int *accept_beta02,
                                       int *accept_zeta2,
                                       int *nClass_DP2,
                                       gsl_rng *rr);


extern void BAFT_DPscr_update_y(gsl_vector *y1L,
                                gsl_vector *y1U,
                                gsl_vector *y2L,
                                gsl_vector *y2U,
                                gsl_vector *y1L_neginf,
                                gsl_vector *y2L_neginf,
                                gsl_vector *y1U_posinf,
                                gsl_vector *y2U_posinf,
                                gsl_vector *c0_neginf,
                                gsl_vector *y1_NA,
                                gsl_matrix *X1,
                                gsl_matrix *X2,
                                gsl_matrix *X3,
                                gsl_vector *y1,
                                gsl_vector *y2,
                                gsl_vector *beta1,
                                gsl_vector *beta2,
                                gsl_vector *beta3,
                                gsl_vector *r1,
                                gsl_vector *r2,
                                gsl_vector *r3,
                                gsl_vector *mu1_all,
                                gsl_vector *mu2_all,
                                gsl_vector *mu3_all,
                                gsl_vector *mu3_vec,
                                gsl_vector *zeta1_all,
                                gsl_vector *zeta2_all,
                                gsl_vector *zeta3_all,
                                gsl_vector *zeta3_vec,
                                gsl_vector *r1Uniq,
                                gsl_vector *r2Uniq,
                                gsl_vector *r3Uniq,
                                gsl_vector *r1Uniq_count,
                                gsl_vector *r2Uniq_count,
                                gsl_vector *r3Uniq_count,
                                gsl_vector *gamma,
                                int *nClass_DP1,
                                int *nClass_DP2,
                                int *nClass_DP3,
                                double mu03,
                                double sigSq03,
                                double a03,
                                double b03,
                                double tau3,
                                gsl_rng *rr);




extern void BAFT_DPscr_update_beta1(gsl_vector *y1_NA,
                                    gsl_vector *c0,
                                    gsl_vector *c0_neginf,
                                    gsl_matrix *X1,
                                    gsl_vector *y1,
                                    gsl_vector *y2,
                                    gsl_vector *beta1,
                                    gsl_vector *gamma,
                                    gsl_vector *r1,
                                    gsl_vector *mu1_all,
                                    gsl_vector *zeta1_all,
                                    gsl_vector *r1Uniq,
                                    gsl_vector *r1Uniq_count,
                                    int *nClass_DP1,
                                    double beta1_prop_var,
                                    gsl_vector *accept_beta1);


extern void BAFT_DPscr_update_beta2(gsl_vector *y1_NA,
                                    gsl_vector *c0,
                                    gsl_vector *c0_neginf,
                                    gsl_matrix *X2,
                                    gsl_vector *y1,
                                    gsl_vector *y2,
                                    gsl_vector *beta2,
                                    gsl_vector *gamma,
                                    gsl_vector *r2,
                                    gsl_vector *mu2_all,
                                    gsl_vector *zeta2_all,
                                    gsl_vector *r2Uniq,
                                    gsl_vector *r2Uniq_count,
                                    int *nClass_DP2,
                                    double beta2_prop_var,
                                    gsl_vector *accept_beta2);

extern void BAFT_DPscr_update_beta3(gsl_vector *y1_NA,
                                    gsl_vector *c0,
                                    gsl_vector *c0_neginf,
                                    gsl_matrix *X3,
                                    gsl_vector *y1,
                                    gsl_vector *y2,
                                    gsl_vector *beta3,
                                    gsl_vector *gamma,
                                    gsl_vector *r3,
                                    gsl_vector *mu3_all,
                                    gsl_vector *zeta3_all,
                                    gsl_vector *r3Uniq,
                                    gsl_vector *r3Uniq_count,
                                    int *nClass_DP3,
                                    double beta3_prop_var,
                                    gsl_vector *accept_beta3);



extern void BAFT_DPscr_update_gamma(gsl_vector *y1_NA,
                                    gsl_vector *c0,
                                    gsl_vector *c0_neginf,
                                    gsl_matrix *X1,
                                    gsl_matrix *X2,
                                    gsl_matrix *X3,
                                    gsl_vector *y1,
                                    gsl_vector *y2,
                                    gsl_vector *beta1,
                                    gsl_vector *beta2,
                                    gsl_vector *beta3,
                                    gsl_vector *gamma,
                                    gsl_vector *r1,
                                    gsl_vector *r2,
                                    gsl_vector *r3,
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
                                    int *nClass_DP1,
                                    int *nClass_DP2,
                                    int *nClass_DP3,
                                    double theta,
                                    double gamma_prop_var,
                                    gsl_vector *accept_gamma);

extern void BAFT_DPscr_update_theta(gsl_vector *gamma,
                                    double *theta,
                                    double a_theta,
                                    double b_theta);

extern void BAFT_DPscr_update_tau(int *n,
                                  double *tau,
                                  double aTau,
                                  double bTau,
                                  int *nClass_DP);



extern void matrixInv(gsl_matrix *X, gsl_matrix *Xinv);

extern void c_quadform_vMv(gsl_vector *v,
                           gsl_matrix *Minv,
                           double     *value);

extern void c_dmvnorm2(gsl_vector *x,
                       gsl_vector *mu,
                       double     sigma,
                       gsl_matrix *AInv,
                       double     *value);

extern void c_riwishart_general(int v,
                                gsl_matrix *X_ori,
                                gsl_matrix *sample);

extern void c_rtnorm(double mean,
                     double sd,
                     double LL,
                     double UL,
                     int LL_neginf,
                     int UL_posinf,
                     double *value);



extern void log_Jpdf_Upper_BAFT_DP(int i,
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
                                   double *value);




extern void log_Jpdf_Lower_BAFT_DP(int i,
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
                                   double *value);




extern void c_rigamma(double *temp,
                      double alpha,
                      double beta);

extern double Qfunc_BAFT_DP(double V,
                            double mu0,
                            double zeta0,
                            double a0,
                            double b0);


extern int c_multinom_sample(gsl_rng *rr,
                             gsl_vector *prob,
                             int length_prob);


extern int MFunc_BAFT_DP(double yL,
                         double yU,
                         int yL_neginf,
                         int yU_posinf,
                         gsl_vector *mu_all,
                         gsl_vector *zeta_all,
                         gsl_vector *rUniq_count,
                         int u,
                         double eta,
                         gsl_rng *rr);

extern void log_fg_BAFT_DP(double y,
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
                           double *logSg);





