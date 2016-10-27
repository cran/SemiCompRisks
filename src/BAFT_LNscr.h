
extern int c_multinom_sample(gsl_rng *rr,
                             gsl_vector *prob,
                             int length_prob);


extern double logistic(double x);



extern void BAFT_LNscr_update_y(gsl_vector *y1L,
                                gsl_vector *y1U,
                                gsl_vector *y2L,
                                gsl_vector *y2U,
                                gsl_vector *y1L_neginf,
                                gsl_vector *y2L_neginf,
                                gsl_vector *y1U_posinf,
                                gsl_vector *y2U_posinf,
                                gsl_vector *y1_NA,
                                gsl_matrix *X1,
                                gsl_matrix *X2,
                                gsl_matrix *X3,
                                gsl_vector *y1,
                                gsl_vector *y2,
                                gsl_vector *beta1,
                                gsl_vector *beta2,
                                gsl_vector *beta3,
                                gsl_vector *gamma,
                                double beta01,
                                double beta02,
                                double beta03,
                                double sigSq1,
                                double sigSq2,
                                double sigSq3);

extern void BAFT_LNscr_update_beta1(gsl_vector *y1_NA,
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
                                    double beta01,
                                    double beta02,
                                    double beta03,
                                    double sigSq1,
                                    double sigSq2,
                                    double sigSq3,
                                    double beta1_prop_var,
                                    gsl_vector *accept_beta1);

extern void BAFT_LNscr_update_beta2(gsl_vector *y1_NA,
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
                                    double beta01,
                                    double beta02,
                                    double beta03,
                                    double sigSq1,
                                    double sigSq2,
                                    double sigSq3,
                                    double beta2_prop_var,
                                    gsl_vector *accept_beta2);

extern void BAFT_LNscr_update_beta3(gsl_vector *y1_NA,
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
                                    double beta01,
                                    double beta02,
                                    double beta03,
                                    double sigSq1,
                                    double sigSq2,
                                    double sigSq3,
                                    double beta3_prop_var,
                                    gsl_vector *accept_beta3);

extern void BAFT_LNscr_update_beta01(gsl_vector *y1_NA,
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
                                     double *beta01,
                                     double beta02,
                                     double beta03,
                                     double sigSq1,
                                     double sigSq2,
                                     double sigSq3,
                                     double beta01_prop_var,
                                     int *accept_beta01);

extern void BAFT_LNscr_update_beta02(gsl_vector *y1_NA,
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
                                     double beta01,
                                     double *beta02,
                                     double beta03,
                                     double sigSq1,
                                     double sigSq2,
                                     double sigSq3,
                                     double beta02_prop_var,
                                     int *accept_beta02);

extern void BAFT_LNscr_update_beta03(gsl_vector *y1_NA,
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
                                     double beta01,
                                     double beta02,
                                     double *beta03,
                                     double sigSq1,
                                     double sigSq2,
                                     double sigSq3,
                                     double beta03_prop_var,
                                     int *accept_beta03);


extern void BAFT_LNscr_update_sigSq1(gsl_vector *y1_NA,
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
                                     double beta01,
                                     double beta02,
                                     double beta03,
                                     double *sigSq1,
                                     double sigSq2,
                                     double sigSq3,
                                     double a_sigSq1,
                                     double b_sigSq1,
                                     double sigSq1_prop_var,
                                     int *accept_sigSq1);

extern void BAFT_LNscr_update_sigSq2(gsl_vector *y1_NA,
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
                                     double beta01,
                                     double beta02,
                                     double beta03,
                                     double sigSq1,
                                     double *sigSq2,
                                     double sigSq3,
                                     double a_sigSq2,
                                     double b_sigSq2,
                                     double sigSq2_prop_var,
                                     int *accept_sigSq2);


extern void BAFT_LNscr_update_sigSq3(gsl_vector *y1_NA,
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
                                     double beta01,
                                     double beta02,
                                     double beta03,
                                     double sigSq1,
                                     double sigSq2,
                                     double *sigSq3,
                                     double a_sigSq3,
                                     double b_sigSq3,
                                     double sigSq3_prop_var,
                                     int *accept_sigSq3);


extern void BAFT_LNscr_update_gamma(gsl_vector *y1_NA,
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
                                    double beta01,
                                    double beta02,
                                    double beta03,
                                    double sigSq1,
                                    double sigSq2,
                                    double sigSq3,
                                    double theta,
                                    double gamma_prop_var,
                                    gsl_vector *accept_gamma);

extern void BAFT_LNscr_update_theta(gsl_vector *gamma,
                                    double *theta,
                                    double a_theta,
                                    double b_theta);


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


extern void log_Jpdf_Upper_BAFT_LN(int i,
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
                                   double *value);



extern void log_Jpdf_Lower_BAFT_LN(int i,
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
                                   double *value);




extern void c_rigamma(double *temp,
                      double alpha,
                      double beta);












