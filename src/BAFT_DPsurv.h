



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

extern double logistic(double x);

extern double fmixTN(double y,
                     double yL,
                     double yU,
                     int yL_neginf,
                     int yU_posinf,
                     gsl_vector *mu_all,
                     gsl_vector *zeta_all,
                     gsl_vector *rUniq_count,
                     int u,
                     double eta);

extern void BAFT_DPsurv_update_mu_zeta(gsl_vector *yL,
                                       gsl_vector *yU,
                                       gsl_vector *yU_posinf,
                                       gsl_vector *c0,
                                       gsl_vector *c0_neginf,
                                       gsl_matrix *X,
                                       gsl_vector *y,
                                       gsl_vector *beta,
                                       gsl_vector *r,
                                       gsl_vector *mu_all,
                                       gsl_vector *zeta_all,
                                       gsl_vector *rUniq,
                                       gsl_vector *rUniq_count,
                                       double tau,
                                       double a0,
                                       double b0,
                                       double mu0,
                                       double sigSq0,
                                       double beta0_prop_var,
                                       double zeta_prop_var,
                                       int *accept_beta0,
                                       int *accept_sigSq,
                                       int *nClass_DP,
                                       gsl_rng *rr);

extern void BAFT_DPsurv_update_y(gsl_vector *yL,
                                 gsl_vector *yU,
                                 gsl_vector *yL_neginf,
                                 gsl_vector *yU_posinf,
                                 gsl_vector *c0,
                                 gsl_matrix *X,
                                 gsl_vector *y,
                                 gsl_vector *beta,
                                 gsl_vector *r,
                                 gsl_vector *mu_all,
                                 gsl_vector *zeta_all,
                                 gsl_vector *rUniq,
                                 gsl_vector *rUniq_count,
                                 int *nClass_DP,
                                 gsl_rng *rr);

extern void BAFT_DPsurv_update_beta(gsl_vector *yL,
                                    gsl_vector *yU,
                                    gsl_vector *yU_posinf,
                                    gsl_vector *c0,
                                    gsl_vector *c0_neginf,
                                    gsl_matrix *X,
                                    gsl_vector *y,
                                    gsl_vector *beta,
                                    gsl_vector *r,
                                    gsl_vector *mu_all,
                                    gsl_vector *zeta_all,
                                    gsl_vector *rUniq,
                                    gsl_vector *rUniq_count,
                                    int *nClass_DP,
                                    double beta_prop_var,
                                    gsl_vector *accept_beta);


extern void BAFT_DPsurv_update_tau(int *n,
                                   double *tau,
                                   double aTau,
                                   double bTau,
                                   int *nClass_DP);


extern void matrixInv(gsl_matrix *X, gsl_matrix *Xinv);

extern void c_quadform_vMv(gsl_vector *v,
                           gsl_matrix *Minv,
                           double     *value);

extern int c_multinom_sample(gsl_rng *rr,
                             gsl_vector *prob,
                             int length_prob);


extern void c_dmvnorm2(gsl_vector *x,
                       gsl_vector *mu,
                       double     sigma,
                       gsl_matrix *AInv,
                       double     *value);



extern void c_rtnorm(double mean,
                     double sd,
                     double LL,
                     double UL,
                     int LL_neginf,
                     int UL_posinf,
                     double *value);


