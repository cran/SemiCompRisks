


extern double logistic(double x);

extern void BAFT_LNsurv_update_y(gsl_vector *yL,
                                 gsl_vector *yU,
                                 gsl_vector *yU_posinf,
                                 gsl_vector *c0,
                                 gsl_matrix *X,
                                 gsl_vector *y,
                                 gsl_vector *beta,
                                 double beta0,
                                 double sigSq);

extern void BAFT_LNsurv_update_beta(gsl_vector *yL,
                                    gsl_vector *yU,
                                    gsl_vector *yU_posinf,
                                    gsl_vector *c0,
                                    gsl_vector *c0_neginf,
                                    gsl_matrix *X,
                                    gsl_vector *y,
                                    gsl_vector *beta,
                                    double beta0,
                                    double sigSq,
                                    double beta_prop_var,
                                    gsl_vector *accept_beta);

extern void BAFT_LNsurv_update_beta0(gsl_vector *yL,
                                     gsl_vector *yU,
                                     gsl_vector *yU_posinf,
                                     gsl_vector *c0,
                                     gsl_vector *c0_neginf,
                                     gsl_matrix *X,
                                     gsl_vector *y,
                                     gsl_vector *beta,
                                     double *beta0,
                                     double sigSq,
                                     double beta0_prop_var,
                                     int *accept_beta0);


extern void BAFT_LNsurv_update_sigSq(gsl_vector *yL,
                                     gsl_vector *yU,
                                     gsl_vector *yU_posinf,
                                     gsl_vector *c0,
                                     gsl_vector *c0_neginf,
                                     gsl_matrix *X,
                                     gsl_vector *y,
                                     gsl_vector *beta,
                                     double beta0,
                                     double *sigSq,
                                     double a_sigSq,
                                     double b_sigSq,
                                     double sigSq_prop_var,
                                     int *accept_sigSq);

extern void matrixInv(gsl_matrix *X, gsl_matrix *Xinv);

extern void c_quadform_vMv(gsl_vector *v,
                           gsl_matrix *Minv,
                           double     *value);

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


extern int c_multinom_sample(gsl_rng *rr,
                             gsl_vector *prob,
                             int length_prob);


