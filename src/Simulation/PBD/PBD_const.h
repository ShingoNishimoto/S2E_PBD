#ifndef PBD_CONST_H_
#define PBD_CONST_H_

// constant parameters
extern const double L1_frequency; //[MHz]
extern const double L2_frequency;
extern const double L1_lambda; //[m];
extern const double L2_lambda; //[m]
extern const double mu_e; //GM_E m^3/s^2
extern const double J2_const; //無次元 重力J2項
extern const double Earth_Radius; //m

// A-priori noise
extern const double sigma_r_ini;    //[m]
extern const double sigma_v_ini;    //[m/s]
extern const double sigma_acc_r_ini;    //[nm/s^2]
extern const double sigma_acc_t_ini;    //[nm/s^2]
extern const double sigma_acc_n_ini;    //[nm/s^2]
extern const double sigma_cdt_ini;    //[m]
extern const double sigma_N_ini;    //[cycle]

// process noise
extern const double sigma_r_process;    //[m]
extern const double sigma_v_process;    //[m/s]
extern const double sigma_acc_r_process;    //[nm/s^2]
extern const double sigma_acc_t_process;    //[nm/s^2]
extern const double sigma_acc_n_process;    //[nm/s^2]
extern const double sigma_cdt_process;    //[m]
extern const double sigma_N_process;    //[cycle]

// gauss markov parameters
extern const double tau_a; //[s]
extern const double tau_cdt; //[s]

extern const double alpha; // forgetting factor
#endif
