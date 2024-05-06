#include"1PIECE_main.cpp"
#include<iostream>
#include<fstream>
#include<stdio.h>
#include<string>

// #include <charconv>

int main(){

float temps;
clock_t t1,t2;


t1 = clock();
std::tuple<double,double,double,double,double,double,double,double,double,double,double> rheology;
std::tuple<double,double,double,double,double,double,double,double> melt;
std::tuple<bool,double,double,double,double,double,double,double,double> melt_param;
std::tuple<double,double,double,double,double,double,double,double,double,double,double,double,double> thermo;
Eigen::MatrixXd RAD = Eigen::MatrixXd::Zero(4,4);

double Tm0,Ts, dTc0, Dl_init, Dc_init , dDl_init, dDc_init, LMBD_cr_Ni, LMBD_cr_Si, t0, tstop, Rp, Rc, f ,Masse_Mars, gl,gc,dt_min,dt_max,dt_init,DTMAX,dr_min,dr_max,t_acc=0.0,t_acc_m=0.0,fmag=0.0,fbase=0.0, HPE_factor = 1;
double X1,X2,Sat_expo,CH2O_p, KH2Ocr, gamma_cr,phi_min_cr;
int N_soliq, N_melt,Nc;
bool unlinear_phi,LMBDdt,RK4,Pressure,Steady, ecrit_time_b, ecrit_tech_b, ecrit_profil_b,URG_STOP,adv_lith_bool,adv_interf_bool;

auto solidus = std::make_tuple(1E2,1E2,1E2,1E2,1E2,1E2);
auto liquidus = std::make_tuple(1E2,1E2,1E2,1E2,1E2,1E2);



// /Users/valentinbonnetgibet/Documents/These/Code Mars/Mars Cpp Dev/
std::string physical_param = "physical_parameters.txt";
std::string numerical_param = "numerical_parameters.txt";
std::string adress_rad = "Radioelement_Wanke.txt";
std::string adress_solidus = "Solidus.txt";


// std::string fichier_time = "/Users/valentinbonnetgibet/Documents/These/Code Mars/Mars Cpp Dev/time.txt";
// std::string fichier_profil = "/Users/valentinbonnetgibet/Documents/These/Code Mars/Mars Cpp Dev/profil.txt";

Lecture_param(physical_param,rheology,melt,melt_param,thermo,Tm0,Ts,dTc0,Dl_init,Dc_init,dDl_init, dDc_init, LMBD_cr_Ni,LMBD_cr_Si,t0,tstop,Rp,Rc,f,Masse_Mars,gl,gc,fmag,fbase,HPE_factor,X1,X2,Sat_expo,CH2O_p,KH2Ocr,gamma_cr,phi_min_cr);
lecture_num(numerical_param,dt_min,dt_max, dt_init,DTMAX,dr_min,dr_max,Nc,N_soliq,N_melt,t_acc,t_acc_m,unlinear_phi,LMBDdt,RK4,Pressure,Steady,ecrit_profil_b,ecrit_time_b,ecrit_tech_b, URG_STOP,adv_lith_bool,adv_interf_bool);
lecture_rad(adress_rad,RAD,tstop,HPE_factor);
lecture_sol(adress_solidus,solidus,liquidus);



double Tm,Tp,Tc,Phicrmax_N = 0,Phicrmax_S = 0,Dl_N,Dl_S,Dc_N,Dc_S,LMBD_ath,LMBD_lith_N,LMBD_lith_S,LMBD_cr_N,LMBD_cr_S,dmax_cr_N,dmax_cr_S;
double  dBLC_N,dBLC_S;
std::tuple<double,double> dTdz_N;
std::tuple<double,double> dTdz_S;
std::tuple<double,double,double> agecrust_N;
std::tuple<double,double,double> agecrust_S;

double Phi_eff_i = 0.1222;

double rho_c = std::get<4>(thermo);
double Dref = 0.2/3.*(pow(Rp,3.)-pow(Rc,3.))/pow(Rp,2.);
double Vc = 4./3.*M_PI*pow(Rc,3.);
double Ac = 4.*M_PI*pow(Rc,2.);
double Vp = 4./3.*M_PI*pow(Rp,3.) - Vc;
double t = 0;
double rho_p = (Masse_Mars - Vc*rho_c)/Vp;

if(unlinear_phi == 0){std::get<2>(melt_param)=1.0;}


char* dossier_time;
char* dossier_tech;
char* dossier_profilN;
char* dossier_profilS;

std:: string tmp1 = "time.txt";
dossier_time = &tmp1[0];
std:: string tmp4 = "tech.txt";
dossier_tech = &tmp4[0];
std::string tmp2 = "profilN.txt";
dossier_profilN = &tmp2[0];
std::string tmp3 = "profilS.txt";
dossier_profilS = &tmp3[0];

std::cout << "Début MODEL" << std::endl;

Parametric_model(thermo,rheology,melt,melt_param,RAD,solidus,liquidus,
Rp,Rc,f,Vc,Ac,Vp,rho_p,gc,gl,dr_min,dr_max,Dref,dt_max,dt_min,dt_init,
Tm0,dTc0,Ts,t0,tstop,t_acc,t_acc_m,DTMAX,fmag,fbase,Dl_init,Dc_init,dDl_init,dDc_init,LMBD_cr_Ni, Phi_eff_i, CH2O_p, Sat_expo, X1, X2,KH2Ocr, gamma_cr,phi_min_cr,
RK4,Steady,LMBDdt,N_melt, Nc,
t,Tm,Tc,Tp,Phicrmax_N,Phicrmax_S,Dl_N,Dl_S,Dc_N,Dc_S,dBLC_N,dBLC_S,agecrust_N,agecrust_S,LMBD_ath,LMBD_lith_N,LMBD_lith_S,LMBD_cr_N,LMBD_cr_S,dmax_cr_N,dmax_cr_S,dTdz_N,dTdz_S,
dossier_time,dossier_tech, dossier_profilN,dossier_profilS,ecrit_profil_b,ecrit_time_b,ecrit_tech_b,URG_STOP,adv_lith_bool,adv_interf_bool
);



t2 =clock();
temps =(float) (t2-t1)/CLOCKS_PER_SEC;
std::cout << "Temps d'execution : " << temps << " s." << std::endl;


}