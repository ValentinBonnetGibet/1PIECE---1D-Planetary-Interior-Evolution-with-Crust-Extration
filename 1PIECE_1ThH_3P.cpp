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
std::tuple<double,double,double,double> agecrust_N;
std::tuple<double,double,double,double> agecrust_S;

double Phi_vis_i = 0.1;
double Phi_obj = 0.3;

double rho_c = std::get<4>(thermo);
double Dref = 0.2/3.*(pow(Rp,3.)-pow(Rc,3.))/pow(Rp,2.);
double Vc = 4./3.*M_PI*pow(Rc,3.);
double Ac = 4.*M_PI*pow(Rc,2.);
double Vp = 4./3.*M_PI*pow(Rp,3.) - Vc;
double t = t0 *std::get<11>(thermo);
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

double eta_0 = std::get<0>(rheology);
double A =std::get<1>(rheology);
double R = std::get<2>(rheology);
double V = std::get<3>(rheology);
double Tref =std::get<5>(rheology);
double Pref =std::get<6>(rheology);
double beta_u = std::get<7>(rheology);
double Ra_crit_u = std::get<9>(rheology);
double C_m=std::get<0>(thermo);
double rho_cr=std::get<3>(thermo);
double k_m=std::get<5>(thermo);
double a_m=std::get<7>(thermo);

double DTsolidus=std::get<5>(melt)*Dc_init/Dref;
double delta_guess = 50E3;
double Phi_guess,Phi_eff_guess,Va_guess,dmadtm_guess,Rl_guess,LMBD_guess,Vm_guess,Ra_guess,eta_guess;
Rl_guess = Rp - (Dl_init-Dc_init + delta_guess) ;
// First Guess using 75E3 for delta_u

double Pm_guess = (Dc_init * gl * rho_cr + rho_p * (Rp-Rl_guess-Dc_init) * gl)/1E9 ;
double Tliq_guess = std::get<0>(liquidus)+std::get<1>(liquidus)*Pm_guess+std::get<2>(liquidus)*Pm_guess*Pm_guess+std::get<3>(liquidus)*Pm_guess*Pm_guess*Pm_guess;
double Tsol_guess = std::get<0>(solidus)+std::get<1>(solidus)*Pm_guess+std::get<2>(solidus)*Pm_guess*Pm_guess+std::get<3>(solidus)*Pm_guess*Pm_guess*Pm_guess+DTsolidus;
double Tm0_guess_new = Phi_obj * (Tliq_guess-Tsol_guess) + Tsol_guess;

double Tm0_guess = 0;

while (std::abs(Tm0_guess - Tm0_guess_new) > 1 )
{

Tm0_guess = Tm0_guess_new;
// Calculating delta_u with Tm0 by calculating Ra at Tm
double Tl_guess = Tm0_guess - std::get<4>(rheology) * R/A * Tm0_guess*Tm0_guess;
double Dm = k_m / rho_p / C_m;

melting(Phi_guess,Phi_eff_guess,Va_guess,dmadtm_guess,Tm0_guess,a_m,C_m,gc,gl,N_melt,Dref,std::get<4>(melt),std::get<5>(melt),rho_cr,rho_p,Dc_init,Dc_init,Rl_guess,Rp,delta_guess,std::get<3>(melt),std::get<6>(melt),
melt_param, LMBD_guess,
solidus,liquidus,rheology);  

Volume(Vm_guess,Rl_guess,Rl_guess,Rc,0);  

Phi_vis_i = Phi_guess * Va_guess / Vm_guess;

eta_guess = eta_0*std::exp(  (A + Pm_guess*1E9 * V )/ (R  * Tm0_guess ) - (A + Pref * V) / ( Tref * R) + std::get<10>(rheology) * Phi_vis_i   ) ;
Ra_guess = a_m*rho_p*gl* (Tm0_guess-Tl_guess) *pow(Rp-Rc,3.)/Dm/eta_guess ;

delta_guess = (Rp-Dl_init-Rc) * pow(Ra_crit_u/Ra_guess,beta_u);

// Second Guess using delta_guess for delta_u
Pm_guess = (Dc_init * gl * rho_cr + rho_p * (Dl_init-Dc_init + delta_guess) * gl)/1E9 ;

Tliq_guess = std::get<0>(liquidus)+std::get<1>(liquidus)*Pm_guess+std::get<2>(liquidus)*Pm_guess*Pm_guess+std::get<3>(liquidus)*Pm_guess*Pm_guess*Pm_guess;
Tsol_guess = std::get<0>(solidus)+std::get<1>(solidus)*Pm_guess+std::get<2>(solidus)*Pm_guess*Pm_guess+std::get<3>(solidus)*Pm_guess*Pm_guess*Pm_guess+DTsolidus;

Tm0_guess_new = Phi_obj * (Tliq_guess-Tsol_guess) + Tsol_guess;
    /* code */
}

Tm0 = Tm0_guess_new;

std::cout << "Début MODEL" << std::endl;

Parametric_model(thermo,rheology,melt,melt_param,RAD,solidus,liquidus,
Rp,Rc,f,Vc,Ac,Vp,rho_p,gc,gl,dr_min,dr_max,Dref,dt_max,dt_min,dt_init,
Tm0,dTc0,Ts,t0,tstop,t_acc,t_acc_m,DTMAX,fmag,fbase,Dl_init,Dc_init,dDl_init,dDc_init,LMBD_cr_Ni, Phi_vis_i, CH2O_p, Sat_expo, X1, X2,KH2Ocr, gamma_cr,phi_min_cr,
RK4,Steady,LMBDdt,N_melt, Nc,
t,Tm,Tc,Tp,Phicrmax_N,Phicrmax_S,Dl_N,Dl_S,Dc_N,Dc_S,dBLC_N,dBLC_S,agecrust_N,agecrust_S,LMBD_ath,LMBD_lith_N,LMBD_lith_S,LMBD_cr_N,LMBD_cr_S,dmax_cr_N,dmax_cr_S,dTdz_N,dTdz_S,
dossier_time,dossier_tech, dossier_profilN,dossier_profilS,ecrit_profil_b,ecrit_time_b,ecrit_tech_b,URG_STOP,adv_lith_bool,adv_interf_bool
);

std::cout << "Initial Phi : " << Phi_vis_i << std::endl;
std::cout << "Temps extraction nord (Gyr) : "  << std::get<3>(agecrust_N)/1E9  << std::endl;
std::cout << "Temps extraction sud (Gyr) :"  << std::get<3>(agecrust_S)/1E9  << std::endl;


t2 =clock();
temps =(float) (t2-t1)/CLOCKS_PER_SEC;
std::cout << "Temps d'execution : " << temps << " s." << std::endl;


}