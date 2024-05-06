#include "1PIECE_main.cpp"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <random>
#include <mpi.h>

bool lecture_MPI_init(std::string const &adresse,int const &rank, 
double const &Dc,double const &Rp, double const &Rc, double const &f, double const &Vp, double const &rhop, double const &Dref,
double const &g, double const &Cm, double const &am, int const &N_melt, double const & Di,double const &km, double const &Dl0_min,double const &Dl0_max, double const &dDl0,
double const &Phi_obj, double &Phi,
std::tuple<bool,double,double,double,double,double,double,double,double> const &melt_param,std::tuple<double,double,double,double,double,double> const &solidus, std::tuple<double,double,double,double,double,double> const &liquidus, std::tuple<double,double,double,double,double,double,double,double> const &melt,
std::tuple<double,double,double,double,double,double,double,double,double,double,double> const &rheology,
double &eta0, double &Tm0, double &k0, double &rhocr, double &kcr, double &A,double &V,double &fmag,double &fbase,double &CH2O,double &Dl0, double &dDl_init, double &dTc0);

bool lecture_MPI_restart(std::string const &adresse, int const &rank,  int &n_start,
double const &Dc,double const &Rp, double const &Rc, double const &f, double const &Vp, double const &rhop, double const &Dref,
double const &g, double const &Cm, double const &am, int const &N_melt, double const & Di, double const &km, double const &Dl0_min, double const &Dl0_max, double const &dDl0,
double const &Phi_obj, double &Phi, 
std::tuple<bool,double,double,double,double,double,double,double,double> const &melt_param,std::tuple<double,double,double,double,double,double> const &solidus, std::tuple<double,double,double,double,double,double> const &liquidus, std::tuple<double,double,double,double,double,double,double,double> const &melt,
std::tuple<double,double,double,double,double,double,double,double,double,double,double> const &rheology,
double &eta0, double &Tm0, double &k0, double &rhocr, double &kcr, double &A,double &V,double &fmag,double &fbase,double &CH2O,double &Dl0, double &dDl_nit, double &dTc0);

bool lecture_random(std::string const &adresse,double &Phi_obj, int &N_MCMC, bool &MCMC_restart,
    double &k0_max, double &k0_min, double &dk0,
    double &eta0_max, double &eta0_min, double &deta0,
    double &fmag_max, double  &fmag_min, double &dfmag,
    double &fbase_max, double  &fbase_min, double &dfbase,
    double &Tc0_max, double  &Tc0_min, double &dTc0,
    double &A_max, double  &A_min, double &dA,
    double &V_max, double &V_min, double &dV,
    double &rhocr_max, double &rhocr_min, double &drhocr,
    double &kcr_max, double &kcr_min, double &dkcr,
    double &dDl_max, double &dDl_min, double &d2Dl,
    double &CH2O_max, double &CH2O_min, double &dCH2O,
    double &Tm0_max, double &Tm0_min, double &dTm,
    double &Dl0_max, double &Dl0_min, double &dDl0
   

);

bool lecture_likelihood(std::string const &adresse,
    double &Dc_insight_mean, double &Dc_insight_sigma, double &Dl_insight_mean, double &Dl_insight_sigma, 
    double &Tp_insight_mean, double &Tp_insight_sigma,
    bool &Dl_likelihood_bool, bool &dDc_likelihood_bool, bool &Dc_likelihood_bool, bool &Tp_likelihood_bool,
    double &thick_slope, double &thick_y0, double &dich_slope, double &dich_y0, double &dich_slope_sigma, double &dich_y0_sigma, double &delta_rhom
);

double random_param(
    double &k0, double const &k0_max, double const &k0_min, double const &dk0,
    double &eta0, double const &eta0_max, double const &eta0_min, double const &deta0,
    double &fmag, double const &fmag_max, double const &fmag_min, double const &dfmag,
    double &fbase, double const &fbase_max, double const &fbase_min, double const &dfbase,
    double &Tc0, double const &Tc0_max, double const &Tc0_min, double const &dTc0,
    double &A, double const &A_max, double const &A_min, double const &dA,
    double &V, double const &V_max, double const &V_min, double const &dV,
    double &rhocr, double const &rhocr_max, double const &rhocr_min, double const &drhocr,
    double &kcr, double const &kcr_max, double const &kcr_min, double const &dkcr,
    double &dDl_init, double const &dDl_max, double const &dDl_min, double const &d2Dl,
    double &CH2O, double const &CH2O_max, double const &CH2O_min, double const &dCH2O,
    double &Tm0, double const &Tm0_max, double const &Tm0_min, double const &dTm,
    double &Dl0,  double const &Dl0_max,  double const &Dl0_min, double const &dDl0,
    double const &Dc,double const &Rp, double const &Rc, double const &f, double const &Vp, double const &rhop, double const &Dref,
    double const &g, double const &km, double const &Cm, double const &am, int const &N_melt, double const & Di,
    double const &Phi_obj, double &Phi,
    std::tuple<bool,double,double,double,double,double,double,double,double> const &melt_param,std::tuple<double,double,double,double,double,double> const &solidus, std::tuple<double,double,double,double,double,double> const &liquidus, std::tuple<double,double,double,double,double,double,double,double> const &melt,
    std::tuple<double,double,double,double,double,double,double,double,double,double,double> const &rheology

);

double rand_01();

double rand_gauss(double const &mean, double const &sigma);

double likelihood(double const &DcN, double const &DcS, double const &DlN, double const &DlS, double const &Tp, double const &Dc_insight_mean, double const &Dc_insight_sigma,
double const &Dl_insight_mean, double const &Dl_insight_sigma, double const &Tp_insight_mean, double const &Tp_insight_sigma,
bool const &Dl_likelihood_bool,bool const &dDc_likelihood_bool, bool const &Dc_likelihood_bool, bool const &Tp_likelihood_bool,
double const &thick_slope, double const &thick_y0, double const &dich_slope, double const &dich_y0, double const &dich_slope_sigma, double const &dich_y0_sigma,
double const &delta_rhom, double const &f, double const &Rp, double const &rhocr, double const &rhom, double const &dBLC_N, double const &dBLC_S

);


int main(int argc,char* argv[]){

int nb_procs,rank;  
MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD,&nb_procs); 
MPI_Comm_rank(MPI_COMM_WORLD,&rank);

float temps, temps2;
clock_t t1,t2,t3;

t1 = clock();
std::tuple<double,double,double,double,double,double,double,double,double,double,double> rheology;
std::tuple<double,double,double,double,double,double,double,double> melt;
std::tuple<bool,double,double,double,double,double,double,double,double> melt_param;
std::tuple<double,double,double,double,double,double,double,double> unlinear_table;
std::tuple<double,double,double,double,double,double,double,double,double,double,double,double,double> thermo;

Eigen::MatrixXd RAD = Eigen::MatrixXd::Zero(4,4);

double Tm0,Ts, dTc0, Dl_init, Dc_init , dDl_init, dDc_init, LMBD_cr_Ni, LMBD_cr_Si, t0, tstop, Rp, Rc, f ,Masse_Mars, gl,gc,dt_min,dt_max,dt_init,DTMAX,dr_min,dr_max,t_acc=0.0,t_acc_m=0.0,fmag=0.0,fbase=0.0, HPE_factor = 1;
double X1,X2,Sat_expo,CH2O_p, KH2Ocr, gamma_cr,phi_min_cr;
int N_soliq, N_melt,Nc;
bool unlinear_phi,LMBDdt,RK4,Pressure,Steady, ecrit_time_b, ecrit_tech_b, ecrit_profil_b,URG_STOP,adv_lith_bool,adv_interf_bool, MCMC_restart;
// bool melt_bol = 1;
auto solidus = std::make_tuple(1E2,1E2,1E2,1E2,1E2,1E2);
auto liquidus = std::make_tuple(1E2,1E2,1E2,1E2,1E2,1E2);
bool dDc_likelihood_bool,Dc_likelihood_bool, Dl_likelihood_bool, Tp_likelihood_bool;
bool test_run = 1;

// /Users/valentinbonnetgibet/Documents/These/Code Mars/Mars Cpp Dev/
std::string physical_param = "physical_parameters.txt";
std::string numerical_param = "numerical_parameters.txt";
std::string adress_rad = "Radioelement_Wanke.txt";
std::string adress_solidus = "Solidus.txt";
std::string adresse_random = "MCMC_random.txt";
std::string adresse_likelihood = "Likehood_param.txt";
std::string adresse_MPI_init = "MPI_MCMC_entry.txt";

char* adresse_ecriture_MCMC;
std:: string adresse_string_ecriture_MCMC = "MCMC_res_" + std::to_string(rank) + ".txt";
adresse_ecriture_MCMC = &adresse_string_ecriture_MCMC[0];
 


double Dc_insight_mean, Dc_insight_sigma, Dl_insight_mean, Dl_insight_sigma;
double Tp_insight_mean, Tp_insight_sigma, thick_slope, thick_y0;
double dich_slope, dich_y0, dich_slope_sigma, dich_y0_sigma, delta_rhom;


double k0_max,k0_min,dk0,eta0_max, eta0_min, deta0;
double fmag_max, fmag_min, dfmag, fbase_max, fbase_min, dfbase;
double DTc0_max, DTc0_min, dDTc0;
double A_max, A_min, dA,V_max, V_min, dV;
double rhocr_max, rhocr_min,drhocr,kcr_max,kcr_min,dkcr;
double dDl_max, dDl_min, d2Dl,CH2O_max, CH2O_min,dCH2O;
double Tm0_max = 1, Tm0_min, dTm, dDl0, Dl0_min, Dl0_max;
double Phi_obj,Phi_init;
int N_MCMC,n_start;
int rejected =0, accepted = 0;
double acceptance = 0.0;

double logproba1,logproba2, prior;

Lecture_param(physical_param,rheology,melt,melt_param,thermo,Tm0,Ts,dTc0,Dl_init,Dc_init,dDl_init, dDc_init, LMBD_cr_Ni,LMBD_cr_Si,t0,tstop,Rp,Rc,f,Masse_Mars,gl,gc,fmag,fbase,HPE_factor,X1,X2,Sat_expo,CH2O_p,KH2Ocr,gamma_cr,phi_min_cr);
lecture_num(numerical_param,dt_min,dt_max, dt_init,DTMAX,dr_min,dr_max,Nc,N_soliq,N_melt,t_acc,t_acc_m,unlinear_phi,LMBDdt,RK4,Pressure,Steady,ecrit_profil_b,ecrit_time_b,ecrit_tech_b, URG_STOP,adv_lith_bool,adv_interf_bool);
lecture_sol(adress_solidus,solidus,liquidus);
lecture_random(adresse_random, Phi_obj,N_MCMC, MCMC_restart,
k0_max,k0_min,dk0,eta0_max, eta0_min, deta0,
fmag_max, fmag_min, dfmag, fbase_max, fbase_min, dfbase,
DTc0_max, DTc0_min, dDTc0,
A_max, A_min, dA,V_max, V_min, dV,
rhocr_max, rhocr_min,drhocr,kcr_max,kcr_min,dkcr,
dDl_max, dDl_min, d2Dl,CH2O_max, CH2O_min,dCH2O,
Tm0_max, Tm0_min, dTm, Dl0_max, Dl0_min, dDl0
 
);

lecture_likelihood(adresse_likelihood,
Dc_insight_mean, Dc_insight_sigma, Dl_insight_mean, Dl_insight_sigma, 
Tp_insight_mean, Tp_insight_sigma,
Dl_likelihood_bool, dDc_likelihood_bool, Dc_likelihood_bool,Tp_likelihood_bool, thick_slope, thick_y0,
dich_slope, dich_y0, dich_slope_sigma, dich_y0_sigma, delta_rhom
);


double Tm,Tc,Phicrmax_N = 0,Phicrmax_S = 0,Dl_N,Dl_S,Dc_N,Dc_S,LMBD_ath,LMBD_lith_N,LMBD_lith_S,LMBD_cr_N,LMBD_cr_S,dmax_cr_N,dmax_cr_S, Tp, dBLC_N,dBLC_S;
std::tuple<double,double> dTdz_N;
std::tuple<double,double> dTdz_S;
std::tuple<double,double,double> agecrust_N;
std::tuple<double,double,double> agecrust_S;

double Phi_eff_i = Phi_obj;

double rho_c=std::get<4>(thermo);
double Dref = 0.2/3.*(pow(Rp,3.)-pow(Rc,3.))/pow(Rp,2.);
double Vc = 4./3.*M_PI*pow(Rc,3.);
double Ac = 4.*M_PI*pow(Rc,2.);
double Vp = 4./3.*M_PI*pow(Rp,3.) - Vc;
double Vm,Vcr, rho_m;
double t = 0;
double rho_p = (Masse_Mars - Vc*rho_c)/Vp;

if(MCMC_restart == 0){
lecture_MPI_init(adresse_MPI_init, rank, 
Dc_init, Rp, Rc, f, Vp, rho_p, Dref,
gl, std::get<0>(thermo), std::get<7>(thermo),  N_melt, std::get<9>(thermo),std::get<5>(thermo),
Dl0_min, Dl0_max, dDl0,
Phi_obj, Phi_init, 
melt_param,solidus, liquidus, melt,rheology,
std::get<0>(rheology), Tm0, std::get<0>(melt), std::get<3>(thermo), std::get<6>(thermo), std::get<1>(rheology), std::get<3>(rheology), fmag, fbase, CH2O_p, Dl_init, dDl_init, dTc0);

std::ofstream fichier_MCMC(adresse_ecriture_MCMC);
fichier_MCMC <<  "HEADER " << std::endl;
fichier_MCMC.close();
n_start = 0;
}

else{

lecture_MPI_restart(adresse_string_ecriture_MCMC, rank, n_start,
Dc_init, Rp, Rc, f, Vp, rho_p, Dref,
gl, std::get<0>(thermo), std::get<7>(thermo),  N_melt, std::get<9>(thermo),std::get<5>(thermo),
Dl0_min, Dl0_max, dDl0,
Phi_obj, Phi_init, 
melt_param,solidus, liquidus, melt,rheology,
std::get<0>(rheology), Tm0, std::get<0>(melt), std::get<3>(thermo), std::get<6>(thermo), std::get<1>(rheology), std::get<3>(rheology), fmag, fbase, CH2O_p, Dl_init, dDl_init, dTc0);
    
}


lecture_rad(adress_rad,RAD,tstop,HPE_factor);

if(unlinear_phi == 0){std::get<1>(unlinear_table)=1;}


char* dossier_time;
char* dossier_tech;
char* dossier_profilN;
char* dossier_profilS;

std:: string tmp1 = "time"+std::to_string(rank)+".txt";
dossier_time = &tmp1[0];
std:: string tmp4 = "tech"+std::to_string(rank)+".txt";
dossier_tech = &tmp4[0];
std::string tmp2 = "profilN.txt";
dossier_profilN = &tmp2[0];
std::string tmp3 = "profilS.txt";
dossier_profilS = &tmp3[0];




t3 = clock();

temps2 = (float) (t3-t1)/CLOCKS_PER_SEC;
std::cout << rank << " : " << "Before init, temps = "<< temps2 <<std::endl;


std::cout << "DEBUGG, km : " << std::get<5>(thermo) << "dDl init : "<< dDl_init << ", Dc init : " << dDc_init << std::endl;

test_run = Parametric_model(thermo,rheology,melt,melt_param,RAD,solidus,liquidus,
Rp,Rc,f,Vc,Ac,Vp,rho_p,gc,gl,dr_min,dr_max,Dref,dt_max,dt_min,dt_init,
Tm0,dTc0,Ts,t0,tstop,t_acc,t_acc_m,DTMAX,fmag,fbase,Dl_init,Dc_init,dDl_init,dDc_init,LMBD_cr_Ni, Phi_eff_i, CH2O_p, Sat_expo, X1, X2,KH2Ocr, gamma_cr,phi_min_cr,
RK4,Steady,LMBDdt,N_melt, Nc,
t,Tm,Tc,Tp,Phicrmax_N,Phicrmax_S,Dl_N,Dl_S,Dc_N,Dc_S,dBLC_N,dBLC_S,agecrust_N,agecrust_S,LMBD_ath,LMBD_lith_N,LMBD_lith_S,LMBD_cr_N,LMBD_cr_S,dmax_cr_N,dmax_cr_S,dTdz_N,dTdz_S,
dossier_time,dossier_tech, dossier_profilN,dossier_profilS,ecrit_profil_b,ecrit_time_b,ecrit_tech_b,URG_STOP,adv_lith_bool,adv_interf_bool
);


Volume(Vm,Rp-Dc_N,Rp-Dc_S,Rc,f);
Vcr=Vp-Vm;
rho_m=(Vp*rho_p-Vcr*std::get<3>(thermo))/Vm;

logproba1 = likelihood(Dc_N,Dc_S,Dl_N,Dl_S,Tp,
Dc_insight_mean, Dc_insight_sigma, Dl_insight_mean, Dl_insight_sigma, 
Tp_insight_mean, Tp_insight_sigma, thick_slope, thick_y0,
Dl_likelihood_bool,dDc_likelihood_bool, Dc_likelihood_bool, Tp_likelihood_bool,
dich_slope, dich_y0, dich_slope_sigma, dich_y0_sigma, delta_rhom,
f,Rp,std::get<3>(thermo),rho_m,dBLC_N,dBLC_S);


double k0_prev = std::get<0>(melt);
double eta0_prev = std::get<0>(rheology);
double fmag_prev = fmag;
double fbase_prev = fbase;
double dTc0_prev = dTc0;
double A_prev = std::get<1>(rheology);
double V_prev = std::get<3>(rheology);
double rhocr_prev =std::get<3>(thermo);
double kcr_prev =std::get<6>(thermo);
double dDl_prev = dDl_init;
double CH2Op_prev = CH2O_p;
double Tm0_prev = Tm0;
double Dl_init_prev = Dl_init;
double Phi_init_prev = Phi_init;
double DcN_prev = Dc_N;
double DcS_prev = Dc_S;
double DlN_prev = Dl_N;
double DlS_prev = Dl_S;
double Tp_prev = Tp;
double Tm_prev = Tm;
double Tc_prev = Tc;
double Phicrmax_N_prev = Phicrmax_N;
double Phicrmax_S_prev = Phicrmax_S;
double dmax_cr_N_prev = dmax_cr_N;
double dmax_cr_S_prev = dmax_cr_S;
std::tuple<double,double> dTdz_N_prev = dTdz_N;
std::tuple<double,double> dTdz_S_prev = dTdz_S;
std::tuple<double,double,double> agecrust_N_prev =  agecrust_N;
std::tuple<double,double,double> agecrust_S_prev = agecrust_S;
double dBLC_N_prev =dBLC_N;
double dBLC_S_prev =dBLC_S;
double rho_m_prev = rho_m;
double LMBD_cr_N_prev = LMBD_cr_N;
double LMBD_cr_S_prev = LMBD_cr_S;




t2 = clock();
temps2 = (float) (t2-t1)/CLOCKS_PER_SEC;
std::cout << rank << " : rank," << " Before star of the chain, temps = "<< temps2 <<std::endl;

// DEBUT DU LA BOUCLE
for(int i = n_start; i< N_MCMC; ++i){


t2 = clock();
t = 0;
prior = random_param(
std::get<0>(melt), k0_max, k0_min, dk0,
std::get<0>(rheology), eta0_max, eta0_min, deta0,
fmag, fmag_max, fmag_min, dfmag,
fbase, fbase_max, fbase_min, dfbase,
dTc0, DTc0_max, DTc0_min, dDTc0,
std::get<1>(rheology) , A_max, A_min, dA,
std::get<3>(rheology), V_max, V_min, dV,
std::get<3>(thermo), rhocr_max, rhocr_min, drhocr,
std::get<6>(thermo), kcr_max, kcr_min, dkcr,
dDl_init, dDl_max, dDl_min, d2Dl,
CH2O_p, CH2O_max, CH2O_min, dCH2O,
Tm0, Tm0_max, Tm0_min, dTm,
Dl_init, Dl0_max, Dl0_min, dDl0,
Dc_init, Rp, Rc, f, Vp, rho_p, Dref,
gl, std::get<5>(thermo), std::get<0>(thermo), std::get<7>(thermo),  N_melt, std::get<9>(thermo),
Phi_obj, Phi_init, melt_param,solidus, liquidus, melt,rheology);



if(prior > 0){
std::cout << rank << " : rank," << "New Model choosen :" << std::endl;
std::cout << rank << " : rank,"  << " k0 : " << std::get<0>(melt) << " eta0 : " << std::get<0>(rheology) << ", fmag : " << fmag 
<< ", fbase : " << fbase << ", A : " << std::get<1>(rheology) << std::endl;
std::cout << "V : " << std::get<3>(rheology)  << " rhocr : " << std::get<3>(thermo) 
<< ", kcr : " <<std::get<6>(thermo) << std::endl;
std::cout << "CH2O : " << CH2O_p << ", Tm0 : " << Tm0 << ", dDl_init : " << dDl_init <<", dTc0 : " << dTc0 << " Dl_init : " << Dl_init 
<< ", Phi_init : " << Phi_init  << std::endl;



test_run = Parametric_model(thermo,rheology,melt,melt_param,RAD,solidus,liquidus,
Rp,Rc,f,Vc,Ac,Vp,rho_p,gc,gl,dr_min,dr_max,Dref,dt_max,dt_min,dt_init,
Tm0,dTc0,Ts,t0,tstop,t_acc,t_acc_m,DTMAX,fmag,fbase,Dl_init,Dc_init,dDl_init,dDc_init,LMBD_cr_Ni, Phi_eff_i, CH2O_p, Sat_expo, X1, X2,KH2Ocr, gamma_cr,phi_min_cr,
RK4,Steady,LMBDdt,N_melt, Nc,
t,Tm,Tc,Tp,Phicrmax_N,Phicrmax_S,Dl_N,Dl_S,Dc_N,Dc_S,dBLC_N,dBLC_S,agecrust_N,agecrust_S,LMBD_ath,LMBD_lith_N,LMBD_lith_S,LMBD_cr_N,LMBD_cr_S,dmax_cr_N,dmax_cr_S,dTdz_N,dTdz_S,
dossier_time,dossier_tech, dossier_profilN,dossier_profilS,ecrit_profil_b,ecrit_time_b,ecrit_tech_b,URG_STOP,adv_lith_bool,adv_interf_bool
);

}
else
{
    std::cout << rank << " : rank," << "Prior null" << std::endl;
    test_run = 1;
}


if(test_run == 0){
Volume(Vm,Rp-Dc_N,Rp-Dc_S,Rc,f);
Vcr=Vp-Vm;
rho_m=(Vp*rho_p-Vcr*std::get<3>(thermo))/Vm;

logproba2 = likelihood(Dc_N,Dc_S,Dl_N,Dl_S,Tp,
Dc_insight_mean, Dc_insight_sigma, Dl_insight_mean, Dl_insight_sigma, Tp_insight_mean, Tp_insight_sigma,
Dl_likelihood_bool,dDc_likelihood_bool, Dc_likelihood_bool, Tp_likelihood_bool,
thick_slope, thick_y0, dich_slope, dich_y0, dich_slope_sigma, dich_y0_sigma, delta_rhom,
f,Rp,std::get<3>(thermo),rho_m,dBLC_N,dBLC_S);
}
else{
    logproba2 = std::log(0);
}


t3 = clock();

temps = (float) (t3-t2)/CLOCKS_PER_SEC;

if(logproba2 - logproba1 > std::log(rand_01())){
logproba1 = logproba2;

accepted++;
acceptance = (double)accepted  / ((double)accepted + (double)rejected ) * 100.0;

FILE * fichier = fopen(adresse_ecriture_MCMC,"a");
std::fprintf(fichier, " %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
std::get<0>(rheology),  Tm0, std::get<0>(melt),std::get<3>(thermo),  std::get<6>(thermo), std::get<1>(rheology),std::get<3>(rheology), fmag, fbase, CH2O_p, dDl_init, dTc0, Dl_init,Phi_init, temps,
Tp, Tm, Dc_N, Dc_S, Dl_N, Dl_S, Tc, Phicrmax_N, Phicrmax_S, dmax_cr_N,dmax_cr_S, std::get<0>(agecrust_N) ,std::get<0>(agecrust_S), std::get<1>(agecrust_N) ,std::get<1>(agecrust_S),std::get<2>(agecrust_N) ,std::get<2>(agecrust_S), rho_m,LMBD_cr_N,LMBD_cr_S,std::get<0>(dTdz_N), std::get<0>(dTdz_S), std::get<1>(dTdz_N), std::get<1>(dTdz_S), dBLC_N,dBLC_S, logproba1, acceptance, 1.0
);



fclose(fichier);
std::cout << rank << " : rank," << "Boucle n = " << i  << " accepting, acceptance : " << acceptance << " %" << std::endl;

k0_prev = std::get<0>(melt);
eta0_prev = std::get<0>(rheology);
fmag_prev = fmag;
fbase_prev = fbase;
dTc0_prev = dTc0;
A_prev = std::get<1>(rheology);
V_prev = std::get<3>(rheology);
rhocr_prev =std::get<3>(thermo);
kcr_prev =std::get<6>(thermo);
dDl_prev =dDl_init;
CH2Op_prev = CH2O_p;
Tm0_prev = Tm0;
Dl_init_prev = Dl_init;
Phi_init_prev = Phi_init;
DcN_prev = Dc_N;
DcS_prev = Dc_S;
DlN_prev = Dl_N;
DlS_prev = Dl_S;
Tp_prev = Tp;
Tm_prev = Tm;
Tc_prev = Tc;
Phicrmax_N_prev = Phicrmax_N;
Phicrmax_S_prev = Phicrmax_S;
dmax_cr_N_prev = dmax_cr_N;
dmax_cr_S_prev = dmax_cr_S;
dTdz_N_prev = dTdz_N;
dTdz_S_prev = dTdz_S;
agecrust_N_prev =  agecrust_N;
agecrust_S_prev = agecrust_S;
dBLC_N_prev =dBLC_N;
dBLC_S_prev =dBLC_S;
rho_m_prev = rho_m;
LMBD_cr_N_prev = LMBD_cr_N;
LMBD_cr_S_prev = LMBD_cr_S;
dTc0_prev = dTc0 ; 

}

else{
    
std::get<0>(melt) = k0_prev;
std::get<0>(rheology) = eta0_prev;
fmag = fmag_prev;
fbase = fbase_prev;
dTc0 = dTc0_prev;
std::get<1>(rheology)= A_prev;
std::get<3>(rheology)= V_prev;
std::get<3>(thermo) = rhocr_prev ;
std::get<6>(thermo) = kcr_prev;
dDl_init = dDl_prev;
CH2O_p = CH2Op_prev;
Tm0 = Tm0_prev;
Dl_init = Dl_init_prev;
Phi_init = Phi_init_prev;
Dc_N = DcN_prev;
Dc_S =  DcS_prev;
Dl_N = DlN_prev;
Dl_S = DlS_prev;
Tp = Tp_prev;
Tm = Tm_prev;
Tc = Tc_prev;
Phicrmax_N = Phicrmax_N_prev;
Phicrmax_S = Phicrmax_S_prev;
dmax_cr_N = dmax_cr_N_prev;
dmax_cr_S = dmax_cr_S_prev;
dTdz_N = dTdz_N_prev;
dTdz_S = dTdz_S_prev;
agecrust_N=  agecrust_N_prev ;
agecrust_S = agecrust_S_prev;
dBLC_N = dBLC_N_prev;
dBLC_S = dBLC_S_prev;
rho_m = rho_m_prev ;
LMBD_cr_N = LMBD_cr_N_prev ;
LMBD_cr_S = LMBD_cr_S_prev ; 
dTc0 = dTc0_prev ; 


rejected++;
acceptance = (double)accepted  / ((double)accepted + (double)rejected ) * 100.0;
std::cout << rank << " : rank," << "Boucle n = " << i  << " rejecting, acceptance : " << acceptance << " %" << std::endl;


FILE * fichier = fopen(adresse_ecriture_MCMC,"a");
std::fprintf(fichier, " %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e \n",
std::get<0>(rheology),  Tm0, std::get<0>(melt),std::get<3>(thermo),  std::get<6>(thermo), std::get<1>(rheology),std::get<3>(rheology), fmag, fbase, CH2O_p,dDl_init, dTc0, Dl_init,Phi_init, temps,
Tp, Tm, Dc_N, Dc_S, Dl_N, Dl_S, Tc, Phicrmax_N, Phicrmax_S, dmax_cr_N,dmax_cr_S, std::get<0>(agecrust_N) ,std::get<0>(agecrust_S), std::get<1>(agecrust_N) ,std::get<1>(agecrust_S),std::get<2>(agecrust_N) ,std::get<2>(agecrust_S), rho_m,LMBD_cr_N,LMBD_cr_S,std::get<0>(dTdz_N), std::get<0>(dTdz_S), std::get<1>(dTdz_N), std::get<1>(dTdz_S), dBLC_N,dBLC_S, logproba1, acceptance, 0.0
);
   

fclose(fichier);


}

}

t2 = clock();;

temps2 = (float) (t3-t1)/CLOCKS_PER_SEC;

std::cout << rank << " : rank," << "After, temps = "<< temps2 << ", rejected  = " << rejected << ", acceptance : " << acceptance << "% ."<< std::endl;

MPI_Finalize();

}


double likelihood(double const &DcN, double const &DcS, double const &DlN, double const &DlS, double const &Tp, double const &Dc_insight_mean, double const &Dc_insight_sigma,
double const &Dl_insight_mean, double const &Dl_insight_sigma, double const &Tp_insight_mean, double const &Tp_insight_sigma,
bool const &Dl_likelihood_bool, bool const &dDc_likelihood_bool,  bool const &Dc_likelihood_bool, bool const &Tp_likelihood_bool,
double const &thick_slope, double const &thick_y0, double const &dich_slope, double const &dich_y0, double const &dich_slope_sigma, double const &dich_y0_sigma,
double const &delta_rhom, double const &f, double const &Rp, double const &rhocr, double const &rhom, double const &dBLC_N, double const &dBLC_S

){

    double logproba, logp1, logp2, logp3, logp4;
    double rho_ratio = rhocr / (rhom-rhocr-delta_rhom);
    double Dc_avg = (Rp-pow((f*pow(Rp-DcN,3.)+(1-f)*pow(Rp-DcS,3.)),1./3.))/1000;
    double Dc_calc = Dc_avg - rho_ratio * thick_slope - thick_y0;
    double dDc_calc = (DcS - DcN)/1000;
    double Dl_calc = (Rp-pow(f*pow(Rp-DlN-dBLC_N,3.)+(1-f)*pow(Rp-DlS-dBLC_S,3.),1./3.))/1000;
    double dDc_insight_mean = dich_slope * rho_ratio + dich_y0;
    double dDc_insight_sigma = (dich_slope_sigma*rho_ratio+dich_y0_sigma) ;
   

    // std::cout << " dDc_calc : " << dDc_calc << ", dDc_insight_mean : " << dDc_insight_mean << ", dDc_insight_sigma : " << dDc_insight_sigma << ", ratio : " << rho_ratio << std::endl;
    if(Dc_likelihood_bool == 1){
        logp1 = -0.5*pow((Dc_calc-Dc_insight_mean)/Dc_insight_sigma,2.) - std::log(Dc_insight_sigma * std::sqrt(2*M_PI));
    
    } 
    else{
        logp1 = 0.0;
       
    }

    if(dDc_likelihood_bool == 1){
        
        logp2 = -0.5*pow((dDc_calc-dDc_insight_mean)/dDc_insight_sigma,2.) - std::log(dDc_insight_sigma * std::sqrt(2*M_PI));
    } 
    else{
       
        logp2 = 0.0;
    }
    if(Dl_likelihood_bool == 1){
         logp3 = -0.5*pow((Dl_calc-Dl_insight_mean)/Dl_insight_sigma,2.) - std::log(Dl_insight_sigma * std::sqrt(2*M_PI));
        } 
    else{
        logp3 = 0.0;
    }
    if(Tp_likelihood_bool == 1){
        logp4 = -0.5*pow((Tp-Tp_insight_mean)/Tp_insight_sigma,2.) - std::log(Tp_insight_sigma * std::sqrt(2*M_PI));
        } 
    else{
        logp4 = 0.0;
    }
    // std::cout << "proba Dc : " << logp1 << ", proba dDc : " << logp2 << ", proba Dl : " << logp3 << ", proba Tp : " << logp4 << std::endl;

    logproba = logp1 + logp2 + logp3 +  logp4;

    return logproba;
}

bool lecture_likelihood(std::string const &adresse,
    double &Dc_insight_mean, double &Dc_insight_sigma, double &Dl_insight_mean, double &Dl_insight_sigma, 
    double &Tp_insight_mean, double &Tp_insight_sigma,
    bool &Dl_likelihood_bool, bool &dDc_likelihood_bool, bool &Dc_likelihood_bool, bool &Tp_likelihood_bool,
    double &thick_slope, double &thick_y0, double &dich_slope, double &dich_y0, double &dich_slope_sigma, double &dich_y0_sigma, double &delta_rhom
)
{

    std::ifstream fichier(adresse);
    bool test = true;

    if(fichier)
    {
        std::string ligne_txt;
        std::string text_ign;
        getline(fichier,ligne_txt);
        fichier.ignore();
        fichier >> text_ign >> Dc_insight_mean >> Dc_insight_sigma >> Dc_likelihood_bool >>
        text_ign >> Dl_insight_mean >> Dl_insight_sigma >>   Dl_likelihood_bool >>
        text_ign >> Tp_insight_mean >> Tp_insight_sigma >>   Tp_likelihood_bool >>
        text_ign >> thick_slope >>  thick_y0 >> 
        text_ign >> dich_slope >>  dich_y0 >> dDc_likelihood_bool >>
        text_ign >> dich_slope_sigma >>  dich_y0_sigma >>
        text_ign >> delta_rhom;

        fichier.close();
    }

    else
    {
        std::cout << "ERREUR: Impossible d'ouvrir le fichier lecture_likelihood en lecture." << std::endl;
        test = false;
    }

return test;

}


double random_param(
    double &k0, double const &k0_max, double const &k0_min, double const &dk0,
    double &eta0, double const &eta0_max, double const &eta0_min, double const &deta0,
    double &fmag, double const &fmag_max, double const &fmag_min, double const &dfmag,
    double &fbase, double const &fbase_max, double const &fbase_min, double const &dfbase,
    double &Tc0, double const &Tc0_max, double const &Tc0_min, double const &dTc0,
    double &A, double const &A_max, double const &A_min, double const &dA,
    double &V, double const &V_max, double const &V_min, double const &dV,
    double &rhocr, double const &rhocr_max, double const &rhocr_min, double const &drhocr,
    double &kcr, double const &kcr_max, double const &kcr_min, double const &dkcr,
    double &dDl_init, double const &dDl_max, double const &dDl_min, double const &d2Dl,
    double &CH2O, double const &CH2O_max, double const &CH2O_min, double const &dCH2O,
    double &Tm0, double const &Tm0_max, double const &Tm0_min, double const &dTm,
    double &Dl0,  double const &Dl0_max,  double const &Dl0_min, double const &dDl0,
    double const &Dc,double const &Rp, double const &Rc, double const &f, double const &Vp, double const &rhop, double const &Dref,
    double const &g, double const &km, double const &Cm, double const &am, int const &N_melt, double const & Di,
    double const &Phi_obj, double &Phi,
    std::tuple<bool,double,double,double,double,double,double,double,double> const &melt_param,std::tuple<double,double,double,double,double,double> const &solidus, std::tuple<double,double,double,double,double,double> const &liquidus, std::tuple<double,double,double,double,double,double,double,double> const &melt,
    std::tuple<double,double,double,double,double,double,double,double,double,double,double> const &rheology

){

    double Dl02 = Dl0;
    double prior = 1;
    
    // k0 = k0 * std::exp(dk0 * (rand_01() - 0.5) * 2.);
    k0 = std::exp(rand_gauss(std::log(k0),dk0));
    if(k0 > k0_max || k0 < k0_min)
    {
        prior = 0;
        std::cout << "k0 out" << std::endl;
    }
    
    // eta0 = eta0 * std::exp(deta0 * (rand_01() - 0.5) * 2);
    eta0 = std::exp(rand_gauss(std::log(eta0),deta0));
    if( eta0 > eta0_max || eta0 < eta0_min)
    {
        prior = 0;
        std::cout << "eta0 out" << std::endl;
    }

    // fmag = fmag + dfmag * (rand_01()- 0.5) * 2;
    fmag = rand_gauss(fmag,dfmag);
    if( fmag > fmag_max || fmag < fmag_min)
    {
        prior = 0;
        std::cout << "fmag out" << std::endl;
    }

 
    // fbase = fbase + dfbase * (rand_01() - 0.5) * 2;
    fbase = rand_gauss(fbase,dfbase);
    if( fbase > fbase_max || fbase < fbase_min)
    {
        prior = 0;
        std::cout << "fbase out" << std::endl;
    }

    //A =  A + dA * (rand_01() - 0.5) * 2;
    A = rand_gauss(A,dA);
    if( A > A_max || A < A_min)
    {
        prior = 0;
        std::cout << "A out,A = " << A << ",A_min = " << A_min << ", A_max = " << A_max <<", dA = " << dA << std::endl;
    }
 
    // V =  V + dV * (rand_01() - 0.5) * 2;
    V = rand_gauss(V,dV);
    if( V > V_max || V < V_min)
    {
        prior = 0;
        std::cout << "V out" << std::endl;
    }
    
    //rhocr =  rhocr + drhocr * (rand_01() - 0.5) * 2;
    rhocr = rand_gauss(rhocr,drhocr);
    if( rhocr > rhocr_max || rhocr < rhocr_min)
    {
        prior = 0;
        std::cout << "rhocr out, rhocr = " << rhocr << ",rhocr_min = " << rhocr_min << ", rhocr_max = " << rhocr_max <<", drhocr = " << drhocr << std::endl;
    
    }

    // kcr =  kcr + dkcr * (rand_01() - 0.5) * 2;
    kcr = rand_gauss(kcr,dkcr);
    if( kcr > kcr_max || kcr < kcr_min)
    {
        prior = 0;
        std::cout << "kcr out" << std::endl;
    }

    // CH2O = CH2O + dCH2O * (rand_01() - 0.5) * 2; 
    CH2O = rand_gauss(CH2O,dCH2O);
    if( CH2O > CH2O_max || CH2O < CH2O_min)
    {
        prior = 0;
        std::cout << "CH2O out" << std::endl;
    }


    // Tm0 = Tm0 + dTm * (rand_01() - 0.5) * 2;
    Tm0 = rand_gauss(Tm0,dTm);
    if(   Tm0 > Tm0_max ||   Tm0 < Tm0_min  )
    {
        prior = 0;
        std::cout << "Tm0 out" << std::endl;
    }
 
    Tc0 = rand_gauss(Tc0,dTc0);
    if(   Tc0 > Tc0_max ||   Tc0 < Tc0_min  )
    {
        prior = 0;
        std::cout << "dTc0 out" << std::endl;
        std::cout << "dTc0 out,  Tc0 = " <<  Tc0 << ",dTc0 = " << Tc0_min << ", dTc0_max = " << Tc0_max <<", d2Tc0 = " << dTc0 << std::endl;
    }
 //  double &dDl, double const &dDl_max, double const &dDl_min, double const &d2Dl,
    dDl_init = rand_gauss(dDl_init,d2Dl);
    if( dDl_init > dDl_max || dDl_init < dDl_min)
    {
        prior = 0;
        std::cout << "dDl out, dDl = " << dDl_init << ",dDl_min = " << dDl_min << ", dDl_max = " << dDl_max <<", d2Dl = " << d2Dl << std::endl;
    }

    if(prior > 0){
    double Tref =std::get<5>(rheology);
    double Pref =std::get<6>(rheology);
    double R = std::get<2>(rheology);
    double Ra_crit_u = std::get<9>(rheology);
    double beta_u = std::get<7>(rheology);

    double Phi_eff, Va, dmadtm;
    double dBLC0 = 80E3;
    Dl02 = Dl0_min;

    double Ra_u,eta_u,dBLC,Tl,DT_u;
    double Rm = Rp - Dc; 
    double Vm;
    Volume(Vm,Rm,Rm,Rc,f);
    double Vcr=Vp-Vm;
    double rhom = (Vp*rhop-Vcr*rhocr)/Vm;
    double Dm = km / Cm / rhom;
    double LMBD_liq = 1 / (Phi_obj+(1-Phi_obj)*Di);
    double Rl, Pm;
    double Tsol,Tliq, DT_unl,phi,Ps;
    Va = 0;

    if(std::get<0>(melt_param) == 1){
    double frac_unlin_max = std::get<1>(melt_param);
    double x0 = std::get<2>(melt_param);
    double x1 = std::get<3>(melt_param);
    double x2 = std::get<4>(melt_param);
    double x3 = std::get<5>(melt_param);
    double x4 = std::get<6>(melt_param);
    double x5 = std::get<7>(melt_param);
    double x6 = std::get<8>(melt_param); 
   
    litho_temp(Tl,Tm0,rheology);
    DT_u = Tm0 - Tl;
    do
    {
   
    Dl02 = Dl02 + dDl0;
    Rl = Rp - Dl0;
    Pm = g*(rhocr*Dc+rhom*dBLC0); 
 
    eta_u=eta0*std::exp( (A + Pm * V )/ (R  * Tm0 ) - (A + Pref * V) / ( Tref * R) );   
    Ra_u=am*rhom*g* DT_u *pow(Rl-Rc,3.)/Dm/eta_u;  
    
    dBLC= (Rl-Rc) * pow(Ra_crit_u/Ra_u,beta_u);
    Pm = g*(rhocr*Dc+rhom*(Dl02+dBLC));
    Ps = Pm/1E9;
    Tsol=std::get<0>(solidus)+std::get<1>(solidus)*Ps+std::get<2>(solidus)*Ps*Ps+std::get<3>(solidus)*Ps*Ps*Ps;
    Tliq=std::get<0>(liquidus)+std::get<1>(liquidus)*Ps+std::get<2>(liquidus)*Ps*Ps+std::get<3>(liquidus)*Ps*Ps*Ps;
     
    DT_unl=(Tm0-Tsol)/((Tliq-Tsol)*frac_unlin_max);
    phi = std::min(x6*pow(DT_unl,6.) +x5*pow(DT_unl,5.) +x4*pow(DT_unl,4.) +x3*pow(DT_unl,3.) +x2*pow(DT_unl,2.) +x1*DT_unl +x0, (Tm0-Tsol)/(Tliq-Tsol));
    


    } while(phi>Phi_obj && Dl02 < Dl0_max);
    }
    else{
    
    double Mcpx = std::get<1>(melt_param);
    double beta1 = std::get<2>(melt_param);
    double beta2 = std::get<3>(melt_param);
    double r0 = std::get<4>(melt_param);
    double r1 = std::get<5>(melt_param);
    double Rcpx; double Fcpx; double Tcpx; double Tprim;
    double Tliq_lherz;

    litho_temp(Tl,Tm0,rheology);
    DT_u = Tm0 - Tl;
    do
    {
   
    Dl02 = Dl02 + dDl0;
    Rl = Rp - Dl0;
    Pm = g*(rhocr*Dc+rhom*dBLC0); 

    eta_u=eta0*std::exp( (A + Pm * V )/ (R  * Tm0 ) - (A + Pref * V) / ( Tref * R) );   
    Ra_u=am*rhom*g* DT_u *pow(Rl-Rc,3.)/Dm/eta_u;  
    
    dBLC= (Rl-Rc) * pow(Ra_crit_u/Ra_u,beta_u);
    Pm = g*(rhocr*Dc+rhom*(dBLC+Dl02));
    Ps = Pm/1E9;
    Tsol=std::get<0>(solidus)+std::get<1>(solidus)*Ps+std::get<2>(solidus)*Ps*Ps+std::get<3>(solidus)*Ps*Ps*Ps;
    Tliq=std::get<0>(liquidus)+std::get<1>(liquidus)*Ps+std::get<2>(liquidus)*Ps*Ps+std::get<3>(liquidus)*Ps*Ps*Ps;
    
    Rcpx = r0 * r1 * Pm;
    Fcpx = Mcpx / Rcpx ; 
    Tliq_lherz = Tliq - (Tliq-Tsol)/2;
    Tcpx = pow(Fcpx,beta1)*(Tliq_lherz - Tsol) + Tsol;
    Tprim = (Tm0 - Tsol) / (Tliq_lherz - Tsol);
    phi = pow(Tprim,beta1);
    if(phi > Fcpx){
        phi = Fcpx + (1-Fcpx) * pow((Tm0 -Tcpx)/(Tliq-Tcpx),beta2);
    }  


    } while(phi>Phi_obj && Dl02 < Dl0_max);

    }
   

    melting(Phi,Phi_eff,Va,dmadtm,Tm0,am,Cm,g,g,N_melt,Dref,std::get<4>(melt),std::get<5>(melt),
    rhocr,rhom,Dc,Dc,Rl,Rp,dBLC,std::get<3>(melt),std::get<6>(melt),
    melt_param, LMBD_liq,
    solidus,liquidus,rheology);

    Dl0 = Dl02 ;
    }

    return prior;
}

bool lecture_random(std::string const &adresse,double &Phi_obj, int &N_MCMC, bool &MCMC_restart,
    double &k0_max, double &k0_min, double &dk0,
    double &eta0_max, double &eta0_min, double &deta0,
    double &fmag_max, double  &fmag_min, double &dfmag,
    double &fbase_max, double  &fbase_min, double &dfbase,
    double &Tc0_max, double  &Tc0_min, double &dTc0,
    double &A_max, double  &A_min, double &dA,
    double &V_max, double &V_min, double &dV,
    double &rhocr_max, double &rhocr_min, double &drhocr,
    double &kcr_max, double &kcr_min, double &dkcr,
    double &dDl_max, double &dDl_min, double &d2Dl,
    double &CH2O_max, double &CH2O_min, double &dCH2O,
    double &Tm0_max, double &Tm0_min, double &dTm,
    double &Dl0_max, double &Dl0_min, double &dDl0
    

){
    std::ifstream fichier(adresse);
    bool test = true;

    if(fichier)
    {
        std::string ligne_txt;
        std::string text_ign;
        getline(fichier,ligne_txt);
        fichier.ignore();
        fichier >> text_ign >> k0_max >> k0_min >> dk0 >> 
        text_ign >> eta0_max >> eta0_min >> deta0 >> 
        text_ign >> fmag_max >> fmag_min >> dfmag >> 
        text_ign >> fbase_max >> fbase_min >> dfbase >> 
        text_ign >> Tc0_max >> Tc0_min >> dTc0 >> 
        text_ign >> A_max >> A_min >> dA >> 
        text_ign >> V_max >> V_min >> dV >> 
        text_ign >> rhocr_max >> rhocr_min >> drhocr >> 
        text_ign >> kcr_max >> kcr_min >> dkcr >> 
        text_ign >> dDl_max >> dDl_min >> d2Dl >> 
        text_ign >> CH2O_max >> CH2O_min >> dCH2O >> 
        text_ign >> Tm0_max >> Tm0_min >> dTm >> 
        text_ign >> Dl0_max >> Dl0_min >> dDl0 >>
        text_ign >> N_MCMC >>
        text_ign >> Phi_obj >>
        text_ign >> MCMC_restart;
        
        fichier.close();

    }

    else
    {
        std::cout << "ERREUR: Impossible d'ouvrir le fichier lecture_random en lecture." << std::endl;
        test = false;
    }

return test;
}

double rand_01(){

    double RANDMAX = 100000;
    std::random_device rd;
    std::mt19937 re(rd());
    std::uniform_int_distribution<int> distrib{0, int(RANDMAX)};

    double nombre = distrib(re)/RANDMAX;


    return nombre;

}

double rand_gauss(double const &mean, double const &sigma){

    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<> distrib{mean, sigma};

    double nombre = distrib(gen);


    return nombre;

}



bool lecture_MPI_init(std::string const &adresse,int const &rank,
double const &Dc,double const &Rp, double const &Rc, double const &f, double const &Vp, double const &rhop, double const &Dref,
double const &g, double const &Cm, double const &am, int const &N_melt, double const & Di, double const &km, double const &Dl0_min, double const &Dl0_max, double const &dDl0,
double const &Phi_obj, double &Phi, 
std::tuple<bool,double,double,double,double,double,double,double,double> const &melt_param,std::tuple<double,double,double,double,double,double> const &solidus, std::tuple<double,double,double,double,double,double> const &liquidus, std::tuple<double,double,double,double,double,double,double,double> const &melt,
std::tuple<double,double,double,double,double,double,double,double,double,double,double> const &rheology,
double &eta0, double &Tm0, double &k0, double &rhocr, double &kcr, double &A,double &V,double &fmag,double &fbase,double &CH2O,double &Dl0, double &dDl_init, double &dTc0){

    std::ifstream fichier(adresse);
    bool test = true;
    int cpt = 0;

    if(fichier)
    {   
        // fmag, fbase, CH2O_p,dDl_init, dTc0, Dl_init,
        std::string ligne_txt;
        while(fichier){ 
        if(cpt == rank){
            fichier >> eta0 >> Tm0 >> k0 >> rhocr >> kcr >> A >> V >> fmag >> fbase >> CH2O >> dDl_init >> dTc0; 
        }
        getline(fichier,ligne_txt);
        cpt++;
        }

        fichier.close();
    
    // CALCUL Dl0
    double Dl02;
    double Tref =std::get<5>(rheology);
    double Pref =std::get<6>(rheology);
    double R = std::get<2>(rheology);
    double Ra_crit_u = std::get<9>(rheology);
    double beta_u = std::get<7>(rheology);

    double Phi_eff, Va, dmadtm;
    double dBLC0 = 80E3;
    Dl02 = Dl0_min;

    double Ra_u,eta_u,dBLC,Tl,DT_u;
    double Rm = Rp - Dc; 
    double Vm;
    Volume(Vm,Rm,Rm,Rc,f);
    double Vcr=Vp-Vm;
    double rhom = (Vp*rhop-Vcr*rhocr)/Vm;
    double Dm = km / Cm / rhom;
    double LMBD_liq = 1 / (Phi_obj+(1-Phi_obj)*Di);
    double Rl, Pm;
    double Tsol,Tliq, DT_unl,phi,Ps;
    Va = 0;

    if(std::get<0>(melt_param) == 1){
    double frac_unlin_max = std::get<1>(melt_param);
    double x0 = std::get<2>(melt_param);
    double x1 = std::get<3>(melt_param);
    double x2 = std::get<4>(melt_param);
    double x3 = std::get<5>(melt_param);
    double x4 = std::get<6>(melt_param);
    double x5 = std::get<7>(melt_param);
    double x6 = std::get<8>(melt_param); 
   
    litho_temp(Tl,Tm0,rheology);
    DT_u = Tm0 - Tl;
    do
    {
   
    Dl02 = Dl02 + dDl0;
    Rl = Rp - Dl0;
    Pm = g*(rhocr*Dc+rhom*dBLC0); 
 
    eta_u=eta0*std::exp( (A + Pm * V )/ (R  * Tm0 ) - (A + Pref * V) / ( Tref * R) );   
    Ra_u=am*rhom*g* DT_u *pow(Rl-Rc,3.)/Dm/eta_u;  
    
    dBLC= (Rl-Rc) * pow(Ra_crit_u/Ra_u,beta_u);
    Pm = g*(rhocr*Dc+rhom*(Dl02+dBLC));
    Ps = Pm/1E9;
    Tsol=std::get<0>(solidus)+std::get<1>(solidus)*Ps+std::get<2>(solidus)*Ps*Ps+std::get<3>(solidus)*Ps*Ps*Ps;
    Tliq=std::get<0>(liquidus)+std::get<1>(liquidus)*Ps+std::get<2>(liquidus)*Ps*Ps+std::get<3>(liquidus)*Ps*Ps*Ps;
     
    DT_unl=(Tm0-Tsol)/((Tliq-Tsol)*frac_unlin_max);
    phi = std::min(x6*pow(DT_unl,6.) +x5*pow(DT_unl,5.) +x4*pow(DT_unl,4.) +x3*pow(DT_unl,3.) +x2*pow(DT_unl,2.) +x1*DT_unl +x0, (Tm0-Tsol)/(Tliq-Tsol));
    


    } while(phi>Phi_obj && Dl02 < Dl0_max);
    }

    else{
    
    double Mcpx = std::get<1>(melt_param);
    double beta1 = std::get<2>(melt_param);
    double beta2 = std::get<3>(melt_param);
    double r0 = std::get<4>(melt_param);
    double r1 = std::get<5>(melt_param);
    double Rcpx; double Fcpx; double Tcpx; double Tprim;
    double Tliq_lherz;

    litho_temp(Tl,Tm0,rheology);
    DT_u = Tm0 - Tl;
    do
    {
   
    Dl02 = Dl02 + dDl0;
    Rl = Rp - Dl0;
    Pm = g*(rhocr*Dc+rhom*dBLC0); 

    eta_u=eta0*std::exp( (A + Pm * V )/ (R  * Tm0 ) - (A + Pref * V) / ( Tref * R) );   
    Ra_u=am*rhom*g* DT_u *pow(Rl-Rc,3.)/Dm/eta_u;  
    
    dBLC= (Rl-Rc) * pow(Ra_crit_u/Ra_u,beta_u);
    Pm = g*(rhocr*Dc+rhom*(dBLC+Dl02));
    Ps = Pm/1E9;
    Tsol=std::get<0>(solidus)+std::get<1>(solidus)*Ps+std::get<2>(solidus)*Ps*Ps+std::get<3>(solidus)*Ps*Ps*Ps;
    Tliq=std::get<0>(liquidus)+std::get<1>(liquidus)*Ps+std::get<2>(liquidus)*Ps*Ps+std::get<3>(liquidus)*Ps*Ps*Ps;
    
    Rcpx = r0 * r1 * Pm;
    Fcpx = Mcpx / Rcpx ; 
    Tliq_lherz = Tliq - (Tliq-Tsol)/2;
    Tcpx = pow(Fcpx,beta1)*(Tliq_lherz - Tsol) + Tsol;
    Tprim = (Tm0 - Tsol) / (Tliq_lherz - Tsol);
    phi = pow(Tprim,beta1);
    if(phi > Fcpx){
        phi = Fcpx + (1-Fcpx) * pow((Tm0 -Tcpx)/(Tliq-Tcpx),beta2);
    }  


    } while(phi>Phi_obj && Dl02 < Dl0_max);

    }

    melting(Phi,Phi_eff,Va,dmadtm,Tm0,am,Cm,g,g,N_melt,Dref,std::get<4>(melt),std::get<5>(melt),
    rhocr,rhom,Dc,Dc,Rl,Rp,dBLC,std::get<3>(melt),std::get<6>(melt),
    melt_param, LMBD_liq,
    solidus,liquidus,rheology);

    Dl0 = Dl02 ;

    }

    else
    {
        std::cout<< rank << " : rank, " << "ERREUR: Impossible d'ouvrir le fichier MCMC_init en lecture." << std::endl;
        test = false;
    }

return test;
}

bool lecture_MPI_restart(std::string const &adresse, int const &rank, int &n_start,
double const &Dc,double const &Rp, double const &Rc, double const &f, double const &Vp, double const &rhop, double const &Dref,
double const &g, double const &Cm, double const &am, int const &N_melt, double const & Di, double const &km, double const &Dl0_min, double const &Dl0_max, double const &dDl0,
double const &Phi_obj, double &Phi, 
std::tuple<bool,double,double,double,double,double,double,double,double> const &melt_param,std::tuple<double,double,double,double,double,double> const &solidus, std::tuple<double,double,double,double,double,double> const &liquidus, std::tuple<double,double,double,double,double,double,double,double> const &melt,
std::tuple<double,double,double,double,double,double,double,double,double,double,double> const &rheology,
double &eta0, double &Tm0, double &k0, double &rhocr, double &kcr, double &A,double &V,double &fmag,double &fbase,double &CH2O,double &Dl0,  double &dDl_init, double &dTc0){

    std::cout << "test 2" << std::endl;
    std::ifstream fichier(adresse);
    bool test = true;
    std::cout << "test 3" << std::endl;
    n_start = 0;
    if(fichier)
    {   
        std::string ligne_txt;
        while(fichier){
        n_start++;      
        getline(fichier,ligne_txt);
        fichier >> eta0 >> Tm0 >> k0 >> rhocr >> kcr >> A >> V >> fmag >> fbase >> CH2O >> dDl_init >> dTc0;
        }
        std::cout<< rank << " : rank, eta0 : " << eta0 << ", Tm0 : " << Tm0 << ", k0 : " << k0 << ", rhocr : " << rhocr << ", kcr : " << kcr << ", A : " << A << ", V : " << V << ", fmag : " << fmag << ", fbase : " << fbase << ", CH2O : " << CH2O << std::endl;
        fichier.close();
    
    // CALCUL Dl0
    double Dl02;
    double Tref =std::get<5>(rheology);
    double Pref =std::get<6>(rheology);
    double R = std::get<2>(rheology);
    double Ra_crit_u = std::get<9>(rheology);
    double beta_u = std::get<7>(rheology);

    double Phi_eff, Va, dmadtm;
    double dBLC0 = 80E3;
    Dl02 = Dl0_min;

    double Ra_u,eta_u,dBLC,Tl,DT_u;
    double Rm = Rp - Dc; 
    double Vm;
    Volume(Vm,Rm,Rm,Rc,f);
    double Vcr=Vp-Vm;
    double rhom = (Vp*rhop-Vcr*rhocr)/Vm;
    double Dm = km / Cm / rhom;
    double LMBD_liq = 1 / (Phi_obj+(1-Phi_obj)*Di);
    double Rl, Pm;
    double Tsol,Tliq, DT_unl,phi,Ps;
    Va = 0;

    if(std::get<0>(melt_param) == 1){
    double frac_unlin_max = std::get<1>(melt_param);
    double x0 = std::get<2>(melt_param);
    double x1 = std::get<3>(melt_param);
    double x2 = std::get<4>(melt_param);
    double x3 = std::get<5>(melt_param);
    double x4 = std::get<6>(melt_param);
    double x5 = std::get<7>(melt_param);
    double x6 = std::get<8>(melt_param); 
   
    litho_temp(Tl,Tm0,rheology);
    DT_u = Tm0 - Tl;
    do
    {
   
    Dl02 = Dl02 + dDl0;
    Rl = Rp - Dl0;
    Pm = g*(rhocr*Dc+rhom*dBLC0); 
 
    eta_u=eta0*std::exp( (A + Pm * V )/ (R  * Tm0 ) - (A + Pref * V) / ( Tref * R) );   
    Ra_u=am*rhom*g* DT_u *pow(Rl-Rc,3.)/Dm/eta_u;  
    
    dBLC= (Rl-Rc) * pow(Ra_crit_u/Ra_u,beta_u);
    Pm = g*(rhocr*Dc+rhom*(Dl02+dBLC));
    Ps = Pm/1E9;
    Tsol=std::get<0>(solidus)+std::get<1>(solidus)*Ps+std::get<2>(solidus)*Ps*Ps+std::get<3>(solidus)*Ps*Ps*Ps;
    Tliq=std::get<0>(liquidus)+std::get<1>(liquidus)*Ps+std::get<2>(liquidus)*Ps*Ps+std::get<3>(liquidus)*Ps*Ps*Ps;
     
    DT_unl=(Tm0-Tsol)/((Tliq-Tsol)*frac_unlin_max);
    phi = std::min(x6*pow(DT_unl,6.) +x5*pow(DT_unl,5.) +x4*pow(DT_unl,4.) +x3*pow(DT_unl,3.) +x2*pow(DT_unl,2.) +x1*DT_unl +x0, (Tm0-Tsol)/(Tliq-Tsol));
    


    } while(phi>Phi_obj && Dl02 < Dl0_max);
    }
    else{
    
    double Mcpx = std::get<1>(melt_param);
    double beta1 = std::get<2>(melt_param);
    double beta2 = std::get<3>(melt_param);
    double r0 = std::get<4>(melt_param);
    double r1 = std::get<5>(melt_param);
    double Rcpx; double Fcpx; double Tcpx; double Tprim;
    double Tliq_lherz;

    litho_temp(Tl,Tm0,rheology);
    DT_u = Tm0 - Tl;
    do
    {
   
    Dl02 = Dl02 + dDl0;
    Rl = Rp - Dl0;
    Pm = g*(rhocr*Dc+rhom*dBLC0); 

    eta_u=eta0*std::exp( (A + Pm * V )/ (R  * Tm0 ) - (A + Pref * V) / ( Tref * R) );   
    Ra_u=am*rhom*g* DT_u *pow(Rl-Rc,3.)/Dm/eta_u;  
    
    dBLC= (Rl-Rc) * pow(Ra_crit_u/Ra_u,beta_u);
    Pm = g*(rhocr*Dc+rhom*(dBLC+Dl02));
    Ps = Pm/1E9;
    Tsol=std::get<0>(solidus)+std::get<1>(solidus)*Ps+std::get<2>(solidus)*Ps*Ps+std::get<3>(solidus)*Ps*Ps*Ps;
    Tliq=std::get<0>(liquidus)+std::get<1>(liquidus)*Ps+std::get<2>(liquidus)*Ps*Ps+std::get<3>(liquidus)*Ps*Ps*Ps;
    
    Rcpx = r0 * r1 * Pm;
    Fcpx = Mcpx / Rcpx ; 
    Tliq_lherz = Tliq - (Tliq-Tsol)/2;
    Tcpx = pow(Fcpx,beta1)*(Tliq_lherz - Tsol) + Tsol;
    Tprim = (Tm0 - Tsol) / (Tliq_lherz - Tsol);
    phi = pow(Tprim,beta1);
    if(phi > Fcpx){
        phi = Fcpx + (1-Fcpx) * pow((Tm0 -Tcpx)/(Tliq-Tcpx),beta2);
    }  


    } while(phi>Phi_obj && Dl02 < Dl0_max);

    }

    melting(Phi,Phi_eff,Va,dmadtm,Tm0,am,Cm,g,g,N_melt,Dref,std::get<4>(melt),std::get<5>(melt),
    rhocr,rhom,Dc,Dc,Rl,Rp,dBLC,std::get<3>(melt),std::get<6>(melt),
    melt_param, LMBD_liq,
    solidus,liquidus,rheology);

    Dl0 = Dl02 ;

    }

    else
    {
        std::cout << rank << " : rank, " << "ERREUR: Impossible d'ouvrir le fichier MCMC_resume en lecture." << std::endl;
        test = false;
    }

return test;
}

    





    



