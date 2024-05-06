#include<iostream>
#include<vector>
#include<array>
#include<cmath>
#include<string>
#include<stdio.h>
#include<time.h>

#include "Eigen/Dense"
#include "1PIECE_lecture.cpp"
#include "1PIECE_HEATPROD.cpp"
#include "1PIECE_Differential_eq.cpp"
#include "1PIECE_Conduction.cpp"
#include "1PIECE_ecriture.cpp"


void Internal_energy(Eigen::VectorXd const &T_N, Eigen::VectorXd const &T_S, Eigen::VectorXd const &PhiN,Eigen::VectorXd const &PhiS,Eigen::VectorXd const &dVN_r,Eigen::VectorXd const &dVS_r,Eigen::VectorXd const &PN_r,Eigen::VectorXd const &PS_r, int const &Nm_N, int const &Nm_S,
double const &rho_c, double const &rho_m, double const &rho_cr, double const &Cc, double const &Ccr,double const &Cm, double const &f,double const &L_m, double const &L_cr, double const &epsi_m,double const &epsi_c,double const &Path,
double const &Tm,double const &Tc,double const &qsurf_N, double const &qsurf_S, double const &phi_avg, double const &Va, double const &dt, double const &Vc, double const &Vcm, double const &Ap, double const &Vp,
double &dU_tot, double &Qs, double &U_old, double &H_tot, double &T_avg);

// Function given one thermal evolution

// thermo is a table with the thermodynics values
// rheology is a table with rheology values
// solidus - liquidus are tables with the 4 factor of a degree 3 polynom for the solidus - liquidus curves
// RAD is a matrix which containts the heat producing constant

// _cr relative to the crust, _m relative to the all mantle, _ath to the convective mantle, _lith to the lithoshperic mantle, _c relative to the core 
// N/S => North / South  hemisphere

bool Parametric_model(std::tuple<double,double,double,double,double,double,double,double,double,double,double,double,double> &thermo, std::tuple<double,double,double,double,double,double,double,double,double,double,double> &rheology, std::tuple<double,double,double,double,double,double,double,double> &melt, std::tuple<bool,double,double,double,double,double,double,double,double> const &melt_param,
Eigen::MatrixXd RAD, std::tuple<double,double,double,double,double,double> const &solidus, std::tuple<double,double,double,double,double,double> const &liquidus,
double const &Rp, double const &Rc, double const &f, double const &Vc, double const &Ac, double const &Vp, double const &rho_p, double const &gc, double const &gl, double const &dr_min, double const dr_max,double const &Dref, double const &dt_max, double const &dt_min, double const &dt_init,
double const &Tm0, double const &dTc0, double const &Ts, double const &tstart, double const &tstop, double const &t_acc, double const &t_acc_m, double const &DTMAX, double const &fmag, double const &fbase,
double const &Dl_init, double const &Dc_init, double const &dDl_init, double const &dDc_init, double const &LMBD_cr_i, double const &Phi_eff_i, double const &CH2O_p, double const &Sat_expo, double const &X1, double const &X2, double const &K_cr, double const &gamma_cr, double const &phi_min_cr,
bool const  &RK4, bool const &steady, bool const &LMBDdt_bol, int const &N_melt, int const &Nc,
double &t, double &Tm, double &Tc,double &Tp, double &Phicrmax_N, double &Phicrmax_S, double &Dl_N,double &Dl_S,double &Dc_N,double &Dc_S, double &dBLC_N, double &dBLC_S,std::tuple<double,double,double> &agecrust_N,std::tuple<double,double,double> &agecrust_S,
double &LMBD_ath,double &LMBD_lith_N,double  &LMBD_lith_S,double &LMBD_cr_N,double  &LMBD_cr_S,double  &dmax_cr_N, double  &dmax_cr_S ,std::tuple<double,double>  &dTdz_N,std::tuple<double,double>  &dTdz_S,
char *dossier_time, char *dossier_tech, char* dossier_profilN, char * dossier_profilS, bool const &ecrit_profil_b, bool const & ecrit_time_b,bool const & ecrit_tech_b, bool const &URG_STOP,bool const &adv_lith_bool, bool const &adv_interf_bool
){

/*
std::cout << "DEBUGG - DEBUT" << std::endl;
std::cout << "THERMO : " << std::get<0>(thermo) << ", " << std::get<1>(thermo) << ", " << std::get<2>(thermo) << ", " << std::get<3>(thermo) << ", " << std::get<4>(thermo) << ", " << std::get<5>(thermo) << ", " << std::get<6>(thermo) << ", " << std::get<7>(thermo) << ", " << std::get<8>(thermo) << ", " << std::get<9>(thermo) << ", " << std::get<10>(thermo) << ", " << std::get<11>(thermo)  << std::endl;
std::cout << "RHEOLOGY : " << std::get<0>(rheology) << ", " << std::get<1>(rheology) << ", " << std::get<2>(rheology) << ", " << std::get<3>(rheology) << ", " << std::get<4>(rheology) << ", " << std::get<5>(rheology) << ", " << std::get<6>(rheology) << ", " << std::get<7>(rheology) << ", " << std::get<8>(rheology) << ", " << std::get<9>(rheology) << ", " << std::get<10>(rheology)  << std::endl;
std::cout << "MELT : " << std::get<0>(melt) << ", " << std::get<1>(melt) << ", " << std::get<2>(melt) << ", " << std::get<3>(melt) << ", " << std::get<4>(melt) << ", " << std::get<5>(melt) << ", " << std::get<6>(melt) << ", " << std::get<7>(melt)   << std::endl;
std::cout << "AUTRE : " << rho_p << ", " << Tm0 << ", " << dTc0 << ", " << Ts << ", " << DTMAX << ", " << fmag << ", " << fbase << ", " << Dl_init << ", " <<  Dc_init << ", " << dDl_init << ", " << dDc_init << ", " << LMBD_cr_i << ", " << CH2O_p << std::endl;
std::cout << "DEBUGG - FIN" << std::endl;
*/

if(ecrit_time_b == 1 || ecrit_profil_b == 1) {ecriture_initial(dossier_time,dossier_tech,dossier_profilN,dossier_profilS,ecrit_time_b,ecrit_tech_b,ecrit_profil_b);}

// The differents thermodynamics values, C heat capacity
double C_m=std::get<0>(thermo);
double C_cr=std::get<1>(thermo);
double C_c=std::get<2>(thermo);
double rho_cr=std::get<3>(thermo);
double rho_c=std::get<4>(thermo);
double k_m=std::get<5>(thermo);
double k_cr=std::get<6>(thermo);
double a_m=std::get<7>(thermo);
double L_m =std::get<8>(thermo);
double Di =std::get<9>(thermo);
double epsi_c=std::get<10>(thermo);
double an_s=std::get<11>(thermo);
double L_cr =std::get<12>(thermo);

double Ap = 4.*M_PI*pow(Rp,2.);

double dtmdt = 0, dtmdt1 = 0, dtmdt2 = 0, dtmdt3 = 0, dtmdt4 = 0; 
double dtcdt = 0, dtcdt1 = 0, dtcdt2 = 0, dtcdt3 = 0, dtcdt4 = 0; 

double dDcdt_N = 0, dDcdt_N1 = 0, dDcdt_N2 = 0, dDcdt_N3 = 0, dDcdt_N4 = 0; 
double dDcdt_S = 0, dDcdt_S1 = 0, dDcdt_S2 = 0, dDcdt_S3 = 0, dDcdt_S4 = 0; 
 
double dDldt_N = 0, dDldt_N1 = 0, dDldt_N2 = 0, dDldt_N3 = 0, dDldt_N4 = 0; 
double dDldt_S = 0, dDldt_S1 = 0, dDldt_S2 = 0, dDldt_S3 = 0, dDldt_S4 = 0;

Dl_N = Dl_init + dDl_init;     // Initial lithosphere thickness, the perturbation is add in the north
Dl_S = Dl_init;


std::get<0>(agecrust_N) = 0;
std::get<0>(agecrust_S) = 0;

std::get<1>(agecrust_N) = 0;
std::get<1>(agecrust_S) = 0;

std::get<2>(agecrust_N) = 0;
std::get<2>(agecrust_S) = 0;

double Rl_N2, Rl_S2, Rm_N2, Rm_S2;
double Rl_N = Rp-Dl_N;
double Rl_S = Rp-Dl_S;
double Rl= pow((f*pow(Rl_N,3)+(1-f)*pow(Rl_S,3)),1./3.);

Dc_N = Dc_init + dDc_init;   // Initial crust thickness, the perturbation is add in the north   => a thinner crust in the north need a negative perturbation then. 
Dc_S = Dc_init;

double Rm_N = Rp-Dc_N;
double Rm_S = Rp-Dc_S;
double Rm = pow ((f*pow(Rm_N,3)+(1-f)*pow(Rm_S,3)),1./3.);  // Volume average crust radius and associate thickness
double Dc= Rp-Rm;


// Corresponding Volume

double Vm = 4./3.*M_PI*pow(Rm,3.) - Vc;
double Vcr = Vp-Vm;
double Vcr_N = f*4./3.*M_PI*(pow(Rp,3.)-pow(Rm_N,3.));
double Vcr_S = (1-f)*4./3.*M_PI*(pow(Rp,3.)-pow(Rm_S,3.));
double Acr_N = f*4.0*M_PI*pow(Rm_N,2.0);
double Acr_S = (1-f)*4.0*M_PI*pow(Rm_S,2.0);

//Average mantle density calculation by Mass balance
double rho_m = (Vp*rho_p - Vcr*rho_cr)/Vm;
// 

// LMBD = enrichment factor

LMBD_cr_N = LMBD_cr_i;
LMBD_cr_S = LMBD_cr_i;
LMBD_ath=(Vp*rho_p-rho_cr*(LMBD_cr_N*Vcr_N+LMBD_cr_S*Vcr_S))/(Vm*rho_m);
LMBD_lith_N=LMBD_ath;
LMBD_lith_S=LMBD_ath;

double LMBD_liq_N = LMBD_cr_i;
double LMBD_liq_S = LMBD_cr_i;

Tm = Tm0;
double Tl;
double Tb;

// Initial boudnary layer thickness fix to calculte Tb and Tc
double dBLH = 5E4;
double dBLC_avg = 1.5E5;

dBLC_N = dBLC_avg;
dBLC_S = dBLC_avg;

//Initial Tc calculate with an adiabatic gradient from Tm to the CMB with a overheatin teperature dTc0
Base_temp(Tc,Tm,a_m,C_m,gc,Rl,Rc,dBLC_avg,0);
Tc = Tc + dTc0;

// Pressure at the boundary layer
double Pm_N = gl*rho_cr*Dc_N + gl*rho_m*(Rp-Rl_N-dBLC_N-Dc);
double Pm_S = gl*rho_cr*Dc_S + gl*rho_m*(Rp-Rl_S-dBLC_S-Dc);
double Pb = gl*rho_cr*Dc+ gc*rho_m*(Rp-Rc+dBLH-Dc);

// Tl, Tb initial calculation
litho_temp(Tl,Tm,rheology);
Base_temp(Tb,Tm,a_m,C_m,gc,Rl,Rc,dBLH,dBLC_avg);


double Ra_avg = 1E9;
auto Pm = std::make_tuple(Pm_N,Pm_S,Pb);   // tuple creation to stock the pressures

std::tuple<int,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,bool> THERMOGEO_N;
std::tuple<int,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,bool> THERMOGEO_S;

// Calculation of the number of point for each lithosphere
int N_N = (Rp-Rl_N)/dr_min;
int N_S = (Rp-Rl_S)/dr_min;
int Nm_N = 0;
int Nm_S = 0;


// Vector of temperature, and melt fraction
Eigen::VectorXd T_N = Eigen::VectorXd::Ones(N_N);
Eigen::VectorXd T_S = Eigen::VectorXd::Ones(N_S);
Eigen::VectorXd rPhi_N = Eigen::VectorXd::Zero(N_N);
Eigen::VectorXd rPhi_S = Eigen::VectorXd::Zero(N_S);
Eigen::VectorXd T1_N = Eigen::VectorXd::Ones(N_N);
Eigen::VectorXd T1_S = Eigen::VectorXd::Ones(N_S);

Eigen::VectorXd Tsol_N = Eigen::VectorXd::Ones(N_N);
Eigen::VectorXd Tsol_S = Eigen::VectorXd::Ones(N_S);

Eigen::VectorXd Tliq_N = Eigen::VectorXd::Ones(N_N);
Eigen::VectorXd Tliq_S = Eigen::VectorXd::Ones(N_S);

Eigen::VectorXd DTsolcr_N = Eigen::VectorXd::Ones(N_N);
Eigen::VectorXd DTsolcr_S = Eigen::VectorXd::Ones(N_S);

double Hb_N = 0, Vmelt_N =0, Tcr_N =0, HaT_N = 0;
double Hb_S = 0, Vmelt_S =0, Tcr_S =0, HaT_S = 0;

double P_cr_N; 
double P_lith_N;

double P_cr_S;
double P_lith_S;

double P_ath;

// Initial heat production in the lithosphere/crust
Internal_heat(P_cr_N,0,LMBD_cr_N,rho_cr,an_s,RAD);
Internal_heat(P_cr_S,0,LMBD_cr_S,rho_cr,an_s,RAD);
Internal_heat(P_lith_N,0,LMBD_lith_N,rho_m,an_s,RAD);
Internal_heat(P_lith_S,0,LMBD_lith_N,rho_m,an_s,RAD);

double ql_N;
double ql_S;

double qsurf_N;
double qsurf_S;

double epsi_m;
double eta_u_N, eta_u_S;
double Phi_avg;
double Phi_N;
double Phi_S;
double Phi_eff_N = Phi_eff_i;
double Phi_eff_S = Phi_eff_i;
double Phi_cr_avg_N;
double Phi_cr_avg_S;
double Phi_lith_avg_N;
double Phi_lith_avg_S;
double Phi_cr_max_N;
double Phi_cr_max_S;
double Vcr_melt_N;
double Vcr_melt_S;
double Vlith_melt_N;
double Vlith_melt_S;
double qmoho_N;
double qmoho_S;
double Va_avg;
double dmadtm;
double qm_N;
double qm_S;
double qcr_N;
double qcr_S;
double qc;
double Vath;
Phicrmax_N = 0;
Phicrmax_S = 0;
dmax_cr_N = 0;
dmax_cr_S = 0;
double dU =0,Qs =0, P_prim =0, H_tot = 0, St=0, U_old =1E30;
double T_avg;




// Initial geometry of the ithosphere/crust 
geothermo(THERMOGEO_N,N_N,Rp,Rl_N,Rm_N,k_cr,k_m,C_cr,C_m,rho_cr,rho_m,dDcdt_N,P_cr_N,P_lith_N,Hb_N,L_m,L_cr,fbase,fmag,gl,Tl,adv_lith_bool);
geothermo(THERMOGEO_S,N_S,Rp,Rl_S,Rm_S,k_cr,k_m,C_cr,C_m,rho_cr,rho_m,dDcdt_S,P_cr_S,P_lith_S,Hb_S,L_m,L_cr,fbase,fmag,gl,Tl,adv_lith_bool);

// The initial temperature profil is obtain with a steady state - analytical solution
analytique(T_N,std::get<1>(THERMOGEO_N),std::get<2>(THERMOGEO_N),ql_N,qsurf_N,N_N,std::get<0>(THERMOGEO_N),Rl_N,Rm_N,Rp,Ts,Tl,P_cr_N,P_lith_N,k_m,k_cr,rho_cr,C_cr,rho_m,C_m);
analytique(T_S,std::get<1>(THERMOGEO_S),std::get<2>(THERMOGEO_N),ql_S,qsurf_S,N_S,std::get<0>(THERMOGEO_S),Rl_S,Rm_S,Rp,Ts,Tl,P_cr_S,P_lith_S,k_m,k_cr,rho_cr,C_cr,rho_m,C_m);


T1_N=T_N;
T1_S=T_S;


Eigen::VectorXd R1_N = std::get<1>(THERMOGEO_N);
Eigen::VectorXd R1_S = std::get<1>(THERMOGEO_S);
// Temperature at the base of the crust
Tcr_N = T_N(std::get<0>(THERMOGEO_N));
Tcr_S = T_S(std::get<0>(THERMOGEO_S));
// Contact True = crust thickness = lithoshere thickness
std::tuple<bool,bool> contact;
std::get<0>(contact) = 0;
std::get<1>(contact) = 0;


bool stop = 0;
bool melt_bol = 1 ;
double dt = 0;
double dt2 = 0;
double dt_old = dt;
double dt_conduct_N = 0;
double dt_conduct_S = 0;
double t_wmelt = 0;

bool interpol_bool =  1;

// clock_t t1,t2;
// float temps;


// Time loop


while( t < (tstop-tstart)*an_s && stop == 0 && Ra_avg > std::get<9>(rheology)*10.0)
{
    
 
    // t1 = clock();
    
    // Heat production at the time t
    Internal_heat(P_cr_N,t,LMBD_cr_N,rho_cr,an_s,RAD);
    Internal_heat(P_cr_S,t,LMBD_cr_S,rho_cr,an_s,RAD);
    Internal_heat(P_lith_N,t,LMBD_lith_N,rho_m,an_s,RAD);
    Internal_heat(P_lith_S,t,LMBD_lith_N,rho_m,an_s,RAD);
    Internal_heat(P_prim,t,1,rho_p,an_s,RAD);
    Internal_heat(P_ath,t,LMBD_ath,rho_m,an_s,RAD);

    
    
    if(RK4 == 1){   // RK4 schema

    //derivation
    
    derivate(dtmdt1,dtcdt1,dDcdt_N1,dDcdt_S1,dDldt_N1,dDldt_S1,
    P_ath,Tm,Tc,f,Rl_N,Rl_S,Dc_N,Dc_S,
    Rp,Rc,Vc,Ac,Ts,gl,gc,Tb,a_m,k_m,C_m,C_c,C_cr,rho_c,rho_m,
    rho_cr,L_m,LMBD_ath,LMBD_liq_N,LMBD_liq_S,Dref,solidus,liquidus,melt,rheology,RAD,Pm,an_s,ql_N,ql_S,epsi_c,contact,melt_bol,adv_interf_bool,melt_param,N_melt,
    Tl, epsi_m, Ra_avg, dBLH, dBLC_N, dBLC_S,dBLC_avg, eta_u_N,eta_u_S, Phi_avg, Phi_N, Phi_S, Phi_eff_N, Phi_eff_S, Va_avg, dmadtm, qm_N, qm_S, qcr_N, qcr_S, qc,St
);  
    if(t == 0){ dt = dt_init*an_s;} // Initial time step in the parameters
    else{   // time step calculation  
    dt2 =dr_min / (4*std::max({std::abs(dDldt_N1),std::abs(dDldt_S1),std::abs(dDcdt_N1),std::abs(dDcdt_S1)}));  // Crust and lithosphere can't growth more than 1 half finite volume
    if(rPhi_N.maxCoeff()==0){dt_conduct_N = dt_max*an_s;}else{dt_conduct_N = dt_max*an_s/20;}
    if(rPhi_S.maxCoeff()==0){dt_conduct_S = dt_max*an_s;}else{dt_conduct_S = dt_max*an_s/20;}
    dt2=std::min({dt2,dt_conduct_N,dt_conduct_S});
    if(std::abs(dt2*dtmdt1) > DTMAX || std::abs(dt2*dtcdt1) > DTMAX ){dt2=dt_min*an_s;}
    dt2 = std::min(dt*2.0,dt2);} // Max temperature jump in 1 time step DTMAX
    if(dt2/an_s>dt_max){dt=dt_max*an_s;} else if(dt2/an_s < dt_min){dt=dt_min*an_s;} else{dt=dt2;}     // min and max dt compute in parameters to control the calculation time. 
    
    
    derivate(dtmdt2,dtcdt2,dDcdt_N2,dDcdt_S2,dDldt_N2,dDldt_S2,
    t+dt/2,Tm+dtmdt1*dt/2.,Tc+dtcdt1*dt/2.,f,Rl_N-dDldt_N1*dt/2.,Rl_S-dDldt_S1*dt/2.,Dc_N+dDcdt_N1*dt/2.,Dc_S+dDcdt_S1*dt/2.,
    Rp,Rc,Vc,Ac,Ts,gl,gc,Tb,a_m,k_m,C_m,C_c,C_cr,rho_c,rho_m,
    rho_cr,L_m,LMBD_ath,LMBD_liq_N,LMBD_liq_S,Dref,solidus,liquidus,melt,rheology,RAD,Pm,an_s,ql_N,ql_S,epsi_c,contact,melt_bol,adv_interf_bool, melt_param, N_melt, Phi_eff_N, Phi_eff_S);
    
    derivate(dtmdt3,dtcdt3,dDcdt_N3,dDcdt_S3,dDldt_N3,dDldt_S3,
    t+dt/2.,Tm+dtmdt2*dt/2.,Tc+dtcdt2*dt/2.,f,Rl_N-dDldt_N2*dt/2.,Rl_S-dDldt_S2*dt/2.,Dc_N+dDcdt_N2*dt/2.,Dc_S+dDcdt_S2*dt/2.,
    Rp,Rc,Vc,Ac,Ts,gl,gc,Tb,a_m,k_m,C_m,C_c,C_cr,rho_c,rho_m,
    rho_cr,L_m,LMBD_ath,LMBD_liq_N,LMBD_liq_S,Dref,solidus,liquidus,melt,rheology,RAD,Pm,an_s,ql_N,ql_S,epsi_c,contact,melt_bol, adv_interf_bool, melt_param,N_melt, Phi_eff_N, Phi_eff_S);
    
    derivate(dtmdt4,dtcdt4,dDcdt_N4,dDcdt_S4,dDldt_N4,dDldt_S4,
    t+dt,Tm+dtmdt3*dt,Tc+dtcdt3*dt,f,Rl_N-dDldt_N3*dt,Rl_S-dDldt_S3*dt,Dc_N+dDcdt_N3*dt,Dc_S+dDcdt_S3*dt,
    Rp,Rc,Vc,Ac,Ts,gl,gc,Tb,a_m,k_m,C_m,C_c,C_cr,rho_c,rho_m,
    rho_cr,L_m,LMBD_ath,LMBD_liq_N,LMBD_liq_S,Dref,solidus,liquidus,melt,rheology,RAD,Pm,an_s,ql_N,ql_S,epsi_c,contact,melt_bol, adv_interf_bool, melt_param, N_melt, Phi_eff_N, Phi_eff_S);
    
    // std::cout << "dDldt_N 1-4 : "  << dDldt_N1 << ", " << dDldt_N2 << ", "  << dDldt_N3 << ", " << dDldt_N4 << std::endl;
    // std::cout << "dTmdt_N 1-4 : "  << dtmdt1 << ", " << dtmdt2 << ", "  << dtmdt3 << ", " <<dtmdt4 << std::endl;
    dtmdt = (dtmdt1 + 2.* dtmdt2 + 2.* dtmdt3 + dtmdt4)/6.;
    dtcdt = (dtcdt1 + 2.* dtcdt2 + 2.* dtcdt3 + dtcdt4)/6.;
    dDldt_N = (dDldt_N1 + 2.* dDldt_N2 + 2.* dDldt_N3 + dDldt_N4)/6.;
    dDldt_S = (dDldt_S1 + 2.* dDldt_S2 + 2.* dDldt_S3 + dDldt_S4)/6.;
    dDcdt_N = (dDcdt_N1 + 2.* dDcdt_N2 + 2.* dDcdt_N3 + dDcdt_N4)/6.;
    dDcdt_S = (dDcdt_S1 + 2.* dDcdt_S2 + 2.* dDcdt_S3 + dDcdt_S4)/6.;

    


    
    } else { // Implicite scheme
    
    derivate(dtmdt,dtcdt,dDcdt_N,dDcdt_S,dDldt_N,dDldt_S,
    P_ath,Tm,Tc,f,Rl_N,Rl_S,Dc_N,Dc_S,
    Rp,Rc,Vc,Ac,Ts,gl,gc,Tb,a_m,k_m,C_m,C_c,C_cr,rho_c,rho_m,
    rho_cr,L_m,LMBD_ath,LMBD_liq_N,LMBD_liq_S,Dref,solidus,liquidus,melt,rheology,RAD,Pm,an_s,ql_N,ql_S,epsi_c,contact,melt_bol,adv_interf_bool,melt_param, N_melt,
    Tl, epsi_m, Ra_avg, dBLH, dBLC_N, dBLC_S,dBLC_avg, eta_u_N,eta_u_S, Phi_avg, Phi_N, Phi_S, Phi_eff_N, Phi_eff_S, Va_avg, dmadtm, qm_N, qm_S, qcr_N, qcr_S, qc, St
);
    if(t == 0){ dt = dt_init*an_s;} else{
    dt2 =dr_min / (4*std::max({std::abs(dDldt_N),std::abs(dDldt_S),std::abs(dDcdt_N),std::abs(dDcdt_S)}));
    if( t/an_s > t_acc && t_wmelt/an_s > t_acc_m){
        dt_conduct_N = dt_max*an_s;
        dt_conduct_S = dt_max*an_s;
    }
    else{
    if(rPhi_N.maxCoeff()==0){dt_conduct_N = dt_max*an_s;}else{dt_conduct_N = dt_max*an_s/20;}
    if(rPhi_S.maxCoeff()==0){dt_conduct_S = dt_max*an_s;}else{dt_conduct_S = dt_max*an_s/20;}
    }
    dt2=std::min({dt2,dt_conduct_N,dt_conduct_S});
    if(std::abs(dt2*dtmdt) > DTMAX ){dt2=DTMAX/dtmdt;}
    dt2 = std::min(dt*2,dt2);
    if(dt2/an_s>dt_max){dt=dt_max*an_s;} else if(dt2/an_s < dt_min){dt=dt_min*an_s;} else{dt=dt2;}
    }
   
    }

    // H_mag = HaT * T(r) + Hb
    if(dDcdt_N  > 0){

        Vmelt_N = dDcdt_N * dt * 4. * M_PI* pow(Rm_N,2.);
        Hb_N  = dDcdt_N * rho_cr  * (Acr_N/Vcr_N) ;
        HaT_N =  -dDcdt_N * C_cr*  rho_cr * (Acr_N/Vcr_N);
       

    } else { Vmelt_N =0; Hb_N = 0; HaT_N = 0;}

    if(dDcdt_S  > 0){

        Vmelt_S = dDcdt_S * dt * 4. * M_PI* pow(Rm_S,2.);
        Hb_S  = dDcdt_S * rho_cr  * (Acr_S/Vcr_S);
        HaT_S = -dDcdt_S * C_cr*  rho_cr * (Acr_S/Vcr_S);
        

    } else { Vmelt_S = 0; Hb_S = 0; HaT_S = 0;}

   
    if(steady == 1){ //in case of steady state
        
        // Number of points N/S calculation
        if(Dl_N/(N_N+1) > dr_min){N_N = N_N + 1;}
        if(Dl_S/(N_S+1) > dr_min){N_S = N_S + 1;}


        
        //Geometry
        geothermo(THERMOGEO_N,N_N,Rp,Rl_N,Rm_N,k_cr,k_m,C_cr,C_m,rho_cr,rho_m,dDcdt_N,P_cr_N,P_lith_N,Hb_N,L_m,L_cr,fbase,fmag,gl,Tl,adv_lith_bool);
        geothermo(THERMOGEO_S,N_S,Rp,Rl_S,Rm_S,k_cr,k_m,C_cr,C_m,rho_cr,rho_m,dDcdt_S,P_cr_S,P_lith_S,Hb_S,L_m,L_cr,fbase,fmag,gl,Tl,adv_lith_bool);
 
        std::get<0>(contact) = std::get<14>(THERMOGEO_N);
        std::get<1>(contact) = std::get<14>(THERMOGEO_S);

    } else {  // Unsteady state

        
        if( t/an_s > t_acc && t_wmelt/an_s > t_acc_m){ // if the time is higher than the time of acceleration (in parameters) and the time without melt (t_wmelt) is higher than the time in paramters
        // Acceleration
        
        geothermo_acc(THERMOGEO_N,N_N,Nc,dr_max,Rp,Rl_N,Rm_N,k_cr,k_m,C_cr,C_m,rho_cr,rho_m,P_cr_N,P_lith_N,gl);
        geothermo_acc(THERMOGEO_S,N_S,Nc,dr_max,Rp,Rl_S,Rm_S,k_cr,k_m,C_cr,C_m,rho_cr,rho_m,P_cr_S,P_lith_S,gl);
    
        
        } else {  // Normal computation

        if(Dl_N/(N_N+1) > dr_min){N_N = N_N + 1;}
        if(Dl_S/(N_S+1) > dr_min){N_S = N_S + 1;}
    
        // std::cout <<" N_N = " << N_N  << std::endl;
        geothermo(THERMOGEO_N,N_N,Rp,Rl_N,Rm_N,k_cr,k_m,C_cr,C_m,rho_cr,rho_m,dDcdt_N,P_cr_N,P_lith_N,Hb_N,L_m,L_cr,fbase,fmag,gl,Tl,adv_lith_bool);
        // std::cout <<" N_S = " << N_S  << std::endl;
        geothermo(THERMOGEO_S,N_S,Rp,Rl_S,Rm_S,k_cr,k_m,C_cr,C_m,rho_cr,rho_m,dDcdt_S,P_cr_S,P_lith_S,Hb_S,L_m,L_cr,fbase,fmag,gl,Tl,adv_lith_bool);

        }

        std::get<0>(contact) = std::get<14>(THERMOGEO_N);
        std::get<1>(contact) = std::get<14>(THERMOGEO_S);

    }

   

    Nm_N = std::get<0>(THERMOGEO_N);
    Nm_S = std::get<0>(THERMOGEO_S);

    T_N.resize(N_N);
    T_S.resize(N_S);

    
    
    if(interpol_bool == 1){
    interpolate(T1_N,R1_N,T_N,std::get<1>(THERMOGEO_N),Tl,Ts);
    interpolate(T1_S,R1_S,T_S,std::get<1>(THERMOGEO_S),Tl,Ts);
    }

    else{
    no_interpolate(T1_N,T_N,Tl,Ts);
    no_interpolate(T1_S,T_S,Tl,Ts); 
    }
   

    
    if( t/an_s > t_acc && t_wmelt/an_s > t_acc_m){

        rPhi_N = Eigen::VectorXd::Zero(N_N);
        rPhi_S = Eigen::VectorXd::Zero(N_S);

        Tsol_N = Eigen::VectorXd::Zero(N_N);
        Tsol_S = Eigen::VectorXd::Zero(N_S);

        Tliq_N = Eigen::VectorXd::Zero(N_N);
        Tliq_S = Eigen::VectorXd::Zero(N_S);

        DTsolcr_N = Eigen::VectorXd::Zero(N_N);
        DTsolcr_S = Eigen::VectorXd::Zero(N_S);


        Vcr_melt_N = 0;  Phi_cr_avg_N = 0; Vlith_melt_N = 0; Phi_lith_avg_N = 0;
        Vcr_melt_S = 0;  Phi_cr_avg_S = 0; Vlith_melt_S = 0; Phi_lith_avg_S = 0;
    }
    
    else{

    Tsol_N.resize(N_N);
    Tsol_S.resize(N_S);
    Tliq_N.resize(N_N);
    Tliq_S.resize(N_S);
    rPhi_N.resize(N_N);
    rPhi_S.resize(N_S);
    DTsolcr_N.resize(N_N);
    DTsolcr_S.resize(N_S);

    lith_melting(Vcr_melt_N, Phi_cr_avg_N,Phi_cr_max_N,Vlith_melt_N, Phi_lith_avg_N,T_N,rPhi_N,Tsol_N,Tliq_N,DTsolcr_N,std::get<9>(THERMOGEO_N),std::get<10>(THERMOGEO_N),solidus,liquidus,Nm_N,N_N,CH2O_p*LMBD_cr_N,X1,X2,Sat_expo,K_cr,gamma_cr,phi_min_cr,Di);
    lith_melting(Vcr_melt_S, Phi_cr_avg_S,Phi_cr_max_S,Vlith_melt_S, Phi_lith_avg_S,T_S,rPhi_S,Tsol_S,Tliq_S,DTsolcr_S,std::get<9>(THERMOGEO_S),std::get<10>(THERMOGEO_S),solidus,liquidus,Nm_S,N_S,CH2O_p*LMBD_cr_S,X1,X2,Sat_expo,K_cr,gamma_cr,phi_min_cr,Di);
   

    }

 

    Volume(Vath,Rl_N,Rl_S,Rc,f);
    
    if(ecrit_tech_b){
    Internal_energy(T_N,T_S,rPhi_N,rPhi_S,std::get<10>(THERMOGEO_N),std::get<10>(THERMOGEO_S),std::get<8>(THERMOGEO_N),std::get<8>(THERMOGEO_S),std::get<0>(THERMOGEO_N),std::get<0>(THERMOGEO_S),rho_c,rho_m,rho_cr,C_c,C_cr,C_m,f,L_m,L_cr,epsi_m,epsi_c,P_ath,Tm,Tc,qsurf_N,qsurf_S,Phi_avg,
    Va_avg,dt_old,Vc,Vath,Ap,Vp,dU,Qs,U_old,H_tot,T_avg);
    dt_old = dt;  
    }

    

    if(Phi_N == 0 && Phi_cr_avg_N == 0 && Phi_S == 0 && Phi_cr_avg_S == 0 && dtmdt < 0 ){ melt_bol = 0; t_wmelt += dt;} else{ melt_bol = 1; t_wmelt = 0;}

    
    if(steady == 1){  // Temperature profil calculation in steady state
        
        analytique(T_N,std::get<1>(THERMOGEO_N),std::get<2>(THERMOGEO_N),ql_N,qsurf_N,N_N,Nm_N,Rl_N,Rm_N,Rp,Ts,Tl,P_cr_N,P_lith_N,k_m,k_cr,rho_cr,C_cr,rho_m,C_m);
        analytique(T_S,std::get<1>(THERMOGEO_S),std::get<2>(THERMOGEO_N),ql_S,qsurf_S,N_S,Nm_S,Rl_N,Rm_N,Rp,Ts,Tl,P_cr_S,P_lith_S,k_m,k_cr,rho_cr,C_cr,rho_m,C_m);
        
    } else {  // Temperature profile calculation in unsteady state
    
        unsteady(N_N,dt,ql_N,qsurf_N,T_N,rPhi_N,DTsolcr_N,THERMOGEO_N,HaT_N,Vcr_N,fmag,fbase,Tl);
        unsteady(N_S,dt,ql_S,qsurf_S,T_S,rPhi_S,DTsolcr_S,THERMOGEO_S,HaT_S,Vcr_S,fmag,fbase,Tl);
          
    }

    // std::cout << "ql_S = " << ql_S << std::endl;

    if(ecrit_time_b == 1){ecriture_time(t/an_s,Tl,Tm,Tc,Tb,Tcr_N,Tcr_S, Dl_N,Dl_S,Dc_N,Dc_S,dBLC_N,dBLC_S,dBLH,Phi_eff_N,Phi_N,Phi_eff_S,Phi_S,rho_m,Phi_cr_avg_N,Phi_cr_avg_S,Phi_cr_max_N,Phi_cr_max_S,Vcr_melt_N,Vcr_melt_S,Phi_lith_avg_N,Phi_lith_avg_S,Vlith_melt_N,Vlith_melt_S,P_ath,LMBD_ath,LMBD_lith_N,LMBD_lith_S,LMBD_cr_N,LMBD_cr_S,Vcr_melt_N,Vcr_melt_S,qmoho_N,qmoho_S,dossier_time);}
    if(ecrit_tech_b == 1){ecriture_tech(t/an_s,dt,qsurf_N,qsurf_S,qc,T_avg,Ra_avg,P_prim,H_tot,Qs,dU,epsi_m,St,dmadtm,ql_N,ql_S,HaT_N,HaT_S,dossier_tech);}
//P_prim,H_tot,Qs,dU,epsi_m,St,dmadtm
    T1_N.resize(N_N);
    T1_S.resize(N_S);

    T1_N = T_N;
    T1_S = T_S; 

    R1_N.resize(N_N);
    R1_S.resize(N_S);

    R1_N = std::get<1>(THERMOGEO_N);
    R1_S = std::get<1>(THERMOGEO_S);

    
    Tcr_N = T_N(Nm_N);
    Tcr_S = T_S(Nm_S);
    
    
  
    if(Nm_N > 0){
        qmoho_N = std::get<7>(THERMOGEO_N)(Nm_N)*(T_N(Nm_N-1)-T_N(Nm_N));
        
    }
    else{
        qmoho_N = std::get<7>(THERMOGEO_N)(Nm_N+1)*(T_N(Nm_N)-T_N(Nm_N+1));
       
    }
 
    if(Nm_S > 0){
    qmoho_S= std::get<7>(THERMOGEO_S)(Nm_S)*(T_S(Nm_S-1)-T_S(Nm_S));

    }
    else{
        
        qmoho_S= std::get<7>(THERMOGEO_S)(Nm_S+1)*(T_S(Nm_S)-T_S(Nm_S+1));
       
    }
   
    // std::cout << qmoho_N << ", "<< qmoho_S << std::endl;

    if(Phi_cr_avg_N > Phicrmax_N){Phicrmax_N=Phi_cr_avg_N;}
    if(Phi_cr_avg_S > Phicrmax_S){Phicrmax_S=Phi_cr_avg_S;}
    if((Vcr_melt_N / Acr_N) > dmax_cr_N){dmax_cr_N = (Vcr_melt_N / Acr_N);}
    if((Vcr_melt_S / Acr_S) > dmax_cr_S){dmax_cr_S = (Vcr_melt_S / Acr_S);}

    if(Phi_cr_avg_N > 0){
        std::get<0>(agecrust_N) = std::get<0>(agecrust_N) + dt/an_s;
        std::get<1>(agecrust_N) = t/an_s;
        if(std::get<2>(agecrust_N) == 0){std::get<2>(agecrust_N) = t/an_s; }
        }
    if(Phi_cr_avg_S > 0){
       std::get<0>(agecrust_S) = std::get<0>(agecrust_S) + dt/an_s;
       std::get<1>(agecrust_S) = t/an_s;
       if(std::get<2>(agecrust_S) == 0){std::get<2>(agecrust_S) = t/an_s; }
        
    }
    
    // New values for the next time step
    Tm = Tm + dtmdt * dt;
    Tc = Tc + dtcdt * dt;
    Rl_N2 = Rl_N - dDldt_N * dt;
    Rl_S2 = Rl_S - dDldt_S * dt;
    Rm_N2 = Rm_N - dDcdt_N * dt;
    Rm_S2 = Rm_S - dDcdt_S * dt;



    if(Rm_N2 < Rl_N2){Rl_N2 = Rm_N2; dDldt_N=(Rl_N2-Rl_N)/dt;}
    if(Rm_S2 < Rl_S2){Rl_S2 = Rm_S2; dDldt_S=(Rl_S2-Rl_S)/dt;}
    
    
    enrichment(LMBD_cr_N,LMBD_cr_S,LMBD_lith_N,LMBD_lith_S,LMBD_ath,LMBD_liq_N,LMBD_liq_S,f,LMBDdt_bol,Di,Phi_N,Phi_S,
    Rp, Rc, Vmelt_N,Vmelt_S,
    Rl_N, Rl_N2, Rl_S, Rl_S2,
    Rm_N, Rm_N2, Rm_S, Rm_S2,
    rho_m, rho_cr);

    Rl_N=Rl_N2;
    Rl_S=Rl_S2;
    Rm_N=Rm_N2;
    Rm_S=Rm_S2;

    // std::cout << "Rm_S-Rl_S : "<< Rm_S-Rl_S  << std::endl;

    Dl_N = Rp - Rl_N;
    Dl_S = Rp - Rl_S;
    Dc_N = Rp - Rm_N;
    Dc_S = Rp - Rm_S;
    Dc = Rp-pow(f*pow(Rm_N,3.)+(1-f)*pow(Rm_S,3.),1./3.);

    std::get<0>(Pm) = rho_cr*gl*Dc_N + rho_m*gl*(Dl_N-Dc_N+dBLC_N);
    std::get<1>(Pm) = rho_cr*gl*Dc_S + rho_m*gl*(Dl_S-Dc_S+dBLC_S);
    std::get<2>(Pm) = rho_cr*gl*Dc + rho_m*gc*(Rp-Rc-dBLH-Dc);
    
    Volume(Vm,Rm_N,Rm_S,Rc,f);
    Vcr=Vp-Vm;
    rho_m=(Vp*rho_p-Vcr*rho_cr)/Vm;

    Vcr_N = f*4./3.*M_PI*(pow(Rp,3.)-pow(Rm_N,3.));
    Vcr_S = (1-f)*4./3.*M_PI*(pow(Rp,3.)-pow(Rm_S,3.));
    Acr_N = f*4.0*M_PI*pow(Rm_N,2.0);
    Acr_S = (1-f)*4.0*M_PI*pow(Rm_S,2.0);

   // std::cout << "FIN DE BOUCLE" << std::endl;
     
    t = t + dt;
    if(ecrit_profil_b == 1){
    ecriture_profil(t/an_s,T1_N,R1_N,std::get<2>(THERMOGEO_N),std::get<9>(THERMOGEO_N),rPhi_N,std::get<11>(THERMOGEO_N),std::get<6>(THERMOGEO_N),Tsol_N,std::get<13>(THERMOGEO_N),N_N,dossier_profilN);
    ecriture_profil(t/an_s,T1_S,R1_S,std::get<2>(THERMOGEO_S),std::get<9>(THERMOGEO_S),rPhi_S,std::get<11>(THERMOGEO_S),std::get<6>(THERMOGEO_S),Tsol_S,std::get<13>(THERMOGEO_S),N_S,dossier_profilS);
    }
   
    if(Rl_N > Rp || Rl_S > Rp){ stop = 1;}
    if(URG_STOP){
        if(std::get<0>(contact) == 1 || std::get<1>(contact) == 1){
        std::cout << "STOP - CONTACT, " << "Dc_S = "<< Dc_S << ", Dl_S = " << Dl_S << ",t = " << t/an_s << std::endl;
        stop = 1;
        }
        if((Dc_N - dr_min) > Dc_S ){
        std::cout << "STOP - N > S : DcrN =  " << Dc_N  << ", DcrS = " << Dc_S << ", t = " << t/an_s << std::endl;
        stop = 1;
        }
    }
}


Rl = pow((f*pow(Rl_N,3)+(1-f)*pow(Rl_S,3)),1./3.);
Tp = Tm - a_m*gl*Tm/C_m*(Rp-Rl+dBLC_avg);

if(N_N > 3){
std::get<0>(dTdz_N) = -(T1_N(2)-T1_N(0))/(R1_N(2)-R1_N(0));}
else{std::get<0>(dTdz_N) =0.0;}
if(N_S > 3){
std::get<0>(dTdz_S) = -(T1_S(2)-T1_S(0))/(R1_S(2)-R1_S(0));}
else{std::get<0>(dTdz_S) =0.0;}

if(Nm_N > 2){
std::get<1>(dTdz_N) =   -(T_N(Nm_N-2)-T_N(Nm_N-1))/(R1_N(Nm_N-2)-R1_N(Nm_N-1));}
else if(Nm_N > 1){
std::get<1>(dTdz_N) =   -(T_N(Nm_N-1)-T_N(Nm_N))/(R1_N(Nm_N-1)-R1_N(Nm_N));}
else{std::get<1>(dTdz_N) =  std::get<0>(dTdz_N); }

if(Nm_S > 2){
std::get<1>(dTdz_S) =  -(T_S(Nm_S-2)-T_S(Nm_S-1))/(R1_S(Nm_S-2)-R1_S(Nm_S-1));}
else if(Nm_S > 1){
std::get<1>(dTdz_S) =  -(T_S(Nm_S-1)-T_S(Nm_S))/(R1_S(Nm_S-1)-R1_S(Nm_S));}
else{std::get<1>(dTdz_S) = std::get<0>(dTdz_S);}

if(ecrit_time_b == 1){ecriture_time(t/an_s,Tl,Tm,Tc,Tb,Tcr_N,Tcr_S, Dl_N,Dl_S,Dc_N,Dc_S,dBLC_N,dBLC_S,dBLH,Phi_eff_N,Phi_N,Phi_eff_S,Phi_S,rho_m,Phi_cr_avg_N,Phi_cr_avg_S,Phi_cr_max_N,Phi_cr_max_S,Vcr_melt_N,Vcr_melt_S,Phi_lith_avg_N,Phi_lith_avg_S,Vlith_melt_N,Vlith_melt_S,P_ath,LMBD_ath,LMBD_lith_N,LMBD_lith_S,LMBD_cr_N,LMBD_cr_S,Vcr_melt_N,Vcr_melt_S,qmoho_N,qmoho_S,dossier_time);}
if(ecrit_tech_b == 1){ecriture_tech(t/an_s,dt,qsurf_N,qsurf_S,qc,T_avg,Ra_avg,P_prim,H_tot,Qs,dU,epsi_m,St,dmadtm,ql_N,ql_S,HaT_N,HaT_S,dossier_tech);}

if(stop == 0){std::cout << "Run termine sans probleme" << std::endl;} else { std::cout << "Run termine avec probleme" << std::endl;}
return stop;
}


void Internal_energy(Eigen::VectorXd const &T_N, Eigen::VectorXd const &T_S, Eigen::VectorXd const &PhiN,Eigen::VectorXd const &PhiS,Eigen::VectorXd const &dVN_r,Eigen::VectorXd const &dVS_r,Eigen::VectorXd const &PN_r,Eigen::VectorXd const &PS_r, int const &Nm_N, int const &Nm_S,
double const &rho_c, double const &rho_m, double const &rho_cr, double const &Cc, double const &Ccr,double const &Cm, double const &f,double const &L_m, double const &L_cr, double const &epsi_m, double const &epsi_c, double const &Path,
double const &Tm,double const &Tc,double const &qsurf_N, double const &qsurf_S, double const &phi_avg, double const &Va, double const &dt, double const &Vc, double const &Vcm, double const &Ap, double const &Vp,
double &dU_tot, double &Qs, double &U_old, double &H_tot, double &T_avg){

    double U_cm= 0, U_c = 0, U_N = 0, U_S = 0, H_N = 0, H_S = 0, H_cm = 0 ,T_cm =  0, T_avg_N =0, T_avg_S =0;
    int N_N = T_N.rows();
    int N_S = T_S.rows();
    U_cm = rho_m * (Vcm * Cm * Tm * epsi_m + L_m * phi_avg * Va);
    T_cm = epsi_m * Tm;
    U_c = rho_c * Vc * Cc * Tc * epsi_c;
    H_cm = Vcm * Path;
    
   for(int i = 0; i < N_N; i++){
    if(i<Nm_N){
        U_N = U_N +  rho_m * dVN_r(i) * ( Cm * T_N(i) +  L_m * PhiN(i));
        
       
    }
    else{
        U_N = U_N + rho_cr * dVN_r(i) * ( Ccr  * T_N(i) +  L_cr * PhiN(i));
        
    }
    T_avg_N = T_N(i)* dVN_r(i);
    H_N = H_N + PN_r(i)*dVN_r(i);
    }


    for(int i = 0; i < N_S; i++){
    if(i<Nm_S){
        U_S = U_S + rho_m * dVS_r(i) * ( Cm * T_S(i) +  L_m * PhiS(i));
    }
    else{
        U_S = U_S + rho_cr * dVS_r(i) * ( Ccr  * T_S(i) +  L_cr * PhiS(i));
    }
    T_avg_S = T_S(i)* dVS_r(i);
    H_S = H_S + PS_r(i)*dVS_r(i);
  
   }

    dU_tot = (U_c + U_cm + U_N * f + U_S * (1.-f) - U_old)/dt;
    U_old = U_c + U_cm + U_N * f + U_S * (1.-f);
    Qs = Ap * (f*qsurf_N + (1.-f) * qsurf_S);
    H_tot = H_cm + f * H_N + (1.-f) * H_S;
    T_avg = (T_cm * Vcm + T_avg_N * f  + T_avg_S * (1.-f) )/Vp;

}
