#include "1PIECE_main.cpp"
#include<iostream>
#include<fstream>
#include <stdio.h>
#include<string>
#include <mpi.h>



int main(int argc,char* argv[]){

int nb_procs,rank;  
MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD,&nb_procs); 
MPI_Comm_rank(MPI_COMM_WORLD,&rank);

float temps;
clock_t t1,t2,t1_interm,t2_interm;

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
bool unlinear_phi,LMBDdt,RK4,Pressure,Steady, ecrit_time_b, ecrit_tech_b, ecrit_profil_b,URG_STOP,adv_lith_bool,adv_interf_bool, test_run =1;
// bool melt_bol = 1;
auto solidus = std::make_tuple(1E2,1E2,1E2,1E2,1E2,1E2);
auto liquidus = std::make_tuple(1E2,1E2,1E2,1E2,1E2,1E2);


double k0_init,k0_end,eta0_init,eta0_end;
int N_k0,N_eta0,N_tot,N_start = 0,n_run;

// /Users/valentinbonnetgibet/Documents/These/Code Mars/Mars Cpp Dev/
std::string physical_param = "physical_parameters.txt";
std::string numerical_param = "numerical_parameters.txt";
std::string adress_rad = "Radioelement_Wanke.txt";
std::string adress_solidus = "Solidus.txt";
std::string explor_param = "explor_param.txt";


// std::string fichier_time = "/Users/valentinbonnetgibet/Documents/These/Code Mars/Mars Cpp Dev/time.txt";
// std::string fichier_profil = "/Users/valentinbonnetgibet/Documents/These/Code Mars/Mars Cpp Dev/profil.txt";

Lecture_param(physical_param,rheology,melt,melt_param,thermo,Tm0,Ts,dTc0,Dl_init,Dc_init,dDl_init, dDc_init, LMBD_cr_Ni,LMBD_cr_Si,t0,tstop,Rp,Rc,f,Masse_Mars,gl,gc,fmag,fbase,HPE_factor,X1,X2,Sat_expo,CH2O_p,KH2Ocr,gamma_cr,phi_min_cr);
lecture_num(numerical_param,dt_min,dt_max, dt_init,DTMAX,dr_min,dr_max,Nc,N_soliq,N_melt,t_acc,t_acc_m,unlinear_phi,LMBDdt,RK4,Pressure,Steady,ecrit_profil_b,ecrit_time_b,ecrit_tech_b, URG_STOP,adv_lith_bool,adv_interf_bool);
lecture_sol(adress_solidus,solidus,liquidus);
lecture_explo(explor_param,k0_init,k0_end,eta0_init,eta0_end,N_k0,N_eta0);
lecture_rad(adress_rad,RAD,tstop,HPE_factor);


double Tm,Tc,Phicrmax_N = 0,Phicrmax_S = 0,Dl_N,Dl_S,Dc_N,Dc_S,LMBD_ath,LMBD_lith_N,LMBD_lith_S,LMBD_cr_N,LMBD_cr_S,dmax_cr_N,dmax_cr_S, Tp, dBLC_N,dBLC_S;
std::tuple<double,double> dTdz_N;
std::tuple<double,double> dTdz_S;
std::tuple<double,double,double,double> agecrust_N;
std::tuple<double,double,double,double> agecrust_S;

double Phi_eff_i = 0.2;

double rho_c=std::get<4>(thermo);
double Dref = 0.2/3.*(pow(Rp,3.)-pow(Rc,3.))/pow(Rp,2.);
double Vc = 4./3.*M_PI*pow(Rc,3.);
double Ac = 4.*M_PI*pow(Rc,2.);
double Vp = 4./3.*M_PI*pow(Rp,3.) - Vc;
double rho_p = (Masse_Mars - Vc*rho_c)/Vp;
double rho_m;
double Vm,Vcr;
double t = 0;
if(unlinear_phi == 0){std::get<1>(unlinear_table)=1;}
N_tot = (N_k0 * N_eta0) / nb_procs;


Eigen::VectorXd k0_exp = Eigen::VectorXd::Ones(N_k0);
Eigen::VectorXd eta0_exp = Eigen::VectorXd::Ones(N_eta0);
std::string adress_resum_string;
char* adress_resum;

if(rank <10){
adress_resum_string = "run_resume_0"+std::to_string(rank)+".txt";
}
else{
adress_resum_string = "run_resume_" + std::to_string(rank)+ ".txt";
}


adress_resum = &adress_resum_string[0];
 

std::ifstream test_resume(adress_resum);
if(test_resume.fail() == 1){
std::cout << rank << " : rank," << "run_resume_rank.txt don't find, start a new calculation" ;
FILE * resume = fopen(adress_resum ,"w");
std::fprintf(resume, "Nbr run   k0  eta0    Tm  Tc  Tcr_max_N   Tcr_max_S   Dl_N    Dl_S    Dc_N    Dc_S    LMBD_ath    LMBD_lith_N LMBD_lith_S LMBD_cr_N   LMBD_cr_S \n");
fclose(resume);
}
else {
    std::string line,tmp ;
    while(std::getline(test_resume, line)){
        
        N_start ++ ;
       
    }
    test_resume.close();
    
    N_start --;
    tmp ="rm -rf run_" + std::to_string(N_start);
    char* cmd_rm= &tmp[0];
    std::system(cmd_rm);
    std::cout << rank << " : rank," << "run_resume_rank.txt find, start the calculation at N = " << N_start << std::endl;
    
}


double logk0_init = std::log10(k0_init);
double logk0_end = std::log10(k0_end);
double logeta0_init = std::log10(eta0_init);
double logeta0_end = std::log10(eta0_end);

double dk0 = (logk0_end-logk0_init)/(N_k0-1);
double deta0 = (logeta0_end-logeta0_init)/(N_eta0-1);


for(int i=0;i < N_k0; i++){

k0_exp(i) = pow(10,logk0_init+i*dk0);
eta0_exp(i) = pow(10,logeta0_init+i*deta0);

}

int k_i, eta_i;

n_run = rank*N_k0/nb_procs*N_eta0;
std::cout << rank << " : rank," <<  "  n_run start :" << n_run << std::endl;
if(N_start == 0){ k_i = rank*N_k0/nb_procs; eta_i = 0;}
else{k_i=rank*N_k0/nb_procs+N_start/N_tot; eta_i=N_start-k_i*N_k0;}


std::string cmd;
std::string adresse;
char* cmd_mk;
char* dossier_time;
char* dossier_tech;
char* dossier_profilN;
char* dossier_profilS;


for(int i = N_start; i< N_tot; ++i){

if(ecrit_time_b == 1 || ecrit_profil_b == 1){

adresse = "r"+std::to_string(rank) + "_" + std::to_string(i);
std:: string tmp1 = adresse +"_time.txt";
dossier_time = &tmp1[0];
std:: string tmp4 = adresse +"_tech.txt";
dossier_tech = &tmp4[0];
std::string tmp2 = adresse +"_profilN.txt";
dossier_profilN = &tmp2[0];
std::string tmp3 = adresse +"_profilS.txt";
dossier_profilS = &tmp3[0];
}


std::get<0>(rheology) = eta0_exp(eta_i);
std::get<0>(melt) = k0_exp(k_i);
std::cout << rank << " : rank," << i <<  " / " << N_tot << " done" << "  eta0 = " << eta0_exp(eta_i) << "  k0 = " << k0_exp(k_i) << "         nom dossier :  " << dossier_time << std::endl;

t1_interm = clock();
t = 0;

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


t2_interm = clock();
std::cout << rank << " : rank," <<  "  time :" << (float) (t2_interm-t1_interm)/CLOCKS_PER_SEC <<", test_run : "<< test_run <<std::endl;

FILE * resume = fopen(adress_resum ,"a");
std::fprintf(resume, "%i %i %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",i,rank,t,k0_exp(k_i),eta0_exp(eta_i),Tp, Tm, Dc_N, Dc_S, Dl_N, Dl_S, Tc, std::get<3>(agecrust_N) ,std::get<3>(agecrust_S), rho_m,LMBD_cr_N,LMBD_cr_S,LMBD_lith_N,LMBD_lith_S,LMBD_ath);
fclose(resume);


n_run = n_run + 1;
++eta_i;
if(eta_i == N_eta0){eta_i=0;++k_i;}


}


t2 =clock();
temps =(float) (t2-t1)/CLOCKS_PER_SEC;
std::cout << rank << " : rank," << "Temps d'execution : " << temps << " s." << std::endl;

MPI_Finalize();

}