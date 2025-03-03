#include<iostream>
#include<fstream>
#include<vector>
#include<array>
#include<cmath>
#include<string>
#include "Eigen/Dense"



bool ecriture_initial(char* dossier_time,char* dossier_tech, char * dossier_profilN, char* dossier_profilS, bool const &ecrit_time_b, bool const &ecrit_tech_b, bool const &ecrit_profil_b){
    
    bool test = 1;
    if(ecrit_time_b == 1){
    std::ofstream fichier_time(dossier_time);
    

    if(fichier_time){

        fichier_time << "t" << " " << "dt" << " " << "Ra" << " " << "Tl" << " " << "Tm" << " " << "Tc"<< " " << "Tb" << " " <<  "Dl_N" << " " <<  "Dl_S" << " "
        <<  "Dc_N" << " " <<  "Dc_S" << " " << "dBLC_N" << " " << "dBLC_S" << " " << "dBLH" << " " <<  "Phi_N" << " " <<  "Phi_S" << " " <<  "Phi_eff_N" << " " <<  "Phi_eff_S" << " " << "rho_m" 
         << " " << "ql_N" << " " << "ql_S" << " " << "qm_N" << " " << "qm_S" << " " << "Path" << " " << "LMBD_ath"
        << " " << "LMBD_lith_N" << " " << "LMBD_lith_S" << " " << "LMBD_cr_N" << " " << "LMBD_cr_S" << std::endl;
        fichier_time.close();

        // ,t,dt,Ra,Tl,Tm,Tb,Tc,Tcr_N,Tcr_S,Dl_N,Dl_S,Dc_N,Dc_S,dBLC_N,dBLC_S,dBLH,eta_u_N,eta_u_S,Phi_N,Phi_S,Phi_eff_N,Phi_eff_S,rho_m,ql_N,ql_S,qm_N,qm_S,Path,LMBD_ath,LMBD_lith_N,LMBD_lith_S,LMBD_cr_N,LMBD_cr_S,P_prim, P_tot,dQs,dU,epsi_m,St,dmadtm
    }
    else
    {
        std::cout << "ERREUR: Impossible de creer le fichier " << "time.txt" << std::endl;
        test = 0;
    }

    }

    if(ecrit_tech_b == 1){
    std::ofstream fichier_time(dossier_tech);
    

    if(fichier_time){

        fichier_time << "t" << " " << "dt" << " " << "Ra" << " " 
        << " " << "Pprim" << " " << "Ptot" << " " << "dQs" << " " << "dU" << " " << "epsim" << " " << "St" << " " << "dmadtm" << std::endl;
        fichier_time.close();

        // ,t,dt,Ra,Tl,Tm,Tb,Tc,Tcr_N,Tcr_S,Dl_N,Dl_S,Dc_N,Dc_S,dBLC_N,dBLC_S,dBLH,eta_u_N,eta_u_S,Phi_N,Phi_S,Phi_eff_N,Phi_eff_S,rho_m,ql_N,ql_S,qm_N,qm_S,Path,LMBD_ath,LMBD_lith_N,LMBD_lith_S,LMBD_cr_N,LMBD_cr_S,P_prim, P_tot,dQs,dU,epsi_m,St,dmadtm
    }
    else
    {
        std::cout << "ERREUR: Impossible de creer le fichier " << "time.txt" << std::endl;
        test = 0;
    }

    }

    if(ecrit_profil_b == 1){

    

    std::ofstream fichier_profilN(dossier_profilN);

    if(fichier_profilN){

        fichier_profilN << "format : 1st line time (yr), 2nd line N point at this time, table : R, T etc., line break" << std::endl;
        fichier_profilN.close();
    }
    else
    {
        std::cout << "ERREUR: Impossible de creer le fichier " << "profil.txt" << std::endl;
        test = 0;
    }
    
    std::ofstream fichier_profilS(dossier_profilS);

    if(fichier_profilS){

        fichier_profilS << "format : 1st line time (yr), 2nd line N point at this time, table : R, T etc., line break" << std::endl;
        fichier_profilS.close();
    }
    else
    {
        std::cout << "ERREUR: Impossible de creer le fichier " << "profil.txt" << std::endl;
        test = 0;
    }

    }

    return test;
}

void ecriture_time(double const &t,double const &Tl, double const &Tm, double const &Tc,double const &Tb,double const &Tcr_N,double const &Tcr_S ,double const &Dl_N,double const &Dl_S, double const &Dc_N, double const &Dc_S, double const &dBLC_N, double const &dBLC_S, double const &dBLH,
double const &Phi_eff_N,double const &Phi_N, double const &Phi_eff_S,double const &Phi_S, double const &rho_m,double const &Phi_cr_avg_N, double const &Phi_cr_avg_S, double const &Phi_cr_max_N, double const &Phi_cr_max_S, double const &Vmelt_cr_N, double const &Vmelt_cr_S,double const &Phi_lith_avg_N, double const &Phi_lith_avg_S,double const &Vmelt_lith_N, double const &Vmelt_lith_S,
double const &Path, double const &LMBD_ath, double const &LMBD_lith_N, double const &LMBD_lith_S, double const &LMBD_cr_N, double const &LMBD_cr_S, double const &LMBD_H2O_N, double const &LMBD_H2O_S,double const &qmoho_N, double const &qmoho_S,
char *dossier
){

    
    FILE * fichier = fopen(dossier,"a");

   
    std::fprintf(fichier, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e \n",t,Tl,Tm,Tb,Tc,Tcr_N,Tcr_S,Dl_N,Dl_S,Dc_N,Dc_S,dBLC_N,dBLC_S,dBLH,Phi_N,Phi_S,Phi_eff_N,Phi_eff_S,rho_m,Phi_cr_avg_N,Phi_cr_avg_S,Phi_cr_max_N,Phi_cr_max_S,Vmelt_cr_N,Vmelt_cr_S,Phi_lith_avg_N,Phi_lith_avg_S,Vmelt_lith_N,Vmelt_lith_S,Path,LMBD_ath,LMBD_lith_N,LMBD_lith_S,LMBD_cr_N,LMBD_cr_S,LMBD_H2O_N,LMBD_H2O_S,qmoho_N,qmoho_S);
    fclose(fichier);

}


void ecriture_tech(double const &t, double const &dt, double const &qsurf_N,double const &qsurf_S, double const &qc,double const &T_avg,double const &Ra,
double const &P_prim, double const &P_tot, double const &dQs, double const &dU, double const &epsi_m, double const &St, double const &dmadtm, double const &ql_N, double const &ql_S, double const &HaT_N, double const &HaT_S, double const &Va_N, double const &Va_S,
char *dossier
){

    
    FILE * fichier = fopen(dossier,"a");

   
    std::fprintf(fichier, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e  \n",t,dt,qsurf_N,qsurf_S,qc,T_avg,Ra,P_prim, P_tot,dQs,dU,epsi_m,St,dmadtm,ql_N,ql_S,HaT_N,HaT_S,Va_N,Va_S);
    fclose(fichier);

}


void ecriture_profil(double const &t, Eigen::VectorXd const &T, Eigen::VectorXd const &R, Eigen::VectorXd const &dr,Eigen::VectorXd const &P_r,Eigen::VectorXd const &Phi_r ,Eigen::VectorXd const &u_r,Eigen::VectorXd const &au_r,Eigen::VectorXd const &Tsol,Eigen::VectorXd const &Hb_r,int const &N,
char *dossier)
{

    FILE * profil = fopen(dossier,"a");
    std::fprintf(profil,"%f\n", t/1E6);
    std::fprintf(profil,"%i\n", N);

    for(int i = 0; i < N; ++i){
    std::fprintf(profil, "%e %e %e %e %e %e %e %e %e \n",R(i),T(i),Phi_r(i),dr(i),P_r(i),u_r(i),au_r(i),Tsol(i),Hb_r(i));
    }

    std::fprintf(profil,"\n");
    fclose(profil);


}


void ecriture_exploration(int const &nbr_run,double const &eta0, double const &k0,double const &Dl_N,double const &Dl_S,double const &Dc_N,double const &Dc_S,double const &Tm,double const &Tc,double const &Tcr_max,
double const &LMBD_ath,double const &LMBD_lith_N,double const &LMBD_lith_S,double const &LMBD_cr_N,double const &LMBD_cr_S
 ){

    
    FILE * fichier = fopen("run_resume.txt" ,"a");

    
    std::fprintf(fichier, "\n");
    fclose(fichier);

}

/*
int main(){

std::string fichier_time = "/Users/valentinbonnetgibet/Documents/These/Code Mars/Mars Cpp Dev/time.txt";
std::string fichier_profil = "/Users/valentinbonnetgibet/Documents/These/Code Mars/Mars Cpp Dev/profil.txt"; 
ecriture_initial(fichier_time,fichier_profil);
int N = 100;
double t,Ra,Tm,Tc;


for(int i = 0; i < N; ++i){

    t=i*500;
    Ra=1E6*i;
    Tm=1.8E3+i;
    Tc=2E3-i;
    ecriture_time(fichier_time,t,Ra,Tm,Tc);
}

}*/