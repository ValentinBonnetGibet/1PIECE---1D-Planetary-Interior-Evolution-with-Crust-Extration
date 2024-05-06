#include<iostream>
#include<vector>
#include<array>
#include<cmath>
#include<string>
#include "Eigen/Dense"


void Internal_heat(double &H,double const &t, double const &LMBD,double const &rho, double const &an_s,Eigen::MatrixXd const &RAD){

    int N = RAD.rows();
    double Q = 0.0;
    double tyr=t/an_s;

    for(int i = 0; i < N; ++i){
    
        Q = Q + RAD(0,i)*RAD(1,i)*RAD(3,i)*std::exp(-tyr*std::log(2.0)/RAD(2,i));
    }

    H = rho * LMBD * Q;

}

void enrichment(double &LMBD_cr_N,double &LMBD_cr_S, double &LMBD_lith_N, double &LMBD_lith_S, double &LMBD_ath,double &LMBD_liq_N,double &LMBD_liq_S, double const &f,bool const &LMBDdt_bol, double const &Di, double const &Phi_N, double const &Phi_S,
double const &Rp, double const &Rc, double const &Vmeltt_N, double const &Vmeltt_S,
double const &Rl1_N, double const &Rl2_N, double const &Rl1_S, double const &Rl2_S,
double const &Rm1_N, double const &Rm2_N, double const &Rm1_S, double const &Rm2_S,
double const &rho_m, double const &rho_cr
){

    double V1_ath=4./3.*M_PI*(f*(pow(Rl1_N,3.)-pow(Rc,3.))+(1.-f)*(pow(Rl1_S,3.)-pow(Rc,3.)));
    double V2_ath=4./3.*M_PI*(f*(pow(Rl2_N,3.)-pow(Rc,3.))+(1.-f)*(pow(Rl2_S,3.)-pow(Rc,3.)));
 
    double V1_cr_N=f*4./3.*M_PI*(pow(Rp,3.)-pow(Rm1_N,3.));
    double V1_cr_S=(1.-f)*4./3.*M_PI*(pow(Rp,3.)-pow(Rm1_S,3.));
    
    double V2_cr_N = f*4./3.*M_PI*(pow(Rp,3.)-pow(Rm2_N,3.));
    double V2_cr_S = (1.-f)*4./3.*M_PI*(pow(Rp,3.)-pow(Rm2_S,3.));
    
    double V1_lith_N=f*4./3.*M_PI*(pow(Rm1_N,3.)-pow(Rl1_N,3.));
    double V1_lith_S=(1.-f)*4./3.*M_PI*(pow(Rm1_S,3.)-pow(Rl1_S,3.));
    
    double V2_lith_N=f*4./3.*M_PI*(pow(Rm2_N,3.)-pow(Rl2_N,3.));
    double V2_lith_S=(1.-f)*4./3.*M_PI*(pow(Rm2_S,3.)-pow(Rl2_S,3.)) ;  
    
    double DV_lith_N=V2_lith_N-V1_lith_N;
    double DV_lith_S=V2_lith_S-V1_lith_S;
    
    double DV_cr_N=V2_cr_N-V1_cr_N;
    double DV_cr_S=V2_cr_S-V1_cr_S;
    
    
    double Vmelt_N=Vmeltt_N*f;
    double Vmelt_S=Vmeltt_S*(1-f);

    LMBD_liq_N= LMBD_cr_N;
    LMBD_liq_S= LMBD_cr_S;

    double V_ero_N = 0;
    double V_ero_S = 0;

    double LMBD_ero_N = 0;
    double LMBD_ero_S = 0;

    // double & test = LMBD_cr_N;
    
    if (LMBDdt_bol == 1) {

        if (Phi_N > 0 ){

            LMBD_liq_N = LMBD_ath/(Phi_N+(1-Phi_N)*Di);
         } 
        else {

            LMBD_liq_N=0;

        }
            
            
        if (Phi_S > 0) {

            LMBD_liq_S = LMBD_ath/(Phi_S+(1-Phi_S)*Di);

        }
        else {

            LMBD_liq_S = 0;
        }
            

    }  
    
    if(DV_cr_N < 0 && DV_cr_S < 0 ){

    
        V_ero_N=DV_cr_N;
        V_ero_S=DV_cr_S;
        
        LMBD_ero_N=LMBD_cr_N;
        LMBD_ero_S=LMBD_cr_S;
   
       
        LMBD_cr_N = (LMBD_cr_N*V1_cr_N+LMBD_liq_N*Vmelt_N+V_ero_N*LMBD_ero_N)/(V2_cr_N);
        LMBD_cr_S = (LMBD_cr_S*V1_cr_S+LMBD_liq_S*Vmelt_S+V_ero_S*LMBD_ero_S)/(V2_cr_S);
        
        LMBD_ath = (LMBD_ath*V1_ath*rho_m-LMBD_liq_N*Vmelt_N*rho_cr+V_ero_N*LMBD_ero_N*rho_cr-LMBD_liq_S*Vmelt_S*rho_cr+V_ero_S*LMBD_ero_S*rho_cr)/(V2_ath*rho_m);
        
        LMBD_lith_N=LMBD_ath ;
        LMBD_lith_S=LMBD_ath ;
       // std::cout << " 1 :"<< test << std::endl; 
    }
    else if (DV_cr_N < 0 && DV_lith_S < 0) {
        
        
        V_ero_N=DV_cr_N;
        V_ero_S=DV_lith_S;
        
        LMBD_ero_N=LMBD_cr_N ;
        LMBD_ero_S=LMBD_lith_S ;
        
        LMBD_ath = (LMBD_ath*V1_ath*rho_m-LMBD_liq_N*Vmelt_N*rho_cr-V_ero_N*LMBD_ero_N*rho_cr-LMBD_liq_S*Vmelt_S*rho_cr-V_ero_S*LMBD_ero_S*rho_m)/(V2_ath*rho_m);
        
        LMBD_cr_N = (LMBD_cr_N*V1_cr_N+LMBD_liq_N*Vmelt_N+V_ero_N*LMBD_ero_N)/(V2_cr_N);
        LMBD_cr_S = (LMBD_cr_S*V1_cr_S+LMBD_liq_S*Vmelt_S)/(V2_cr_S);
        
        LMBD_lith_N=LMBD_ath;
        if (V2_lith_S > 0 ){LMBD_lith_S=(LMBD_lith_N*V1_lith_S+V_ero_S*LMBD_ero_S)/(V2_lith_S);} else{LMBD_lith_S=LMBD_ath;}
        // std::cout << " 2 :"<< test << std::endl;
       
    }
    else if (DV_cr_S < 0 && DV_lith_N > 0) {
        
        
        V_ero_N=0;
        V_ero_S=DV_cr_S;
        
        LMBD_ero_N=0;
        LMBD_ero_S=LMBD_cr_S ;
        
        LMBD_ath = (LMBD_ath*V1_ath*rho_m-LMBD_liq_N*Vmelt_N*rho_cr-LMBD_liq_S*Vmelt_S*rho_cr-V_ero_S*LMBD_ero_S*rho_cr)/(V2_ath*rho_m);
        
        LMBD_cr_N = (LMBD_cr_N*V1_cr_N+LMBD_liq_N*Vmelt_N)/(V2_cr_N);
        LMBD_cr_S = (LMBD_cr_S*V1_cr_S+LMBD_liq_S*Vmelt_S+V_ero_S*LMBD_ero_S)/(V2_cr_S);
        
        LMBD_lith_S=LMBD_ath;
        if (V2_lith_N > 0 ){LMBD_lith_N=(LMBD_lith_N*V1_lith_N+DV_lith_N*LMBD_ath)/(V2_lith_N);} else{LMBD_lith_N=LMBD_ath;}
        // std::cout << " 3 :"<< test << std::endl;
       
    }
    else if (DV_cr_N < 0 && DV_lith_S > 0) {
        
        
        V_ero_S=0;
        V_ero_N=DV_cr_N;
        
        LMBD_ero_S=0;
        LMBD_ero_N=LMBD_cr_N ;
        
        LMBD_ath = (LMBD_ath*V1_ath*rho_m-LMBD_liq_S*Vmelt_S*rho_cr-LMBD_liq_N*Vmelt_N*rho_cr-V_ero_N*LMBD_ero_N*rho_cr)/(V2_ath*rho_m);
        
        LMBD_cr_S = (LMBD_cr_S*V1_cr_S+LMBD_liq_S*Vmelt_S)/(V2_cr_S);
        LMBD_cr_N = (LMBD_cr_N*V1_cr_N+LMBD_liq_N*Vmelt_N+V_ero_N*LMBD_ero_N)/(V2_cr_N);
        
        LMBD_lith_N=LMBD_ath;
        if (V2_lith_S > 0 ){LMBD_lith_S=(LMBD_lith_S*V1_lith_S+DV_lith_S*LMBD_ath)/(V2_lith_S);} else{LMBD_lith_S=LMBD_ath;}
        // std::cout << " 4 :"<< test << std::endl;
       
    }
    else if (DV_cr_S < 0 && DV_lith_N < 0) {

        
        V_ero_S=DV_cr_S;
        V_ero_N=DV_lith_N;
        
        LMBD_ero_S=LMBD_cr_S;
        LMBD_ero_N=LMBD_lith_N;
        
        
        LMBD_ath = (LMBD_ath*V1_ath*rho_m-LMBD_liq_S*Vmelt_S*rho_cr-V_ero_S*LMBD_ero_S*rho_cr-LMBD_liq_N*Vmelt_N*rho_cr-V_ero_N*LMBD_ero_N*rho_m)/(V2_ath*rho_m);
        
        
        LMBD_cr_S = (LMBD_cr_S*V1_cr_S+LMBD_liq_S*Vmelt_S+V_ero_S*LMBD_ero_S)/(V2_cr_S);
        LMBD_cr_N = (LMBD_cr_N*V1_cr_N+LMBD_liq_N*Vmelt_N)/(V2_cr_N);
        
        LMBD_lith_S=LMBD_ath;
        if(V2_lith_N > 0) {LMBD_lith_N=(LMBD_lith_N*V1_lith_N+V_ero_N*LMBD_ero_N)/(V2_lith_N);} else { LMBD_lith_N=LMBD_ath ;}
       // std::cout << " 5 :"<< test << std::endl;
        
    }

    else if(DV_lith_N < 0 && DV_lith_S < 0){

        V_ero_N=DV_lith_N;
        V_ero_S=DV_lith_S;
        
        LMBD_ero_N=LMBD_lith_N;
        LMBD_ero_S=LMBD_lith_S;
        
        
        LMBD_cr_N = (LMBD_cr_N*V1_cr_N+LMBD_liq_N*Vmelt_N)/(V2_cr_N);
        LMBD_cr_S = (LMBD_cr_S*V1_cr_S+LMBD_liq_S*Vmelt_S)/(V2_cr_S);
        
        LMBD_ath= (LMBD_ath*V1_ath*rho_m-LMBD_liq_N*Vmelt_N*rho_cr-LMBD_liq_S*Vmelt_S*rho_cr-LMBD_ero_N*rho_m*DV_lith_N-LMBD_ero_S*rho_m*DV_lith_S)/(V2_ath*rho_m);
        
        if(V2_lith_N > 0 ){LMBD_lith_N=(LMBD_lith_N*V1_lith_N*rho_m+LMBD_lith_N*rho_m*DV_lith_N)/(V2_lith_N*rho_m);} else{LMBD_lith_N=LMBD_ath; } 
        if(V2_lith_S > 0 ){LMBD_lith_S=(LMBD_lith_S*V1_lith_S*rho_m+LMBD_lith_S*rho_m*DV_lith_S)/(V2_lith_S*rho_m);} else{LMBD_lith_S=LMBD_ath; }
       // std::cout << " 6 :"<< test << std::endl;

    }
    else if(DV_lith_N > 0 && DV_lith_S < 0) {


        V_ero_N=DV_lith_N;
        V_ero_S=DV_lith_S;
        
        LMBD_ero_N=LMBD_ath;
        LMBD_ero_S=LMBD_lith_S;
        
        
        LMBD_cr_N = (LMBD_cr_N*V1_cr_N+LMBD_liq_N*Vmelt_N)/(V2_cr_N);
        LMBD_cr_S = (LMBD_cr_S*V1_cr_S+LMBD_liq_S*Vmelt_S)/(V2_cr_S);
        
        LMBD_ath = (LMBD_ath*V1_ath*rho_m-LMBD_liq_N*Vmelt_N*rho_cr-LMBD_liq_S*Vmelt_S*rho_cr-LMBD_ero_N*rho_m*DV_lith_N-LMBD_ero_S*rho_m*DV_lith_S)/(V2_ath*rho_m);
        
        if(V2_lith_N > 0 ){LMBD_lith_N=(LMBD_lith_N*V1_lith_N*rho_m+LMBD_lith_N*rho_m*DV_lith_N)/(V2_lith_N*rho_m);} else{LMBD_lith_N=LMBD_ath; } 
        if(V2_lith_S > 0 ){LMBD_lith_S=(LMBD_lith_S*V1_lith_S*rho_m+LMBD_lith_S*rho_m*DV_lith_S)/(V2_lith_S*rho_m);} else{LMBD_lith_S=LMBD_ath; }
       // std::cout << " 7 :"<< test << std::endl;
        
    } 
    else if(DV_lith_N < 0 and DV_lith_S > 0){
        
        V_ero_S=DV_lith_S;
        V_ero_N=DV_lith_N;
        
        LMBD_ero_S=LMBD_ath;
        LMBD_ero_N=LMBD_lith_N;
        
        
        LMBD_cr_S = (LMBD_cr_S*V1_cr_S+LMBD_liq_S*Vmelt_S)/(V2_cr_S);
        LMBD_cr_N = (LMBD_cr_N*V1_cr_N+LMBD_liq_N*Vmelt_N)/(V2_cr_N);
        
        LMBD_ath = (LMBD_ath*V1_ath*rho_m-LMBD_liq_N*Vmelt_N*rho_cr-LMBD_liq_S*Vmelt_S*rho_cr-LMBD_ero_N*rho_m*DV_lith_N-LMBD_ero_S*rho_m*DV_lith_S)/(V2_ath*rho_m);
        
        if(V2_lith_S > 0 ){LMBD_lith_S=(LMBD_lith_S*V1_lith_S*rho_m+LMBD_lith_S*rho_m*DV_lith_S)/(V2_lith_S*rho_m);} else{LMBD_lith_S=LMBD_ath; } 
        if(V2_lith_N > 0 ){LMBD_lith_N=(LMBD_lith_N*V1_lith_N*rho_m+LMBD_lith_N*rho_m*DV_lith_N)/(V2_lith_N*rho_m);} else{LMBD_lith_N=LMBD_ath; }
      //  std::cout << " 8 :"<< test << std::endl;
    }
        
    else {
    
        LMBD_cr_N = (LMBD_cr_N*V1_cr_N+LMBD_liq_N*Vmelt_N)/(V2_cr_N);
        LMBD_cr_S = (LMBD_cr_S*V1_cr_S+LMBD_liq_S*Vmelt_S)/(V2_cr_S);
        
        LMBD_ath = (LMBD_ath*V1_ath*rho_m-LMBD_liq_N*Vmelt_N*rho_cr-LMBD_ath*rho_m*DV_lith_N-LMBD_liq_S*Vmelt_S*rho_cr-LMBD_ath*rho_m*DV_lith_S)/(V2_ath*rho_m);
        
        if(V2_lith_N > 0 ){LMBD_lith_N=(LMBD_lith_N*V1_lith_N*rho_m+LMBD_ath*rho_m*DV_lith_N)/(V2_lith_N*rho_m);} else{LMBD_lith_N=LMBD_ath; }
        if(V2_lith_S > 0 ){LMBD_lith_S=(LMBD_lith_S*V1_lith_S*rho_m+LMBD_ath*rho_m*DV_lith_S)/(V2_lith_S*rho_m);} else{LMBD_lith_S=LMBD_ath; }
       // std::cout << " 9 :"<< test << std::endl;
    }
}