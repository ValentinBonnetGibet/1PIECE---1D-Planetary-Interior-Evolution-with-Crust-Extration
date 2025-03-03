#include<iostream>
#include<vector>
#include<array>
#include<cmath>
#include<string>
#include "Eigen/Dense"
// #include "Radioactive_heat.cpp"


// Volume calculation (V) of a double spherical shell  between two radius R_N and R_S and a same radius at the base R_bot with a volume ratio of N of f
void Volume(double &V, double const & R_N,double const &R_S,double const &Rbot,double const &f){

    V = 4./3.*M_PI*(f*(R_N*R_N*R_N-Rbot*Rbot*Rbot)+(1-f)*(R_S*R_S*R_S-Rbot*Rbot*Rbot));
}

// Area calculation (A) of a double spherical shell geometry for two radius R_N, R_S and a area ratio of R_n of f
void Area(double &A, double const & R_N,double const &R_S,double const &f){

    A = 4.*M_PI*(f*(R_N*R_N)+(1-f)*(R_S*R_S));
}

// Temperature at the base of the lid, following the decription of Davaille & Jaupart 1993
void litho_temp(double &Tl, double const &Tm,std::tuple<double,double,double,double,double,double,double,double,double,double,double> const &rheology){

    Tl = Tm - std::get<4>(rheology) *(std::get<2>(rheology)/std::get<1>(rheology)*Tm*Tm);

}

// Temperature at the top of the lower thermal boudnary layer
void Base_temp(double &Tb,double const &Tm,double const &a_m, double const &Cm,double const &gc,double const& Rl_avg,double const &Rc,double const &dBLC, double const &dBLH){

    Tb= Tm + a_m*gc*Tm/Cm*(Rl_avg-Rc-dBLH-dBLC);

}

// Calculation of the average melt fraction (ma) and effective melt fraction (ma_eff) over the melt Volume (Va) and the melt production at Tm (dmadtm)
void melting(double &ma,double &ma_eff,double &Va,double &dmadtm,double const &Tm,double am,double Cm, double gc, double gl, int N,double const &Dref,
double Pmax, double DTsol,double const &rho_cr,double const &rho_m,double Dc,double Dc_avg,double Rl,double Rp, double const &dBLC, double phi_crit, double const &dTm,
std::tuple<bool,double,double,double,double,double,double,double,double> const &melt_param, double const &LMBD_liq,
std::tuple<double,double,double,double,double,double> const &solidus, std::tuple<double,double,double,double,double,double> const &liquidus, std::tuple<double,double,double,double,double,double,double,double,double,double,double> const &rheology)
{
    double Vab = 0,Vmelt = 0,Vmeltb = 0, Vmelt_eff = 0, mab=0; //
    
    double P0 = rho_cr*gl*Dc; //Pressure at the base of the crust
    
    double Rm = Rp-Dc;        // moho radius
    double Rmax = Rm - (Pmax-P0)/(rho_m*gl);  // Radius of the pressure Pmax where the flotability is negative
    double dr = (Rl-Rmax)/(N-1);              // space step for the discretisation
    double DTsolidus=DTsol*Dc_avg/Dref;        // Increase of the solidus temperature in function of the crust extracted
    
    double Tmb= Tm+dTm;                        // Temperature to calcule dTmdT (ma(Tm) - ma(Tm+dTm)) / dTm
    double Rdblc = Rl-dBLC;                    // Upper boundary layer radius

    if(Rl > Rmax){                              
    double Tl =0, Tlb = 0;

    litho_temp(Tl,Tm,rheology);             // Tl for Tm
    litho_temp(Tlb,Tmb,rheology);         // Tl for Tm + dTm

    double R;                         
    double P;                         
    double T;                          //Temperature
    double Tb;                        // Tb Temperature for ma(eff+dTm)
    double Tsol;                        
    double Tliq;                        
    double phi=0, phib=0;               // melt fraction at T,P for ma(Tm), ma(Tm+dTm) calculation
    double phi_eff;                      // effective melt fraction
    double V;                               // Volume at R
    double DT =0, DT_unl =0, DT_unlb = 0;   

    if(std::get<0>(melt_param) == 1){

    
    double frac_unlin_max = std::get<1>(melt_param);
    double x0 = std::get<2>(melt_param);
    double x1 = std::get<3>(melt_param);
    double x2 = std::get<4>(melt_param);
    double x3 = std::get<5>(melt_param);
    double x4 = std::get<6>(melt_param);
    double x5 = std::get<7>(melt_param);
    double x6 = std::get<8>(melt_param); 
   
    if(x0 == 1){    // constant melt production
        
        for(int i = 0; i < N-1; ++i){

        R = Rmax+i*dr;       // Radius of the Volume i
        P = (P0  + rho_m*gl*(Rm-R))/1E9 ;  // Pressure at R
        Tsol=std::get<0>(solidus)+std::get<1>(solidus)*P+std::get<2>(solidus)*P*P+std::get<3>(solidus)*P*P*P+DTsolidus;  // Solidus temperature
        Tliq=std::get<0>(liquidus)+std::get<1>(liquidus)*P+std::get<2>(liquidus)*P*P+std::get<3>(liquidus)*P*P*P;  // liquidus temperature
        if(R <= Rdblc){ T=Tm+am*gc/Cm*Tm*(Rdblc-R); Tb=Tmb+am*gc/Cm*Tmb*(Rdblc-R);} else if(R > Rdblc && R < Rl) {T=Tm+(Tl-Tm)/dBLC*(R-Rdblc); Tb=Tmb+(Tlb-Tmb)/dBLC*(R-Rdblc);} else {T=0;}
        DT=Tliq-Tsol;   
        V=4*M_PI*R*R*dr;
        
        if(T > Tsol){
            
            DT_unl=(T-Tsol)/(DT*frac_unlin_max);
            phi = (T-Tsol)/DT;
            phi = std::min(phi,1.0);
            phi_eff= std::max(phi-phi_crit,0.0); 
            
            Va = Va + V;
            Vmelt = Vmelt + V * phi;
            Vmelt_eff = Vmelt_eff + V * phi_eff;

            if(Tb > Tsol){
            
            DT_unlb=(Tb-Tsol)/(DT*frac_unlin_max);
            phib = (Tb-Tsol)/DT;
            phib = std::min(phib,1.0);
            Vab = Vab + V;
            Vmeltb = Vmeltb + V * phib;   
        }
       
            
        }


        // std::cout << "melt : " << T <<  ", " << Tsol << ", " << P << ", " << Tm << ", " << Vmelt << ", " << phi <<  ", " << R <<  ", " << Rdblc <<  ", " <<  Rmax <<  ", " << Rl << std::endl;
    }


    } else{ // non constant melt production

    
    for(int i = 0; i < N-1; ++i){

        R = Rmax+i*dr;
        P = (P0  + rho_m*gl*(Rm-R))/1E9 ;
        Tsol=std::get<0>(solidus)+std::get<1>(solidus)*P+std::get<2>(solidus)*P*P+std::get<3>(solidus)*P*P*P+DTsolidus;
        Tliq=std::get<0>(liquidus)+std::get<1>(liquidus)*P+std::get<2>(liquidus)*P*P+std::get<3>(liquidus)*P*P*P;
        if(R <= Rdblc){ T=Tm+am*gc/Cm*Tm*(Rdblc-R); Tb=Tmb+am*gc/Cm*Tmb*(Rdblc-R);} else if(R > Rdblc && R < Rl) {T=Tm+(Tl-Tm)/dBLC*(R-Rdblc); Tb=Tmb+(Tlb-Tmb)/dBLC*(R-Rdblc);} else {T=0;}
        DT=Tliq-Tsol;
        V=4*M_PI*R*R*dr;
        

        if(T > Tsol){
            
            DT_unl=(T-Tsol)/(DT*frac_unlin_max);
            phi = std::min(x6*pow(DT_unl,6.) +x5*pow(DT_unl,5.) +x4*pow(DT_unl,4.) +x3*pow(DT_unl,3.) +x2*pow(DT_unl,2.) +x1*DT_unl +x0, (T-Tsol)/DT);
            phi = std::min(phi,1.0);
            phi_eff= std::max(phi-phi_crit,0.0);

            Va = Va + V;
            Vmelt = Vmelt + V * phi;
            Vmelt_eff = Vmelt_eff + V * phi_eff;

            if(Tb > Tsol){
            
            DT_unlb=(Tb-Tsol)/(DT*frac_unlin_max);
            phib = std::min(x6*pow(DT_unlb,6.) +x5*pow(DT_unlb,5.) +x4*pow(DT_unlb,4.) +x3*pow(DT_unlb,3.) +x2*pow(DT_unlb,2.) +x1*DT_unlb +x0, (Tb-Tsol)/DT);
            phib = std::min(phib,1.0);
            
            Vab = Vab + V;
            Vmeltb = Vmeltb + V * phib;   
            
        }
        }

        
        
        // std::cout << "melt : " << T <<  ", " << Tsol << ", " << P << ", " << Tm << ", " << Vmelt << ", " << phi <<  ", " << R <<  ", " << Rdblc <<  ", " <<  Rmax <<  ", " << Rl << std::endl;
    }
    }
    }
    else{
        //KHATZ et al 2003

            double Mcpx = std::get<1>(melt_param);
            double beta1 = std::get<2>(melt_param);
            double beta2 = std::get<3>(melt_param);
            double r0 = std::get<4>(melt_param);
            double r1 = std::get<5>(melt_param);
            double gamma = std::get<6>(melt_param);
            double K = std::get<7>(melt_param);
            double XH2O = std::get<8>(melt_param) * LMBD_liq; 
            double DTH2O = K * pow(XH2O,gamma);

            double Rcpx; double Fcpx; double Tcpx; double Tprim;
            double Tliq_lherz;
           

           for(int i = 0; i < N-1; ++i){

            R = Rmax+i*dr;       // Radius of the Volume i
            P = (P0  + rho_m*gl*(Rm-R))/1E9 ;  // Pressure at R
            Tsol=std::get<0>(solidus)+std::get<1>(solidus)*P+std::get<2>(solidus)*P*P+std::get<3>(solidus)*P*P*P + DTsolidus - DTH2O;  // Solidus temperature
            Tliq=std::get<0>(liquidus)+std::get<1>(liquidus)*P+std::get<2>(liquidus)*P*P+std::get<3>(liquidus)*P*P*P;  // liquidus temperature
            if(R <= Rdblc){ T=Tm+am*gc/Cm*Tm*(Rdblc-R); Tb=Tmb+am*gc/Cm*Tmb*(Rdblc-R);} else if(R > Rdblc && R < Rl) {T=Tm+(Tl-Tm)/dBLC*(R-Rdblc); Tb=Tmb+(Tlb-Tmb)/dBLC*(R-Rdblc);} else {T=0;}
            DT=Tliq-Tsol;   
            V=4*M_PI*R*R*dr;
            if(T > Tsol){

                Rcpx = r0 * r1 * P;
                Fcpx = Mcpx / Rcpx ; 
                Tliq_lherz = Tliq - DT/2;
                Tcpx = pow(Fcpx,beta1)*(Tliq_lherz - Tsol) + Tsol;
                Tprim = (T - Tsol) / (Tliq_lherz - Tsol);
                phi = pow(Tprim,beta1);
                if(phi > Fcpx){
                    phi = Fcpx + (1-Fcpx) * pow((T-Tcpx)/(Tliq-Tcpx),beta2);
                }
                phi = std::min(phi,1.0);
                phi_eff= std::max(phi-phi_crit,0.0); 
                
                Va = Va + V;
                Vmelt = Vmelt + V * phi;
                Vmelt_eff = Vmelt_eff + V * phi_eff;

                if(Tb > Tsol){
                
                Rcpx = r0 * r1 * P;
                Fcpx = Mcpx / Rcpx ; 
                Tliq_lherz = Tliq - DT/2;
                Tcpx = pow(Fcpx,beta1)*(Tliq_lherz - Tsol) + Tsol;
                Tprim = (Tb - Tsol) / (Tliq_lherz - Tsol);
                phib = pow(Tprim,beta1);
                if(phi > Fcpx){
                    phib = Fcpx + (1-Fcpx) * pow((Tb-Tcpx)/(Tliq-Tcpx),beta2);
                }
                phib = std::min(phib,1.0);
                Vab = Vab + V;
                Vmeltb = Vmeltb + V * phib;
            }
            }


    }

    }

    if(Va > 0)
        {
            ma = Vmelt/Va ; 
            ma_eff = Vmelt_eff/Va;
            
        } 
    else
    {
            ma = 0. ;
            ma_eff = 0. ;
        }

     
       if(Vab > 0.)
       {
            mab = Vmeltb/Vab ; 
            
        } 
        else{
            mab = 0. ;  
        }
    }
    

    else {
        ma = 0.;
        Va = 0.;
        ma_eff = 0.;
        mab= 0.;
    }
   
    dmadtm = (mab-ma)/dTm;
    


}


// epsilon = Tmoy / Tm : analityc solution
double epsilon(double const &Tm, double const&Tc, double const &Tl, double const &Tb,double const &Rc, double const &Rl,double const &dBLH, double const &dBLC){

    double Rm = Rl - dBLC;
    double Rb = Rc + dBLH;

    double a1 =(Tb-Tc)/dBLH;
    double b1 = Tc - a1*Rc;

    double a2 =(Tm-Tb)/(Rm-Rb);
    double b2 = Tb - a2*Rb;

    double a3 =(Tl-Tm)/dBLC;
    double b3 = Tl - a3*Rl;
    
    double Tmoy = 3./((Rl*Rl*Rl-Rc*Rc*Rc)) * (a1/4.*(pow(Rb,4.)-pow(Rc,4.))+b1/3*(pow(Rb,3.)-pow(Rc,3.)) +a2/4.*(pow(Rm,4.)-pow(Rb,4.))+b2/3.*(pow(Rm,3.)-pow(Rb,3.)) + a3/4.*(pow(Rl,4.)-pow(Rm,4.))+b3/3.*(pow(Rl,3.)-pow(Rm,3.)));
    double epsi = Tmoy/Tm;

    return epsi;

}


// calculation of the derivative dXdT for k1 and Calculation at t of some variable (Tl,Tb, dBLC ...)

void derivate(double &dtmdt,double &dtcdt, double &dDcdt_N, double &dDcdt_S, double &dDldt_N, double &dDldt_S,
double const &t0,double const &Tm, double const &Tc,double const &f,double const &Rl_N,double const &Rl_S,double const &Dc_N,double const &Dc_S,
double const &Rp, double const &Rc, double const &Vc, double const &Ac, double const &Ts,double const &gl,double const &gc,double const &Tb0,double const &a_m,double const &k_m,double const &C_m,double const &C_c,
double const &C_cr,double const &rho_c,double const &rho_m,double const &rho_cr,double const &L,double const &LMBD_ath, double const &LMBD_liq_N, double const &LMBD_liq_S, double const &Dref,
std::tuple<double,double,double,double,double,double> const &solidus, std::tuple<double,double,double,double,double,double> const &liquidus, std::tuple<double,double,double,double,double,double,double,double> const &melt,
std::tuple<double,double,double,double,double,double,double,double,double,double,double> const &rheology, Eigen::MatrixXd const &RAD,
std::tuple<double,double,double> const &Pm, double const &an_s, double const &ql_N, double const &ql_S, double const &epsi_c, std::tuple<bool,bool> contact, bool const &melt_bol, double const &adv_interf_bool,
std::tuple<bool,double,double,double,double,double,double,double,double> const &melt_param, int const &N_melt, double const &Phi_vis_N, double const &Phi_vis_S
){

    double D_m=k_m/rho_m/C_m;     //Mantle diffusivity   thermal conductivity / mantle density / mantle heat capacity
    double Qm =0;                 // Mantle Heat production
    double Tl=0;                   // Lid base temperature
    double Tb=Tb0;                 // previous Tb
    
    double Vm;                      // Mantle volume
    Internal_heat(Qm,t0,LMBD_ath,rho_m,an_s,RAD);        // RAD : matrice with the heat production data, LMBD_ath : mantle enrichment factor
    Volume(Vm,Rl_N,Rl_S,Rc,f);                            // Rl_S/N lithosphere radius North - South
    double Am_N=4*M_PI*Rl_N*Rl_N*f;
    double Am_S=4*M_PI*Rl_S*Rl_S*(1-f);
    double Rl = pow((f*pow(Rl_N,3.)+(1-f)*pow(Rl_S,3.)),1./3.);  // Volume average of the lithosphere radius
    double Dc = Rp-pow(f*pow(Rp-Dc_N,3.)+(1-f)*pow(Rp-Dc_S,3.),1./3.);  // Volume average of the crust thickness
    
    double P_mean=rho_cr*gl*Dc+rho_m*gl*(Rp-Rl);          //  Average Pressure at the lithosphere depth
    litho_temp(Tl,Tm,rheology);                           // rheology -> table with the rheology constants

    double DT_u;                                     // temperature jump in the upper boundary layer
    if(Tc > Tb) {DT_u = Tm-Tl+Tc-Tb;} else{ DT_u = Tm-Tl;}             

    double eta_0 = std::get<0>(rheology);
    double A =std::get<1>(rheology);
    double R = std::get<2>(rheology);
    double V = std::get<3>(rheology);
    double Tref =std::get<5>(rheology);
    double Pref =std::get<6>(rheology);
    double beta_u = std::get<7>(rheology);
    double beta_c = std::get<8>(rheology);
    double Ra_crit_u = std::get<9>(rheology);

    // viscosity for the three boundary layers (up N/S and average)
    double eta_u_N=eta_0*std::exp(  (A + std::get<0>(Pm) * V )/ (R  * Tm ) - (A + Pref * V) / ( Tref * R) + std::get<10>(rheology) * Phi_vis_N  ) ;
    double eta_u_S=eta_0*std::exp(  (A + std::get<1>(Pm) * V ) / (R  * Tm ) - (A + Pref * V) / ( Tref * R) + std::get<10>(rheology) * Phi_vis_S ) ;
    double eta_u_avg=eta_0*std::exp(  (A + P_mean * V ) / (R  * Tm ) - (A + Pref * V) / ( Tref * R) ) ;
   
    // Rayleigh number
    double Ra_u_N=a_m*rho_m*gl* DT_u *pow(Rl_N-Rc,3.)/D_m/eta_u_N;
    double Ra_u_S=a_m*rho_m*gl* DT_u *pow(Rl_S-Rc,3.)/D_m/eta_u_S;
   
   // Boundary layer thickness
    double dBLC_N= std::min((Rl-Rc) * pow(Ra_crit_u/Ra_u_N,beta_u),(Rl-Rc)/2.);
    
    double dBLC_S=std::min((Rl-Rc) * pow(Ra_crit_u/Ra_u_S,beta_u),(Rl-Rc)/2.);
    double dBLC_avg=Rl - pow(f*pow((Rl_N-dBLC_N),3.)+(1-f)*pow((Rl_S-dBLC_S),3.),1./3.);
    

    
    double DT_i;   // lid temperature jump 
    if(Tc > Tb  ){ DT_i= (Tm-Ts) + (Tc-Tb);} else {DT_i= (Tm-Ts);}
    double Ra_i = a_m*rho_m*gl* DT_i *pow(Rp-Rc,3.)/D_m/eta_u_avg; // intern Rayeligh
    double Ra_crit_c=0.28*pow(Ra_i,0.21);  // Rayleigh critic at the base of the mantle

    double DT_c = std::abs( Tc-Tm);  // Temperature jump in the lower boundary layer
    double eta_c = eta_0*std::exp(  (A + std::get<2>(Pm) * V ) / (R  * (Tb+Tc)/2 ) - (A + Pref * V) / ( Tref * R)  ) ; // viscosity at the base of the mantle
   
    double Ra_c = a_m* rho_m * gc * DT_c * pow((Rl-Rc),3.)/D_m/eta_c;   // Rayleigh number at the base of the mantle
    
    double dBLH=std::min((Rl-Rc)*pow(Ra_crit_c/Ra_c,beta_c),(Rl-Rc)/3.);  // Hot boundary layer thickness


    Base_temp(Tb,Tm,a_m,C_m,gc,Rl,Rc,dBLC_avg,dBLH);  // New Tb
    

    double St = 0; 

    double qcr_N = 0 ;
    double qcr_S = 0 ;
    double qm_N = 0;
    double qm_S = 0;
    double qc = 0;

    dDcdt_N = 0;
    dDcdt_S = 0;

    double w_N = 0, w_S = 0;

    double Phi_N=0,Phi_eff_N =0 ,Va_N =0 ,dmadtm_N = 0;
    double Phi_S=0,Phi_eff_S =0 ,Va_S=0 ,dmadtm_S = 0;

    double epsi_m = epsilon(Tm,Tc,Tl,Tb,Rc,Rl,dBLH,dBLC_avg);   // ratio Tmoy/Tm
    
    
    if(dBLC_N > 0){qm_N = k_m*(Tm-Tl)/dBLC_N;}              // Convective Heat flux N/S
    if(dBLC_S > 0){qm_S = k_m*(Tm-Tl)/dBLC_S;} 
    if(dBLH > 0){qc = k_m*(Tc-Tb)/dBLH;}                     // CMB heat flux

    
    if(melt_bol == 1) {
        
        melting(Phi_N,Phi_eff_N,Va_N,dmadtm_N,Tm,a_m,C_m,gc,gl,N_melt,Dref,std::get<4>(melt),std::get<5>(melt),rho_cr,rho_m,Dc_N,Dc,Rl_N,Rp,dBLC_N,std::get<3>(melt),std::get<6>(melt), melt_param, LMBD_liq_N,
        solidus,liquidus,rheology);                         // Melt fraction under the northern lithosphere
        
        melting(Phi_S,Phi_eff_S,Va_S,dmadtm_S,Tm,a_m,C_m,gc,gl,N_melt,Dref,std::get<4>(melt),std::get<5>(melt),rho_cr,rho_m,Dc_S,Dc,Rl_S,Rp,dBLC_S,std::get<3>(melt),std::get<6>(melt),
        melt_param, LMBD_liq_S,
        solidus,liquidus,rheology);                             // Melt fraction under the southern lithosphere
        
        double Va_avg=f*Va_N+(1-f)*Va_S;                         // Average melt volume
        double dmadtm=0;               
        if(Va_avg > 0){ dmadtm = (Va_N*dmadtm_N*f+Va_S*dmadtm_S*(1-f))/Va_avg;}  // Average melt production
        St = L*Va_avg*dmadtm/C_m/Vm;                              // Stefan number
       
        double k_phi_N=std::get<0>(melt) * (1 - Phi_eff_N) * pow(Phi_eff_N,3.);      // Mantle permeability N/S
        double k_phi_S=std::get<0>(melt) * (1 - Phi_eff_S) * pow(Phi_eff_S,3.); 

        double rho_diff_N = rho_m - std::get<7>(melt)  * (1+std::get<0>(Pm)/1E9*std::get<1>(melt));    // Delta rho = f(P) N/S
        double rho_diff_S = rho_m - std::get<7>(melt)  * (1+std::get<1>(Pm)/1E9*std::get<1>(melt));
        
        if(rho_diff_N < 0){rho_diff_N = 0;}
        if(rho_diff_S < 0){rho_diff_S = 0;}
       
        w_N = k_phi_N / std::get<2>(melt) * rho_diff_N * gl;                     // Darcy velocity N/S
        w_S = k_phi_S / std::get<2>(melt) * rho_diff_S * gl;                     //

        double u_N = D_m / (Rl-Rc) * pow((Ra_u_N / Ra_crit_u), 2*beta_u) ;             // Mantle velocity N/S
        double u_S = D_m / (Rl-Rc) * pow((Ra_u_S / Ra_crit_u), 2*beta_u) ; 
        
        w_N = std::min(w_N, u_N * Phi_eff_N);                                           
        w_S = std::min(w_S, u_S * Phi_eff_S);
        
        qcr_N = rho_cr*w_N*(L+C_cr*(Tm-Tl)) ;                                               // magma Heat flux N/s
        qcr_S = rho_cr*w_S*(L+C_cr*(Tm-Tl)) ;
        dDcdt_N = w_N*pow(Rl_N,2.)/pow(Rp-Dc_N,2.);    // Crust growth N
        dDcdt_S = w_S*pow(Rl_S,2.)/pow(Rp-Dc_S,2.);  // Crust growth S
        
        if(std::get<0>(contact) == 0){         // The crust growth controll the lid growth in case of contact

            dDldt_N = (-qm_N-ql_N)/(rho_m*(C_m*(Tm-Tl)+L*Phi_N)) + w_N * adv_interf_bool;

        } else {
                                           
            dDldt_N = std::max((-qm_N-ql_N)/(rho_m*(C_m*(Tm-Tl)+L*Phi_N)) + w_N * adv_interf_bool,w_N);
        }

        if(std::get<1>(contact) == 0){
            
            
            dDldt_S = (-qm_S-ql_S)/(rho_m*(C_m*(Tm-Tl)+L*Phi_S)) + w_S * adv_interf_bool;

        } else {

            dDldt_S = std::max((-qm_S-ql_S)/(rho_m*(C_m*(Tm-Tl)+L*Phi_S)) + w_S * adv_interf_bool,w_S);
        
        }
        
    }

    dtmdt = (-(qm_N+qcr_N)*Am_N-(qm_S+qcr_S)*Am_S+Qm*Vm+qc*Ac)/(C_m*rho_m*Vm*epsi_m*(St+1));    
    dtcdt = (-qc * Ac)/(C_c*rho_c*Vc*epsi_c);


}

// Same function but with less return and calculation for k2,k3,k4 

void derivate(double &dtmdt,double &dtcdt, double &dDcdt_N, double &dDcdt_S, double &dDldt_N, double &dDldt_S,
double const &Qm,double const &Tm, double const &Tc,double const &f,double const &Rl_N,double const &Rl_S,double const &Dc_N,double const &Dc_S,
double const &Rp, double const &Rc, double const &Vc, double const &Ac, double const &Ts,double const &gl,double const &gc,double &Tb,double const &a_m,double const &k_m,double const &C_m,double const &C_c,
double const &C_cr,double const &rho_c,double const &rho_m,double const &rho_cr,double const &L,double const &LMBD_ath, double const &LMBD_liq_N, double const &LMBD_liq_S, double const &Dref,
std::tuple<double,double,double,double,double,double> const &solidus, std::tuple<double,double,double,double,double,double> const &liquidus, std::tuple<double,double,double,double,double,double,double,double> const &melt,
std::tuple<double,double,double,double,double,double,double,double,double,double,double> const &rheology, Eigen::MatrixXd const &RAD,
std::tuple<double,double,double> const &Pm, double const &an_s, double const &ql_N, double const &ql_S, double const &epsi_c, std::tuple<bool,bool> contact, bool const &melt_bol,  double const &adv_interf_bool, std::tuple<bool,double,double,double,double,double,double,double,double> const &melt_param, int const &N_melt, double &Phi_vis_N, double &Phi_vis_S,
double &Tl, double &epsi_m,double &Ra_avg, double &dBLH, double &dBLC_N,  double &dBLC_S, double &dBLC_avg, double &eta_u_N,double &eta_u_S, double &Phi_avg,double &Phi_N,double &Phi_S, double &Phi_eff_N, double &Phi_eff_S, double &Va_avg, double &Va_N,double &Va_S, double &dmadtm, double &qm_N, double &qm_S, double &qcr_N, double &qcr_S, double &qc, double &St
){

    double D_m=k_m/rho_m/C_m;
    
    
    double Vm;
    double Vm_N, Vm_S;
   
    
    Volume(Vm,Rl_N,Rl_S,Rc,f);
    Vm_N = f* Vm ;
    Vm_S = (1-f) * Vm;
    double Am_N=4*M_PI*Rl_N*Rl_N*f;
    double Am_S=4*M_PI*Rl_S*Rl_S*(1-f);
    double Rl = pow((f*pow(Rl_N,3.)+(1-f)*pow(Rl_S,3.)),1./3.);
    double Dc = Rp-pow(f*pow(Rp-Dc_N,3.)+(1-f)*pow(Rp-Dc_S,3.),1./3.);
    
    double P_mean=rho_cr*gl*Dc+rho_m*gl*(Rp-Rl-Dc);
    litho_temp(Tl,Tm,rheology);

    double DT_u = Tm-Tl+Tc-Tb ;
    if(Tc < Tb) { DT_u= Tm-Tl;}

    double A =std::get<1>(rheology);
    double V = std::get<3>(rheology);
    double Tref =std::get<5>(rheology);
    double Pref =std::get<6>(rheology);
    double R = std::get<2>(rheology);
    double eta_0 = std::get<0>(rheology);
    double Ra_crit_u = std::get<9>(rheology);
    double beta_u = std::get<7>(rheology);
    double beta_c = std::get<8>(rheology);

  
    eta_u_N=eta_0*std::exp(  (A + std::get<0>(Pm) * V )/ (R  * Tm ) - (A + Pref * V) / ( Tref * R) + std::get<10>(rheology) * Phi_vis_N ) ;
    eta_u_S=eta_0*std::exp(  (A + std::get<1>(Pm) * V ) / (R  * Tm ) - (A + Pref * V) / ( Tref * R) + std::get<10>(rheology) * Phi_vis_S ) ;
    double eta_u_avg=eta_0*std::exp(  (A + P_mean * V ) / (R  * Tm ) - (A + Pref * V) / ( Tref * R)  ) ;
   
    double Ra_u_N=a_m*rho_m*gl* DT_u *pow(Rl_N-Rc,3.)/D_m/eta_u_N;
    double Ra_u_S=a_m*rho_m*gl* DT_u *pow(Rl_S-Rc,3.)/D_m/eta_u_S;
   
    Ra_avg = a_m*rho_m*gl* DT_u *pow(Rl-Rc,3.)/D_m/eta_u_avg;
    
    dBLC_N= std::min((Rl-Rc) * pow(Ra_crit_u/Ra_u_N,beta_u),(Rl-Rc)/2.);
    dBLC_S=std::min((Rl-Rc) * pow(Ra_crit_u/Ra_u_S,beta_u),(Rl-Rc)/2.);
    dBLC_avg= Rl - pow(f*pow((Rl_N-dBLC_N),3.)+(1-f)*pow((Rl_S-dBLC_S),3.),1./3.);
    
    double DT_i;
    if(Tc > Tb ){ DT_i= (Tm-Ts) + (Tc-Tb);} else {DT_i= (Tm-Ts);}
    double Ra_i = a_m*rho_m*gl* DT_i *pow(Rp-Rc,3.)/D_m/eta_u_avg;
    double Ra_crit_c=0.28*pow(Ra_i,0.21);

    double DT_c = std::abs(Tc-Tm);
    
    double eta_c = eta_0*std::exp(  (A + std::get<2>(Pm) * V ) / (R  * (Tb+Tc)/2 ) - (A + Pref * V) / ( Tref * R)  ) ;
    double Ra_c = a_m* rho_m * gc * DT_c * pow((Rl-Rc),3.)/D_m/eta_c;
    
    dBLH=std::min((Rl-Rc)*pow(Ra_crit_c/Ra_c,beta_c),(Rl-Rc)/3.);
 
    Base_temp(Tb,Tm,a_m,C_m,gc,Rl,Rc,dBLC_avg,dBLH);


    qcr_N = 0 ;
    qcr_S = 0 ;

    dDcdt_N = 0;
    dDcdt_S = 0;

    double w_N = 0, w_S = 0;

    Va_avg = 0;
    dmadtm = 0;
    Phi_avg = 0;
    

    Phi_N=0, Phi_eff_N =0; double dmadtm_N = 0;
    Phi_S=0, Phi_eff_S =0; double dmadtm_S = 0;

    if(melt_bol == 1) {
        
        melting(Phi_N,Phi_eff_N,Va_N,dmadtm_N,Tm,a_m,C_m,gc,gl,N_melt,Dref,std::get<4>(melt),std::get<5>(melt),rho_cr,rho_m,Dc_N,Dc,Rl_N,Rp,dBLC_N,std::get<3>(melt),std::get<6>(melt),
        melt_param, LMBD_liq_N,
        solidus,liquidus,rheology);
        
        melting(Phi_S,Phi_eff_S,Va_S,dmadtm_S,Tm,a_m,C_m,gc,gl,N_melt,Dref,std::get<4>(melt),std::get<5>(melt),rho_cr,rho_m,Dc_S,Dc,Rl_S,Rp,dBLC_S,std::get<3>(melt),std::get<6>(melt),
        melt_param, LMBD_liq_S,
        solidus,liquidus,rheology);

        Va_avg=f*Va_N+(1-f)*Va_S;
        dmadtm=0;
        if(Va_avg > 0){ dmadtm = (Va_N*dmadtm_N*f+Va_S*dmadtm_S*(1-f))/Va_avg;}
        if(Va_avg > 0){ Phi_avg = (Va_N*Phi_N*f+Va_S*Phi_S*(1-f))/Va_avg;}
        St = L*Va_avg*dmadtm/C_m/Vm;

        double k_phi_N=std::get<0>(melt) * (1 - Phi_eff_N) * pow(Phi_eff_N,3.);
        double k_phi_S=std::get<0>(melt) * (1 - Phi_eff_S) * pow(Phi_eff_S,3.);

        double rho_diff_N = rho_m - std::get<7>(melt)  * (1+std::get<0>(Pm)/1E9*std::get<1>(melt));    // Delta rho = f(P) N/S
        double rho_diff_S = rho_m - std::get<7>(melt)  * (1+std::get<1>(Pm)/1E9*std::get<1>(melt));
        
        if(rho_diff_N < 0){rho_diff_N = 0;}
        if(rho_diff_S < 0){rho_diff_S = 0;}
        
        
        w_N = k_phi_N / std::get<2>(melt) * rho_diff_N * gl;
        w_S = k_phi_S / std::get<2>(melt) * rho_diff_S * gl;

        double u_N = D_m / (Rl-Rc) * pow((Ra_u_N / Ra_crit_u), 2*beta_u) * Phi_eff_N ;
        double u_S = D_m / (Rl-Rc) * pow((Ra_u_S / Ra_crit_u), 2*beta_u) * Phi_eff_S; 

        w_N = std::min(w_N,u_N);
        w_S = std::min(w_S,u_S);
        
        qcr_N = rho_cr*w_N*(L+C_cr*(Tm-Tl)) ;
        qcr_S = rho_cr*w_S*(L+C_cr*(Tm-Tl)) ;

    }

        qm_N = 0;
        qm_S = 0;
        qc = 0;
    
    
        if(dBLC_N > 0){qm_N = k_m*(Tm-Tl)/dBLC_N;} 
        if(dBLC_S > 0){qm_S = k_m*(Tm-Tl)/dBLC_S;} 
        // std::cout << "qmN : " << qm_N << ", qmS : " << qm_S << std::endl;
        if(dBLH > 0){qc = k_m*(Tc-Tb)/dBLH;}

        dDcdt_N = w_N*pow(Rl_N,2.)/pow(Rp-Dc_N,2.);                     // Crust growth N
        dDcdt_S = w_S*pow(Rl_S,2.)/pow(Rp-Dc_S,2.);                     // Crust growth S

       
       
        if(std::get<0>(contact) == 0){         // The crust growth growth is the lid thickness is thicker than the crust
            
                                           
            dDldt_N = (-qm_N-ql_N)/(rho_m*(C_m*(Tm-Tl)+L*Phi_N)) + w_N * adv_interf_bool;

        } else {
            
            dDldt_N = std::max((-qm_N-ql_N)/(rho_m*(C_m*(Tm-Tl)+L*Phi_N)) + w_N * adv_interf_bool,w_N);
        }

        if(std::get<1>(contact) == 0){
            
            dDldt_S = (-qm_S-ql_S)/(rho_m*(C_m*(Tm-Tl)+L*Phi_S)) + w_S * adv_interf_bool;

        } else {
            dDldt_S = std::max((-qm_S-ql_S)/(rho_m*(C_m*(Tm-Tl)+L*Phi_S)) + w_S * adv_interf_bool,w_S);
        
        }
    
    Phi_vis_N = Phi_N*Va_N / Vm_N;
    Phi_vis_S = Phi_S*Va_S / Vm_S;
    // std::cout <<" Phi vis N = " <<Phi_vis_N << ", phi vis S = " << Phi_vis_S << std::endl;
    epsi_m = epsilon(Tm,Tc,Tl,Tb,Rc,Rl,dBLH,dBLC_avg);
    dtmdt = (-(qm_N+qcr_N)*Am_N-(qm_S+qcr_S)*Am_S+Qm*Vm+qc*Ac)/(C_m*rho_m*Vm*epsi_m*(St+1));
    dtcdt = (-qc * Ac)/(C_c*rho_c*Vc*epsi_c);

}