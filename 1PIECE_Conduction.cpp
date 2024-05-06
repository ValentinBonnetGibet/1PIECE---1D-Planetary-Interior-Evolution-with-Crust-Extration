#include<iostream>
#include<vector>
#include<array>
#include<cmath>
#include<string>
#include "Eigen/Dense"
#include "unsupported/Eigen/Splines"

class SplineFunction {
public:
  SplineFunction(Eigen::VectorXd const &x_vec,
                 Eigen::VectorXd const &y_vec)
    : x_min(x_vec.minCoeff()),
      x_max(x_vec.maxCoeff()),
      // Spline fitting here. X values are scaled down to [0, 1] for this.
      spline_(Eigen::SplineFitting<Eigen::Spline<double, 1>>::Interpolate(
                y_vec.transpose(),
                 // No more than cubic spline, but accept short vectors.

                std::min<int>(x_vec.rows() - 1, 3),
                scaled_values(x_vec)))
  { }

  double operator()(double x) const {
    // x values need to be scaled down in extraction as well.
    return spline_(scaled_value(x))(0);
  }

private:
  // Helpers to scale X values down to [0, 1]
  double scaled_value(double x) const {
    return (x - x_min) / (x_max - x_min);
  }

  Eigen::RowVectorXd scaled_values(Eigen::VectorXd const &x_vec) const {
    return x_vec.unaryExpr([this](double x) { return scaled_value(x); }).transpose();
  }

  double x_min;
  double x_max;

  // Spline of one-dimensional "points."
  Eigen::Spline<double, 1> spline_;
};

// Analytical solution 
void analytique(Eigen::VectorXd &T, Eigen::VectorXd const &R,Eigen::VectorXd const &dr_r, double &ql, double &qsurf, int const &N, int const &Nm, double const &Rl,double const &Rm, double const &Rp, double const &Tsurf, double const &Tl, double const &P_cr, double const &P_m, double const &k_m, double const &k_cr, double const &rho_cr, double const &Cp_cr,double const &rho_m,double const &Cp_m)
{
    double H_cr = P_cr/k_cr;
    double H_m = P_m/k_m;
    
    double a = -1*Rl;
    double b = -1*Rl*Rl/6*H_m-Tl;
    double c = -1*Rp;
    double d = -1*Rp*Rp/6*H_cr-Tsurf;
    double e = k_cr/k_m;
    double f = Rm*Rm*Rm/3*(H_m-H_cr*k_cr/k_m);
    double g = Rm;
    double h = -Rm*Rm/6*H_m;
    double k = Rm*Rm/6*H_cr;

    double B2 = 1/(1/e/c-1/a+1/e/g-1/g)*(b-d+f/e/c+f/e/g-h-k);
    double B1 = (B2-f)/e;
    double A1 = (f-B2)/(e*c)-d;
    double A2 = -1*B2/a-b;

    
    for (int i = 0; i < N; ++i){


    if (R(i) < Rm )
    {

    T(i)=A2- (B2/R(i))-(R(i)*R(i))*H_m/6;
    
    }

    else 
    {

    T(i)=A1-B1/R(i)-(R(i)*R(i))*H_cr/6.;

    }
    
    }


 ql = -(B2*k_m/(Rl*Rl)-Rl/3.*P_m);
 qsurf = -(B1*k_cr/(Rp*Rp)-Rp/3.*P_cr);

}



// Calculation of the geometry for a Finite volume calculation
void geothermo(std::tuple<int,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,bool> &table, int const &N,
double const &Rp,double const &Rl,double const &Rm,double const &k_cr,double const &k_m,double const &C_cr,double const &C_m,double const &rho_cr,double const &rho_m,double const &dDcrdt,double const &H_cr_rad,double const &H_m, double const &Hb, double const &L_m, double const &L_cr, double const &fbase, double const &fmag, double const &g, double const &Tl, bool const &advection_bool)
{
    double dr_i=(Rp-Rl)/(N-1);
    
    Eigen::VectorXd R(N,1);
    Eigen::VectorXd k_r(N,1);
    Eigen::VectorXd Cp_r(N,1);
    Eigen::VectorXd rho_r(N,1);
    Eigen::VectorXd au_r(N,1);
    Eigen::VectorXd ad_r(N,1);
    Eigen::VectorXd dr_r(N,1);
    Eigen::VectorXd dV_r(N,1);
    Eigen::VectorXd H_r(N,1);
    Eigen::VectorXd Pr(N,1);
    Eigen::VectorXd L_r(N,1);
    Eigen::VectorXd Hb_r(N,1);
    Eigen::VectorXd ur = Eigen::VectorXd::Zero(N);
    
    
    bool contact = 0;
    int Nm = 0;
    double dr_cr = 0;
    double dr_m = 0;
    double Vcr = 4/3*M_PI * (pow(Rp,3.) - pow(Rm,3.));
    double P0 = rho_cr*g*(Rp-Rm);
   

    if(Rm-Rl < dr_i/4)     // if lithospheric mantle thinner than 1/4 of the space step : only crust
    {
        
        Nm=0;   
        dr_cr=(Rp-Rm)/(N-1);
        dr_r(0)=dr_cr/2;
        dr_r(N-1)=dr_cr/2;

        k_r(0)=k_m;
        Cp_r(0)=C_cr;
        rho_r(0)=rho_cr;
        H_r(0)=H_cr_rad;
        Hb_r(0)=0;
        L_r(0) = L_cr;
        contact = 1;
    
    } else if(Rm-Rl < dr_i/2){      // Else if is it's thicker than 1/2 space step :
        
        Nm=1;
        dr_m=Rm-Rl;
        dr_cr=(Rp-Rm)/(N-1.5);
        dr_r(0)=dr_m;
        dr_r(N-1)=dr_cr/2;

        k_r(0)=k_m;
        Cp_r(0)=C_m;
        rho_r(0)=rho_m;
        H_r(0)=H_m;
        Hb_r(0)=0;
        L_r(0) = L_m;
       
   
    } else {   
        
        Nm=std::round((Rm-Rl)/dr_i-0.5) + 1;   // Number of point for the lithospheric mantle
        dr_m=(Rm-Rl)/(Nm-0.5);                  // space step lithospheric mantle
        dr_cr=(Rp-Rm)/(N-Nm-0.5);           // space step crust
        dr_r(0)=dr_m/2;                      // space step of the first and last volume
        dr_r(N-1)=dr_cr/2;

        k_r(0)=k_m;
        Cp_r(0)=C_m;
        rho_r(0)=rho_m;
        H_r(0)=H_m;
        Hb_r(0)=0;
        L_r(0) = L_m;
        
    }
    

    for(int i = 1; i < N-1; ++i){
    if(i < Nm){ 
        dr_r(i)=dr_m;
        k_r(i)=k_m;
        Cp_r(i)=C_m;
        rho_r(i)=rho_m;
        H_r(i)=H_m;
        Hb_r(i) = 0.0;
        L_r(i) = L_m;
    } else {
        dr_r(i)=dr_cr;
        k_r(i)=k_cr;
        Cp_r(i)=C_cr;
        rho_r(i)=rho_cr;
        H_r(i)=H_cr_rad;
        Hb_r(i) = Hb;
        L_r(i) = L_cr;
    }}
    
    dr_r(N-1)=dr_cr/2;
    k_r(N-1)=k_cr;
    Cp_r(N-1)=C_cr;
    rho_r(N-1)=rho_cr;
    H_r(N-1)=H_cr_rad;
    Hb_r(N-1) = Hb;
    L_r(N-1) = L_cr;
  
  
    if(Nm > 0){R(0)=Rl; Pr(0) = P0 + rho_m * g * (Rm-Rl); } else {R(0)=Rm; Pr(0) = P0;}
   
    R(1)=R(0)+dr_r(0)+0.5*dr_r(1);
   
    if(R(1) < Rm){ Pr(1) = P0 + rho_m * g * (Rm - R(1)); }
    else{Pr(1) = rho_cr * g * (Rp-R(1)); }
    
    double dru = 0; // dr+
    double drd = 0; // dr-
    double fu = 0;  
    double fd = 0;
    double ku = 0;  // k+
    double kd = 0;  // k-
    
    dru=dr_r(0)+0.5*dr_r(1);
    ku=1/((1-fu)/k_r(0)+fu/k_r(1));
    au_r(0)=ku/dru;
    ad_r(0) = 0;
    dru=0.5*dr_r(1)+0.5*dr_r(2);
    drd=0.5*dr_r(1)+dr_r(0);
    fu=dr_r(2)/2/dru;
    fd=dr_r(0)/2/drd;
    ku=1/((1-fu)/k_r(1)+fu/k_r(2));
    kd=1/((1-fd)/k_r(1)+fd/k_r(0));
    au_r(1)=ku/dru;   
    ad_r(1)=kd/drd;

   
    for(int i = 2; i < N-2; ++i){

    R(i)=R(i-1)+0.5*dr_r(i-1)+0.5*dr_r(i);
    dru=0.5*dr_r(i)+0.5*dr_r(i+1);
    drd=0.5*dr_r(i)+0.5*dr_r(i-1);
    fu=dr_r(i+1)/2/dru;
    fd=dr_r(i-1)/2/drd;
    ku=1/((1-fu)/k_r(i)+fu/k_r(i+1));
    kd=1/((1-fd)/k_r(i)+fd/k_r(i-1));
    au_r(i)=ku/dru;
    ad_r(i)=kd/drd;
    if(R(i) < Rm){ Pr(i) = P0 + rho_m * g * (Rm - R(i)); }
    else{Pr(i) = rho_cr * g * (Rp-R(i)); }

    }

   
    R(N-2)=R(N-3)+0.5*dr_r(N-3)+0.5*dr_r(N-2);
    if(R(N-2) < Rm){ Pr(N-2) = P0 + rho_m * g * (Rm-R(N-2)); }
    else{Pr(N-2) = rho_cr * g * (Rp-R(N-2)); }
    dru=0.5*dr_r(N-2)+dr_r(N-1);
    drd=0.5*dr_r(N-2)+0.5*dr_r(N-3);
    fu=dr_r(N-1)/2/dru;
    fd=dr_r(N-3)/2/drd;
    ku=1/((1-fu)/k_r(N-2)+fu/k_r(N-1));
    kd=1/((1-fd)/k_r(N-2)+fd/k_r(N-3));
    au_r(N-2)=ku/dru;
    ad_r(N-2)=kd/drd;
    

    
    drd=dr_r(N-1)+0.5*dr_r(N-2);
    fd=dr_r(N-2)/2/drd;
    kd=1/((1-fd)/k_r(N-1)+fd/k_r(N-2));

    au_r(N-1) = 0;
    ad_r(N-1) = kd/drd; 

    R(N-1)=Rp;
    Pr(N-1) = 0.0;
    Pr = Pr / 1E9; 


    
    dV_r(0) = 4./3.*M_PI*(pow(R(0)+dr_r(0),3.0)-pow(R(0),3.0));
    for(int i = 1; i < N-1; i++){
        dV_r(i) = 4./3.*M_PI*(pow(R(i)+dr_r(i)/2.,3.0)-pow(R(i)-dr_r(i)/2.,3.0));
    }
    
    dV_r(N-1) = 4./3.*M_PI*(pow(R(N-1),3.0)-pow(R(N-1)-dr_r(N-1),3.0));


    if(advection_bool){
    for(int i = 0; i < N; i++){
    if(i < Nm) {
        ur(i) = dDcrdt * pow(Rm/(R(i)),2.0);
        }

    else if(i == Nm){
        ur(i) =  dDcrdt;
    }

    else{
        ur(i) = (1-fmag) * dDcrdt *  pow(Rm/(R(i)),2.0) + (1-fbase) * fmag * dDcrdt * pow(Rm/R(i),2.0) * (pow(Rp,3.0)-pow(R(i),3.0))/(pow(Rp,3.0)-pow(Rm,3.0)) ;
           }
    }
    ur(N-1) = (1-fmag) * dDcrdt *  pow(Rm/Rp,2.0);
    }
    
    
    std::get<0>(table)=Nm;
    std::get<1>(table)=R;
    std::get<2>(table)=dr_r;
    std::get<3>(table)=Cp_r;
    std::get<4>(table)=rho_r;
    std::get<5>(table)=k_r;
    std::get<6>(table)=au_r;
    std::get<7>(table)=ad_r;
    std::get<8>(table)=H_r;
    std::get<9>(table)=Pr;
    std::get<10>(table)=dV_r;
    std::get<11>(table) = ur;
    std::get<12>(table) = L_r;
    std::get<13>(table) = Hb_r;
    std::get<14>(table)=contact;
    
  
   
}


void lith_melting(double &Vmeltcr, double &Phicr_avg, double &Phicr_max,double &Vmeltlith, double &Philith_avg, Eigen::VectorXd const &T, Eigen::VectorXd &Phi_r, Eigen::VectorXd &Tsol,Eigen::VectorXd &Tliq, Eigen::VectorXd &DT_sol_cr, Eigen::VectorXd const &Pr,Eigen::VectorXd const &dV_r,
std::tuple<double,double,double,double,double,double> const &solidus, std::tuple<double,double,double,double,double,double> const &liquidus, int const &Nm, int const &N, double const &CH2O_cr, double const &X1, double const X2, double const &lmbd, double const &KH2O_cr, double const &gamma_cr, double const &phi_min_cr, double const &Di)
{
    

    Phi_r = Eigen::VectorXd::Zero(N) ;
    Vmeltcr = 0;
    double VxPhicr_cum  = 0 ;
    double VxPhilith_cum  = 0 ;
    double sat = 0, DTH2O = 0;
    double CH2O_cr_crit = CH2O_cr /  (phi_min_cr+(1-phi_min_cr)*Di);

    for(int i = 0; i < N; ++i)
    {
        if(i < Nm){
        Tsol(i)=std::get<0>(solidus)+std::get<1>(solidus)*Pr(i)+std::get<2>(solidus)*Pr(i)*Pr(i)+std::get<3>(solidus)*Pr(i)*Pr(i)*Pr(i);  // Solidus temperature  without DTsol !!!!!
        Tliq(i)=std::get<0>(liquidus)+std::get<1>(liquidus)*Pr(i)+std::get<2>(liquidus)*Pr(i)*Pr(i)+std::get<3>(liquidus)*Pr(i)*Pr(i)*Pr(i);  // liquidus temperature

        if(T(i)<=Tsol(i))
        {
            Phi_r(i) = 0.0;
        }
        else if(T(i) < Tliq(i))
        {
            Phi_r(i) =  (T(i) - Tsol(i))/(Tliq(i)-Tsol(i));
            Vmeltlith = Vmeltlith + dV_r(i) ;
            VxPhilith_cum = VxPhilith_cum + dV_r(i)* Phi_r(i) ;

        }
        else if(T(i) > Tliq(i))
        {
            Phi_r(i)=1.0;
            Vmeltlith = Vmeltlith + dV_r(i) ;
            VxPhilith_cum = VxPhilith_cum + dV_r(i) ;
        
        }
    

        }

        else{

            sat = X1 * pow(Pr(i),lmbd) * X2 * Pr(i);
            if(sat > CH2O_cr_crit){
                DTH2O =KH2O_cr*pow(CH2O_cr,gamma_cr);
                Tsol(i) = std::get<4>(solidus) + std::get<5>(liquidus) * Pr(i) - DTH2O;
                Tliq(i) = std::get<4>(liquidus) + std::get<5>(liquidus) * Pr(i);

            }
            else{
                DTH2O =KH2O_cr*pow(sat,gamma_cr);
                Tsol(i) = std::get<4>(solidus) + std::get<5>(liquidus) * Pr(i) - DTH2O;
                Tliq(i) = std::get<4>(liquidus) + std::get<5>(liquidus) * Pr(i);


            }

       
            
        if(T(i)<=Tsol(i))
         {
            Phi_r(i) = 0.0;
         }
         else if(T(i) < Tliq(i))
         {
            Phi_r(i) = (T(i) - Tsol(i))/(Tliq(i)-Tsol(i)) + phi_min_cr;
            Vmeltcr = Vmeltcr + dV_r(i) ;
            VxPhicr_cum = VxPhicr_cum + dV_r(i)* Phi_r(i) ;
         }
         else if(T(i) > Tliq(i))
            {
            Phi_r(i)=1.0;
            Vmeltcr = Vmeltcr + dV_r(i) ;
            VxPhicr_cum = VxPhicr_cum + dV_r(i);
            }
        
        }
    DT_sol_cr(i) = Tliq(i) - Tsol(i);
    
     
    }

        if(Vmeltcr > 0){ Phicr_avg = VxPhicr_cum / Vmeltcr;} else {Phicr_avg=0;}
        if(Vmeltlith > 0){Philith_avg = VxPhilith_cum / Vmeltlith;} else{Philith_avg = 0;}
        Phicr_max = Phi_r((Eigen::lastN(N-Nm))).maxCoeff();
    
    
}


// Unsteady conduction
void unsteady(int const &N ,double const &dt,double &ql,double &qsurf,Eigen::VectorXd &T_r,Eigen::VectorXd const &Phi_r, Eigen::VectorXd &DT_sol_cr, std::tuple<int,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,bool> &table, double const &Pcrmag,double const &Vcr, double const &fmag, double const &fbase,double const &Tl)
{
    int Nm = std::get<0>(table);
    Eigen::VectorXd R=std::get<1>(table);
    Eigen::VectorXd dr_r=std::get<2>(table);
    Eigen::VectorXd Cp_r=std::get<3>(table);
    Eigen::VectorXd rho_r=std::get<4>(table);
    Eigen::VectorXd au_r=std::get<6>(table);
    Eigen::VectorXd ad_r=std::get<7>(table);
    Eigen::VectorXd P_r=std::get<8>(table);
    Eigen::VectorXd dV_r=std::get<10>(table);
    Eigen::VectorXd u_r=std::get<11>(table);
    Eigen::VectorXd L_r=std::get<12>(table);
    Eigen::VectorXd Hb_r=std::get<13>(table);
    
    R = R - 0.5 * dr_r;
    double HaT =0, Hb = 0;
    double gamma = 0;
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(N,N);
    Eigen::VectorXd H_r = Eigen::VectorXd::Zero(N);
    
    // Implicite scheme with finite volume and spherical geometry
    // Matrix creation

    for(int i = 1; i < N-1; ++i){
       
        if(i < Nm){
        HaT = 0.0;
        Hb_r(i) = 0;}

        else if(i == Nm){
        HaT = Pcrmag * fmag * ((1-fbase)+ fbase*Vcr/dV_r(Nm));
        Hb_r(i) =  Hb_r(i) *  fmag  *(Cp_r(i)*Tl+L_r(i)-L_r(i)*Phi_r(i))*((1-fbase)+ fbase*Vcr/dV_r(Nm));
        }

        else{ 
        HaT = Pcrmag * fmag * (1-fbase);
        Hb_r(i) = Hb_r(i) * fmag * (1-fbase) * (Cp_r(i)*Tl+L_r(i)-L_r(i)*Phi_r(i));
        }
        
        if(Phi_r(i) < 0 && Phi_r(i) > 1){gamma = rho_r(i)*Cp_r(i);} else {gamma = rho_r(i)*Cp_r(i)*(1+L_r(i)/Cp_r(i)/DT_sol_cr(i));}
   
        M(i,i-1)=-dt/(gamma)*((ad_r(i)-u_r(i)*dr_r(i)/2.0/(dr_r(i+1)+dr_r(i))*rho_r(i)*Cp_r(i))/dr_r(i));
        M(i,i)=1+dt/(gamma)*((ad_r(i)+u_r(i)*dr_r(i)/2.0/(dr_r(i+1)+dr_r(i))*rho_r(i)*Cp_r(i))/dr_r(i)+ (R(i)+dr_r(i))*(R(i)+dr_r(i))/R(i)/R(i)/dr_r(i) * (au_r(i)-u_r(i+1)*dr_r(i+1)/2.0/(dr_r(i+1)+dr_r(i))*rho_r(i+1)*Cp_r(i+1)) - HaT );
        M(i,i+1)=-dt/(gamma)*((R(i)+dr_r(i))*(R(i)+dr_r(i))/R(i)/R(i)/dr_r(i) * (au_r(i)+u_r(i+1)*dr_r(i+1)/2.0/(dr_r(i+1)+dr_r(i))*rho_r(i+1)*Cp_r(i+1)));
        
        H_r(i) = (P_r(i) + Hb_r(i) - u_r(i)*rho_r(i)*L_r(i)/dr_r(i)*Phi_r(i) + (R(i)+dr_r(i))*(R(i)+dr_r(i))/R(i)/R(i)/dr_r(i) * u_r(i)*rho_r(i+1)*L_r(i)*Phi_r(i+1))*dt/(gamma);

    } 

    M(0,0)=1;
    M(N-1,N-1)=1;
    
    qsurf = (T_r(N-2)-T_r(N-1))* ad_r(N-1); // surface Heat flux

    T_r=M.colPivHouseholderQr().solve(T_r+H_r);    // Matrix inversion

    ql = (T_r(1)-T_r(0))*au_r(0); // Heat flux at the base
    
    // std::cout << "T_r(1) = " << T_r(1)<< ", T_r(0) = " << T_r(0) << ", Tl = " << Tl << std::endl;


}


// Interpolation function
bool interpolate(Eigen::VectorXd const &T1, Eigen::VectorXd const &R1, Eigen::VectorXd &T2, Eigen::VectorXd const &R2, double const &Tl, double const &Ts ){
    bool test = 1 ;

    int N1 = R1.rows(); 
    int N2 = R2.rows();
    
    
    
    if (R1(0) < R2(0)){
        
        SplineFunction interpol(R1, T1);
        for(int i = 1; i < N1; ++i){

        
        T2(i) = interpol(R2(i));

            }

        


    } else if(N1 == N2 ) {
        SplineFunction interpol(R1, T1);
        for(int i = 0; i < N1; ++i){

    
         T2(i) = interpol(R2(i));

            }
         
        test = 0;

    } else if(N1 > N2 ) {

        SplineFunction interpol(R1, T1);
        for(int i = 0; i < N2; ++i){

        
         T2(i) = interpol(R2(i));

            }
         
        test = 0;
    } 

    else {
        
        SplineFunction interpol(R1, T1);
        for(int i = 0; i < N1; ++i){

            
         // T2(i+1)    = T1(i);
            T2(i) = interpol(R2(i));
            }

        
        test = 0;
    }

    T2(0) = Tl;   
    T2(N2-1) = Ts;
    
    return test;
}


bool no_interpolate(Eigen::VectorXd const &T1, Eigen::VectorXd &T2, double const &Tl, double const &Ts ){

    bool test = 1 ;

    int N1 = T1.rows(); 
    int N2 = T2.rows();

    if(N1 == N2)
    {

        T2 = T1;

    
    }
    else if(N1 < N2)
    {
       for(int i = 1; i < N2; ++i){

           T2(i) = T1(i-1);

       }
    }
    else
    {
        test = 0;
        std::cout << "interpolation problems" << std::endl;

    }

    T2(0) = Tl;   
    T2(N2-1) = Ts;
    
    return test;
}

// Same function as geothermo but in the case of an acceleration with a higher dr and a constant number of point in the crust
void geothermo_acc(std::tuple<int,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd,bool> &table,int &N,int const &Nc,double const &dr_max,double const &Rp,double const &Rl,double const &Rm,double const &k_cr,double const &k_m,double const C_cr,double const &C_m,double const &rho_cr,double const &rho_m,double const &H_cr,double const &H_m,double const &g)
{

    
    int Nm =std::round((Rm-Rl)/dr_max - 0.5);
    N = Nc + Nm;    

    Eigen::VectorXd R(N,1);
    Eigen::VectorXd k_r(N,1);
    Eigen::VectorXd Cp_r(N,1);
    Eigen::VectorXd rho_r(N,1);
    Eigen::VectorXd au_r(N,1);
    Eigen::VectorXd ad_r(N,1);
    Eigen::VectorXd dr_r(N,1);
    Eigen::VectorXd H_r(N,1);
    Eigen::VectorXd Pr(N,1);
    Eigen::VectorXd dV_r(N,1);
    Eigen::VectorXd L_r(N,1);
    Eigen::VectorXd Hb_r(N,1);
    Eigen::VectorXd ur = Eigen::VectorXd::Zero(N);
    bool contact = 0;

    
    double dr_m=(Rm-Rl)/(Nm-0.5);
    double dr_cr = (Rp-Rm)/(Nc-0.5);
    double P0 = rho_cr*g*(Rp-Rm);

    dr_r(0)=dr_m/2;
    k_r(0)=k_m;
    Cp_r(0)=C_m;
    rho_r(0)=rho_m;
    H_r(0)=H_m;
    

    for(int i = 1; i < N-1; ++i){
    if(i < Nm){
        dr_r(i)=dr_m;
        k_r(i)=k_m;
        Cp_r(i)=C_m;
        rho_r(i)=rho_m;
        H_r(i)=H_m;
    } else {
        dr_r(i)=dr_cr;
        k_r(i)=k_cr;
        Cp_r(i)=C_cr;
        rho_r(i)=rho_cr;
        H_r(i)=H_cr;
    }}
    
    dr_r(N-1)=dr_cr/2;
    k_r(N-1)=k_cr;
    Cp_r(N-1)=C_cr;
    rho_r(N-1)=rho_cr;
    H_r(N-1)=H_cr;
    
    if(Nm > 0){R(0)=Rl; Pr(0) = P0 + rho_m * g * (Rm-Rl); } else {R(0)=Rm; Pr(0) = P0;}

    R(1)=R(0)+dr_r(0)+0.5*dr_r(1);
    if(R(1) < Rm){ Pr(1) = P0 + rho_m * g * (Rm - R(1)); }
    else{Pr(1) = rho_cr * g * (Rp-R(1)); }

    double dru = 0;
    double drd = 0;
    double fu = 0;
    double fd = 0;
    double ku = 0;
    double kd = 0;

    dru=dr_r(0)+0.5*dr_r(1);
    ku=1/((1-fu)/k_r(0)+fu/k_r(1));
    au_r(0)=ku/dru;
    ad_r(0)=0;
    dru=0.5*dr_r(1)+0.5*dr_r(2);
    drd=0.5*dr_r(1)+dr_r(0);
    fu=dr_r(2)/2/dru;
    fd=dr_r(0)/2/drd;
    ku=1/((1-fu)/k_r(1)+fu/k_r(2));
    kd=1/((1-fd)/k_r(1)+fd/k_r(0));
    au_r(1)=ku/dru;
    ad_r(1)=kd/drd;
    

    for(int i = 2; i < N-2; ++i){

    R(i)=R(i-1)+0.5*dr_r(i-1)+0.5*dr_r(i);
    dru=0.5*dr_r(i)+0.5*dr_r(i+1);
    drd=0.5*dr_r(i)+0.5*dr_r(i-1);
    fu=dr_r(i+1)/2/dru;
    fd=dr_r(i-1)/2/drd;
    ku=1/((1-fu)/k_r(i)+fu/k_r(i+1));
    kd=1/((1-fd)/k_r(i)+fd/k_r(i-1));
    au_r(i)=ku/dru;
    ad_r(i)=kd/drd;
    if(R(i) < Rm){ Pr(i) = P0 + rho_m * g * (Rm - R(i)); }
    else{Pr(i) = rho_cr * g * (Rp-R(i)); }

    }
   
    R(N-2)=R(N-3)+0.5*dr_r(N-3)+0.5*dr_r(N-2);
    if(R(N-2) < Rm){ Pr(N-2) = P0 + rho_m * g * (Rm-R(N-2)); }
    else{Pr(N-2) = rho_cr * g * (Rp-R(N-2)); }
    dru=0.5*dr_r(N-2)+dr_r(N-1);
    drd=0.5*dr_r(N-2)+0.5*dr_r(N-3);
    fu=dr_r(N-1)/2/dru;
    fd=dr_r(N-3)/2/drd;
    ku=1/((1-fu)/k_r(N-2)+fu/k_r(N-1));
    kd=1/((1-fd)/k_r(N-2)+fd/k_r(N-3));
    au_r(N-2)=ku/dru;
    ad_r(N-2)=kd/drd;

    drd=dr_r(N-1)+0.5*dr_r(N-2);
    fd=dr_r(N-2)/2/drd;
    kd=1/((1-fd)/k_r(N-1)+fd/k_r(N-2));

    au_r(N-1) = 0;
    ad_r(N-1) = kd/drd; 

    R(N-1)=Rp;
    Pr(N-1) = 0.0;
    Pr = Pr / 1E9; 

    dV_r(0) = 4./3.*M_PI*(pow(R(0)+dr_r(0),3.0)-pow(R(0),3.0));
    for(int i = 1; i < N-1; i++){
        dV_r(i) = 4./3.*M_PI*(pow(R(i)+dr_r(i)/2.,3.0)-pow(R(i)-dr_r(i)/2.,3.0));
    }
    dV_r(N-1) = 4./3.*M_PI*(pow(R(N-1),3.0)-pow(R(N-1)-dr_r(N-1),3.0));


    std::get<0>(table)=Nm;
    std::get<1>(table)=R;
    std::get<2>(table)=dr_r;
    std::get<3>(table)=Cp_r;
    std::get<4>(table)=rho_r;
    std::get<5>(table)=k_r;
    std::get<6>(table)=au_r;
    std::get<7>(table)=ad_r;
    std::get<8>(table)=H_r;
    std::get<9>(table)=Pr;
    std::get<10>(table)=dV_r;
    std::get<11>(table)=ur;
    std::get<12>(table)=L_r;
    std::get<13>(table)=Hb_r;
    std::get<14>(table)=contact;
   
    
}

    
