#include<iostream>
#include<fstream>
#include<vector>
#include<array>
#include<cmath>
#include<string>
#include "Eigen/Dense"

bool Lecture_param(std::string const & adresse,std::tuple<double,double,double,double,double,double,double,double,double,double,double> &rheology,
                    std::tuple<double,double,double,double,double,double,double,double> &melt,
                    std::tuple<bool,double,double,double,double,double,double,double,double> &melt_param,  
                    std::tuple<double,double,double,double,double,double,double,double,double,double,double,double,double> &thermo,
 double &Tm0,double &Ts,double &dTc0,double &Dl_i,double &Dc_i,double &dDl_i,double &dDc_i,double &LMBD_cr_N,double &LMBD_cr_S,
                   double &t0,double &tstop,double &Rp,double &Rc,double &f,double &Masse_Mars,double &gl,double &gc, double &fmag, double &fbase, double &HPE_factor,
                   double &X1, double &X2, double &Sat_expo, double &CH2O_p, double &K_cr, double &gamma_cr, double &phi_min_cr)
{
    std::ifstream fichier(adresse);
    bool test = true;
    
    
    if(fichier)
    {
        std::string ligne_txt;
        std::string text_ign;
        getline(fichier,ligne_txt);

        fichier.ignore();

        fichier >> text_ign >> std::get<0>(rheology) >> text_ign >> std::get<1>(rheology) >> text_ign >> std::get<2>(rheology) >> text_ign >>
        std::get<3>(rheology) >> text_ign >> std::get<4>(rheology) >> text_ign >> std::get<5>(rheology) >> text_ign >> std::get<6>(rheology) >> text_ign >>
        std::get<7>(rheology) >> text_ign >> std::get<8>(rheology) >> text_ign >> std::get<9>(rheology) >> text_ign >> std::get<10>(rheology);

        fichier.ignore();

        getline(fichier,ligne_txt);

        fichier.ignore();

        fichier >> text_ign >> std::get<0>(melt) >> text_ign >> std::get<1>(melt) >> text_ign >> std::get<2>(melt) >> text_ign >> std::get<3>(melt) >> text_ign >> std::get<4>(melt) >> text_ign >> std::get<5>(melt) >> text_ign >> std::get<6>(melt) >> text_ign >> fmag >>  text_ign >> fbase  >> text_ign >> std::get<7>(melt);
        
        fichier.ignore();

        getline(fichier,ligne_txt);

        fichier.ignore();

        fichier >> text_ign >> std::get<0>(melt_param);

        fichier.ignore();

        getline(fichier,ligne_txt);

        fichier.ignore();

        if (std::get<0>(melt_param) == 1){

        fichier >> text_ign >> std::get<1>(melt_param) >> text_ign >> std::get<2>(melt_param) >> text_ign >> std::get<3>(melt_param) >> text_ign >> std::get<4>(melt_param) >> text_ign >> std::get<5>(melt_param) >> text_ign >> std::get<6>(melt_param) >> text_ign >> std::get<7>(melt_param) >> text_ign >> std::get<8>(melt_param);
        
        fichier.ignore();

        getline(fichier,ligne_txt);

        fichier.ignore();

        fichier >> text_ign >> text_ign >> text_ign >> text_ign >> text_ign >> text_ign >> text_ign >> text_ign >> text_ign >> text_ign >> text_ign >> text_ign >> text_ign >> text_ign;

        }
        else
        {
        
        fichier >> text_ign >> text_ign >> text_ign >> text_ign >> text_ign >> text_ign >> text_ign >> text_ign >> text_ign >> text_ign >> text_ign >> text_ign >> text_ign >> text_ign >> text_ign >> text_ign;

        fichier.ignore();

        getline(fichier,ligne_txt);

        fichier.ignore();

        fichier >> text_ign >> std::get<1>(melt_param) >> text_ign >> std::get<2>(melt_param) >> text_ign >> std::get<3>(melt_param) >> text_ign >> std::get<4>(melt_param) >> text_ign >> std::get<5>(melt_param) >> text_ign >> std::get<6>(melt_param) >> text_ign >> std::get<7>(melt_param);
        std::get<8>(melt_param) = 0;

        }

        fichier.ignore();

        getline(fichier,ligne_txt);

        fichier.ignore();

        fichier >> text_ign >> Rp >> text_ign >> Rc >> text_ign >> f >> text_ign >> Masse_Mars >> text_ign >> gl >> text_ign >> gc;

        fichier.ignore();

        getline(fichier,ligne_txt);

        fichier.ignore();

        fichier >> text_ign >> Tm0 >> text_ign >> dTc0  >> text_ign >> Ts >> text_ign >> LMBD_cr_N >> text_ign >> LMBD_cr_S >> text_ign >> Dl_i 
        >> text_ign >> dDl_i >> text_ign >> Dc_i >> text_ign >> dDc_i >> text_ign >> t0 >> text_ign >> tstop;

        fichier.ignore();

        getline(fichier,ligne_txt);

        fichier.ignore();

        fichier >> text_ign >> std::get<0>(thermo) >> text_ign >> std::get<1>(thermo) >> text_ign >> std::get<2>(thermo)
        >> text_ign >> std::get<3>(thermo) >> text_ign >> std::get<4>(thermo) >> text_ign >> std::get<5>(thermo) >> text_ign >> std::get<6>(thermo)
        >> text_ign >> std::get<7>(thermo) >> text_ign >> std::get<8>(thermo) >> text_ign >> std::get<9>(thermo) >> text_ign >> std::get<10>(thermo) >> text_ign >> std::get<11>(thermo) >> text_ign >> std::get<12>(thermo) >> text_ign >> HPE_factor;

        fichier.ignore();

        getline(fichier,ligne_txt);

        fichier.ignore();

        fichier >> text_ign >> X1 >> text_ign >> X2 >> text_ign >> Sat_expo >> text_ign >> CH2O_p  >> text_ign >> K_cr >> text_ign >> gamma_cr >> text_ign >> phi_min_cr;

        fichier.close();


        if(std::get<0>(melt_param)==0){
        std::get<8>(melt_param)=CH2O_p;
                                        }               


    }

    else
    {
        std::cout << "ERREUR: Impossible d'ouvrir le fichier physical_parameters en lecture." << std::endl;
        test = false;
    }

    // std::cout << std::get<0>(melt_param) << ", " << std::get<1>(melt_param) << ", " << std::get<8>(melt_param) << ", " << Rp << std::endl;
    return test;
}

bool lecture_num(std::string const &adresse,double &dt_min,double &dt_max,double &dt_init,double &DTMAX,double &dr_min, double &dr_max, int &NC,int &N_soliq,int &N_melt,double &t_acc,double &t_acc_m,bool &unlinear_phi,bool &LMBDdt,bool &RK4,bool &Pressure,bool &Steady, bool &ecriture_profil, bool &ecriture_time,bool &ecriture_tech, bool &URG_STOP, bool &adv_lith_bool, bool &adv_interf_bool){

std::ifstream fichier(adresse);
bool test = true;

if(fichier)
    {
        std::string ligne_txt;
        std::string text_ign;
        getline(fichier,ligne_txt);

        fichier.ignore();

        fichier >> text_ign >> dt_min >> text_ign >> dt_max >> text_ign >> dt_init >> text_ign >>
        DTMAX >> text_ign >> dr_min >> text_ign >> dr_max >> text_ign >> NC >> text_ign >> N_soliq >> text_ign >> N_melt >> text_ign >> t_acc >> text_ign >> t_acc_m;

        fichier.ignore();

        getline(fichier,ligne_txt);

        fichier.ignore();

        fichier >> text_ign >> unlinear_phi >> text_ign >> LMBDdt >> text_ign >> RK4 >> text_ign >> Pressure >> text_ign >> Steady >> text_ign >> ecriture_profil >> text_ign >> ecriture_time >> text_ign >> ecriture_tech >> text_ign >> URG_STOP >> text_ign >> adv_lith_bool >> text_ign >> adv_interf_bool;

        
        fichier.close();

    }

    else
    {
        std::cout << "ERREUR: Impossible d'ouvrir le fichier numerical_parameters en lecture." << std::endl;
        test = false;
    }



return test;

}

bool lecture_rad(std::string const &adresse,Eigen::MatrixXd &RAD, double const &t0, double const &HPE_factor){

    std::ifstream fichier(adresse);
    bool test = true;

    if(fichier)
    {
        std::string ligne_txt;
        std::string text_ign;
        getline(fichier,ligne_txt);
        fichier.ignore();
        fichier >> text_ign >> RAD(0,0) >> RAD(1,0) >> RAD(2,0) >> RAD(3,0) >>
        text_ign >> RAD(0,1) >> RAD(1,1) >> RAD(2,1) >> RAD(3,1) >>
        text_ign >> RAD(0,2) >> RAD(1,2) >> RAD(2,2) >> RAD(3,2) >>
        text_ign >> RAD(0,3) >> RAD(1,3) >> RAD(2,3) >> RAD(3,3);
        
        RAD(0,0)=HPE_factor * RAD(0,0)/std::exp(-t0*std::log(2.0)/RAD(2,0));
        RAD(0,1)=HPE_factor * RAD(0,1)/std::exp(-t0*std::log(2.0)/RAD(2,1));
        RAD(0,2)=HPE_factor * RAD(0,2)/std::exp(-t0*std::log(2.0)/RAD(2,2));
        RAD(0,3)=HPE_factor * RAD(0,3)/std::exp(-t0*std::log(2.0)/RAD(2,3));


        
        fichier.close();

    }

    else
    {
        std::cout << "ERREUR: Impossible d'ouvrir le fichier rad en lecture." << std::endl;
        test = false;
    }

return test;
}

bool lecture_sol(std::string const &adresse, std::tuple<double,double,double,double,double,double> &solidus, std::tuple<double,double,double,double,double,double> &liquidus){

    std::ifstream fichier(adresse);
    bool test = true;

    if(fichier)
    {
        std::string ligne_txt;
        std::string text_ign;
        getline(fichier,ligne_txt);
        fichier.ignore();
        fichier >> 
        std::get<0>(solidus) >> std::get<0>(liquidus) >> text_ign >>
        std::get<1>(solidus) >> std::get<1>(liquidus) >> text_ign >>
        std::get<2>(solidus) >> std::get<2>(liquidus) >> text_ign >>
        std::get<3>(solidus) >> std::get<3>(liquidus) >> text_ign >>
        std::get<4>(solidus) >> std::get<4>(liquidus) >> text_ign >>
        std::get<5>(solidus) >> std::get<5>(liquidus);
        
        fichier.close();

    }

    else
    {
        std::cout << "ERREUR: Impossible d'ouvrir le fichier solidus en lecture." << std::endl;
        test = false;
    }

return test;

}



bool lecture_explo(std::string const &adresse,double &k0_init,double &k0_end,double &eta0_init,double &eta0_end, int &N_k0, int &N_eta0){

std::ifstream fichier(adresse);
bool test = true;

if(fichier)
    {
        std::string ligne_txt;
        std::string text_ign;
        getline(fichier,ligne_txt);

        fichier.ignore();

        fichier >> text_ign >> k0_init >> text_ign >> k0_end >> text_ign >> eta0_init >> text_ign >>
        eta0_end >> text_ign >> N_k0 >> text_ign >> N_eta0;
        
        fichier.close();
        
    }

    else
    {
        std::cout << "ERREUR: Impossible d'ouvrir le fichier explo_parameters en lecture." << std::endl;
        test = false;
    }



return test;

}

