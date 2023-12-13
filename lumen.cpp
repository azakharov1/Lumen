#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <cstdio>
#include <math.h>
#include <cstdlib>
#include <time.h>
#include <vector>
#include <iomanip>
// #include <tr1/random>
//#include "windows.h"
//#include <unistd.h>
#include <algorithm>
#include <iterator>
#include <ctime>
#include </usr/local/opt/libomp/include/omp.h>

double phi1=2*M_PI/6,phi2=2*M_PI/6,phi0=2*M_PI/6; //angle occupied by a cell in lumen
double R_cell_ini=10.0*pow(10,-6); //ini cell size [m]
double cleft_h=20*pow(10,-9); //cleft width
double RT=8.3145*300;

double c_ext=300;//external concentration of ions [mM]
double p_ext=101325.0; // external pressure [Pa]

double cortex_h = 0.6*pow(10,-6); //cortex thickness [m]
double sigma_act1 =100,sigma_act2 =100; //active cortex contraction  [Pa]
double Kapical1=10000,Kbasal1=10000,  Kapical2=10000,Kbasal2=10000; //effective stiffness of cortex  [Pa]

double dOsm_cr = 40*pow(10,9); //critical osmotic pressure  [Pa]

double A_lum1=pow(cleft_h*2,2)*sin(phi1)/2, A_lum1_ini=A_lum1, A_lum2=pow(cleft_h*2,2)*sin(phi2)/2, A_lum2_ini=A_lum2; // Lumen area
double A_cell1=phi1/2*pow(R_cell_ini+cleft_h,2)-A_lum1, A_cell1_ini=A_cell1, A_cell2=phi2/2*pow(R_cell_ini+cleft_h,2)-A_lum2, A_cell2_ini=A_cell2; // Cell area

double beta1,beta2;//apical side curvature angle

double L_a_ini1,L_a_ini2;//ini apical length
double L_cleft_ini=R_cell_ini;//ini cleft length
double L_b_ini1=phi0*(R_cell_ini+cleft_h),L_b_ini2=phi0*(R_cell_ini+cleft_h);//ini basal length
double apical_growth1=1.0, apical_growth2=1.0;

double Osm_ext=RT*c_ext;// external osmotic pressure

double cleft_perm=1.00; //cleft permeability factor

double Diff1=cleft_perm*5.00*pow(10,-9),Diff2=5.00*pow(10,-9);;//ion diffusion through cleft
double omega1=1.0*pow(10,-9),omega2=1.0*pow(10,-9);//cell membrane ion permeability

double L_pc1=cleft_perm*2.0*pow(10,-14),L_pc2=2.0*pow(10,-14); //cleft water permeability
double L_pm1=5.0*pow(10,-12),L_pm2=5.0*pow(10,-12); //cell membrane water permeability

//active pumping rates
double active_basal1=5.0*pow(10,-6), active_basal2=5.00*pow(10,-6), active_bas1=0.5*pow(10,-6), active_bas2=0.5*pow(10,-6);
double active_apical1=1.0*pow(10,-6),active_apical2=1.0*pow(10,-6), active_apic1=0.1*pow(10,-6), active_apic2=0.1*pow(10,-6);

double E_gel1=200,E_gel2=200; //hydrogel elastic modulus [Pa]

double Lumen1_radius()
{
    double R = sqrt(A_lum1 / (sin(phi1)/2 - pow(sin(phi1/2)/sin(beta1),2)*(beta1-sin(2*beta1)/2)));
	return R;
}
double Lumen2_radius()
{
    double R = sqrt(A_lum2 / (sin(phi2)/2 - pow(sin(phi2/2)/sin(beta2),2)*(beta2-sin(2*beta2)/2)));
	return R;
}

double Basal1_radius()
{
    double R = sqrt(2*(A_lum1+A_cell1)/phi1);
	return R;
}
double Basal2_radius()
{
    double R = sqrt(2*(A_lum2+A_cell2)/phi2);
	return R;
}
double Apical_radius1()
{
    double R = Lumen1_radius()*sin(phi1/2)/sin(beta1);
	return R;
}
double Apical_radius2()
{
    double R = Lumen2_radius()*sin(phi2/2)/sin(beta2);
	return R;
}

double Cleft_length1()
{
    double L = Basal1_radius() - Lumen1_radius();
	return L;
}
double Cleft_length2()
{
    double L = Basal2_radius() - Lumen2_radius();
	return L;
}

double Apical1_length()
{
    double L = 2*beta1*Lumen1_radius()*sin(phi1/2)/sin(beta1);
	return L;
}
double Apical2_length()
{
    double L = 2*beta2*Lumen2_radius()*sin(phi2/2)/sin(beta2);
	return L;
}

double Basal1_length()
{
    double L = phi1*Basal1_radius();
	return L;
}
double Basal2_length()
{
    double L = phi2*Basal2_radius();
	return L;
}

double IonFluxBasal1(double c_cell1)
{
    double dOsm=RT*(c_ext-c_cell1);
    double F = Basal1_length()*(omega1*dOsm + active_basal1);
    return F;
}
double IonFluxBasal2(double c_cell2)
{
    double dOsm=RT*(c_ext-c_cell2);
    double F = Basal2_length()*(omega2*dOsm + active_basal2);
    return F;
}

double IonFluxApical1(double c_cell1, double c_lum1)
{
    double dOsm=RT*(c_cell1 - c_lum1);
    double F = Apical1_length()*(omega1*dOsm + active_apical1);
    return F;
}
double IonFluxApical2(double c_cell2, double c_lum2)
{
    double dOsm=RT*(c_cell2 - c_lum2);
    double F = Apical2_length()*(omega2*dOsm + active_apical2);
    return F;
}

double IonCleftLeak1(double c_lum1)
{
    double F = Diff1*cleft_h*(c_lum1-c_ext)/Cleft_length1();
    return F;
}
double IonCleftLeak2(double c_lum2)
{
    double F = Diff2*cleft_h*(c_lum2-c_ext)/Cleft_length2();
    return F;
}

double IonFluxCleft1(double c_cell1, double c_lum1)
{
    double c_cleft=(c_lum1+c_ext)/2.0;
    double F = Cleft_length1()*omega1*RT*(c_cell1-c_cleft);
    return F;
}
double IonFluxCleft2(double c_cell2, double c_lum2)
{
    double c_cleft=(c_lum2+c_ext)/2.0;
    double F = Cleft_length2()*omega2*RT*(c_cell2-c_cleft);
    return F;
}

double P_cell1()
{
    double strain=Basal1_length()/(L_b_ini1)-1.0;
    double strain_g=Basal1_radius()/(R_cell_ini+cleft_h)-1.0;
    double p = p_ext + 2*cortex_h/Basal1_radius()*(Kbasal1/2*strain+sigma_act1)+E_gel1*strain_g;
    return p;
}
double P_cell2()
{
    double strain=Basal2_length()/(L_b_ini2)-1.0;
    double strain_g=Basal2_radius()/(R_cell_ini+cleft_h)-1.0;
    double p = p_ext + 2*cortex_h/Basal2_radius()*(Kbasal2/2*strain+sigma_act2)+E_gel2*strain_g;
    return p;
}

double WaterFluxBasal1(double c_cell1)
{
    double F = Basal1_length()*L_pm1*((p_ext-P_cell1()) - RT*(c_ext-c_cell1));
    return F;
}
double WaterFluxBasal2(double c_cell2)
{
    double F = Basal2_length()*L_pm2*((p_ext-P_cell2()) - RT*(c_ext-c_cell2));
    return F;
}

double P_lum1()
{
    double strain=Apical1_length()/(apical_growth1*L_a_ini1)-1.0;
    double p = P_cell1() - 2*cortex_h/Apical_radius1()*(Kapical1/2*strain+sigma_act1);
    return p;
}
double P_lum2()
{
    double strain=Apical2_length()/(apical_growth2*L_a_ini2)-1.0;
    double p = P_cell2() - 2*cortex_h/Apical_radius2()*(Kapical2/2*strain+sigma_act2);
    return p;
}

double WaterFluxApical1(double c_cell1, double c_lum1)
{
    double F = Apical1_length()*L_pm1*((P_cell1()-P_lum1()) - RT*(c_cell1-c_lum1));
    return F;
}
double WaterFluxApical2(double c_cell2, double c_lum2)
{
    double F = Apical2_length()*L_pm2*((P_cell2()-P_lum2()) - RT*(c_cell2-c_lum2));
    return F;
}

double WaterCleftLeak1(double c_lum1)
{
    double F = L_pc1*cleft_h*((P_lum1()-p_ext)-RT*(c_lum1-c_ext))/Cleft_length1();
    return F;
}
double WaterCleftLeak2(double c_lum2)
{
    double F = L_pc2*cleft_h*((P_lum2()-p_ext)-RT*(c_lum2-c_ext))/Cleft_length2();
    return F;
}

double WaterFluxCleft1(double c_cell1, double c_lum1)
{
    double p_cleft=(P_lum1()+p_ext)/2.0;
    double c_cleft=(c_lum1+c_ext)/2.0;
    double F = Cleft_length1()*L_pm1*((P_cell1()-p_cleft) - RT*(c_cell1-c_cleft));
    return F;
}
double WaterFluxCleft2(double c_cell2, double c_lum2)
{
    double p_cleft=(P_lum2()+p_ext)/2.0;
    double c_cleft=(c_lum2+c_ext)/2.0;
    double F = Cleft_length2()*L_pm2*((P_cell2()-p_cleft) - RT*(c_cell2-c_cleft));
    return F;
}

double f(double x)
{
    double f= sqrt(2*(A_lum2+A_cell2)/x) - sqrt(A_lum2 / (0.5*sin(x) - pow(sin(x/2)/sin(beta2),2)*(beta2-0.5*sin(2*beta2)))) - Cleft_length1();
	return f;
}

double secant(double x1, double x2, double E)
{
	double n = 0, xm, x0, c;
	if (f(x1) * f(x2) < 0) {
		do {
			x0 = (x1 * f(x2) - x2 * f(x1)) / (f(x2) - f(x1));
			c = f(x1) * f(x0);
			x1 = x2; x2 = x0;
			n++;
			if (c == 0) break;
			xm = (x1 * f(x2) - x2 * f(x1)) / (f(x2) - f(x1));
		} while (fabs(xm - x0) >= E);
        return x0;
	} else
    {
        std::cout << "Can not find a root in the given interval"<<std::endl;
        std::cout << " A_lum2:"<<A_lum2<< " A_cell2:"<<A_cell2<< " beta2:"<<beta2<< " Lc:"<<Cleft_length1()<<std::endl;
        std::cin.get();
        return 0;

    }
}

    
int main()
{

    FILE *fp = NULL;
    char dfname[80];
        sprintf(dfname, "Data_c1.dat");
        fp = fopen(dfname,"w"); 
        // fprintf(fp, "#time #c_cell1/c_ext #c_lum1/c_ext #A_cell1/A_cell1_ini #A_lum1/A_lum1_ini #A_lum1/A_cell1_ini #Lumen1_radius #Cleft_length1 #P_cell1 #P_lum1 #Apical1_length #Basal1_length #new_beta1\n");
        fclose(fp);


        sprintf(dfname, "Data_c2.dat");
        fp = fopen(dfname,"w");
        fclose(fp);

int Nsteps=0;

    beta1=M_PI*0.01/180;
    beta2=M_PI*0.01/180;

double new_beta1=0.0, adhesion=1.0*pow(10,-10);// M_PI/6;//

double dy=0.01, totA;


std::cout<<" phi1: "<<phi1<<" ("<<phi1/M_PI*180<<" grad)"<<" beta1: "<<beta1<<" ("<<beta1/M_PI*180<<" grad)"<<std::endl;

double c_cell1=1.0*c_ext,c_cell2=1.0*c_ext,c_lum1=c_ext,c_lum2=1.0*c_ext,ions_cell1,ions_cell2, ions_lum1,ions_lum2, d_ion_cell1, d_ion_cell2, d_ion_lum1,d_ion_lum2,  d_w_cell1, d_w_cell2, d_w_lum1, d_w_lum2;


double time=0.0, dt=0.00002, time_max=100*3600;//*dt;

    ions_cell1=c_cell1*A_cell1; ions_cell2=c_cell2*A_cell2;
    ions_lum1=c_lum1*A_lum1; ions_lum2=c_lum2*A_lum2;
     
     Nsteps=0;
     int ouput =1000000, ouputFile =100000;
     double  prev_Lum_R=0.0, Lum_R=0.0, fc=0, new_phi=phi1;
     bool exit=false, check1=true, check2=true;

     L_a_ini1=2*beta1*cleft_h*sin(phi0/2)/sin(beta1);
     L_a_ini2=2*beta1*cleft_h*sin(phi0/2)/sin(beta1);

while(time<=time_max and !exit){
    
             if(time>=0.5*3600 and check1)//abrupt change in values
                {
                    // Diff*=1.25;
                    // L_pc*=2.0;
                    // L_pm2*=1.1;
                     // Diff2*=1.5;
                     // active_basal2*=1.5;
                      // apical_growth2 = 1.2;
                    // apical_growth2 = 1.03;
                      // phi1*=0.9;
                    // omega2 *= 0.9;
                    check1=false;

                }

                if(time>=3*3600 )
                    {
                    // Diff2 = Diff1*(1+0.5*time/3600);
                    // Diff2 = Diff2 + Diff2*0.0001;
                        // Diff2=Diff1*0.5;
                        // L_pc2=0.5*L_pc1;
                    }
                if(time>=0.0*3600 and time<1.0*3600)
                    {
                    active_basal1 = active_bas1*time/3600 /*+ active_bas1*0.01*/;
                    active_apical1 = active_apic1*time/3600;
                    }
                if(time>=0.0*3600 and time<1.0*3600)
                    {
                     active_basal2 = active_bas2*time/3600 /*+ active_bas2*0.01*/;
                     active_apical2 = active_apic2*time/3600;
                       apical_growth2 = 1.0 + 0.2*time/3600;
                    }

            if(time>=1*3600 and check2)//abrupt change in values
                {
                // new_phi=secant(pow(10,-20), 2*M_PI/3, pow(10,-10));
                // phi2=new_phi;
                    // Diff2*=1.5;
                    // active_basal2*=1.5;
                    // apical_growth2 = 1.1;
                check2=false;
                }

    if(Nsteps%ouput==0){
        std::cout<<time/3600<<" dt:"<<dt<<" ("<<Nsteps<<")   A_cell1:"<<A_cell1<<"   A_lum1:"<<A_lum1<<"   P_cell1:"<<P_cell1()<<"   P_lum1:"<<P_lum1()<<"   L_a:"<<Apical1_length()<<"   L_b:"<<Basal1_length()<<"   L_c:"<<Cleft_length1()<<std::endl;

        std::cout<<time/3600<<" dt:"<<dt<<" ("<<Nsteps<<")   A_cell2:"<<A_cell2<<"   A_lum2:"<<A_lum2<<"   P_cell2:"<<P_cell2()<<"   P_lum2:"<<P_lum2()<<"   L_a:"<<Apical2_length()<<"   L_b:"<<Basal2_length()<<"   L_c:"<<Cleft_length2()<<std::endl;
        }


    d_ion_cell1 = IonFluxBasal1(c_cell1) - IonFluxApical1(c_cell1, c_lum1) - 2.0*IonFluxCleft1(c_cell1, c_lum1);
    d_ion_cell2 = IonFluxBasal2(c_cell2) - IonFluxApical2(c_cell2, c_lum2) - 2.0*IonFluxCleft2(c_cell2, c_lum2);


    d_ion_lum1 = IonFluxApical1(c_cell1, c_lum1) - IonCleftLeak1(c_lum1);
    d_ion_lum2 = IonFluxApical2(c_cell2, c_lum2) - IonCleftLeak2(c_lum2);
     
    d_w_cell1 = WaterFluxBasal1(c_cell1) - WaterFluxApical1(c_cell1, c_lum1) - 2.0*WaterFluxCleft1(c_cell1, c_lum1);
    d_w_cell2 = WaterFluxBasal2(c_cell2) - WaterFluxApical2(c_cell2, c_lum2) - 2.0*WaterFluxCleft2(c_cell2, c_lum2);
    
    d_w_lum1 = WaterFluxApical1(c_cell1, c_lum1) - WaterCleftLeak1(c_lum1);
    d_w_lum2 = WaterFluxApical2(c_cell2, c_lum2) - WaterCleftLeak2(c_lum2);
    
    // update number of ions
    ions_cell1 += dt*d_ion_cell1;
    ions_cell2 += dt*d_ion_cell2;

    ions_lum1 += dt*d_ion_lum1;
    ions_lum2 += dt*d_ion_lum2;
    
    if(ions_cell1<=0){ions_cell1=0; exit=true; std::cout<<"C_cell=0"<<std::endl;}
    if(ions_lum1<=0){ions_lum1=0; exit=true;std::cout<<"C_lum=0"<<std::endl;}
    
    // update area
    A_cell1+=dt*d_w_cell1;
    A_cell2+=dt*d_w_cell2;

    A_lum1+=dt*d_w_lum1;
    A_lum2+=dt*d_w_lum2;
    
    if(A_cell1<=0){A_cell1=0; exit=true;std::cout<<"A_cell=0"<<std::endl;}
    if(A_lum1<=0){A_lum1=0; exit=true;std::cout<<"A_lum=0"<<std::endl;}
    
    c_cell1 = ions_cell1/A_cell1;
    c_cell2 = ions_cell2/A_cell2;

    c_lum1 = ions_lum1/A_lum1;
    c_lum2 = ions_lum2/A_lum2;

    
    if(Nsteps%ouput==0)
        {
        std::cout<<"     Ncells:"<<2*M_PI/phi1<<"     c_cell1:"<<c_cell1/c_ext<<"   c_lum1:"<<c_lum1/c_ext<<"   A_cell1:"<<A_cell1/A_cell1_ini<<"   A_lum1:"<<A_lum1/A_cell1_ini<<"   beta1:"<<beta1<<"   R_lum1/R_cell:"<<Lumen1_radius()/R_cell_ini<<"     phi1:"<<phi1<<"     Lc:"<<Cleft_length1()/R_cell_ini<<std::endl;
        std::cout<<"     Ncells:"<<2*M_PI/phi2<<"     c_cell2:"<<c_cell2/c_ext<<"   c_lum2:"<<c_lum2/c_ext<<"   A_cell2:"<<A_cell2/A_cell1_ini<<"   A_lum2:"<<A_lum2/A_cell1_ini<<"   beta2:"<<beta2<<"   R_lum2/R_cell:"<<Lumen2_radius()/R_cell_ini<<"     phi2:"<<phi2<<"     Lc:"<<Cleft_length2()/R_cell_ini<<" Diff2:"<<Diff2<<std::endl;
        // std::cout<<"    cleft:"<<Cleft_length1()/R_cell_ini<<"    tension_c:"<<tension_c<<"   tension_a:"<<tension_a<<"   ratio:"<<tension_c/tension_a<<"   new_beta1:"<<new_beta1<<std::endl;
        std::cout<<std::endl;

        system("gnuplot plots.plt");
        }
    
   if(Nsteps%ouputFile==0){
                sprintf(dfname, "Data_c1.dat");
                fp = fopen(dfname,"a"); 
                fprintf(fp, "%.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %f %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e\n",/*1*/time,  c_cell1/c_ext, c_lum1/c_ext, A_cell1/A_cell1_ini, /*5*/A_lum1/A_lum1_ini, A_lum1/A_cell1_ini, Lumen1_radius()/R_cell_ini, Cleft_length1()/R_cell_ini, P_cell1()/p_ext, /*10*/P_lum1()/p_ext, Apical1_length()/(apical_growth1*L_a_ini1), Basal1_length()/L_b_ini1, beta1*180/M_PI, d_w_cell1/A_cell1, /*15*/d_w_lum1/A_lum1, WaterFluxBasal1(c_cell1),WaterFluxApical1(c_cell1, c_lum1),WaterCleftLeak1(c_lum1),d_ion_cell1,/*20*/d_ion_lum1,phi1,IonFluxBasal1(c_cell1),IonFluxApical1(c_cell1, c_lum1),IonCleftLeak1(c_lum1),/*25*/Basal1_radius()/(R_cell_ini+cleft_h), Apical1_length()/R_cell_ini,WaterFluxCleft1(c_cell1,c_lum1),IonFluxCleft1(c_cell1,c_lum1));
                fclose(fp);

                sprintf(dfname, "Data_c2.dat");
                fp = fopen(dfname,"a");
                // fprintf(fp, "%.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %f %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e\n",/*1*/time,  c_cell2/c_ext, c_lum2/c_ext, A_cell2/A_cell2_ini, /*5*/A_lum2/A_lum2_ini, A_lum2/A_cell2_ini, Lumen2_radius()/R_cell_ini, Cleft_length2()/R_cell_ini, P_cell2()/p_ext, /*10*/P_lum2()/p_ext, Apical2_length()/R_cell_ini, Basal2_length()/R_cell_ini, new_beta1*180/M_PI, d_w_cell2/A_cell2, /*15*/d_w_lum2/A_lum2,WaterFluxBasal2(c_cell2),WaterFluxApical2(c_cell2, c_lum2),WaterCleftLeak1(c_lum2),d_ion_cell2,d_ion_lum2,new_phi);
                fprintf(fp, "%.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %f %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e\n",/*1*/time,  c_cell2/c_ext, c_lum2/c_ext, A_cell2/A_cell2_ini, /*5*/A_lum2/A_lum2_ini, A_lum2/A_cell2_ini, Lumen2_radius()/R_cell_ini, Cleft_length2()/R_cell_ini, P_cell2()/p_ext, /*10*/P_lum2()/p_ext, Apical2_length()/(apical_growth2*L_a_ini2), Basal2_length()/L_b_ini2, beta2*180/M_PI, d_w_cell2/A_cell2, /*15*/d_w_lum2/A_lum2, WaterFluxBasal2(c_cell2),WaterFluxApical2(c_cell2, c_lum2),WaterCleftLeak2(c_lum2),d_ion_cell2,/*20*/d_ion_lum2,phi2,IonFluxBasal2(c_cell2),IonFluxApical2(c_cell2, c_lum2),IonCleftLeak2(c_lum2),/*25*/Basal2_radius()/(R_cell_ini+cleft_h), Apical2_length()/R_cell_ini,WaterFluxCleft2(c_cell2,c_lum2),IonFluxCleft2(c_cell2,c_lum2));
                fclose(fp);

                Lum_R = Lumen1_radius()/R_cell_ini;
                  // if(fabs(prev_Lum_R-Lum_R)/prev_Lum_R<=pow(10,-7)){exit=true; std::cout<<"    dR small -> exit"<<std::endl;}
                prev_Lum_R = Lum_R;

                }
 Nsteps++;   
 time+=dt;
      if(Nsteps==1000000){dt = dt*10;}//time step increase in steady state
      if(Nsteps==11000000){dt = dt*10;}//time step increase in steady state

}

return 0;
}
