// ==============================================================================
// This file is part of THOR.
//
//     THOR is free software : you can redistribute it and / or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
//
//     THOR is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
//     GNU General Public License for more details.
//
//     You find a copy of the GNU General Public License in the main
//     THOR directory under <license.txt>.If not, see
//     <http://www.gnu.org/licenses/>.
// ==============================================================================
//
//
//
//
// Description: Hydrogen recombination/dissociation source term
//
//
// Method:
//
// Known limitations: None.
//
// Known issues: None.
//
// If you use this code please cite the following references:
//
//
// Current Code Owner: Russell Deitrick (EEG), russell.deitrick@csh.unibe.ch
//
// History:
// Version Date       Comment
// ======= ====       =======
//
//
////////////////////////////////////////////////////////////////////////

#define Runiv 8.3144621
#define qbond 2.14e8  //hydrogen dissociation energy from Bell & Cowan 2018

// __global__ void dGibbs(double temp, double *dG) {
//   //calculates change in Gibbs free energy for H (polyfit to Heng's Appdx D values)
//   *dG = 2.1370867596206315e-17*temp*temp*temp*temp*temp +
//          -3.8689132818241159e-13*temp*temp*temp*temp +
//          2.7275438366298867e-09*temp*temp*temp +
//          -9.6170574202103724e-06*temp*temp +
//          -0.043948876890469453*temp +
//          216.81259827590887;
// }

__global__ void ComputeMixH(double *temperature_d,
                                     double *pt_d         ,
                                     double *pressure_d   ,
                                     double *Rho_d        ,
                                     double *mixH_d       ,
                                     double  P_Ref        ,
                                     double  Rd           ,
                                     double  Cp           ,
                                     int     num          ){

// calculate mixing ratio of atomic H via formulae in Heng book
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    int nv = gridDim.y;
    int lev = blockIdx.y;

    double dG, Kprime, temp;

    if (id < num){
      // dGibbs<<<1,1>>>(temperature_d[id*nv+lev], dG);
      temp = temperature_d[id*nv+lev];
      dG = 2.1370867596206315e-17*temp*temp*temp*temp*temp +
             -3.8689132818241159e-13*temp*temp*temp*temp +
             2.7275438366298867e-09*temp*temp*temp +
             -9.6170574202103724e-06*temp*temp +
             -0.043948876890469453*temp +
             216.81259827590887;
      Kprime = exp(2000*dG/Runiv/temperature_d[id*nv+lev])*pressure_d[id*nv+lev]/100000;
      mixH_d[id*nv+lev] = (-1.0+sqrt(1.0+8*Kprime))/(4*Kprime);
    }
}

__device__ double chemical_time_H(double *Mh_d){

  return 1e-3;
}

__device__ double maxT_lev(double *temperature_d, int lev, int num, int nv) {
  int i;
  double maxT = 0.0;

  for (i=0; i<num; i++) {
    if (temperature_d[i*nv+lev] > maxT) maxT = temperature_d[i*nv+lev];
  }

  return maxT;
}

__global__ void recomb_H(double *Mh_d         ,
                            double *pressure_d   ,
                            double *Rho_d        ,
                            double *temperature_d,
                            double *mixH_d       ,
                            double  Gravit       ,
                            double  Cp           ,
                            double  Rd           ,
                            double *Altitude_d   ,
                            double *Altitudeh_d  ,
                            double *lonlat_d     ,
                            double  time_step    ,
                            double *areasT       ,
                            int     num          ){

    int id = blockIdx.x * blockDim.x + threadIdx.x;
    int nv = gridDim.y;
    int lev = blockIdx.y;

    if (id < num){
      double Tau_dyn, Tau_chem, dbar, vh, Tau;
      //calculate dynamical time scale
      dbar = sqrt(areasT[id]);    //approximate "width" of cell
      vh = sqrt(Mh_d[id*3*nv+lev*3+0]*Mh_d[id*3*nv+lev*3+0] +
                Mh_d[id*3*nv+lev*3+1]*Mh_d[id*3*nv+lev*3+1] +
                Mh_d[id*3*nv+lev*3+2]*Mh_d[id*3*nv+lev*3+2])/Rho_d[id*nv+lev]; //horizontal speed
      Tau_dyn = dbar/vh;

      //calculate chemical time scale
      Tau_chem = chemical_time_H(Mh_d);

      double mixH_tmp, dT;

      if (Tau_dyn > Tau_chem) {
        mixH_tmp = mixH_d[id*nv+lev];
        Tau = Tau_dyn;
      } else {
        double dG, Kprime, temp;
        temp = maxT_lev(temperature_d,lev,num,nv);
        dG = ((((2.1370867596206315e-17*temp+-3.8689132818241159e-13)*temp +\
             2.7275438366298867e-09)*temp - 9.6170574202103724e-06)*temp +\
             -0.043948876890469453)*temp + 216.81259827590887;
        Kprime = exp(2000*dG/Runiv/temperature_d[id*nv+lev])*pressure_d[id*nv+lev]/100000;
        mixH_tmp = (-1.0+sqrt(1.0+8*Kprime))/(4*Kprime);
        Tau = Tau_chem;
      }

      dT = qbond/Cp*(mixH_tmp - mixH_d[id*nv+lev])/Tau;
      temperature_d[id*nv+lev] += dT;
    }
}
