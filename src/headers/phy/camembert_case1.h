
#pragma once

__device__ double camembert_k2_18b_interp(double P, double *P_IC, double *T_IC, int n_pressures) {
  double T = 0;
  int i_lower = -1;
  //interpolation here!
  //first search for nearest pressures in P_IC_h
  if (P <= P_IC[0]) {
    T = T_IC[0]; //padding sides with isotherms
  } else if (P >= P_IC[n_pressures - 1]) {
    T = T_IC[n_pressures - 1]; //padding sides with isotherms
  } else {
    for (int i = 0; i < n_pressures - 1; i++) {
      if (P >= P_IC[i] && P < P_IC[i+1]){
        i_lower = i;
        break;
      }
    }
    T = T_IC[i_lower] + (P - P_IC[i_lower]) *
          (T_IC[i_lower+1] - T_IC[i_lower]) /
          (P_IC[i_lower+1] - P_IC[i_lower]);
  }
  return T;
}

__global__ void camembert_k2_18b_force(double *pressure_d,
                                       double *Rho_d,
                                       double *temperature_d,
                                       double *TtendencyTF_d,
                                       double Gravit,
                                       double Cp,
                                       double Rd,
                                       double *lonlat_d,
                                       double *P_IC_d,
                                       double *T_IC_d,
                                       double dTeqmax,
                                       double time_step,
                                       int n_pressures,
                                       int num) {

   int id  = blockIdx.x * blockDim.x + threadIdx.x;
   int nv  = gridDim.y;
   int lev = blockIdx.y;

   if (id < num) {
      double p = pressure_d[id*nv+lev], p_low_tau = 1e2, p_low_dT = 10, p_hi = 1e6;
      double tau_rad_H2, tau_rad;
      double mu = 8.31446261815324*1000/Rd;
      double T0, Teq, dTeq, slope;
      double lat = lonlat_d[id * 2 + 1];
      double lon = lonlat_d[id * 2];

      //radiative time-scales
      if (p <= p_low_tau) {
        tau_rad_H2 = 1e4;
      } else if ((p > p_low_tau) && (p < p_hi)) {
        tau_rad_H2 = pow(10.0,2.5) * pow(p,0.75);
      } else {
        tau_rad_H2 = 1e7;
      }

      tau_rad = tau_rad_H2 * Cp / (3.5*Rd) * 2 / mu;

      //set up Delta T
      if (p <= p_low_dT) {
        dTeq = dTeqmax;
      } else if (p >= p_hi) {
        dTeq = 0.0;
      } else {
        slope = -dTeqmax / (log10(p_hi)-log10(p_low_dT));
        dTeq = slope*log10(p) - slope*log10(p_low_dT) + dTeqmax;
      }

      //interpolate to get init T at this p
      T0 = camembert_k2_18b_interp(p, P_IC_d, T_IC_d, n_pressures);

      //get Teq
      if ((lon <= 0.5 * M_PI) || (lon > 1.5 * M_PI)) {
        //dayside
        Teq = T0 + dTeq * (fabs(cos(lat)*cos(lon)) - 0.5);
      } else {
        //nightside
        Teq = T0 - 0.5*dTeq;
      }

      double ttmp = temperature_d[id*nv+lev];

      temperature_d[id*nv+lev] = (Teq/tau_rad + temperature_d[id*nv+lev]/time_step) /
                                 (1.0/tau_rad + 1.0/time_step);
      
      //compute temperature tendency
      TtendencyTF_d[id*nv+lev] = (temperature_d[id*nv+lev] - ttmp) / time_step;
   }

}
