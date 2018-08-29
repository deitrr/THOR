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
// Description: Computes temperature and pressure.
//
//
// Method: We use the assumption that the atmosphere can be treated as an ideal gas.
//
// Known limitations: None
//
// Known issues: None.
//
//
// If you use this code please cite the following reference:
//
//       [1] Mendonca, J.M., Grimm, S.L., Grosheintz, L., & Heng, K., ApJ, 829, 115, 2016
//
// Current Code Owner: Joao Mendonca, EEG. joao.mendonca@csh.unibe.ch
//
// History:
// Version Date       Comment
// ======= ====       =======
//
// 1.0     16/08/2017 Released version  (JM)
//
////////////////////////////////////////////////////////////////////////

__global__ void Compute_temperature_only (double *temperature_d,
                                     double *pressure_d   ,
                                     double *Rho_d        ,
                                     double  Rd           ,
                                     int     num          ){

    int id = blockIdx.x * blockDim.x + threadIdx.x;
    int nv = gridDim.y;
    int lev = blockIdx.y;

    // Computes absolute temperature
    if(id < num) temperature_d[id*nv + lev] = pressure_d[id*nv + lev]/(Rd*Rho_d[id*nv + lev]);
}

__global__ void Compute_temperature (double *temperature_d,
                                     double *pt_d         ,
                                     double *pressure_d   ,
                                     double *Rho_d        ,
                                     double  P_Ref        ,
                                     double  Rd           ,
                                     double *CpT_d        ,
                                     int     num          ,
                                     bool    CpTemp       ){

    int id = blockIdx.x * blockDim.x + threadIdx.x;
    int nv = gridDim.y;
    int lev = blockIdx.y;

    // Computes absolute and potential temperature
    if(id < num) {
      double Cv = CpT_d[id*nv+lev] - Rd;
      double CvoCp = Cv / CpT_d[id*nv+lev];

      double Cp0 = 18567.9235606, T0 = 3150.95, nu = 0.183215136;
      double k0 = Rd/Cp0;

      temperature_d[id*nv + lev] = pressure_d[id*nv + lev]/(Rd*Rho_d[id*nv + lev]);
      if (CpTemp)  {
        pt_d[id * nv + lev] = pow( pow(temperature_d[id*nv+lev],nu) + k0*nu*pow(T0,nu)*log(P_Ref/(pressure_d[id*nv+lev])), 1.0/nu);
      } else {
        pt_d[id * nv + lev] = (P_Ref / (Rd * Rho_d[id*nv + lev]))*pow(pressure_d[id*nv + lev] / P_Ref, CvoCp);
      }
    }
}

__global__ void Compute_pressure    (double *pressure_d   ,
                                     double *temperature_d,
                                     double *Rho_d        ,
                                     double  Rd           ,
                                     int     num          ){


    int id = blockIdx.x * blockDim.x + threadIdx.x;
    int nv = gridDim.y;
    int lev = blockIdx.y;

    // Computes absolute pressure
    if (id < num) {
      pressure_d[id*nv + lev] = temperature_d[id*nv + lev] * Rd*Rho_d[id*nv + lev];
      // if (isnan(pressure_d[id*nv+lev])){
      //   printf("%d, %d\n",id,lev);
      // }
    }
}


__global__ void isnan_check(double *array, int width, int height, bool *check){

//
//  Description: Check for nan in the temperature array.
//

  int idx = threadIdx.x+blockDim.x*blockIdx.x;

  while (idx < width){
    for (int i = 0; i < height; i++)
      if (isnan(array[(i*width) + idx])) {
        *check = true;
      }
      idx += gridDim.x+blockDim.x;
    }
}

__global__ void isnan_loop(double *t, int num, int nv) {
  int i;

  for (i=0;i<num;i++) {
    if (isnan(t[i*nv])) {
      printf("id = %d, lev = %d, t = %f\n",i,0,t[i*nv]);
    }
  }
}
