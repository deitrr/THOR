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
// Description: Physics modules.
//
//
// Method: This version just includes the held-suarez test.
//
// Known limitations: None
//
// Known issues: None
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

#include "../headers/esp.h"
#include "../headers/phy/profx_auxiliary.h"
#include "../headers/phy/profx_held_suarez.h"
#include "../headers/phy/profx_shallowHJ_hs.h"
#include "../headers/phy/profx_tidalearth_hs.h"
#include "../headers/phy/profx_H_recomb.h"

__host__ void ESP::ProfX(int    planetnumber, // Planet ID
                         int    nstep       , // Step number
                         int    hstest      , // Held-Suarez test option
                         double time_step   , // Time-step [s]
                         double Omega       , // Rotation rate [1/s]
                         double Cp          , // Specific heat capacity [J/kg/K]
                         double Rd          , // Gas constant [J/kg/K]
                         double Mmol        , // Mean molecular mass of dry air [kg]
                         double mu          , // Atomic mass unit [kg]
                         double kb          , // Boltzmann constant [J/K]
                         double P_Ref       , // Reference pressure [Pa]
                         double Gravit      , // Gravity [m/s^2]
                         double A           , // Planet radius [m]
                         bool hh2recomb     ){// option of atomic H<->H2


//
//  Number of threads per block.
    const int NTH = 256;

//  Specify the block sizes.
    dim3 NB((point_num / NTH) + 1, nv, 1);

//  Computes the initial temperature.
    Compute_temperature <<< NB, NTH >>> (temperature_d,
                                         pt_d         ,
                                         pressure_d   ,
                                         Rho_d        ,
                                         P_Ref        ,
                                         Rd           ,
                                         Cp           ,
                                         point_num    );
//  Check for nan.
    check_h = false;
    cudaMemcpy(check_d, &check_h, sizeof(bool), cudaMemcpyHostToDevice);
    isnan_check<<< 16, NTH >>>(temperature_d, nv, point_num, check_d);
    cudaMemcpy(&check_h, check_d, sizeof(bool), cudaMemcpyDeviceToHost);
    if(check_h){
       printf("\n\n Error in NAN check!\n");
       exit(EXIT_FAILURE);
    }

///////////////////////
// HELD SUAREZ TEST  //
///////////////////////
//
    if (planetnumber == 1) {
      if (hstest == 1) {
        cudaDeviceSynchronize();
        held_suarez<<< NB, NTH >>> (Mh_d         ,
                                    pressure_d   ,
                                    Rho_d        ,
                                    temperature_d,
                                    Gravit       ,
                                    Cp           ,
                                    Rd           ,
                                    Altitude_d   ,
                                    Altitudeh_d  ,
                                    lonlat_d     ,
                                    time_step    ,
                                    point_num    );
      } else if (hstest == 2) {
        cudaDeviceSynchronize();
        tidalearth_hs<<< NB, NTH >>> (Mh_d         ,
                                    pressure_d   ,
                                    Rho_d        ,
                                    temperature_d,
                                    Gravit       ,
                                    Cp           ,
                                    Rd           ,
                                    Altitude_d   ,
                                    Altitudeh_d  ,
                                    lonlat_d     ,
                                    time_step    ,
                                    point_num    );
      } else if (hstest == 3) {
        cudaDeviceSynchronize();
        shallowHJ_hs<<< NB, NTH >>> (Mh_d         ,
                                    pressure_d   ,
                                    Rho_d        ,
                                    temperature_d,
                                    Gravit       ,
                                    Cp           ,
                                    Rd           ,
                                    Altitude_d   ,
                                    Altitudeh_d  ,
                                    lonlat_d     ,
                                    time_step    ,
                                    point_num    );
      }
      if (hh2recomb) {
        cudaDeviceSynchronize();
        recomb_H<<< NB, NTH >>> (Mh_d         ,
                                    pressure_d   ,
                                    Rho_d        ,
                                    temperature_d,
                                    mixH_d       ,
                                    Gravit       ,
                                    Cp           ,
                                    Rd           ,
                                    Altitude_d   ,
                                    Altitudeh_d  ,
                                    lonlat_d     ,
                                    time_step    ,
                                    point_num    );
        cudaDeviceSynchronize();
        ComputeMixH<<< NB, NTH >>> (temperature_d,
                                    pt_d         ,
                                    pressure_d   ,
                                    Rho_d        ,
                                    mixH_d       ,
                                    P_Ref        ,
                                    Rd           ,
                                    Cp           ,
                                    point_num          );
      }
    }
//
////////////////////////

    if(planetnumber != 1){
        printf("Planet value incorrect! (see in file planet.h)");
        exit(EXIT_FAILURE);
    }

//  Computes the new pressures.
    cudaDeviceSynchronize();
    Compute_pressure <<< NB, NTH >>> (pressure_d   ,
                                      temperature_d,
                                      Rho_d        ,
                                      Rd           ,
                                      point_num    );
//
//END OF INTEGRATION
//
}
