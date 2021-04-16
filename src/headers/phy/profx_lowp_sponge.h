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
// Description:
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
//       [1] Mendonca, J.M., Grimm, S.L., Grosheintz, L., & Heng, K., ApJ, 829, 115, 2016
//
//       [2] Held, I. M., & Suarez, M. J. 1994, Bullentin of the American
//           Meteorological Society
//
// Current Code Owners: Joao Mendonca (joao.mendonca@space.dtu.dk)
//                      Russell Deitrick (russell.deitrick@csh.unibe.ch)
//                      Urs Schroffenegger (urs.schroffenegger@csh.unibe.ch)
//
// History:
// Version Date       Comment
// ======= ====       =======
// 2.0     30/11/2018 Released version (RD & US)
// 1.0     16/08/2017 Released version  (JM)
//
////////////////////////////////////////////////////////////////////////


__global__ void lowp_sponge(double *Mh_d,
                            double *Wh_d,
                            double *pressure_d,
                            double *Rho_d,
                            double *Altitude_d,
                            double *Altitudeh_d,
                            double  Pup,
                            double  Pdown,
                            double  kmax,
                            double  Rd,
                            double  time_step,
                            int     num) {

    int id  = blockIdx.x * blockDim.x + threadIdx.x;
    int nv  = gridDim.y;
    int lev = blockIdx.y;

    if (id < num) {
        //double Pdown  = 1e-3 * 1e5; // lowest pressure to apply sponge
        // double Pup    = 1e-6 * 1e5; // sponge reaches maximum at this pressure
        double eta    = log10(pressure_d[id * nv + lev]) / log10(Pup);
        double eta_sp = log10(Pdown) / log10(Pup);
        // double kmax   = 1e-3; //max value of sponge
        double ksp;
        double Tlim = 1000;

        //calculate strength of sponge
        if (pressure_d[id * nv + lev] > Pdown) {
            ksp = 0;
        }
        else if (pressure_d[id * nv + lev] < Pup) {
            ksp = kmax;
        }
        else {
            ksp = kmax * pow(sin(M_PI / 2 * (eta - eta_sp) / (1 - eta_sp)), 2);
        }

        //update horizontal momentum implicitly
        for (int k = 0; k < 3; k++) {
            Mh_d[id * 3 * nv + lev * 3 + k] =
                Mh_d[id * 3 * nv + lev * 3 + k] / (1.0 + ksp * time_step);
        }

        //vertical needs to be done a bit differently
        //interpolate pressure to interfaces
        //no need to damp top interface so indexing stays the same

        if (lev != 0) {
            double pint = pressure_d[id * nv + lev - 1]
                          + (pressure_d[id * nv + lev] - pressure_d[id * nv + lev - 1])
                                * (Altitudeh_d[lev] - Altitude_d[lev - 1])
                                / (Altitude_d[lev] - Altitude_d[lev - 1]);
            double eta_int, ksp_int;
            eta_int = log10(pint) / log10(Pup);
            //calculate strength of sponge
            if (pint > Pdown) {
                ksp_int = 0;
            }
            else if (pint < Pup) {
                ksp_int = kmax;
            }
            else {
                ksp_int = kmax * pow(sin(M_PI / 2 * (eta_int - eta_sp) / (1 - eta_sp)), 2);
            }
            Wh_d[id * (nv + 1) + lev] = Wh_d[id * (nv + 1) + lev] / (1.0 + ksp_int * time_step);
        }

        //also temperature
        double T = pressure_d[id * nv + lev] / Rd / Rho_d[id * nv + lev];
        pressure_d[id * nv + lev] =
            Rd * Rho_d[id * nv + lev] * (T / time_step + ksp * Tlim) / (1.0 / time_step + ksp);
    }
}
