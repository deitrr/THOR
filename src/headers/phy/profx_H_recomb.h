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
                            int     num          ){

    int id = blockIdx.x * blockDim.x + threadIdx.x;
    int nv = gridDim.y;
    int lev = blockIdx.y;

    if (id < num){


    }


}
