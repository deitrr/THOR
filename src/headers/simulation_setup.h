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
// Description: Defines simulation parameters and switches that are used throughout the sim
//
//
// Method: -
//
//
// Known limitations: None
//
//
// Known issues: None
//
//
// If you use this code please cite the following reference:
//
//       [1] Mendonca, J.M., Grimm, S.L., Grosheintz, L., & Heng, K., ApJ, 829, 115, 2016
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
#pragma once

// Physical Constants
#define kb_constant 1.38e-23  // Boltzmann constant [J/K]
#define mu_constant 1.660e-27 // Atomic mass unit   [kg]


class SimulationSetup
{

public:
    char simulation_ID[160];

    //////////////
    // BULK     //
    //////////////

    double A;
    double Omega;
    double Gravit;

    ////////////////
    // ATMOSPHERE //
    ////////////////

    double Rd;

    double Cp;
    double Tmean;
    double P_Ref;
    double Top_altitude;
    double Diffc;
    double Diffc_v;
    double DivDampc;

    // Sim
    bool DeepModel;
    bool HyDiff; // Turn on/off hyper-diffusion.
    int  HyDiffOrder;
    bool DivDampP;   // Turn on/off divergence damping.
    bool VertHyDiff; // vertical hyper diffusion
    int  VertHyDiffOrder;
    bool NonHydro; // Turn on/off non-hydrostatic.
    bool globdiag; // calc/output globdiag quantities

    // top sponge layer master switches
    bool RayleighSponge;  // Use sponge layer (rayleigh drag)?
    bool RayleighSpongeT; // include thermal term in sponge layer?
    bool DiffSponge;      // Diffusive sponge

    bool output_mean; //whether or not to output time mean quantities
    bool
         out_interm_momentum; //output intermediate momentum values (start of time step & after profx)
    bool output_diffusion;    //output hyperdiffusion operators, etc

    bool conv_adj;
    int  conv_adj_iter;
    bool gcm_off;
    bool single_column;

    bool rest;

    int n_out;

    SimulationSetup();
};
