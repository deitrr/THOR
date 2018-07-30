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
// Description: Computes the fast modes (see equations 33 to 35 from Mendonca et al. 2016)
//
//
// Method: -
//
// Known limitations: None.
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

template <int NX, int NY>
__global__ void Momentum_Eq (double *M_d        ,
                             double *pressure_d ,
                             double *SlowMh_d   ,
                             double *grad_d     ,
                             double *Altitude_d ,
                             double *DivM_d     ,
                             double  A          ,
                             double *func_r_d   ,
                             double  dt         ,
                             int    *maps_d     ,
                             int     nl_region  ,
                             bool    DeepModel  ){

    int x = threadIdx.x;
    int y = threadIdx.y;
    int ib = blockIdx.x;
    int nv = gridDim.y;
    int lev = blockIdx.y;

    int pt1, pt2, pt3, pt4, pt5, pt6;
    int nhl = nl_region + 2;
    int nhl2 = nhl*nhl;

    int ir = (y + 1)*nhl + x + 1;   // Region index
    int iri;

    __shared__ double nflxv_s[3 * NX*NY];
    __shared__ double pressure_s[(NX + 2)*(NY + 2)];

    int pent_ind = 0;
    int ig, ir2, id, twot;
    double vr;
    double alt, rscale;
    double Mx, My, Mz;
    double funcx, funcy, funcz;

    // Load shared memory
    ig = maps_d[ib*nhl2 + ir];
    id = ig;
    if (x == 0 && y == 0) if (maps_d[ib * nhl2] == -1) pent_ind = 1;
    pressure_s[ir] = pressure_d[ig * nv + lev];

    ///////////////////////////////
    //////////// Halo /////////////
    ///////////////////////////////
    if (x == 0) {
        ir2 = (y + 1) * nhl + x;
        ig = maps_d[ib * nhl2 + ir2];
        pressure_s[ir2] = pressure_d[ig * nv + lev];
    }
    if (x == nhl - 3){
        ir2 = (y + 1) * nhl + x + 2;
        ig = maps_d[ib * nhl2 + ir2];
        pressure_s[ir2] = pressure_d[ig * nv + lev];
    }
    if (y == 0){
        twot = 1;
        ir2 = y * nhl + (x + 1);
        if (x == 0) twot = 2;

        for (int k = 0; k < twot; k++){
            if (k == 1) ir2 = y * nhl + x;
            ig = maps_d[ib * nhl2 + ir2];
            if (ig >= 0) pressure_s[ir2] = pressure_d[ig * nv + lev];
            else         pressure_s[ir2] = 0.0;
        }
    }
    if (y == nhl - 3) {
        twot = 1;
        ir2 = (y + 2) * nhl + (x + 1);
        if (x == nhl - 3) twot = 2;
        for (int k = 0; k < twot; k++){
            if (k == 1) ir2 = (y + 2) * nhl + (x + 2);
            ig = maps_d[ib * nhl2 + ir2];
            pressure_s[ir2] = pressure_d[ig * nv + lev];
        }
    }
    __syncthreads();
    //////////////////////////////////////////////

    iri = y * nl_region + x;

    funcx = func_r_d[id * 3 + 0];
    funcy = func_r_d[id * 3 + 1];
    funcz = func_r_d[id * 3 + 2];

    pt1 = (y + 2)*nhl + x + 1;
    pt2 = (y + 2)*nhl + x + 2;
    pt3 = (y + 1)*nhl + x + 2;
    pt4 = (y    )*nhl + x + 1;
     pt5 = (pent_ind)*((y + 1)*nhl + x) + (!pent_ind)*((y   )*nhl + x);
    pt6 = (y + 1)*nhl + x;

    if (DeepModel){
        alt = Altitude_d[lev];
        rscale = A / (alt + A);
    }
    else rscale = 1.0;

    for (int k = 0; k < 3; k++){
            nflxv_s[iri * 3 + k] = rscale*(grad_d[id * 7 * 3 + 3 * 0 + k] * pressure_s[ir] +
                                           grad_d[id * 7 * 3 + 3 * 1 + k] * pressure_s[pt1] +
                                           grad_d[id * 7 * 3 + 3 * 2 + k] * pressure_s[pt2] +
                                           grad_d[id * 7 * 3 + 3 * 3 + k] * pressure_s[pt3] +
                                           grad_d[id * 7 * 3 + 3 * 4 + k] * pressure_s[pt4] +
                                           grad_d[id * 7 * 3 + 3 * 5 + k] * pressure_s[pt5] +
                                           grad_d[id * 7 * 3 + 3 * 6 + k] * pressure_s[pt6]);
    }

    Mx = (-nflxv_s[iri * 3 + 0] + SlowMh_d[id * 3 * nv + lev * 3 + 0] + DivM_d[id * 3 * nv + lev * 3 + 0])*dt;
    My = (-nflxv_s[iri * 3 + 1] + SlowMh_d[id * 3 * nv + lev * 3 + 1] + DivM_d[id * 3 * nv + lev * 3 + 1])*dt;
    Mz = (-nflxv_s[iri * 3 + 2] + SlowMh_d[id * 3 * nv + lev * 3 + 2] + DivM_d[id * 3 * nv + lev * 3 + 2])*dt;

    vr = Mx * funcx + My * funcy + Mz * funcz;

    Mx += -vr*funcx;
    My += -vr*funcy;
    Mz += -vr*funcz;

    // Updates momenta
    M_d[id * nv * 3 + lev * 3 + 0] += Mx;
    M_d[id * nv * 3 + lev * 3 + 1] += My;
    M_d[id * nv * 3 + lev * 3 + 2] += Mz;

}

template <int NN>
__global__ void Momentum_Eq_Poles (double * M_d        ,
                                   double * pressure_d ,
                                   double * SlowMh_d   ,
                                   double * grad_d     ,
                                   double * Altitude_d ,
                                   double * DivM_d     ,
                                   double   A          ,
                                   double * func_r_d   ,
                                   double   dt         ,
                                   int * point_local_d ,
                                   int  nv             ,
                                   int  num            ,
                                   bool  DeepModel      ){

    int id = blockIdx.x * blockDim.x + threadIdx.x;
    id += num - 2;     // Poles

    __shared__ double grad_p[3 * 7] ;
    __shared__ double func_r_p[3]   ;
    __shared__ double nflxv_p[3]    ;
    __shared__ double pressure_p[NN];
    __shared__ int local_p[NN]      ;

    double vr;
    double alt, rscale;
    double Mx, My, Mz;

    if (id < num){
        for (int i = 0; i < 5; i++) local_p[i] = point_local_d[id * 6 + i];
        func_r_p[0] = func_r_d[id * 3 + 0];
        func_r_p[1] = func_r_d[id * 3 + 1];
        func_r_p[2] = func_r_d[id * 3 + 2];
        for (int i = 0; i < 7; i++) for (int k = 0; k < 3; k++) grad_p[i * 3 + k] = grad_d[id * 7 * 3 + i * 3 + k];

        for (int lev = 0; lev < nv; lev++){

            pressure_p[0] = pressure_d[id * nv + lev];
            for (int i = 1; i < 6; i++)    pressure_p[i] = pressure_d[local_p[i - 1] * nv + lev];

            alt = Altitude_d[lev];

            if (DeepModel)    rscale = A / (alt + A);
            else rscale = 1.0;

            for (int k = 0; k < 3; k++){
                nflxv_p[k] = rscale*(grad_p[3 * 0 + k] * pressure_p[0] +
                                     grad_p[3 * 1 + k] * pressure_p[1] +
                                     grad_p[3 * 2 + k] * pressure_p[2] +
                                     grad_p[3 * 3 + k] * pressure_p[3] +
                                     grad_p[3 * 4 + k] * pressure_p[4] +
                                     grad_p[3 * 5 + k] * pressure_p[5]);
            }

            Mx = (- nflxv_p[0] + SlowMh_d[id * 3 * nv + lev * 3 + 0] + DivM_d[id * 3 * nv + lev * 3 + 0])*dt;
            My = (- nflxv_p[1] + SlowMh_d[id * 3 * nv + lev * 3 + 1] + DivM_d[id * 3 * nv + lev * 3 + 1])*dt;
            Mz = (- nflxv_p[2] + SlowMh_d[id * 3 * nv + lev * 3 + 2] + DivM_d[id * 3 * nv + lev * 3 + 2])*dt;

            vr = Mx * func_r_p[0] + My * func_r_p[1] + Mz * func_r_p[2];

            Mx += -vr*func_r_p[0];
            My += -vr*func_r_p[1];
            Mz += -vr*func_r_p[2];

            // Updates momenta
            M_d[id * nv * 3 + lev * 3 + 0] += Mx;
            M_d[id * nv * 3 + lev * 3 + 1] += My;
            M_d[id * nv * 3 + lev * 3 + 2] += Mz;
        }
    }
}

template <int NX, int NY>
__global__ void Density_Pressure_Eqs(double *pressure_d ,
                                     double *pressurek_d,
                                     double *Rho_d      ,
                                     double *Rhok_d     ,
                                     double *Mh_d       ,
                                     double *Mhk_d      ,
                                     double *Wh_d       ,
                                     double *Whk_d      ,
                                     double *pt_d       ,
                                     double *pth_d      ,
                                     double *SlowRho_d  ,
                                     double *diffpr_d   ,
                                     double *div_d      ,
                                     double *Altitude_d ,
                                     double *Altitudeh_d,
                                     double  Cp         ,
                                     double  Rd         ,
                                     double  A          ,
                                     double  P_Ref      ,
                                     double  dt         ,
                                     int *maps_d        ,
                                     int nl_region      ,
                                     bool DeepModel     ){

    int x = threadIdx.x;
    int y = threadIdx.y;
    int ib = blockIdx.x;
    int nv = gridDim.y;
    int lev = blockIdx.y;

    int pt1, pt2, pt3, pt4, pt5, pt6;
    double div0, div1, div2, div3, div4, div5, div6;
    int nhl = nl_region + 2;
    int nhl2 = nhl*nhl;

    int ir = (y + 1)*nhl + x + 1;   // Region index
    int iri, ir2, twot;

    __shared__ double nflxr_s[NX*NY];
    __shared__ double nflxpt_s[NX*NY];
    __shared__ double v_s[3 * (NX + 2)*(NY + 2)];
    __shared__ double v1_s[3 * (NX + 2)*(NY + 2)];
    __shared__ double pt_s[(NX + 2)*(NY + 2)];

    int pent_ind = 0;
    int ig, id;

    double alt, r2p, r2m, r2l, rscale;
    double wht, whl;
    double wht2,whl2;
    double pht, phl;
    double dz, dwdz;
    double dwptdz;
    double aux, r, p;
    double altht, althl;
    double Cv = Cp - Rd;

    // Load shared memory
    ig = maps_d[ib*nhl2 + ir];
    id = ig;
    if (x == 0 && y == 0) if (maps_d[ib * nhl2] == -1) pent_ind = 1;

    v_s[ir * 3 + 0] = Mh_d[ig * 3 * nv + lev * 3 + 0];
    v_s[ir * 3 + 1] = Mh_d[ig * 3 * nv + lev * 3 + 1];
    v_s[ir * 3 + 2] = Mh_d[ig * 3 * nv + lev * 3 + 2];

    v1_s[ir * 3 + 0] = v_s[ir * 3 + 0] + Mhk_d[ig * 3 * nv + lev * 3 + 0];
    v1_s[ir * 3 + 1] = v_s[ir * 3 + 1] + Mhk_d[ig * 3 * nv + lev * 3 + 1];
    v1_s[ir * 3 + 2] = v_s[ir * 3 + 2] + Mhk_d[ig * 3 * nv + lev * 3 + 2];

    pt_s[ir] = pt_d[ig * nv + lev];

    ///////////////////////////////
    //////////// Halo /////////////
    ///////////////////////////////
    if (x == 0) {
        ir2 = (y + 1) * nhl + x;
        ig = maps_d[ib * nhl2 + ir2];
        v_s[ir2 * 3 + 0] = Mh_d[ig * 3 * nv + lev * 3 + 0];
        v_s[ir2 * 3 + 1] = Mh_d[ig * 3 * nv + lev * 3 + 1];
        v_s[ir2 * 3 + 2] = Mh_d[ig * 3 * nv + lev * 3 + 2];
        v1_s[ir2 * 3 + 0] = v_s[ir2 * 3 + 0] + Mhk_d[ig * 3 * nv + lev * 3 + 0];
        v1_s[ir2 * 3 + 1] = v_s[ir2 * 3 + 1] + Mhk_d[ig * 3 * nv + lev * 3 + 1];
        v1_s[ir2 * 3 + 2] = v_s[ir2 * 3 + 2] + Mhk_d[ig * 3 * nv + lev * 3 + 2];
        pt_s[ir2] = pt_d[ig * nv + lev];
    }
    if (x == nhl - 3){
        ir2 = (y + 1) * nhl + x + 2;
        ig = maps_d[ib * nhl2 + ir2];
        v_s[ir2 * 3 + 0] = Mh_d[ig * 3 * nv + lev * 3 + 0];
        v_s[ir2 * 3 + 1] = Mh_d[ig * 3 * nv + lev * 3 + 1];
        v_s[ir2 * 3 + 2] = Mh_d[ig * 3 * nv + lev * 3 + 2];
        v1_s[ir2 * 3 + 0] = v_s[ir2 * 3 + 0] + Mhk_d[ig * 3 * nv + lev * 3 + 0];
        v1_s[ir2 * 3 + 1] = v_s[ir2 * 3 + 1] + Mhk_d[ig * 3 * nv + lev * 3 + 1];
        v1_s[ir2 * 3 + 2] = v_s[ir2 * 3 + 2] + Mhk_d[ig * 3 * nv + lev * 3 + 2];
        pt_s[ir2] = pt_d[ig * nv + lev];
    }
    if (y == 0){
        twot = 1;
        ir2 = y * nhl + (x + 1);
        if (x == 0) twot = 2;

        for (int k = 0; k < twot; k++){
            if (k == 1) ir2 = y * nhl + x;
            ig = maps_d[ib * nhl2 + ir2];

            if (ig >= 0){
                v_s[ir2 * 3 + 0] = Mh_d[ig * 3 * nv + lev * 3 + 0];
                v_s[ir2 * 3 + 1] = Mh_d[ig * 3 * nv + lev * 3 + 1];
                v_s[ir2 * 3 + 2] = Mh_d[ig * 3 * nv + lev * 3 + 2];
                v1_s[ir2 * 3 + 0] = v_s[ir2 * 3 + 0] + Mhk_d[ig * 3 * nv + lev * 3 + 0];
                v1_s[ir2 * 3 + 1] = v_s[ir2 * 3 + 1] + Mhk_d[ig * 3 * nv + lev * 3 + 1];
                v1_s[ir2 * 3 + 2] = v_s[ir2 * 3 + 2] + Mhk_d[ig * 3 * nv + lev * 3 + 2];
                pt_s[ir2] = pt_d[ig * nv + lev];
            }
            else{
                v_s[ir2 * 3 + 0]  = 0.0;
                v_s[ir2 * 3 + 1]  = 0.0;
                v_s[ir2 * 3 + 2]  = 0.0;
                v1_s[ir2 * 3 + 0] = 0.0;
                v1_s[ir2 * 3 + 1] = 0.0;
                v1_s[ir2 * 3 + 2] = 0.0;
                pt_s[ir2] = 0.0;
            }
        }
    }
    if (y == nhl - 3) {
        twot = 1;
        ir2 = (y + 2) * nhl + (x + 1);
        if (x == nhl - 3) twot = 2;
        for (int k = 0; k < twot; k++){
            if (k == 1) ir2 = (y + 2) * nhl + (x + 2);
            ig = maps_d[ib * nhl2 + ir2];
            v_s[ir2 * 3 + 0] = Mh_d[ig * 3 * nv + lev * 3 + 0];
            v_s[ir2 * 3 + 1] = Mh_d[ig * 3 * nv + lev * 3 + 1];
            v_s[ir2 * 3 + 2] = Mh_d[ig * 3 * nv + lev * 3 + 2];
            v1_s[ir2 * 3 + 0] = v_s[ir2 * 3 + 0] + Mhk_d[ig * 3 * nv + lev * 3 + 0];
            v1_s[ir2 * 3 + 1] = v_s[ir2 * 3 + 1] + Mhk_d[ig * 3 * nv + lev * 3 + 1];
            v1_s[ir2 * 3 + 2] = v_s[ir2 * 3 + 2] + Mhk_d[ig * 3 * nv + lev * 3 + 2];
            pt_s[ir2] = pt_d[ig * nv + lev];
        }
    }
    __syncthreads();
    //////////////////////////////////////////////

    iri = (y  )*nl_region + x;

    pt1 = (y + 2)*nhl + x + 1;
    pt2 = (y + 2)*nhl + x + 2;
    pt3 = (y + 1)*nhl + x + 2;
    pt4 = (y    )*nhl + x + 1;
     pt5 = (pent_ind)*((y + 1)*nhl + x) + (!pent_ind)*((y   )*nhl + x);
    pt6 = (y + 1)*nhl + x;

    altht = Altitudeh_d[lev + 1];
    althl = Altitudeh_d[lev];

    if (DeepModel){
        alt = Altitude_d[lev];
        r2p = pow(altht + A, 2.0);
        r2m = pow(alt + A, 2.0);
        r2l = pow(althl + A, 2.0);
        rscale = A / (alt + A);
    }
    else{
        r2p = 1.0;
        r2m = 1.0;
        r2l = 1.0;
        rscale = 1.0;
    }

    nflxr_s[iri]  = 0.0;
    nflxpt_s[iri] = 0.0;

    for (int k = 0; k < 3; k++){

        div0 = div_d[id * 7 * 3 + 3 * 0 + k];
        div1 = div_d[id * 7 * 3 + 3 * 1 + k];
        div2 = div_d[id * 7 * 3 + 3 * 2 + k];
        div3 = div_d[id * 7 * 3 + 3 * 3 + k];
        div4 = div_d[id * 7 * 3 + 3 * 4 + k];
        div5 = div_d[id * 7 * 3 + 3 * 5 + k];
        div6 = div_d[id * 7 * 3 + 3 * 6 + k];

        nflxr_s[iri] += rscale*(div0 * v_s[ir  * 3 + k] +
                                div1 * v_s[pt1 * 3 + k] +
                                div2 * v_s[pt2 * 3 + k] +
                                div3 * v_s[pt3 * 3 + k] +
                                div4 * v_s[pt4 * 3 + k] +
                                div5 * v_s[pt5 * 3 + k] +
                                div6 * v_s[pt6 * 3 + k]);

        nflxpt_s[iri] += rscale*(div0 * v1_s[ir * 3 + k] * pt_s[ir] +
                                 div1 * v1_s[pt1 * 3 + k]* pt_s[pt1]+
                                 div2 * v1_s[pt2 * 3 + k]* pt_s[pt2] +
                                 div3 * v1_s[pt3 * 3 + k]* pt_s[pt3] +
                                 div4 * v1_s[pt4 * 3 + k]* pt_s[pt4] +
                                 div5 * v1_s[pt5 * 3 + k]* pt_s[pt5] +
                                 div6 * v1_s[pt6 * 3 + k]* pt_s[pt6]);
    }

    if (lev == 0){
        whl = 0.0;
        wht = Wh_d[id*(nv + 1) + lev + 1];
        whl2= 0.0;
        wht2= Wh_d[id*(nv + 1) + lev + 1] + Whk_d[id*(nv + 1) + lev + 1];
        phl = 0.0;
        pht = pth_d[id*(nv + 1) + lev + 1];
    }
    else{
        whl = Wh_d[id*(nv + 1) + lev];
        wht = Wh_d[id*(nv + 1) + lev + 1];
        whl2= Wh_d[id*(nv + 1) + lev    ] + Whk_d[id*(nv + 1) + lev    ];
        wht2= Wh_d[id*(nv + 1) + lev + 1] + Whk_d[id*(nv + 1) + lev + 1];
        phl = pth_d[id*(nv + 1) + lev    ];
        pht = pth_d[id*(nv + 1) + lev + 1];
    }

    dz = altht - althl;
    dwdz = (wht*r2p - whl*r2l) / (dz*r2m);
    dwptdz = (wht2*pht*r2p - whl2*phl*r2l) / (dz*r2m);

    aux = -(nflxpt_s[iri] + dwptdz)*dt;
    r   = Rhok_d[id * nv + lev] + Rho_d[id * nv + lev];
    aux += pt_s[ir]*r;

    // Updates pressure
    p = P_Ref*pow(Rd*aux / P_Ref, Cp / Cv);
    pressure_d[id * nv + lev] = p - pressurek_d[id * nv + lev] + diffpr_d[id * nv + lev] * dt;

    if(isnan(pressure_d[id*nv+lev])){
      printf("Vertical wind gradient too large at %d, %d\n",id,lev);
    }
    if(pressurek_d[id*nv+lev]<-pressure_d[id*nv+lev]){
      printf("Pressure negative at %d, %d\n",id,lev);
    }

    // Updates density
    nflxr_s[iri] += dwdz;
    Rho_d[id * nv + lev] += (SlowRho_d[id * nv + lev] - nflxr_s[iri])*dt;
}

template <int NN>
__global__ void Density_Pressure_Eqs_Poles(double *pressure_d  ,
                                           double *pressurek_d ,
                                           double *Rho_d       ,
                                           double *Rhok_d      ,
                                           double *Mh_d        ,
                                           double *Mhk_d       ,
                                           double *Wh_d        ,
                                           double *Whk_d       ,
                                           double *pt_d        ,
                                           double *pth_d       ,
                                           double *SlowRho_d   ,
                                           double *diffpr_d    ,
                                           double *div_d       ,
                                           double *Altitude_d  ,
                                           double *Altitudeh_d ,
                                           double  Cp          ,
                                           double  Rd          ,
                                           double  A           ,
                                           double  P_Ref       ,
                                           double  dt          ,
                                           int    *point_local_d,
                                           int     num          ,
                                           int     nv           ,
                                           bool    DeepModel    ){

    int id = blockIdx.x * blockDim.x + threadIdx.x;
    id += num - 2;     // Poles

    __shared__ double div_p[3 * 7];
    __shared__ double v_p[3 * NN];
    __shared__ double v1_p[3 * NN];
    __shared__ double pt_p[NN];
    __shared__ int local_p[NN];

    double nflxr_p;
    double nflxpt_p;

    double alt, r2p, r2m, r2l, rscale;
    double wht, whl;
    double wht2, whl2;
    double pht, phl;
    double dz, dwdz;
    double dwptdz;
    double aux, r, p;
    double altht, althl;
    double Cv = Cp - Rd;

    if (id < num){
        for (int i = 0; i < 5; i++)local_p[i] = point_local_d[id * 6 + i];
        for (int i = 0; i < 7; i++) for (int k = 0; k < 3; k++) div_p[i * 3 + k] = div_d[id * 7 * 3 + i * 3 + k];

        for (int lev = 0; lev < nv; lev++){
            v_p[0] = Mh_d[id * 3 * nv + lev * 3 + 0];
            v_p[1] = Mh_d[id * 3 * nv + lev * 3 + 1];
            v_p[2] = Mh_d[id * 3 * nv + lev * 3 + 2];
            v1_p[0] = Mh_d[id * 3 * nv + lev * 3 + 0] + Mhk_d[id * 3 * nv + lev * 3 + 0];
            v1_p[1] = Mh_d[id * 3 * nv + lev * 3 + 1] + Mhk_d[id * 3 * nv + lev * 3 + 1];
            v1_p[2] = Mh_d[id * 3 * nv + lev * 3 + 2] + Mhk_d[id * 3 * nv + lev * 3 + 2];
            pt_p[0] = pt_d[id * nv + lev];
            for (int i = 1; i < 6; i++){
                v_p[i * 3 + 0] = Mh_d[local_p[i - 1] * 3 * nv + lev * 3 + 0];
                v_p[i * 3 + 1] = Mh_d[local_p[i - 1] * 3 * nv + lev * 3 + 1];
                v_p[i * 3 + 2] = Mh_d[local_p[i - 1] * 3 * nv + lev * 3 + 2];
                v1_p[i * 3 + 0] = Mh_d[local_p[i - 1] * 3 * nv + lev * 3 + 0] + Mhk_d[local_p[i - 1] * 3 * nv + lev * 3 + 0];
                v1_p[i * 3 + 1] = Mh_d[local_p[i - 1] * 3 * nv + lev * 3 + 1] + Mhk_d[local_p[i - 1] * 3 * nv + lev * 3 + 1];
                v1_p[i * 3 + 2] = Mh_d[local_p[i - 1] * 3 * nv + lev * 3 + 2] + Mhk_d[local_p[i - 1] * 3 * nv + lev * 3 + 2];
                pt_p[i] = pt_d[local_p[i - 1] * nv + lev];
            }

            if (lev == 0){
                altht = Altitudeh_d[lev + 1];
                althl = Altitudeh_d[lev];
            }

            alt = Altitude_d[lev];

            if (DeepModel){
                r2p = pow(altht + A, 2.0);
                r2m = pow(alt + A, 2.0);
                r2l = pow(althl + A, 2.0);
                rscale = A / (alt + A);
            }
            else{
                r2p = 1.0;
                r2m = 1.0;
                r2l = 1.0;
                rscale = 1.0;
            }

            nflxr_p = 0.0;
            nflxpt_p = 0.0;

            for (int k = 0; k < 3; k++){
                nflxr_p += rscale*(div_p[3 * 0 + k] * v_p[0 * 3 + k] +
                                    div_p[3 * 1 + k] * v_p[1 * 3 + k] +
                                   div_p[3 * 2 + k] * v_p[2 * 3 + k] +
                                   div_p[3 * 3 + k] * v_p[3 * 3 + k] +
                                   div_p[3 * 4 + k] * v_p[4 * 3 + k] +
                                   div_p[3 * 5 + k] * v_p[5 * 3 + k]);

                nflxpt_p += rscale*(div_p[3 * 0 + k] * v1_p[0 * 3 + k] * pt_p[0] +
                                    div_p[3 * 1 + k] * v1_p[1 * 3 + k] * pt_p[1] +
                                    div_p[3 * 2 + k] * v1_p[2 * 3 + k] * pt_p[2] +
                                    div_p[3 * 3 + k] * v1_p[3 * 3 + k] * pt_p[3] +
                                    div_p[3 * 4 + k] * v1_p[4 * 3 + k] * pt_p[4] +
                                    div_p[3 * 5 + k] * v1_p[5 * 3 + k] * pt_p[5]);
            }

            if (lev == 0){
                whl = 0.0;
                wht = Wh_d[id*(nv + 1) + lev + 1];
                whl2 = 0.0;
                wht2 = Wh_d[id*(nv + 1) + lev + 1] + Whk_d[id*(nv + 1) + lev + 1];
                phl = 0.0;
                pht = pth_d[id*(nv + 1) + lev + 1];
            }

            dz = altht - althl;
            dwdz = (wht*r2p - whl*r2l) / (dz*r2m);
            dwptdz = (wht2*pht*r2p - whl2*phl*r2l) / (dz*r2m);

            aux = -(nflxpt_p + dwptdz)*dt;
            r = Rhok_d[id*nv + lev] + Rho_d[id*nv + lev];
            aux += pt_d[id*nv + lev]*r;

            // Updates pressure
            p = P_Ref*pow(Rd*aux / P_Ref, Cp / Cv);
            pressure_d[id*nv + lev] = p - pressurek_d[id*nv + lev] + diffpr_d[id*nv + lev] * dt;

            // Updates density
            nflxr_p += dwdz;
            Rho_d[id*nv + lev] += (SlowRho_d[id*nv + lev] - nflxr_p)*dt;

            if (lev != nv - 1){
                althl = altht;
                altht = Altitudeh_d[lev + 2];
                whl = wht;
                wht = Wh_d[id*(nv + 1) + lev + 2];
                whl2 = wht2;
                wht2 = Wh_d[id*(nv + 1) + lev + 2] + Whk_d[id*(nv + 1) + lev + 2];
                phl = pht;
                pht = pth_d[id*(nv + 1) + lev + 2];
            }
        }
    }
}
