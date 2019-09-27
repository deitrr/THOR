


__global__ void kelt9b(double *Mh_d,
                       double *pressure_d,
                       double *Rho_d,
                       double *temperature_d,
                       double *profx_dP_d,
                       double  Gravit,
                       double  Cp,
                       double  Rd,
                       double *Altitude_d,
                       double *Altitudeh_d,
                       double *lonlat_d,
                       double  time_step,
                       int     num) {

    int id  = blockIdx.x * blockDim.x + threadIdx.x;
    int nv  = gridDim.y;
    int lev = blockIdx.y;

    if (id < num) {
        double trad = 1000.0 * pressure_d[id * nv + lev] / 1e4;
        double Teq  = 4000.0;

        profx_dP_d[id * nv + lev] += Cp / (Cp - Rd) * Rho_d[id * nv + lev] * Rd
                                     * (temperature_d[id * nv + lev] - Teq) / trad;
    }
}
