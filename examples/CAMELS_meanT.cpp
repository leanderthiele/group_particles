/*
 * Define either FOR_ILLUSTRIS or FOR_SIMBA in the compilation script.
 *
 * Call with the command line arguments
 *      [1] group hdf5 file
 *      [2] particle hdf5 file
 *      [3] output directory (excluding final slash)
 */

#include <memory>
#include <string>
#include <cstdio>

#include "group_particles.hpp"
#include "common_fields.hpp"

#if !defined(FOR_ILLUSTRIS) && !defined(FOR_SIMBA)
#   error "Either FOR_SIMBA or FOR_ILLUSTRIS must be defined."
#endif

#ifdef FOR_ILLUSTRIS
#   define SimFields IllustrisFields
#else
#   define SimFields SIMBAFields
#endif

namespace Y_Delta
{
    constexpr const size_t PartType = 0; // gas

    typedef GrpFields<IllustrisFields::GroupPos,
                      IllustrisFields::Group_M_Crit200,
                      // NOTE using R500 here as maybe more useful for X-ray
                      IllustrisFields::Group_R_Crit500> GrpF;
    typedef PrtFields<IllustrisFields::Coordinates,
                      IllustrisFields::InternalEnergy,
                      IllustrisFields::ElectronAbundance> PrtF;
    typedef AllFields<GrpF, PrtF> AF;

    typedef std::pair<double, size_t> grp_Y_t;
    typedef double grp_M_t;

    typedef CallbackUtils::chunk::Single<AF>
        chunk;
    typedef CallbackUtils::name::Illustris<AF, PartType>
        name;
    typedef CallbackUtils::meta::Illustris<AF, PartType>
        meta;
    typedef CallbackUtils::select::LowCutoff<AF, IllustrisFields::Group_M_Crit200>
        grp_select;
    typedef CallbackUtils::radius::Simple<AF, IllustrisFields::Group_R_Crit500>
        grp_radius;
    typedef CallbackUtils::grp_action::StoreGrpProperty<AF, IllustrisFields::Group_M_Crit200, grp_M_t>
        grp_store_M;
    typedef CallbackUtils::prt_action::StorePrtHomogeneous<AF, grp_Y_t>
        prt_compute_Y;
} // namespace Y_Delta

struct Y_Delta_Callback :
    virtual public Callback<Y_Delta::AF>,
    public Y_Delta::chunk, public Y_Delta::name, public Y_Delta::meta,
    public Y_Delta::grp_select, public Y_Delta::grp_radius,
    public Y_Delta::grp_store_M, public Y_Delta::prt_compute_Y
{// {{{
    Y_Delta_Callback (const char *fgrp, const char *fprt) :
        Y_Delta::chunk { fgrp, fprt },
        Y_Delta::grp_select { Mmin },
        Y_Delta::grp_store_M { grp_M },
        Y_Delta::prt_compute_Y { grp_Y }
    { }

    // data (public so user can do something with them once they are assembled)
    std::vector<Y_Delta::grp_M_t> grp_M;
    std::vector<Y_Delta::grp_Y_t> grp_Y;

private :
    using GrpProperties = typename Callback<Y_Delta::AF>::GrpProperties;
    using PrtProperties = typename Callback<Y_Delta::AF>::PrtProperties;
    
    // for calculation of electron pressure
    static constexpr const float gamma = 5.0F/3.0F;
    static constexpr const float XH    = 0.76F;

    // group mass cutoff
    static constexpr const float Mmin = 1e12F;

    // computation of temperature for a single gas particle
    void prt_insert (size_t grp_idx, const GrpProperties &grp, const PrtProperties &prt,
                     float Rsq, Y_Delta::grp_Y_t &data_item) override
    {
        auto e = prt.get<IllustrisFields::InternalEnergy>();
        auto x = prt.get<IllustrisFields::ElectronAbundance>();

        double mu = 4.0/(1.0+3.0*XH+4.0*XH*x);
        double T = (gamma-1.0) * e * mu;

        // get to Kelvin
        T *= 1.211475e-10; // this is 1e6 * mproton[SI] / kB[SI]

        if (T > 1e5)
        {
            data_item.first += T;
            ++data_item.second;
        }
    }
};// }}}

template<typename T>
void vec_to_f (const std::vector<T> &v, const char *s)
{
    std::FILE *f = std::fopen(s, "wb");
    std::fwrite(v.data(), sizeof(T), v.size(), f);
    std::fclose(f);
}

int main (int argc, char **argv)
{
    char *fgrp = argv[1];
    char *fprt = argv[2];
    char *dout = argv[3];

    Y_Delta_Callback y { fgrp, fprt };
    
    group_particles<> ( y );

    char fname_buffer[512];

    // save data to files
    std::sprintf(fname_buffer, "%s/grp_M200c.bin", dout);
    vec_to_f<>(y.grp_M, fname_buffer);
    {
        std::sprintf(fname_buffer, "%s/grp_meanT.bin", dout);
        std::FILE *f = std::fopen(fname_buffer, "w");
        for (auto &d : y.grp_Y)
        {
            double tmp = d.first / (double)d.second;
            std::fwrite(&tmp, sizeof(double), 1, f);
        }
    }

    return 0;
};
