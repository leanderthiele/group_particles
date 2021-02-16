#include <memory>
#include <string>
#include <cstdio>

#include "halo_particles.hpp"
#include "illustris_fields.hpp"

namespace Y_Delta
{
    constexpr const size_t PartType = 0; // gas

    typedef GrpFields<IllustrisFields::GroupPos,
                      IllustrisFields::Group_M_Crit200,
                      IllustrisFields::Group_R_Crit200> GrpF;
    typedef PrtFields<IllustrisFields::Coordinates,
                      IllustrisFields::Masses,
                      IllustrisFields::InternalEnergy,
                      IllustrisFields::ElectronAbundance> PrtF;
    typedef AllFields<GrpF, PrtF> AF;

    typedef double grp_Y_t;
    typedef double grp_M_t;

    typedef CallbackUtils::chunk::Multi<AF>
        chunk;
    typedef CallbackUtils::name::Illustris<AF, PartType>
        name;
    typedef CallbackUtils::meta::Illustris<AF, PartType>
        meta;
    typedef CallbackUtils::select::GrpMassLowCutoff<AF, IllustrisFields::Group_M_Crit200>
        grp_select;
    typedef CallbackUtils::select::PrtAll<AF>
        prt_select;
    typedef CallbackUtils::radius::Simple<AF, IllustrisFields::Group_R_Crit200>
        grp_radius;
    typedef CallbackUtils::action::StoreGrpProperty<AF, IllustrisFields::Group_M_Crit200, grp_M_t>
        grp_store_M;
    typedef CallbackUtils::action::StorePrtHomogeneous<AF, grp_Y_t>
        prt_compute_Y;
} // namespace Y_Delta

struct Y_Delta_Callback :
    virtual public Callback<Y_Delta::AF>,
    public Y_Delta::chunk, public Y_Delta::name, public Y_Delta::meta,
    public Y_Delta::grp_select, public Y_Delta::prt_select,
    public Y_Delta::grp_radius,
    public Y_Delta::grp_store_M, public Y_Delta::prt_compute_Y
{// {{{
    Y_Delta_Callback () :
        Y_Delta::chunk { fgrp, grp_max_idx, fprt, prt_max_idx },
        Y_Delta::grp_select { Mmin },
        Y_Delta::grp_store_M { &grp_M },
        Y_Delta::prt_compute_Y { &grp_Y }
    { }

    // data (public so user can do something with them once they are assembled)
    std::vector<Y_Delta::grp_M_t> grp_M;
    std::vector<Y_Delta::grp_Y_t> grp_Y;

private :
    // for calculation of electron pressure
    static constexpr const float gamma = 5.0F/3.0F;
    static constexpr const float XH    = 0.76F;

    // group mass cutoff
    static constexpr const float Mmin = 1e3F;

    // files
    #define ROOT "/tigress/lthiele/Illustris_300-1_TNG/output/"
    static constexpr const char fgrp[]        = ROOT"groups_099/fof_subhalo_tab_099.%d.hdf5";
    static constexpr const size_t grp_max_idx = 599;
    static constexpr const char fprt[]        = ROOT"snapdir_099/snap_099.%d.hdf5";
    static constexpr const size_t prt_max_idx = 599;
    #undef ROOT

    // The next three functions are required -- they define the functionality of this class
    
    // computation of Compton-Y for a single gas particle
    void prt_insert (size_t grp_idx, const GrpProperties &grp, const PrtProperties &prt,
                     float Rsq, Y_Delta::grp_Y_t &data_item) override
    {
        auto m = prt.get<IllustrisFields::Masses>();
        auto e = prt.get<IllustrisFields::InternalEnergy>();
        auto x = prt.get<IllustrisFields::ElectronAbundance>();

        data_item += 2.0F * (1.0F+XH) / (1.0F+3.0F*XH+4.0F*XH*x)
                     * (gamma-1.0F) * m * e;
    }
};// }}}

template<typename T>
void vec_to_f (const std::vector<T> &v, const char *s)
{
    std::FILE *f = std::fopen(s, "wb");
    std::fwrite(v.data(), sizeof(T), v.size(), f);
    std::fclose(f);
}

int main ()
{
    Y_Delta_Callback y;
    
    halo_particles<> ( y );

    // save data to files
    vec_to_f<>(y.grp_M, "grp_M.bin");
    vec_to_f<>(y.grp_Y, "grp_Y.bin");

    return 0;
};
