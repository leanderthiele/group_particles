#include <memory>
#include <string>
#include <cstdio>
#include <cmath>

#include "group_particles.hpp"
#include "common_fields.hpp"

namespace y_prof
{// {{{
    constexpr const size_t PartType = 0; // gas

    typedef GrpFields<IllustrisFields::GroupCM,
                      IllustrisFields::Group_M_Crit200,
                      IllustrisFields::Group_R_Crit200> GrpF;
    typedef PrtFields<IllustrisFields::Coordinates,
                      IllustrisFields::Masses,
                      IllustrisFields::InternalEnergy,
                      IllustrisFields::ElectronAbundance> PrtF;
    typedef AllFields<GrpF, PrtF> AF;

    // forward declare these classes here, we'll implement them later
    class YProfile;
    class grp_store_P;

    typedef double grp_M_t;
    typedef double grp_R_t;
    typedef double grp_P_t;
    typedef YProfile grp_Y_t;

    // maximum radius (in units of R_200c) -- from Battaglia's plots
    constexpr IllustrisFields::Group_R_Crit200::value_type Rscale = 2.5;

    typedef CallbackUtils::chunk::Multi<AF>
        chunk;
    typedef CallbackUtils::name::Illustris<AF, PartType>
        name;
    typedef CallbackUtils::meta::Illustris<AF, PartType>
        meta;
    typedef CallbackUtils::select::LowCutoff<AF, IllustrisFields::Group_M_Crit200>
        grp_select_M;
    typedef CallbackUtils::select::LowCutoff<AF, IllustrisFields::Group_R_Crit200>
        grp_select_R;
    typedef CallbackUtils::radius::Simple<AF, IllustrisFields::Group_R_Crit200>
        grp_radius;
    typedef CallbackUtils::grp_action::StoreGrpProperty<AF, IllustrisFields::Group_M_Crit200, grp_M_t>
        grp_store_M;
    typedef CallbackUtils::grp_action::StoreGrpProperty<AF, IllustrisFields::Group_R_Crit200, grp_R_t>
        grp_store_R;
    typedef CallbackUtils::prt_action::StorePrtHomogeneous<AF, grp_Y_t>
        prt_compute_Y;
} // namespace y_prof }}}

// implementation of the YProfile class
// We inherit from the Callback template specialization so we can use the
// GrpProperties and PrtProperties types
class y_prof::YProfile
{// {{{
    typedef Callback<y_prof::AF>::GrpProperties GrpProperties;
    typedef Callback<y_prof::AF>::PrtProperties PrtProperties;

    // all internal calculations regarding Y in this type
    typedef double value_type;

    // for calculation of electron pressure
    static constexpr const value_type gamma = 5.0F/3.0F;
    static constexpr const value_type XH    = 0.76F;

    // number of sample points
    static constexpr size_t N = 128;
    std::vector<value_type> pressure;
    std::vector<size_t> num_part;
    coord_t logRmin, logRmax, dlogR;

    // volume of spherical shells, used for save
    value_type shell_vol (size_t idx) const
    {// {{{
        auto logR1 = (value_type)logRmin + (value_type)dlogR * (value_type)idx;
        auto logR2 = (value_type)logRmin + (value_type)dlogR * (value_type)(idx+1UL);
        return 4.0 * M_PI / 3.0 * ( std::exp(3.0*logR2) - std::exp(3.0*logR1) );
    }// }}}
public :
    // we can only construct from a Group
    YProfile () = delete;
    // constructor from a group -- this is the required signature
    YProfile (const GrpProperties &grp)
    {// {{{
        pressure.resize(N, 0.0);
        num_part.resize(N, 0UL);
        logRmin = std::log(0.03F * grp.get<IllustrisFields::Group_R_Crit200>());
        logRmax = std::log(2.50F * grp.get<IllustrisFields::Group_R_Crit200>());
        dlogR = (logRmax - logRmin) / (coord_t)N;
    }// }}}
    // add a particle
    void add (const PrtProperties &prt, coord_t Rsq)
    {// {{{
        coord_t logR = 0.5 * std::log(Rsq);
        if (logR < logRmin || logR > logRmax)
            return;

        // figure out the index of the spherical shell the particle falls into
        size_t idx = (size_t)((logR - logRmin) / dlogR);
        if (idx >= N) // safety check against numerical issues
            return;

        // load the required properties of this particle
        auto m = (value_type)prt.get<IllustrisFields::Masses>();
        auto e = (value_type)prt.get<IllustrisFields::InternalEnergy>();
        auto x = (value_type)prt.get<IllustrisFields::ElectronAbundance>();

        // this is electron pressure * particle volume
        auto Y = 2.0 * (1.0+XH) / (1.0+3.0*XH+4.0*XH*x)
                     * (gamma-1.0) * m * e;

        pressure[idx] += Y;
        ++num_part[idx];
    }// }}}
    // it is assumed that this function is called after all data has been added,
    // because it performs the normalization by volume
    void save (std::FILE *fpressure, std::FILE *fnum_part)
    {// {{{
        for (size_t ii=0; ii != N; ++ii)
            pressure[ii] /= shell_vol(ii);

        std::fwrite(pressure.data(), sizeof(value_type), N, fpressure);
        std::fwrite(num_part.data(), sizeof(size_t), N, fnum_part);
    }// }}}
};// }}}

class y_prof::grp_store_P :
    virtual public Callback<y_prof::AF>,
    public CallbackUtils::meta_init::IllustrisCosmology<AF>,
    public CallbackUtils::grp_action::StoreGrpHomogeneous<AF, y_prof::grp_P_t>
{// {{{
    static constexpr y_prof::grp_P_t GNewton    = 4.30091e4, // (kpc/1e10 Msun) (km/s)^2
                                     rho_crit_0 = 2.775e-8;  // 1e10 Msun/h / (kpc/h)^3

    y_prof::grp_P_t rho_crit () const
    {
        y_prof::grp_P_t a3 = Time * Time * Time;
        return rho_crit_0 * (Omega0/a3 + OmegaLambda);
    }
    y_prof::grp_P_t grp_reduce (const GrpProperties &grp) const override
    {
        auto m = (y_prof::grp_P_t)grp.get<IllustrisFields::Group_M_Crit200>();
        auto r = (y_prof::grp_P_t)grp.get<IllustrisFields::Group_R_Crit200>();
        return 100.0 * GNewton * m * rho_crit() * OmegaBaryon / Omega0 / r;
    }
public :
    grp_store_P (std::vector<y_prof::grp_P_t> &data) :
        CallbackUtils::grp_action::StoreGrpHomogeneous<AF, y_prof::grp_P_t> { data }
    { }
};// }}}

struct y_prof_callback :
    virtual public Callback<y_prof::AF>,
    public y_prof::chunk, public y_prof::name, public y_prof::meta,
    public y_prof::grp_select_M, public y_prof::grp_select_R,
    public y_prof::grp_radius,
    public y_prof::grp_store_M, public y_prof::grp_store_R, public y_prof::grp_store_P,
    public y_prof::prt_compute_Y
{// {{{
    y_prof_callback () :
        y_prof::chunk { fgrp, grp_max_idx, fprt, prt_max_idx },
        y_prof::grp_select_M { Mmin },
        y_prof::grp_select_R { Rmin },
        y_prof::grp_radius { y_prof::Rscale },
        y_prof::grp_store_M { grp_M },
        y_prof::grp_store_R { grp_R },
        y_prof::grp_store_P { grp_P },
        y_prof::prt_compute_Y { grp_Y }
    { }

    // data (public so user can do something with them once they are assembled)
    std::vector<y_prof::grp_M_t> grp_M;
    std::vector<y_prof::grp_R_t> grp_R;
    std::vector<y_prof::grp_P_t> grp_P;
    std::vector<y_prof::grp_Y_t> grp_Y;

private :

    // group mass cutoff
    static constexpr const IllustrisFields::Group_M_Crit200::value_type Mmin = 1e3F;

    // group radius cutoff
    static constexpr const IllustrisFields::Group_R_Crit200::value_type Rmin = 0.0F;
    
    // files
    #define ROOT "/tigress/lthiele/Illustris_300-1_TNG/output/"
    static constexpr const char fgrp[]        = ROOT"groups_099/fof_subhalo_tab_099.%d.hdf5";
    static constexpr const size_t grp_max_idx = 599;
    static constexpr const char fprt[]        = ROOT"snapdir_099/snap_099.%d.hdf5";
    static constexpr const size_t prt_max_idx = 599;
    #undef ROOT
    
    // computation of Compton-Y for a single gas particle
    void prt_insert (size_t grp_idx, const GrpProperties &grp, const PrtProperties &prt,
                     coord_t Rsq, y_prof::grp_Y_t &data_item) override
    {// {{{
        data_item.add(prt, Rsq);
    }// }}}
};// }}}

template<typename T>
void vec_to_f (const std::vector<T> &v, const char *s)
{
    auto f = std::fopen(s, "wb");
    std::fwrite(v.data(), sizeof(T), v.size(), f);
    std::fclose(f);
}

int main ()
{
    y_prof_callback y;
    
    group_particles<> ( y );

    // save data to files
    #define ROOT "y_prof_results_Feb23"
    vec_to_f<>(y.grp_M, ROOT"/grp_M200c.bin");
    vec_to_f<>(y.grp_R, ROOT"/grp_R200c.bin");
    vec_to_f<>(y.grp_P, ROOT"/grp_P200c.bin");
    {
        auto fpressure = std::fopen(ROOT"/grp_pressure_prof.bin", "wb");
        auto fnum_part = std::fopen(ROOT"/grp_num_part_prof.bin", "wb");
        for (auto &prof : y.grp_Y)
            prof.save(fpressure, fnum_part);
        std::fclose(fpressure);
        std::fclose(fnum_part);
    }
    #undef ROOT

    return 0;
};
