#include <memory>
#include <string>
#include <cstdio>
#include <cmath>

#include "halo_particles.hpp"
#include "illustris_fields.hpp"

namespace y_prof
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

    // forward declare these classes here, we'll implement them later
    class YProfile;

    typedef double grp_M_t;
    typedef double grp_R_t;
    typedef YProfile grp_Y_t;

    // maximum radius (in units of R_200c) -- from Battaglia's plots
    constexpr IllustrisFields::Group_R_Crit200::value_type Rscale = 2.5;

    typedef CallbackUtils::chunk::Multi<AF>
        chunk;
    typedef CallbackUtils::name::Illustris<AF, PartType>
        name;
    typedef CallbackUtils::meta::Illustris<AF, PartType>
        meta;
    typedef CallbackUtils::select::GrpMassLowCutoff<AF, IllustrisFields::Group_M_Crit200>
        grp_select;
    typedef CallbackUtils::radius::Simple<AF, IllustrisFields::Group_R_Crit200>
        grp_radius;
    typedef CallbackUtils::action::StoreGrpProperty<AF, IllustrisFields::Group_M_Crit200, grp_M_t>
        grp_store_M;
    typedef CallbackUtils::action::StoreGrpProperty<AF, IllustrisFields::Group_R_Crit200, grp_R_t>
        grp_store_R;
    typedef CallbackUtils::action::StorePrtHomogeneous<AF, grp_Y_t>
        prt_compute_Y;
} // namespace Y_Delta

// implementation of the YProfile class
// We inherit from the Callback template specialization so we can use the
// GrpProperties and PrtProperties types
class y_prof::YProfile :
    virtual public Callback<y_prof::AF>
{// {{{
    // all internal calculations regarding Y in this type
    typedef double value_type;

    // for calculation of electron pressure
    static constexpr const value_type gamma = 5.0F/3.0F;
    static constexpr const value_type XH    = 0.76F;

    // number of sample points
    static constexpr size_t N = 128;
    std::vector<value_type> data;
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
        data.resize(N, 0.0);
        logRmin = std::log(0.03F * grp.get<IllustrisFields::Group_R_Crit200>());
        logRmax = std::log(2.50F * grp.get<IllustrisFields::Group_R_Crit200>());
        dlogR = (logRmin - logRmax) / (coord_t)N;
    }// }}}
    // add a particle
    void add (const PrtProperties &prt, coord_t Rsq)
    {// {{{
        // figure out the index of the spherical shell the particle falls into
        size_t idx = (size_t)((0.5 * std::log(Rsq) - logRmin) / dlogR);
        if (idx >= N) // use wrapping of unsigned here
            return;

        // load the required properties of this particle
        auto m = (value_type)prt.get<IllustrisFields::Masses>();
        auto e = (value_type)prt.get<IllustrisFields::InternalEnergy>();
        auto x = (value_type)prt.get<IllustrisFields::ElectronAbundance>();

        // this is electron pressure * particle volume
        auto Y = 2.0 * (1.0+XH) / (1.0+3.0*XH+4.0*XH*x)
                     * (gamma-1.0) * m * e;

        data[idx] += Y;
    }// }}}
    // it is assumed that this function is called after all data has been added,
    // because it performs the normalization by volume
    void save (std::FILE *f) const
    {// {{{
        for (size_t ii=0; ii != N; ++ii)
            data[ii] /= shell_vol(ii);

        std::fwrite(data.data(), sizeof(value_type), N, f);
    }// }}}
};// }}}

struct y_prof_callback :
    virtual public Callback<y_prof::AF>,
    public y_prof::chunk, public y_prof::name, public y_prof::meta,
    public y_prof::grp_select, public y_prof::grp_radius,
    public y_prof::grp_store_M,
    public y_prof::grp_store_R,
    public y_prof::prt_compute_Y
{// {{{
    y_prof_callback () :
        y_prof::chunk { fgrp, grp_max_idx, fprt, prt_max_idx },
        y_prof::grp_select { Mmin },
        y_prof::grp_radius { y_prof::Rscale },
        y_prof::grp_store_M { grp_M },
        y_prof::grp_store_R { grp_R },
        y_prof::prt_compute_Y { grp_Y }
    { }

    // data (public so user can do something with them once they are assembled)
    std::vector<y_prof::grp_M_t> grp_M;
    std::vector<y_prof::grp_R_t> grp_R;
    std::vector<y_prof::grp_Y_t> grp_Y;

private :

    // group mass cutoff
    static constexpr const IllustrisFields::Group_M_Crit200::value_type Mmin = 1e3F;

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
    
    halo_particles<> ( y );

    // save data to files
    #define ROOT "y_prof_results_Feb23"
    vec_to_f<>(y.grp_M, ROOT"/grp_M200c.bin");
    vec_to_f<>(y.grp_R, ROOT"/grp_R200c.bin");
    {
        auto f = std::fopen(ROOT"/grp_yprof.bin", "wb");
        for (const auto &prof : y.grp_Y)
            prof.save(f);
        std::fclose(f);
    }
    #undef ROOT

    return 0;
};
