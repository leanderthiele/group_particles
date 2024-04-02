/* To specify which particle type to work on,
 * define one (and only one) of :
 *      GAS
 *      DM
 *      STARS
 *      BH
 * in the compilation script.
 */

#include <string>

#include "group_particles.hpp"
#include "common_fields.hpp"

namespace dens_prof
{
    constexpr const uint8_t PartType =
        #if defined GAS 
        0
        #elif defined DM
        1
        #elif defined STARS
        4
        #elif defined BH
        5
        #else
        #   error("one of GAS, DM, STARS, BH must be defined.")
        #endif
        ;

    using GrpF = GrpFields<RockstarFields::pos,
                           RockstarFields::M200b,
                           RockstarFields::R200c>;
    using PrtF = PrtFields<IllustrisFields::Coordinates
                           #ifndef DM
                           , IllustrisFields::Masses
                           #endif // DM
                           >;

    using AF = AllFields<GrpF, PrtF>;

    class dens_prof_data
    {
        using GrpProperties = typename Callback<AF>::GrpProperties;
        using PrtProperties = typename Callback<AF>::PrtProperties;
        using value_type = double;

        static constexpr size_t N = 128;

        std::vector<value_type> dens;
        coord_t logRmin, logRmax, dlogR;

        // volume of spherical shells, used for save
        value_type shell_vol (size_t idx) const
        {// {{{
            // innermost shell goes right to the origin
            if (idx == 0UL)
                return 4.0 * M_PI / 3.0 * std::exp(3.0 * logRmin);
            
            // other shells are between two spheres
            value_type logR1 = (value_type)logRmin + (value_type)dlogR * (value_type)(idx-1UL);
            value_type logR2 = (value_type)logRmin + (value_type)dlogR * (value_type)(idx);
            return 4.0 * M_PI / 3.0 * ( std::exp(3.0*logR2) - std::exp(3.0*logR1) );
        }// }}}

    public :
        dens_prof_data () = delete;

        dens_prof_data(const GrpProperties &grp)
        {
            dens.resize(N, 0.0);
            logRmin = std::log(0.03F * grp.get<RockstarFields::R200c>());
            logRmax = std::log(2.50F * grp.get<RockstarFields::R200c>());
            dlogR = (logRmax - logRmin) / (coord_t)(N-1);
        }

        void prt_insert (size_t, const GrpProperties,
                         const PrtProperties &prt, coord_t Rsq)
        {
            coord_t logR = (coord_t)0.5 * std::log(Rsq);
            if (logR > logRmax)
                return;

            // figure out the index of the spherical shell the particle falls into
            size_t idx;

            if (logR < logRmin)
                idx = 0UL;
            else
                idx = 1UL + (size_t)((logR - logRmin) / dlogR);

            // safety check against numerical issues
            if (idx >= N)
                return;

            #ifndef DM
            const auto m = prt.get<IllustrisFields::Masses>();
            #else // DM
            const value_type m = 1.0;
            #endif // DM

            dens[idx] += m;
        }

        #ifdef DM
        void normalize_mass (value_type unit_mass)
        {
            for (size_t ii=0; ii != N; ++ii)
                dens[ii] *= unit_mass;
        }
        #endif

        void save (std::FILE *f)
        {
            for (size_t ii=0; ii != N; ++ii)
                dens[ii] /= shell_vol(ii);

            std::fwrite(dens.data(), sizeof(value_type), N, f);
        }
    };

    constexpr RockstarFields::R200c::value_type Rscale = 2.5;
    using grp_chunk = CallbackUtils::chunk::SingleGrp<AF>;
    using prt_chunk = CallbackUtils::chunk::MultiPrt<AF>;
    using name = CallbackUtils::name::Illustris<AF, PartType>;
    using meta = CallbackUtils::meta::Illustris<AF, PartType>;
    #ifdef DM
    using masstab = CallbackUtils::meta_init::IllustrisMassTable<AF>;
    #endif // DM
    using grp_radius = CallbackUtils::radius::Simple<AF, RockstarFields::R200c>;
    using prt_compute_dens_prof = CallbackUtils::prt_action::StorePrtHomogeneous<AF, dens_prof_data>;
} // namespace dens_prof

struct dens_prof_callback :
    virtual public Callback<dens_prof::AF>,
    public dens_prof::grp_chunk, public dens_prof::prt_chunk,
    public dens_prof::name, public dens_prof::meta,
    #ifdef DM
    public dens_prof::masstab,
    #endif // DM
    public dens_prof::grp_radius,
    public dens_prof::prt_compute_dens_prof
{
    dens_prof_callback () :
        dens_prof::grp_chunk { fgrp },
        dens_prof::prt_chunk { fprt, prt_max_idx },
        dens_prof::grp_radius { dens_prof::Rscale },
        dens_prof::prt_compute_dens_prof { grp_dens_prof }
    { }

    std::vector<dens_prof::dens_prof_data> grp_dens_prof;

    #ifdef DM
    void normalize_mass ()
    {
        for (auto &dens_prof_item : grp_dens_prof)
            dens_prof_item.normalize_mass(MassTable[dens_prof::PartType]);
    }
    #endif

private :
    using Callback<dens_prof::AF>::GrpProperties;
    using Callback<dens_prof::AF>::PrtProperties;

    #define ROOT "/tigress/lthiele/Illustris_300-1_TNG/"
    static constexpr char fgrp[] = ROOT "rockstar/out_99.hdf5";
    static constexpr char fprt[] = ROOT "output/snapdir_099/snap_099.%d.hdf5";
    static constexpr size_t prt_max_idx = 599;
    #undef ROOT
};

int main ()
{
    dens_prof_callback d;

    group_particles<> ( d );

    #ifdef DM
    d.normalize_mass();
    #endif

    #if defined GAS
    #   define TYPE "GAS"
    #elif defined DM
    #   define TYPE "DM"
    #elif defined STARS
    #   define TYPE "STARS"
    #elif defined BH
    #   define TYPE "BH"
    #endif

    static constexpr const char outfile[] = "dens_prof_for_pdfbaryon" TYPE ".bin";
    auto f = std::fopen(outfile, "w");
    for (auto &dens_prof_item : d.grp_dens_prof)
        dens_prof_item.save(f);
    std::fclose(f);

    #undef TYPE

    return 0;
}
