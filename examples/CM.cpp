/* To specify which particle type to work on,
 * define one (and only one) of :
 *      GAS
 *      DM
 *      STARS
 *      BH
 * in the compilation script.
 */

#include "group_particles.hpp"
#include "common_fields.hpp"

namespace CM
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

    using GrpF = GrpFields<IllustrisFields::GroupPos,
                           IllustrisFields::Group_M_Crit200,
                           IllustrisFields::Group_R_Crit200>;
    using PrtF = PrtFields<IllustrisFields::Coordinates
                           #ifndef DM
                           , IllustrisFields::Masses
                           #endif // DM
                           >;

    using AF = AllFields<GrpF, PrtF>;

    class CMdata
    {
        using GrpProperties = typename Callback<AF>::GrpProperties;
        using PrtProperties = typename Callback<AF>::PrtProperties;
        using value_type = double;

        value_type sumM;
        value_type sumMr[3];

    public :
        CMdata () = delete;

        CMdata(const GrpProperties &)
        {
            sumM = 0.0;
            for (size_t ii=0; ii != 3; ++ii)
                sumMr[ii] = 0.0;
        }

        void prt_insert (size_t, const GrpProperties,
                         const PrtProperties &prt, coord_t)
        {
            const auto x = prt.coord();
            #ifndef DM
            const auto m = prt.get<IllustrisFields::Masses>();
            #else // DM
            const value_type m = 1.0;
            #endif // DM

            sumM += m;
            for (size_t ii=0; ii != 3; ++ii)
                sumMr[ii] += m * x[ii];
        }

        #ifdef DM
        void normalize_mass (value_type unit_mass)
        {
            sumM *= unit_mass;
            for (size_t ii=0; ii != 3; ++ii)
                sumMr[ii] *= unit_mass;
        }
        #endif

        void save (std::FILE *f) const
        {
            std::fwrite(&sumM, sizeof(value_type), 1, f);
            std::fwrite(sumMr, sizeof(value_type), 3, f);
        }
    };

    constexpr IllustrisFields::Group_R_Crit200::value_type Rscale = 1.0;
    using chunk = CallbackUtils::chunk::Multi<AF>;
    using name = CallbackUtils::name::Illustris<AF, PartType>;
    using meta = CallbackUtils::meta::Illustris<AF, PartType>;
    #ifdef DM
    using masstab = CallbackUtils::meta_init::IllustrisMassTable<AF>;
    #endif // DM
    using grp_select_M = CallbackUtils::select::LowCutoff<AF, IllustrisFields::Group_M_Crit200>;
    using grp_radius = CallbackUtils::radius::Simple<AF, IllustrisFields::Group_R_Crit200>;
    using prt_compute_CM = CallbackUtils::prt_action::StorePrtHomogeneous<AF, CMdata>;
} // namespace CM

struct CM_callback :
    virtual public Callback<CM::AF>,
    public CM::chunk, public CM::name, public CM::meta,
    #ifdef DM
    public CM::masstab,
    #endif // DM
    public CM::grp_select_M, public CM::grp_radius,
    public CM::prt_compute_CM
{
    CM_callback () :
        CM::chunk { fgrp, grp_max_idx, fprt, prt_max_idx },
        CM::grp_select_M { Mmin },
        CM::grp_radius { CM::Rscale },
        CM::prt_compute_CM { grp_CM }
    { }

    std::vector<CM::CMdata> grp_CM;

    #ifdef DM
    void normalize_mass ()
    {
        for (auto &CMitem : grp_CM)
            CMitem.normalize_mass(MassTable[CM::PartType]);
    }
    #endif

private :
    using Callback<CM::AF>::GrpProperties;
    using Callback<CM::AF>::PrtProperties;

    static constexpr const IllustrisFields::Group_M_Crit200::value_type Mmin = 1e3F;
    
    #define ROOT "/tigress/lthiele/Illustris_300-1_TNG/output"
    static constexpr const char fgrp[] = ROOT"groups_099/fof_subhalo_tab_099.%d.hdf5";
    static constexpr const char fprt[] = ROOT"snapdir_099/snap_099.%d.hdf5";
    static constexpr const size_t grp_max_idx = 599;
    static constexpr const size_t prt_max_idx = 599;
    #undef ROOT
};

int main ()
{
    CM_callback c;

    group_particles<> ( c );

    #ifdef DM
    c.normalize_mass();
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

    #define ROOT "CM_results"

    auto f = std::fopen(ROOT "/CM_" TYPE ".bin", "wb");
    for (const auto &CMitem : c.grp_CM)
        CMitem.save(f);
    std::fclose(f);

    #undef TYPE
    #undef ROOT


    return 0;
}
