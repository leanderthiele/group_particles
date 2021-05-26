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

// some file scope stuff where we store the previously computed center of mass results
double **CMdata;
size_t Ntypes;
size_t Ngroups;
size_t grp_idx = 0UL;

namespace Mtens
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

    class Mtensdata
    {
        using GrpProperties = typename Callback<AF>::GrpProperties;
        using PrtProperties = typename Callback<AF>::PrtProperties;
        using value_type = double;

        value_type tens[9];
        coord_t grp_coord[3];

    public :
        Mtensdata () = delete;

        Mtensdata(const GrpProperties &)
        {
            value_type sumM = 0.0;
            for (size_t ii=0; ii != 3; ++ii)
                grp_coord[ii] = 0.0;
            for (size_t ii=0; ii != 9; ++ii)
                tens[ii] = 0.0;

            for (size_t ii=0; ii != Ntypes; ++ii)
            {
                sumM += CMdata[ii][grp_idx * 4UL];
                for (size_t jj=0; jj != 3; ++jj)
                    grp_coord[jj] += (coord_t)CMdata[ii][grp_idx * 4UL + 1UL + jj];
            }

            for (size_t ii=0; ii != 3; ++ii)
                grp_coord[ii] /= (coord_t)sumM;

            // now advance the group index to the next group
            ++grp_idx;
        }

        void prt_insert (size_t, const GrpProperties,
                         const PrtProperties &prt, coord_t)
        {
            // position relative to the previously computed center of mass,
            //     consistently taking into account periodic boundary conditions
            const auto x = prt.coord(grp_coord);

            #ifndef DM
            const auto m = prt.get<IllustrisFields::Masses>();
            #else // DM
            const value_type m = 1.0;
            #endif // DM

            const auto Rsq = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];

            for (size_t ii=0; ii != 3; ++ii)
                for (size_t jj=0; jj !=3; ++jj)
                    tens[ii*3+jj] += (value_type)(m * x[ii] * x[jj] / Rsq);

        }

        #ifdef DM
        void normalize_mass (value_type unit_mass)
        {
            for (size_t ii=0; ii != 9; ++ii)
                tens[ii] *= unit_mass;
        }
        #endif

        void save (std::FILE *f) const
        {
            std::fwrite(tens, sizeof(value_type), 9, f);
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
    using prt_compute_Mtens = CallbackUtils::prt_action::StorePrtHomogeneous<AF, Mtensdata>;
} // namespace Mtens

struct Mtens_callback :
    virtual public Callback<Mtens::AF>,
    public Mtens::chunk, public Mtens::name, public Mtens::meta,
    #ifdef DM
    public Mtens::masstab,
    #endif // DM
    public Mtens::grp_select_M, public Mtens::grp_radius,
    public Mtens::prt_compute_Mtens
{
    Mtens_callback () :
        Mtens::chunk { fgrp, grp_max_idx, fprt, prt_max_idx },
        Mtens::grp_select_M { Mmin },
        Mtens::grp_radius { Mtens::Rscale },
        Mtens::prt_compute_Mtens { grp_Mtens }
    { }

    std::vector<Mtens::Mtensdata> grp_Mtens;

    #ifdef DM
    void normalize_mass ()
    {
        for (auto &Mtensitem : grp_Mtens)
            Mtensitem.normalize_mass(MassTable[Mtens::PartType]);
    }
    #endif

private :
    using Callback<Mtens::AF>::GrpProperties;
    using Callback<Mtens::AF>::PrtProperties;

    static constexpr const IllustrisFields::Group_M_Crit200::value_type Mmin = 1e3F;
    
    #define ROOT "/tigress/lthiele/Illustris_300-1_TNG/output"
    static constexpr const char fgrp[] = ROOT"groups_099/fof_subhalo_tab_099.%d.hdf5";
    static constexpr const char fprt[] = ROOT"snapdir_099/snap_099.%d.hdf5";
    static constexpr const size_t grp_max_idx = 599;
    static constexpr const size_t prt_max_idx = 599;
    #undef ROOT
};

int main (int argc, char **argv)
{
    // figure out which center of mass files to use
    Ntypes = (size_t)(argc - 2);
    Ngroups = (size_t)std::atoi(argv[1]);
    CMdata = (double **)std::malloc(Ntypes * sizeof(double *));
    for (size_t ii=0; ii != Ntypes; ++ii)
    {
        CMdata[ii] = (double *)std::malloc(Ngroups * 4UL * sizeof(double));
        auto f = std::fopen(argv[2+ii], "rb");
        auto count = std::fread(CMdata[ii], sizeof(double), 4UL * Ngroups, f);
        assert(count == 4UL * Ngroups);
        std::fclose(f);
    }

    Mtens_callback c;

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

    #define ROOT "Mtens_results"

    auto f = std::fopen(ROOT "/Mtens_" TYPE ".bin", "wb");
    for (const auto &Mtensitem : c.grp_Mtens)
        Mtensitem.save(f);
    std::fclose(f);

    #undef TYPE
    #undef ROOT


    return 0;
}
