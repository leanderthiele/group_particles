/*! @file ne_prof.cpp
 *
 * @brief Example for using the group_particles code,
 *        computing the electron density profiles of halos in an Illustris TNG simulation.
 */

#include <memory>
#include <string>
#include <cstdio>
#include <cmath>

#include "group_particles.hpp"
#include "common_fields.hpp"

/*! @brief assembles some types to use for the electron density profile calculation.
 */
namespace ne_prof
{// {{{
    constexpr const uint8_t PartType = 0; /*!< @brief gas particles */

    /*! @brief for the groups, we want to read their position and 200c masses and radii. */
    using GrpF = GrpFields<IllustrisFields::GroupPos,
                           IllustrisFields::Group_M_Crit200,
                           IllustrisFields::Group_R_Crit200>;

    /*! @brief for the particles, we want to read their position and the quanitities
     *         required to compute the electron density.
     */
    using PrtF = PrtFields<IllustrisFields::Coordinates,
                           IllustrisFields::Masses,
                           IllustrisFields::ElectronAbundance>;

    using AF = AllFields<GrpF, PrtF>; /*!< @brief bundle the group and particle fields */

    /*! @brief stores a single electron electron density profile */
    class NeProfile
    {// {{{
        using GrpProperties = typename Callback<AF>::GrpProperties;
        using PrtProperties = typename Callback<AF>::PrtProperties;

        // all internal calculations regarding ne in this type
        using value_type = double;

        // number of sample points
        static constexpr size_t N = 128;

        std::vector<value_type> e_density;
        std::vector<size_t> num_part;
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
        // we can only construct from a Group
        NeProfile () = delete;

        /*! @brief constructs from group properties (we need to read the radius)
         *
         * As explained in the documentation for #CallbackUtils::prt_action::StorePrtHomogeneous,
         * if the data type stored for each group has a constructor from a
         * `const GrpProperties &`, this constructor will be called during
         * the calls to #Callback::grp_action.
         * We use this fact to read the group radius.
         */
        NeProfile (const GrpProperties &grp)
        {// {{{
            e_density.resize(N, 0.0);
            num_part.resize(N, 0UL);
            logRmin = std::log(0.03F * grp.get<IllustrisFields::Group_R_Crit200>());
            logRmax = std::log(2.50F * grp.get<IllustrisFields::Group_R_Crit200>());
            dlogR = (logRmax - logRmin) / (coord_t)(N-1);
        }// }}}

        /*! @brief adds a particle to the profile
         *
         * This function computes the electron density associated with a particle
         * and adds it (weighted with its volume) to the appropriate profile bin.
         *
         * Note that this exact signature and function name is required,
         * as explained in the documentation for #CallbackUtils::prt_action::StorePrtHomogeneous.
         */
        void prt_insert (size_t, const GrpProperties &,
                         const PrtProperties &prt, coord_t Rsq)
        {// {{{
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

            // load the required properties of this particle
            auto m = (value_type)prt.get<IllustrisFields::Masses>();
            auto x = (value_type)prt.get<IllustrisFields::ElectronAbundance>();

            // proton mass is in 1e10 Msun/h -- note that h is hardcoded here, but ok
            static constexpr const value_type XH = 0.76,
                                              mprot = 0.6774*1.6726219/1.98847 * 1e-67;

            // this is electron density times particle volume
            auto ne = x * XH * m / mprot;

            e_density[idx] += ne;
            num_part[idx]++;
        }// }}}

        /*! @brief append this electron density profile to file.
         *
         * @attention it is assumed that this instance is "dead" after this
         *            function is called!
         */
        void save (std::FILE *fne, std::FILE *fnum_part)
        {// {{{
            for (size_t ii=0; ii != N; ++ii)
                e_density[ii] /= shell_vol(ii);

            std::fwrite(e_density.data(), sizeof(value_type), N, fne);
            std::fwrite(num_part.data(), sizeof(size_t), N, fnum_part);
        }// }}}
    };// }}}
    
    using grp_M_t = double; /*!< @brief type to which M200c will be converted */
    using grp_R_t = double; /*!< @brief type to which R200c will be converted */
    using grp_Ne_t = NeProfile; /*!< @brief type the electron density profiles will be stored in */

    // maximum radius (in units of R_200c) -- from Battaglia's plots
    constexpr IllustrisFields::Group_R_Crit200::value_type Rscale = 2.5;

    /*! @brief group and particle files are in multiple chunks */
    using chunk = CallbackUtils::chunk::Multi<AF>;

    /*! @brief internal file layout is according to Illustris naming conventions */
    using name = CallbackUtils::name::Illustris<AF, PartType>;

    /*! @brief meta-data required by code are according to Illustris conventions */
    using meta = CallbackUtils::meta::Illustris<AF, PartType>;

    /*! @brief apply a minimum-mass selection on groups */
    using grp_select_M =  CallbackUtils::select::LowCutoff<AF, IllustrisFields::Group_M_Crit200>;

    /*! @brief apply a minimum-radius selection on groups (exclude those with negative radius) */
    using grp_select_R = CallbackUtils::select::LowCutoff<AF, IllustrisFields::Group_R_Crit200>;

    /*! @brief compute group radius as a multiple of R200c */
    using grp_radius =  CallbackUtils::radius::Simple<AF, IllustrisFields::Group_R_Crit200>;

    /*! @brief store M200c in a vector */
    using grp_store_M = CallbackUtils::grp_action::StoreGrpProperty<AF, IllustrisFields::Group_M_Crit200, grp_M_t>;

    /*! @brief store R200c in a vector */
    using grp_store_R = CallbackUtils::grp_action::StoreGrpProperty<AF, IllustrisFields::Group_R_Crit200, grp_R_t>;

    /*! @brief store a #ne_prof::NeProfile instance for each group */
    using prt_compute_Ne = CallbackUtils::prt_action::StorePrtHomogeneous<AF, grp_Ne_t>;


} // namespace ne_prof }}}


/*! @brief implementation of the #Callback base for our application
 */
struct ne_prof_callback :
    virtual public Callback<ne_prof::AF>,
    public ne_prof::chunk, public ne_prof::name, public ne_prof::meta,
    public ne_prof::grp_select_M, public ne_prof::grp_select_R,
    public ne_prof::grp_radius,
    public ne_prof::grp_store_M, public ne_prof::grp_store_R,
    public ne_prof::prt_compute_Ne
{// {{{
    ne_prof_callback () :
        ne_prof::chunk { fgrp, grp_max_idx, fprt, prt_max_idx },
        ne_prof::grp_select_M { Mmin },
        ne_prof::grp_select_R { Rmin },
        ne_prof::grp_radius { ne_prof::Rscale },
        ne_prof::grp_store_M { grp_M },
        ne_prof::grp_store_R { grp_R },
        ne_prof::prt_compute_Ne { grp_Ne }
    { }

    std::vector<ne_prof::grp_M_t> grp_M; /*!< @brief M200c data vector */
    std::vector<ne_prof::grp_R_t> grp_R; /*!< @brief R200c data vector */
    std::vector<ne_prof::grp_Ne_t> grp_Ne; /*!< @brief electron density profile data vector */

private :
    using Callback<ne_prof::AF>::GrpProperties;
    using Callback<ne_prof::AF>::PrtProperties;

    // group mass cutoff
    static constexpr const IllustrisFields::Group_M_Crit200::value_type Mmin = 1e4F;

    // group radius cutoff
    static constexpr const IllustrisFields::Group_R_Crit200::value_type Rmin = 0.0F;
    
    // files
    #define ROOT "/tigress/lthiele/Illustris_300-1_TNG/output/"
    static constexpr const char fgrp[]        = ROOT"groups_099/fof_subhalo_tab_099.%d.hdf5";
    static constexpr const size_t grp_max_idx = 599;
    static constexpr const char fprt[]        = ROOT"snapdir_099/snap_099.%d.hdf5";
    static constexpr const size_t prt_max_idx = 599;
    #undef ROOT
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
    ne_prof_callback n;
    
    group_particles<> ( n );

    // save data to files
    #define ROOT "ne_prof_results_Jun22"
    vec_to_f<>(n.grp_M, ROOT"/grp_M200c.bin");
    vec_to_f<>(n.grp_R, ROOT"/grp_R200c.bin");
    {
        auto fne = std::fopen(ROOT"/grp_ne_prof.bin", "wb");
        auto fnum_part = std::fopen(ROOT"/grp_num_part_prof.bin", "wb");
        for (auto &prof : n.grp_Ne)
            prof.save(fne, fnum_part);
        std::fclose(fne);
        std::fclose(fnum_part);
    }
    #undef ROOT

    return 0;
};
