/*! @file mass_prof.cpp
 *
 * @brief Example for using the group_particles code,
 *        computing enclosed mass profiles for various particle types.
 *
 * Define either FOR_ILLUSTRIS or FOR_SIMBA in the compilation script.
 *
 * Define one of GAS, DM, STARS, BH in the compilation script.
 *
 * Call with the command line arguments
 *      [1] group hdf5 file
 *      [2] particle hdf5 file
 *      [3] output directory (excluding final slash)
 */

#include <memory>
#include <string>
#include <cstdio>
#include <cmath>

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

// Illustris DM doesn't have a Masses field
#if defined FOR_SIMBA || !defined DM
#   define MASSES_AVAIL
#endif

/*! @brief assembles some types to use for the profile calculations.
 */
namespace mass_prof
{// {{{
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

    /*! @brief for the groups, we want to read their position and 200c masses and radii. */
    using GrpF = GrpFields<SimFields::GroupPos,
                           SimFields::Group_M_Crit200,
                           SimFields::Group_R_Crit200>;
    /*! @brief for the particles, we want to read their position and the quantities
     *         required to compute the electron pressure, electron density, and temperature.
     */
    using PrtF = PrtFields<SimFields::Coordinates
                           #ifdef MASSES_AVAIL
                           , SimFields::Masses
                           #endif // MASSES_AVAIL
                           >;

    using AF = AllFields<GrpF, PrtF>; /*!< @brief bundle the group and particle fields */

    /*! @brief stores profiles for a single halo */
    class Profile
    {// {{{
        using GrpProperties = typename Callback<AF>::GrpProperties;
        using PrtProperties = typename Callback<AF>::PrtProperties;

        // all internal calculations in this type
        using value_type = double;

        // number of sample points
        static constexpr size_t N = 128;

        std::vector<value_type> encl_mass;
        coord_t logRmin, logRmax, dlogR;

    public :
        // we can only construct from a Group
        Profile () = delete;

        /*! @brief constructs from group properties (we need to read the radius)
         *
         * As explained in the documentation for #CallbackUtils::prt_action::StorePrtHomogeneous,
         * if the data type stored for each group has a constructor from a
         * `const GrpProperties &`, this constructor will be called during
         * the calls to #Callback::grp_action.
         * We use this fact to read the group radius.
         */
        Profile (const GrpProperties &grp)
        {// {{{
            encl_mass.resize(N, 0.0);

            logRmin = std::log(0.03F * grp.get<SimFields::Group_R_Crit200>());
            logRmax = std::log(2.50F * grp.get<SimFields::Group_R_Crit200>());

            dlogR = (logRmax - logRmin) / (coord_t)(N-1);
        }// }}}

        /*! @brief adds a particle to the profiles
         *
         * This function computes the electron pressure, electron density, and temperature
         * associated with a particle and adds it (weighted with its volume)
         * to the appropriate profile bin.
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
            #ifdef MASSES_AVAIL
            auto m = (value_type)prt.get<SimFields::Masses>();
            #else // MASSES_AVAIL
            value_type m = 1.0;
            #endif // MASSES_AVAIL

            for (size_t ii=idx; ii < N; ++ii)
                encl_mass[ii] += m;
        }// }}}

        #ifndef MASSES_AVAIL
        void normalize_mass (value_type unit_mass)
        {// {{{
            for (size_t ii=0; ii != N; ++ii)
                encl_mass[ii] *= unit_mass;
        }// }}}
        #endif

        /*! @brief append these profiles to files.
         *
         * @attention it is assumed that this instance is "dead" after this
         *            function is called!
         */
        void save (std::FILE *f) const
        {// {{{
            std::fwrite(encl_mass.data(), sizeof(value_type), N, f);
        }// }}}
    };// }}}
    
    using grp_M_t = double; /*!< @brief type to which M200c will be converted */
    using grp_R_t = double; /*!< @brief type to which R200c will be converted */
    using grp_P_t = double; /*!< @brief type to which P200c will be converted */
    using grp_prof_t = Profile; /*!< @brief type the pressure profiles will be stored in */

    // maximum radius (in units of R_200c) -- from Battaglia's plots
    constexpr SimFields::Group_R_Crit200::value_type Rscale = 2.5;

    /*! @brief group and particle files are in multiple chunks */
    using chunk = CallbackUtils::chunk::Single<AF>;
    
    /*! @brief internal file layout is according to Illustris naming conventions
     *
     * @note the SIMBA conventions are the same
     */
    using name = CallbackUtils::name::Illustris<AF, PartType>;

    #ifdef FOR_ILLUSTRIS
    /*! @brief meta-data required by code are according to Illustris conventions */
    using meta = CallbackUtils::meta::Illustris<AF, PartType>;
    #else
    /*! @brief meta-data required by code are according to SIMBA conventions */
    using meta = CallbackUtils::meta::SIMBA<AF, PartType>;
    #endif

    #ifndef MASSES_AVAIL
    using masstab = CallbackUtils::meta_init::IllustrisMassTable<AF>;
    #endif // MASSES_AVAIL

    /*! @brief apply a minimum-mass selection on groups */
    using grp_select_M =  CallbackUtils::select::LowCutoff<AF, SimFields::Group_M_Crit200>;

    /*! @brief apply a minimum-radius selection on groups (exclude those with negative radius) */
    using grp_select_R = CallbackUtils::select::LowCutoff<AF, SimFields::Group_R_Crit200>;

    /*! @brief compute group radius as a multiple of R200c */
    using grp_radius =  CallbackUtils::radius::Simple<AF, SimFields::Group_R_Crit200>;

    /*! @brief store M200c in a vector */
    using grp_store_M = CallbackUtils::grp_action::StoreGrpProperty<AF, SimFields::Group_M_Crit200, grp_M_t>;

    /*! @brief store R200c in a vector */
    using grp_store_R = CallbackUtils::grp_action::StoreGrpProperty<AF, SimFields::Group_R_Crit200, grp_R_t>;

    /*! @brief store a #mass_prof::Profile instance for each group */
    using prt_compute_prof = CallbackUtils::prt_action::StorePrtHomogeneous<AF, grp_prof_t>;


} // namespace mass_prof }}}


/*! @brief implementation of the #Callback base for our application
 */
struct mass_prof_callback :
    virtual public Callback<mass_prof::AF>,
    public mass_prof::chunk, public mass_prof::name, public mass_prof::meta,
    #ifndef MASSES_AVAIL
    public mass_prof::masstab,
    #endif
    public mass_prof::grp_select_M, public mass_prof::grp_select_R,
    public mass_prof::grp_radius,
    public mass_prof::grp_store_M, public mass_prof::grp_store_R,
    public mass_prof::prt_compute_prof
{// {{{
    mass_prof_callback (char *fgrp, char *fprt) :
        mass_prof::chunk { fgrp, fprt },
        mass_prof::grp_select_M { Mmin },
        mass_prof::grp_select_R { Rmin },
        mass_prof::grp_radius { mass_prof::Rscale },
        mass_prof::grp_store_M { grp_M },
        mass_prof::grp_store_R { grp_R },
        mass_prof::prt_compute_prof { grp_prof }
    { }

    std::vector<mass_prof::grp_M_t> grp_M; /*!< @brief M200c data vector */
    std::vector<mass_prof::grp_R_t> grp_R; /*!< @brief R200c data vector */
    std::vector<mass_prof::grp_prof_t> grp_prof; /*!< @brief profiles data vector */

    #ifndef MASSES_AVAIL
    void normalize_mass ()
    {
        for (auto &prof_item : grp_prof)
            prof_item.normalize_mass(MassTable[mass_prof::PartType]);
    }
    #endif

private :
    using Callback<mass_prof::AF>::GrpProperties;
    using Callback<mass_prof::AF>::PrtProperties;

    // group mass cutoff -- choose 10^12.5
    static constexpr const SimFields::Group_M_Crit200::value_type Mmin = 316.23F;

    // group radius cutoff
    static constexpr const SimFields::Group_R_Crit200::value_type Rmin = 0.0F;
};// }}}

template<typename T>
void vec_to_f (const std::vector<T> &v, const char *s)
{
    auto f = std::fopen(s, "wb");
    std::fwrite(v.data(), sizeof(T), v.size(), f);
    std::fclose(f);
}

int main (int argc, char **argv)
{
    char *fgrp = argv[1];
    char *fprt = argv[2];
    char *dout = argv[3];

    mass_prof_callback C { fgrp, fprt };
    
    group_particles<> ( C );

    #ifndef MASSES_AVAIL
    C.normalize_mass();
    #endif // MASSES_AVAIL

    char fname_buffer[512];

    // save data to files
    std::sprintf(fname_buffer, "%s/grp_M200c.bin", dout);
    vec_to_f<>(C.grp_M, fname_buffer);

    std::sprintf(fname_buffer, "%s/grp_R200c.bin", dout);
    vec_to_f<>(C.grp_R, fname_buffer);

    #if defined GAS
    #   define TYPE "GAS"
    #elif defined DM
    #   define TYPE "DM"
    #elif defined STARS
    #   define TYPE "STARS"
    #elif defined BH
    #   define TYPE "BH"
    #endif

    {
        std::sprintf(fname_buffer, "%s/grp_mass_encl_prof_" TYPE ".bin", dout);
        auto f = std::fopen(fname_buffer, "wb");

        for (auto &prof : C.grp_prof)
            prof.save(f);

        std::fclose(f);
    }

    #undef TYPE

    return 0;
};
