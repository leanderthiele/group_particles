/*! @file CAMELS_profiles.cpp
 *
 * @brief Example for using the group_particles code,
 *        computing the electron pressure, electron density, and temperature profiles
 *        for halos in the CAMELS simulations.
 *
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

/*! @brief assembles some types to use for the profile calculations.
 */
namespace CAMELS_prof
{// {{{
    constexpr const uint8_t PartType = 0; /*!< @brief gas particles */

    /*! @brief for the groups, we want to read their position and 200c masses and radii. */
    using GrpF = GrpFields<SimFields::GroupPos,
                           SimFields::Group_M_Crit200,
                           SimFields::Group_R_Crit200>;
    /*! @brief for the particles, we want to read their position and the quantities
     *         required to compute the electron pressure, electron density, and temperature.
     */
    using PrtF = PrtFields<SimFields::Coordinates,
                           SimFields::Masses,
                           SimFields::Density,
                           SimFields::InternalEnergy,
                           SimFields::ElectronAbundance>;

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

        std::vector<value_type> pressure;
        std::vector<value_type> electron_density;
        std::vector<value_type> temperature;
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
            pressure.resize(N, 0.0);
            electron_density.resize(N, 0.0);
            temperature.resize(N, 0.0);

            num_part.resize(N, 0UL);
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
            auto m = (value_type)prt.get<SimFields::Masses>();
            auto d = (value_type)prt.get<SimFields::Density>();
            auto e = (value_type)prt.get<SimFields::InternalEnergy>();
            auto x = (value_type)prt.get<SimFields::ElectronAbundance>();

            static constexpr const value_type gamma = 5.0/3.0, XH = 0.76;

            // this is electron pressure * particle volume
            auto P = 4.0 * x * XH / (1.0+3.0*XH+4.0*XH*x)
                         * (gamma-1.0) * m * e;

            // this is electron density * particle volume
            auto rho_e = x * XH * m;

            // this is temperature * particle volume
            auto T = (gamma-1.0) * e * 4.0 / (1.0+3.0*XH+4.0*XH*x) * m / d;

            pressure[idx] += P;
            electron_density[idx] += rho_e;
            temperature[idx] += T;
            ++num_part[idx];
        }// }}}

        /*! @brief append these profiles to files.
         *
         * @attention it is assumed that this instance is "dead" after this
         *            function is called!
         */
        void save (std::FILE *fpressure, std::FILE *felectron_density, std::FILE *ftemperature, std::FILE *fnum_part)
        {// {{{
            for (size_t ii=0; ii != N; ++ii)
            {
                pressure[ii] /= shell_vol(ii);
                electron_density[ii] /= shell_vol(ii);
                temperature[ii] /= shell_vol(ii);
            }

            std::fwrite(pressure.data(), sizeof(value_type), N, fpressure);
            std::fwrite(electron_density.data(), sizeof(value_type), N, felectron_density);
            std::fwrite(temperature.data(), sizeof(value_type), N, ftemperature);
            std::fwrite(num_part.data(), sizeof(size_t), N, fnum_part);
        }// }}}
    };// }}}
    
    using grp_M_t = double; /*!< @brief type to which M200c will be converted */
    using grp_R_t = double; /*!< @brief type to which R200c will be converted */
    using grp_P_t = double; /*!< @brief type to which P200c will be converted */
    using grp_prof_t = Profile; /*!< @brief type the pressure profiles will be stored in */

    /*! @brief stores P200c for each group */
    class grp_store_P :
        virtual public Callback<AF>,
        #ifdef FOR_ILLUSTRIS
        public CallbackUtils::meta_init::IllustrisCosmology<AF>,
        #else
        public CallbackUtils::meta_init::SIMBACosmology<AF>,
        #endif
        public CallbackUtils::grp_action::StoreGrpHomogeneous<AF, grp_P_t>
    {// {{{
        using Callback<AF>::GrpProperties;

        static constexpr grp_P_t GNewton    = 4.30091e4, // (kpc/1e10 Msun) (km/s)^2
                                 rho_crit_0 = 2.775e-8;  // 1e10 Msun/h / (kpc/h)^3

        /*! @brief helper function to compute critical density at this redshift
         *
         * By letting this class inherit from #CallbackUtils::meta_init::IllustrisCosmology
         * or #CallbackUtils::meta_init::SIMBACosmology,
         * we can use the cosmology information in the simulation headers in this function.
         */
        grp_P_t rho_crit () const
        {
            grp_P_t a3 = Time * Time * Time;
            return rho_crit_0 * (Omega0/a3 + OmegaLambda);
        }

        /*! @brief computes P200c from a group's properties
         *
         * This is the required override of the pure virtual in
         * #CallbackUtils::grp_action::StoreGrpHomogeneous.
         */
        grp_P_t grp_reduce (const GrpProperties &grp) const override
        {
            auto m = (grp_P_t)grp.get<SimFields::Group_M_Crit200>();
            auto r = (grp_P_t)grp.get<SimFields::Group_R_Crit200>();
            return 100.0 * GNewton * m * rho_crit() * OmegaBaryon / Omega0 / r;
        }
    public :
        /*! @param data     a zero-length vector into which the P200c values will be written
         */
        grp_store_P (std::vector<grp_P_t> &data) :
            CallbackUtils::grp_action::StoreGrpHomogeneous<AF, grp_P_t> { data }
        { }
    };// }}}

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

    /*! @brief store a #CAMELS_prof::Profile instance for each group */
    using prt_compute_prof = CallbackUtils::prt_action::StorePrtHomogeneous<AF, grp_prof_t>;


} // namespace CAMELS_prof }}}


/*! @brief implementation of the #Callback base for our application
 */
struct CAMELS_prof_callback :
    virtual public Callback<CAMELS_prof::AF>,
    public CAMELS_prof::chunk, public CAMELS_prof::name, public CAMELS_prof::meta,
    public CAMELS_prof::grp_select_M, public CAMELS_prof::grp_select_R,
    public CAMELS_prof::grp_radius,
    public CAMELS_prof::grp_store_M, public CAMELS_prof::grp_store_R, public CAMELS_prof::grp_store_P,
    public CAMELS_prof::prt_compute_prof
{// {{{
    CAMELS_prof_callback (char *fgrp, char *fprt) :
        CAMELS_prof::chunk { fgrp, fprt },
        CAMELS_prof::grp_select_M { Mmin },
        CAMELS_prof::grp_select_R { Rmin },
        CAMELS_prof::grp_radius { CAMELS_prof::Rscale },
        CAMELS_prof::grp_store_M { grp_M },
        CAMELS_prof::grp_store_R { grp_R },
        CAMELS_prof::grp_store_P { grp_P },
        CAMELS_prof::prt_compute_prof { grp_prof }
    { }

    std::vector<CAMELS_prof::grp_M_t> grp_M; /*!< @brief M200c data vector */
    std::vector<CAMELS_prof::grp_R_t> grp_R; /*!< @brief R200c data vector */
    std::vector<CAMELS_prof::grp_P_t> grp_P; /*!< @brief P200c data vector */
    std::vector<CAMELS_prof::grp_prof_t> grp_prof; /*!< @brief profiles data vector */

private :
    using Callback<CAMELS_prof::AF>::GrpProperties;
    using Callback<CAMELS_prof::AF>::PrtProperties;

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

    CAMELS_prof_callback C { fgrp, fprt };
    
    group_particles<> ( C );

    char fname_buffer[512];

    // save data to files
    std::sprintf(fname_buffer, "%s/grp_M200c.bin", dout);
    vec_to_f<>(C.grp_M, fname_buffer);

    std::sprintf(fname_buffer, "%s/grp_R200c.bin", dout);
    vec_to_f<>(C.grp_R, fname_buffer);

    std::sprintf(fname_buffer, "%s/grp_P200c.bin", dout);
    vec_to_f<>(C.grp_P, fname_buffer);

    {
        std::sprintf(fname_buffer, "%s/grp_pressure_prof.bin", dout);
        auto fpressure = std::fopen(fname_buffer, "wb");

        std::sprintf(fname_buffer, "%s/grp_electron_density_prof.bin", dout);
        auto felectron_density = std::fopen(fname_buffer, "wb");

        std::sprintf(fname_buffer, "%s/grp_temperature_prof.bin", dout);
        auto ftemperature = std::fopen(fname_buffer, "wb");

        std::sprintf(fname_buffer, "%s/grp_num_part_prof.bin", dout);
        auto fnum_part = std::fopen(fname_buffer, "wb");

        for (auto &prof : C.grp_prof)
            prof.save(fpressure, felectron_density, ftemperature, fnum_part);

        std::fclose(fpressure);
        std::fclose(felectron_density);
        std::fclose(ftemperature);
        std::fclose(fnum_part);
    }

    return 0;
};
