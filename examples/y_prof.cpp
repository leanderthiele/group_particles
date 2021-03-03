/*! @file y_prof.cpp
 *
 * @brief Example for using the group_particles code,
 *        computing the electron pressure (Compton-y) profiles of halos in an Illustris TNG simulation.
 */

#include <memory>
#include <string>
#include <cstdio>
#include <cmath>

#include "group_particles.hpp"
#include "common_fields.hpp"

/*! @brief assembles some types to use for the pressure profile calculation.
 */
namespace y_prof
{// {{{
    constexpr const uint8_t PartType = 0; /*!< @brief gas particles */

    /*! @brief for the groups, we want to read their position and 200c masses and radii. */
    using GrpF = GrpFields<IllustrisFields::GroupCM,
                           IllustrisFields::Group_M_Crit200,
                           IllustrisFields::Group_R_Crit200>;

    /*! @brief for the particles, we want to read their position and the quanitities
     *         required to compute the electron pressure.
     */
    using PrtF = PrtFields<IllustrisFields::Coordinates,
                           IllustrisFields::Masses,
                           IllustrisFields::InternalEnergy,
                           IllustrisFields::ElectronAbundance>;

    using AF = AllFields<GrpF, PrtF>; /*!< @brief bundle the group and particle fields */

    /*! @brief stores a single electron pressure profile */
    class YProfile
    {// {{{
        using GrpProperties = typename Callback<AF>::GrpProperties;
        using PrtProperties = typename Callback<AF>::PrtProperties;

        // all internal calculations regarding Y in this type
        using value_type = double;

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

        /*! @brief constructs from group properties (we need to read the radius)
         *
         * As explained in the documentation for #CallbackUtils::prt_action::StorePrtHomogeneous,
         * if the data type stored for each group has a constructor from a
         * `const GrpProperties &`, this constructor will be called during
         * the calls to #Callback::grp_action.
         * We use this fact to read the group radius.
         */
        YProfile (const GrpProperties &grp)
        {// {{{
            pressure.resize(N, 0.0);
            num_part.resize(N, 0UL);
            logRmin = std::log(0.03F * grp.get<IllustrisFields::Group_R_Crit200>());
            logRmax = std::log(2.50F * grp.get<IllustrisFields::Group_R_Crit200>());
            dlogR = (logRmax - logRmin) / (coord_t)N;
        }// }}}

        /*! @brief adds a particle to the profile
         *
         * This function computes the electron pressure associated with a particle
         * and adds it (weighted with its volume) to the appropriate profile bin.
         *
         * Note that this exact signature and function name is required,
         * as explained in the documentation for #CallbackUtils::prt_action::StorePrtHomogeneous.
         */
        void prt_insert (size_t, const GrpProperties &,
                         const PrtProperties &prt, coord_t Rsq)
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

            static constexpr const value_type gamma = 5.0/3.0, XH = 0.76;

            // this is electron pressure * particle volume
            auto Y = 2.0 * (1.0+XH) / (1.0+3.0*XH+4.0*XH*x)
                         * (gamma-1.0) * m * e;

            pressure[idx] += Y;
            ++num_part[idx];
        }// }}}

        /*! @brief append this electron pressure profile to file.
         *
         * @attention it is assumed that this instance is "dead" after this
         *            function is called!
         */
        void save (std::FILE *fpressure, std::FILE *fnum_part)
        {// {{{
            for (size_t ii=0; ii != N; ++ii)
                pressure[ii] /= shell_vol(ii);

            std::fwrite(pressure.data(), sizeof(value_type), N, fpressure);
            std::fwrite(num_part.data(), sizeof(size_t), N, fnum_part);
        }// }}}
    };// }}}
    
    using grp_M_t = double; /*!< @brief type to which M200c will be converted */
    using grp_R_t = double; /*!< @brief type to which R200c will be converted */
    using grp_P_t = double; /*!< @brief type to which P200c will be converted */
    using grp_Y_t = YProfile; /*!< @brief type the pressure profiles will be stored in */

    /*! @brief stores P200c for each group */
    class grp_store_P :
        virtual public Callback<AF>,
        public CallbackUtils::meta_init::IllustrisCosmology<AF>,
        public CallbackUtils::grp_action::StoreGrpHomogeneous<AF, grp_P_t>
    {// {{{
        using Callback<AF>::GrpProperties;

        static constexpr grp_P_t GNewton    = 4.30091e4, // (kpc/1e10 Msun) (km/s)^2
                                 rho_crit_0 = 2.775e-8;  // 1e10 Msun/h / (kpc/h)^3

        /*! @brief helper function to compute critical density at this redshift
         *
         * By letting this class inherit from #CallbackUtils::meta_init::IllustrisCosmology,
         * we can use the cosmology information in the Illustris header in this function.
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
            auto m = (grp_P_t)grp.get<IllustrisFields::Group_M_Crit200>();
            auto r = (grp_P_t)grp.get<IllustrisFields::Group_R_Crit200>();
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

    /*! @brief store a #y_prof::YProfile instance for each group */
    using prt_compute_Y = CallbackUtils::prt_action::StorePrtHomogeneous<AF, grp_Y_t>;


} // namespace y_prof }}}


/*! @brief implementation of the #Callback base for our application
 */
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

    std::vector<y_prof::grp_M_t> grp_M; /*!< @brief M200c data vector */
    std::vector<y_prof::grp_R_t> grp_R; /*!< @brief R200c data vector */
    std::vector<y_prof::grp_P_t> grp_P; /*!< @brief P200c data vector */
    std::vector<y_prof::grp_Y_t> grp_Y; /*!< @brief pressure profile data vector */

private :
    using Callback<y_prof::AF>::GrpProperties;
    using Callback<y_prof::AF>::PrtProperties;

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
