/*! @file y_prof.cpp
 *
 * @brief Example for using the group_particles code,
 *        computing the integrated electron pressure (Compton-y) profiles of halos in an Illustris TNG simulation.
 *
 * This code is for debugging purposes only. It will work on only one group and has the optimization that only
 * particle chunks that have particles falling within R200c of this group will be loaded.
 * Thus, in particular, it is not possible to use this code for larger radii or other groups.
 */

#include <memory>
#include <string>
#include <vector>
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

        // this should be equal to the result from Y_Delta_debug
        // -- but it isn't!
        value_type total_Y;

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
            logRmax = std::log(1.00F * grp.get<IllustrisFields::Group_R_Crit200>());
            dlogR = (logRmax - logRmin) / (coord_t)(N-1);
            total_Y = 0.0;
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
            // load the required properties of this particle
            auto m = (value_type)prt.get<IllustrisFields::Masses>();
            auto e = (value_type)prt.get<IllustrisFields::InternalEnergy>();
            auto x = (value_type)prt.get<IllustrisFields::ElectronAbundance>();

            static constexpr const value_type gamma = 5.0/3.0, XH = 0.76;

            // this is electron pressure * particle volume
            auto Y = 2.0 * (1.0+XH) / (1.0+3.0*XH+4.0*XH*x)
                         * (gamma-1.0) * m * e;

            total_Y += Y;
            
            coord_t logR = 0.5 * std::log(Rsq);
            if (logR > logRmax)
                return;

            size_t idx;

            // figure out the index of the spherical shell the particle falls into
            if (logR < logRmin)
                idx = 0UL;
            else
                idx = 1UL + (size_t)((logR - logRmin) / dlogR);
            
            // safety check against numerical issues
            if (idx >= N)
                return;

            for (size_t ii=idx; ii < N; ++ii)
            {
                pressure[ii] += Y;
                ++num_part[ii];
            }
        }// }}}

        /*! @brief append this electron pressure profile to file.
         *
         * @attention it is assumed that this instance is "dead" after this
         *            function is called!
         */
        void save (std::FILE *fpressure, std::FILE *fnum_part)
        {// {{{
            std::fwrite(pressure.data(), sizeof(value_type), N, fpressure);
            std::fwrite(num_part.data(), sizeof(size_t), N, fnum_part);
            std::fprintf(stderr, "total_Y = %.8e\n", total_Y);
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
    constexpr IllustrisFields::Group_R_Crit200::value_type Rscale = 1.0F;

    /*! @brief group and particle files are in multiple chunks */
    using grp_chunk = CallbackUtils::chunk::MultiGrp<AF>;

    /*! @brief internal file layout is according to Illustris naming conventions */
    using name = CallbackUtils::name::Illustris<AF, PartType>;

    /*! @brief meta-data required by code are according to Illustris conventions */
    using meta = CallbackUtils::meta::Illustris<AF, PartType>;

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
    public y_prof::grp_chunk, public y_prof::name, public y_prof::meta,
    public y_prof::grp_radius,
    public y_prof::grp_store_M, public y_prof::grp_store_R, public y_prof::grp_store_P,
    public y_prof::prt_compute_Y
{// {{{
    y_prof_callback () :
        y_prof::grp_chunk { fgrp, grp_max_idx },
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

    // for this specific group!
    bool prt_chunk (size_t chunk_idx, std::string &fname) const override final
    {
        static std::vector<int> prt_chunk_indices = { 34, 151, 153, 154, 155, 156, 157, 158,
                                                     159, 160, 161, 162, 163, 164, 295, 296,
                                                     297, 298, 299, 300, 301, 302, 303, 304,
                                                     305, 306, 307, 309, 310, 311, 312, 313,
                                                     314, 315, 318, 320, 322, 323, 324, 325,
                                                     326, 327, 496, 497, 498, 499, 500, 501,
                                                     502, 503, 504, 505, 506, 507, 508, 509,
                                                     510, 511, 512, 513, 514, 515, 516, 517,
                                                     518, 519, 520, 512, 522, 523, 524, 525,
                                                     526, 527, 528, 529, 530, 531, 532, 533,
                                                     534, 535, 536, 537, 538, 539, 540, 541,
                                                     542, 543, 544, 545, 546, 547, 548, 549,
                                                     550, 551, 552, 553, 554, 555, 556, 557,
                                                     558, 559, 560, 561, 562, 563, 564, 565,
                                                     566, 567, 568, 569, 570, 571, 572, 573,
                                                     574, 575, 576, 577, 578, 579, 580, 581,
                                                     582, 583, 584, 585, 586, 587, 588, 589,
                                                     590, 591, 592, 593, 594, 595, 596, 597,
                                                     598, 599 };

        static char buf[512];

        if (chunk_idx >= prt_chunk_indices.size())
            return false;

        std::sprintf(buf, fprt, prt_chunk_indices[chunk_idx]);
        fname = std::string(buf);
        return true;
    }

    // select only one specific item
    bool grp_select (const GrpProperties &grp) const override final
    {
        const auto m = grp.get<IllustrisFields::Group_M_Crit200>();

        return (m > 4600.0) && (m < 4601.0);
    }
    
    // files
    #define ROOT "/tigress/lthiele/Illustris_300-1_TNG/output/"
    static constexpr const char fgrp[]        = ROOT"groups_099/fof_subhalo_tab_099.%d.hdf5";
    static constexpr const size_t grp_max_idx = 150;
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
    #define ROOT "yintegrated_prof_debug_Apr11"
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
