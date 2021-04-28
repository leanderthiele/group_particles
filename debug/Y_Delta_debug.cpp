#include <memory>
#include <string>
#include <cstdio>

#include "group_particles.hpp"
#include "common_fields.hpp"

namespace Y_Delta
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

    typedef double grp_Y_t;
    typedef double grp_M_t;
    typedef double grp_R_t;

    typedef CallbackUtils::chunk::MultiGrp<AF>
        grp_chunk;
    typedef CallbackUtils::name::Illustris<AF, PartType>
        name;
    typedef CallbackUtils::meta::Illustris<AF, PartType>
        meta;
    typedef CallbackUtils::radius::Simple<AF, IllustrisFields::Group_R_Crit200>
        grp_radius;
    typedef CallbackUtils::grp_action::StoreGrpProperty<AF, IllustrisFields::Group_M_Crit200, grp_M_t>
        grp_store_M;
    typedef CallbackUtils::grp_action::StoreGrpProperty<AF, IllustrisFields::Group_R_Crit200, grp_R_t>
        grp_store_R;
    typedef CallbackUtils::prt_action::StorePrtHomogeneous<AF, grp_Y_t>
        prt_compute_Y;
} // namespace Y_Delta

struct Y_Delta_Callback :
    virtual public Callback<Y_Delta::AF>,
    public Y_Delta::grp_chunk, public Y_Delta::name, public Y_Delta::meta,
    public Y_Delta::grp_radius,
    public Y_Delta::grp_store_M,
    public Y_Delta::grp_store_R,
    public Y_Delta::prt_compute_Y
{// {{{
    Y_Delta_Callback () :
        Y_Delta::grp_chunk { fgrp, grp_max_idx },
        Y_Delta::grp_store_M { grp_M },
        Y_Delta::grp_store_R { grp_R },
        Y_Delta::prt_compute_Y { grp_Y }
    { }

    // data (public so user can do something with them once they are assembled)
    std::vector<Y_Delta::grp_M_t> grp_M;
    std::vector<Y_Delta::grp_R_t> grp_R;
    std::vector<Y_Delta::grp_Y_t> grp_Y;

private :
    // for calculation of electron pressure
    static constexpr const float gamma = 5.0F/3.0F;
    static constexpr const float XH    = 0.76F;

    // group mass cutoff
    static constexpr const float Mmin = 1e3F;

    // files
    #define ROOT "/tigress/lthiele/Illustris_300-1_TNG/output/"
    static constexpr const char fgrp[]        = ROOT"groups_099/fof_subhalo_tab_099.%d.hdf5";
    static constexpr const size_t grp_max_idx = 150;
    static constexpr const char fprt[]        = ROOT"snapdir_099/snap_099.%d.hdf5";
    static constexpr const size_t prt_max_idx = 599;
    #undef ROOT

    // select only one specific item
    bool grp_select (const GrpProperties &grp) const override final
    {
        const auto m = grp.get<IllustrisFields::Group_M_Crit200>();

        return (m > 4600.0) && (m < 4601.0);
    }

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

    // computation of Compton-Y for a single gas particle
    void prt_insert (size_t grp_idx, const GrpProperties &grp, const PrtProperties &prt,
                     float Rsq, Y_Delta::grp_Y_t &data_item) override
    {
        auto m = prt.get<IllustrisFields::Masses>();
        auto e = prt.get<IllustrisFields::InternalEnergy>();
        auto x = prt.get<IllustrisFields::ElectronAbundance>();

        data_item += 2.0F * (1.0F+XH) / (1.0F+3.0F*XH+4.0F*XH*x)
                     * (gamma-1.0F) * m * e;
    }
};// }}}

template<typename T>
void vec_to_f (const std::vector<T> &v, const char *s)
{
    std::FILE *f = std::fopen(s, "wb");
    std::fwrite(v.data(), sizeof(T), v.size(), f);
    std::fclose(f);
}

int main ()
{
    Y_Delta_Callback y;
    
    group_particles<> ( y );

    // save data to files
    #define ROOT "Y_Delta_debug_Apr11"
    vec_to_f<>(y.grp_M, ROOOT"/grp_M.bin");
    vec_to_f<>(y.grp_Y, ROOT"/grp_Y.bin");
    #undef ROOT

    return 0;
};
