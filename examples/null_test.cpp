// This code is a useful test whether the group_particles stuff is working correctly.
#include "group_particles.hpp"
#include "common_fields.hpp"

namespace NullTest
{
    constexpr const size_t PartType = 1;

    typedef GrpFields<IllustrisFields::GroupPos,
                      IllustrisFields::Group_M_Crit200,
                      IllustrisFields::Group_R_Crit200> GrpF;
    typedef PrtFields<IllustrisFields::Coordinates> PrtF;
    typedef AllFields<GrpF, PrtF> AF;

    typedef float grp_M_t;
    typedef size_t grp_N_t;

    typedef CallbackUtils::chunk::Multi<AF>
        chunk;
    typedef CallbackUtils::name::Illustris<AF, PartType>
        name;
    typedef CallbackUtils::meta::Illustris<AF, PartType>
        meta;
    typedef CallbackUtils::select::LowCutoff<AF, IllustrisFields::Group_M_Crit200>
        grp_select;
    typedef CallbackUtils::radius::Simple<AF, IllustrisFields::Group_R_Crit200>
        grp_radius;
    typedef CallbackUtils::grp_action::StoreGrpProperty<AF, IllustrisFields::Group_M_Crit200, grp_M_t>
        grp_store_M;
    typedef CallbackUtils::prt_action::StorePrtHomogeneous<AF, grp_N_t>
        prt_count_prt;
} // namespace NullTest

struct NullTest_Callback :
    virtual public Callback<NullTest::AF>,
    public NullTest::chunk, public NullTest::name, public NullTest::meta,
    public NullTest::grp_select, public NullTest::grp_radius,
    public NullTest::grp_store_M, public NullTest::prt_count_prt
{
    NullTest_Callback () :
        NullTest::chunk { fgrp, grp_max_idx, fprt, prt_max_idx },
        NullTest::grp_select { Mmin },
        NullTest::grp_store_M { grp_M },
        NullTest::prt_count_prt { grp_N }
    { }

    std::vector<NullTest::grp_M_t> grp_M;
    std::vector<NullTest::grp_N_t> grp_N;

private :
    static constexpr const float Mmin = 1e3F;

    // files
    #define ROOT "/tigress/lthiele/Illustris_300-1_Dark/output/"
    static constexpr const char fgrp[]        = ROOT"groups_099/fof_subhalo_tab_099.%d.hdf5";
    static constexpr const size_t grp_max_idx = 74;
    static constexpr const char fprt[]        = ROOT"snapdir_099/snap_099.%d.hdf5";
    static constexpr const size_t prt_max_idx = 74;
    #undef ROOT

    // the (fairly trivial) functionality of this class : count particles
    void prt_insert (size_t grp_idx, const GrpProperties &grp,
                     const PrtProperties &prt, float R, NullTest::grp_N_t &data_item) override
    {
        ++data_item;
    }

};

template<typename T>
void vec_to_f (const std::vector<T> &v, const char *s)
{
    std::FILE *f = std::fopen(s, "wb");
    std::fwrite(v.data(), sizeof(T), v.size(), f);
    std::fclose(f);
}

int main ()
{
    NullTest_Callback n;
    
    group_particles<> ( n );

    // save data to files
    vec_to_f<>(n.grp_M, "./null_test_result_Feb13/grp_M.bin");
    vec_to_f<>(n.grp_N, "./null_test_result_Feb13/grp_N.bin");

    return 0;
};
