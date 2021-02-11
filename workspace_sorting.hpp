#ifndef WORKSPACE_SORTING_HPP
#define WORKSPACE_SORTING_HPP

#include <vector>
#include <utility>
#include <cstdlib>
#include <algorithm>
#include <cstring>
#include <cassert>

template<typename GroupFields, typename ParticleFields>
class Workspace<GroupFields,ParticleFields>::Sorting
{
    typename ParticleFields::coord_t Bsize;
    const size_t Nprt;

    static constexpr const size_t Ncells_side = 16UL;
    static constexpr const size_t Ncells_tot  = Ncells_side * Ncells_side * Ncells_side;
    const typename ParticleFields::coord_t acell;
    
    // stores particle index in original particle order (first) and index in cells (second)
    std::vector<std::pair<size_t, size_t>> prt_indices;

    // stores where particles belonging to a specific cell are in the array
    // special value : Nprt+1 (no particles in cell)
    size_t offsets[Ncells_tot+1UL];

    // this is given by constructor, no memory allocation necessary
    void *tmp_prt_properties[ParticleFields::Nfields];

    // stuff that happens during construction
    void compute_prt_indices ();
    void sort_prt_indices ();
    void reorder_prt_properties ();
    void compute_offsets ();

public :
    Sorting (size_t Nprt_, typename ParticleFields::coord_t Bsize_, void **tmp_prt_properties_);
    ~Sorting ();

    // store the sorted properties here (instance must allocate memory for this!)
    // user can access these
    void *tmp_prt_properties_sorted[ParticleFields::Nfields];

    // user can use this function to find the indices of all particles that may
    // belong to a given group
    void prt_idx_ranges (typename GroupFields::coord_t grp_coord[3], float R,
                         std::vector<std::pair<size_t, size_t>> &out) const;
};

// ----- Implementation -----

template<typename GroupFields, typename ParticleFields>
Workspace<GroupFields,ParticleFields>::Sorting::Sorting (size_t Nprt_,
                                                         typename ParticleFields::coord_t Bsize_,
                                                         void **tmp_prt_properties_) :
    Nprt { Nprt_ }, Bsize { Bsize_ }, tmp_prt_properties { tmp_prt_properties_},
    acell { Bsize_ / Ncells_side }
{// {{{
    compute_prt_indices();
    sort_prt_indices();

    // allocate memory where we can store the sorted particle properties
    for (size_t ii=0; ii != ParticleFields::Nfields; ++ii)
        tmp_prt_properties_sorted[ii] = std::malloc(Nprt * ParticleFields::strides[ii]);

    reorder_prt_properties();
    compute_offsets();
}// }}}

template<typename GroupFields, typename ParticleFields>
Workspace<GroupFields,ParticleFields>::Sorting::~Sorting ()
{// {{{
    for (size_t ii=0; ii != ParticleFields::Nfields; ++ii)
        std::free(tmp_prt_properties_sorted[ii]);
}// }}}

template<typename GroupFields, typename ParticleFields>
Workspace<GroupFields,ParticleFields>::Sorting::compute_prt_indices ()
{// {{{
    prt_indices.reserve(Nprt);
    typename ParticleFields::coord_t *prt_coord = tmp_prt_properties[0];

    #define GRID(x, dir) ((size_t)(x[dir] / acell))

    for (size_t prt_idx=0; prt_idx != Nprt; ++prt_idx, prt_coord+=ParticleFields::dims[0])
        prt_indices.emplace_back(prt_idx, Ncells_side * Ncells_side * GRID(prt_coord, 0)
                                          +             Ncells_side * GRID(prt_coord, 1)
                                          +                           GRID(prt_coord, 2));

    #undef GRID
}// }}}

template<typename GroupFields, typename ParticleFields>
Workspace<GroupFields,ParticleFields>::Sorting::sort_prt_indices ()
{// {{{
    std::sort(prt_indices.begin(), prt_indices.end(),
              [](std::pair<size_t,size_t> a,
                 std::pair<size_t,size_t> b)
              { return a.second < b.second; } );
}// }}}

template<typename GroupFields, typename ParticleFields>
Workspace<GroupFields,ParticleFields>::Sorting::reorder_prt_properties ()
{// {{{
    assert(prt_indices.size() == Nprt);

    for (size_t ii=0; ii != ParticleFields::Nfields; ++ii)
    {
        char *dest = (char *)tmp_prt_properties_sorted[ii];

        for (size_t prt_idx=0; prt_idx != Nprt; ++prt_idx, dest += ParticleFields::strides[ii])
            std::memcpy(dest, (char *)(tmp_prt_properties[ii])
                              + prt_indices[prt_idx].first * ParticleFields::strides[ii],
                        ParticleFields::strides[ii]);
    }
}// }}}

template<typename GroupFields, typename ParticleFields>
Workspace<GroupFields,ParticleFields>::Sorting::compute_offsets ()
{

}

#endif // WORKSPACE_SORTING_HPP
