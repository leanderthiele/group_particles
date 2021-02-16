#ifndef WORKSPACE_SORTING_HPP
#define WORKSPACE_SORTING_HPP

#include <vector>
#include <utility>
#include <cstdlib>
#include <algorithm>
#include <cstring>
#include <cassert>
#include <cmath>

#include "workspace.hpp"
#include "timing.hpp"

// TODO
// put parallel execution policy into std::sort
// (depending on preprocessor flag)

template<typename AFields>
class Workspace<AFields>::Sorting
{// {{{
    typename AFields::ParticleFields::coord_t Bsize;
    const size_t Nprt;

    static constexpr const size_t Ncells_side = 16UL;
    static constexpr const size_t Ncells_tot  = Ncells_side * Ncells_side * Ncells_side;
    const typename AFields::ParticleFields::coord_t acell;
    
    // stores particle index in original particle order (first) and index in cells (second)
    std::vector<std::pair<size_t, size_t>> prt_indices;

    // stores where particles belonging to a specific cell are in the array
    // special value : Nprt (no particles in cell)
    size_t offsets[Ncells_tot+1UL];

    // this is given by constructor, no memory allocation necessary
    void **tmp_prt_properties;

    // stuff that happens during construction
    void compute_prt_indices ();
    void sort_prt_indices ();
    void reorder_prt_properties ();
    void compute_offsets ();

    class Geometry
    {
        // these two helper functions bring the geometric setup into a canonical form
        static void mod_translations (const typename AFields::GroupFields::coord_t grp_coord[AFields::GroupFields::dims[0]],
                                      typename AFields::GroupFields::coord_t cub_coord[AFields::GroupFields::dims[0]],
                                      typename AFields::GroupFields::coord_t periodicity);
        static void mod_reflections (typename AFields::GroupFields::coord_t cub_coord[AFields::GroupFields::dims[0]]);
    public :
        // assumes that the parameters passed have units such that the cube sidelength is unity
        // The cub_coord array will be modified!
        static bool sph_cub_intersect (const typename AFields::GroupFields::coord_t grp_coord[AFields::GroupFields::dims[0]],
                                       typename AFields::GroupFields::coord_t cub_coord[AFields::GroupFields::dims[0]],
                                       typename AFields::GroupFields::coord_t grp_R,
                                       typename AFields::GroupFields::coord_t periodicity);
        Geometry () = delete;
    };

public :
    Sorting (size_t Nprt_, typename AFields::ParticleFields::coord_t Bsize_, void **tmp_prt_properties_);
    Sorting () = delete;
    ~Sorting ();

    // store the sorted properties here (instance must allocate memory for this!)
    // user can access these
    void *tmp_prt_properties_sorted[AFields::ParticleFields::Nfields];

    // user can use this function to find the indices of all particles that may
    // belong to a given group
    std::vector<std::pair<size_t, size_t>> prt_idx_ranges
        (const typename AFields::GroupFields::coord_t grp_coord[AFields::GroupFields::dims[0]],
         const typename AFields::GroupFields::coord_t R) const;
};// }}}

// ----- Implementation -----

template<typename AFields>
Workspace<AFields>::Sorting::Sorting (size_t Nprt_,
                                      typename AFields::ParticleFields::coord_t Bsize_,
                                      void **tmp_prt_properties_) :
    Nprt { Nprt_ }, Bsize { Bsize_ }, tmp_prt_properties { tmp_prt_properties_ },
    acell { Bsize_ / Ncells_side },
    prt_indices { }
{// {{{

    #ifndef NDEBUG
    TIME_PT(t1);
    #endif // NDEBUG
    compute_prt_indices();
    #ifndef NDEBUG
    TIME_MSG(t1,"Sorting::compute_prt_indices");
    #endif // NDEBUG

    #ifndef NDEBUG
    TIME_PT(t2);
    #endif // NDEBUG
    sort_prt_indices();
    #ifndef NDEBUG
    TIME_MSG(t2, "Sorting::sort_prt_indices");
    #endif // NDEBUG

    #ifndef NDEBUG
    TIME_PT(t3);
    #endif // NDEBUG
    // allocate memory where we can store the sorted particle properties
    for (size_t ii=0; ii != AFields::ParticleFields::Nfields; ++ii)
        tmp_prt_properties_sorted[ii] = std::malloc(Nprt * AFields::ParticleFields::strides[ii]);
    #ifndef NDEBUG
    TIME_MSG(t3, "Sorting memory allocation");
    #endif // NDEBUG

    #ifndef NDEBUG
    TIME_PT(t4);
    #endif // NDEBUG
    reorder_prt_properties();
    #ifndef NDEBUG
    TIME_MSG(t4, "Sorting::reorder_prt_properties");
    #endif // NDEBUG

    #ifndef NDEBUG
    TIME_PT(t5);
    #endif // NDEBUG
    compute_offsets();
    #ifndef NDEBUG
    TIME_MSG(t5, "Sorting::compute_offsets");
    #endif
}// }}}

template<typename AFields>
Workspace<AFields>::Sorting::~Sorting ()
{// {{{
    for (size_t ii=0; ii != AFields::ParticleFields::Nfields; ++ii)
        std::free(tmp_prt_properties_sorted[ii]);
}// }}}

template<typename AFields>
void
Workspace<AFields>::Sorting::compute_prt_indices ()
{// {{{
    prt_indices.reserve(Nprt);
    auto *prt_coord = (typename AFields::ParticleFields::coord_t *)tmp_prt_properties[0];

    #define GRID(x, dir) ((size_t)(x[dir] / acell))

    for (size_t prt_idx=0; prt_idx != Nprt;
         ++prt_idx, prt_coord += AFields::ParticleFields::dims[0])
        prt_indices.emplace_back(prt_idx, Ncells_side * Ncells_side * GRID(prt_coord, 0)
                                          +             Ncells_side * GRID(prt_coord, 1)
                                          +                           GRID(prt_coord, 2));

    #undef GRID
}// }}}

template<typename AFields>
void
Workspace<AFields>::Sorting::sort_prt_indices ()
{// {{{
    std::sort(prt_indices.begin(), prt_indices.end(),
              [](std::pair<size_t,size_t> a,
                 std::pair<size_t,size_t> b)
              { return a.second < b.second; } );
}// }}}

template<typename AFields>
void
Workspace<AFields>::Sorting::reorder_prt_properties ()
{// {{{
    assert(prt_indices.size() == Nprt);

    for (size_t ii=0; ii != AFields::ParticleFields::Nfields; ++ii)
    {
        char *dest = (char *)tmp_prt_properties_sorted[ii];

        for (size_t prt_idx=0; prt_idx != Nprt; ++prt_idx, dest += AFields::ParticleFields::strides[ii])
            std::memcpy(dest, (char *)(tmp_prt_properties[ii])
                              + prt_indices[prt_idx].first * AFields::ParticleFields::strides[ii],
                        AFields::ParticleFields::strides[ii]);
    }
}// }}}

template<typename AFields>
void
Workspace<AFields>::Sorting::compute_offsets ()
{// {{{
    // initialize all offsets to default value
    for (size_t ii=0; ii != Ncells_tot+1UL; ++ii)
        offsets[ii] = Nprt;

    offsets[prt_indices[0].second] = 0UL;

    for (size_t jj=1UL; jj != Nprt; ++jj)
        if (prt_indices[jj].second != prt_indices[jj-1UL].second)
            // we have found the beginning of a new cell
            offsets[prt_indices[jj].second] = jj;
}// }}}

template<typename AFields>
std::vector<std::pair<size_t, size_t>>
Workspace<AFields>::Sorting::prt_idx_ranges
    (const typename AFields::GroupFields::coord_t grp_coord[AFields::GroupFields::dims[0]],
     typename AFields::GroupFields::coord_t R) const
{// {{{
    std::vector<std::pair<size_t, size_t>> out;

    typename AFields::GroupFields::coord_t grp_coord_normalized[AFields::GroupFields::dims[0]];
    for (size_t ii=0; ii != AFields::GroupFields::dims[0]; ++ii)
        grp_coord_normalized[ii] = grp_coord[ii] / acell;

    const typename AFields::GroupFields::coord_t R_normalized = R / acell;

    // loop over all cells -- there's a more efficient way of doing this,
    //                        by some preselection based on R,
    //                        but this is fast enough
    for (size_t ii=0; ii != Ncells_tot; ++ii)
    {
        typename AFields::GroupFields::coord_t cub_coord[]
            = { ii / (Ncells_side*Ncells_side),
                (ii/Ncells_side) % Ncells_side,
                ii % Ncells_side };

        if (// check if group has overlap with this cell
            Geometry::sph_cub_intersect(grp_coord_normalized, cub_coord,
                                        R_normalized,
                                        (typename AFields::GroupFields::coord_t)Bsize)
            // check if this cell has any particles
            && offsets[ii] < Nprt)
        {
            std::pair<size_t, size_t> idx_range { offsets[ii], Nprt };

            // find the last one -- taking possibly empty cells that follow into account!
            for (size_t jj=ii+1UL; jj != Ncells_tot; ++jj)
                if (offsets[jj] < Nprt)
                {
                    idx_range.second = offsets[jj];
                    break;
                }
            
            out.push_back(idx_range);
        }
    }

    return out;
}// }}}

template<typename AFields>
inline bool
Workspace<AFields>::Sorting::Geometry::sph_cub_intersect
    (const typename AFields::GroupFields::coord_t grp_coord[AFields::GroupFields::dims[0]],
     typename AFields::GroupFields::coord_t cub_coord[AFields::GroupFields::dims[0]],
     typename AFields::GroupFields::coord_t grp_R,
     typename AFields::GroupFields::coord_t periodicity)
{// {{{
    mod_translations(grp_coord, cub_coord, periodicity);
    mod_reflections(cub_coord);

    return std::hypot(std::max((typename AFields::GroupFields::coord_t)0.0, cub_coord[0]),
                      std::max((typename AFields::GroupFields::coord_t)0.0, cub_coord[1]),
                      std::max((typename AFields::GroupFields::coord_t)0.0, cub_coord[2]) ) < grp_R;
}// }}}

template<typename AFields>
inline void
Workspace<AFields>::Sorting::Geometry::mod_translations
    (const typename AFields::GroupFields::coord_t grp_coord[AFields::GroupFields::dims[0]],
     typename AFields::GroupFields::coord_t cub_coord[AFields::GroupFields::dims[0]],
     typename AFields::GroupFields::coord_t periodicity)
{// {{{
    for (size_t ii=0; ii != AFields::GroupFields::dims[0]; ++ii)
    {
        cub_coord[ii] -= grp_coord[ii];
        if (cub_coord[ii] > 0.5 * periodicity)
            cub_coord[ii] -= periodicity;
        else if (cub_coord[ii] < -0.5 * periodicity)
            cub_coord[ii] += periodicity;
    }
}// }}}

template<typename AFields>
inline void
Workspace<AFields>::Sorting::Geometry::mod_reflections
    (typename AFields::GroupFields::coord_t cub_coord[AFields::GroupFields::dims[0]])
{// {{{
    for (size_t ii=0; ii != AFields::GroupFields::dims[0]; ++ii)
        if (cub_coord[ii] < -0.5)
            cub_coord[ii] = - ( cub_coord[ii] + 1.0 );
}// }}}

#endif // WORKSPACE_SORTING_HPP
