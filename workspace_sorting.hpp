#ifndef WORKSPACE_SORTING_HPP
#define WORKSPACE_SORTING_HPP

#include <vector>
#include <utility>
#include <cstdlib>
#include <algorithm>
#include <cstring>
#include <cassert>
#include <cmath>

template<typename GroupFields, typename ParticleFields>
class Workspace<GroupFields,ParticleFields>::Sorting
{// {{{
    typename ParticleFields::coord_t Bsize;
    const size_t Nprt;

    static constexpr const size_t Ncells_side = 16UL;
    static constexpr const size_t Ncells_tot  = Ncells_side * Ncells_side * Ncells_side;
    const typename ParticleFields::coord_t acell;
    
    // stores particle index in original particle order (first) and index in cells (second)
    std::vector<std::pair<size_t, size_t>> prt_indices;

    // stores where particles belonging to a specific cell are in the array
    // special value : Nprt (no particles in cell)
    size_t offsets[Ncells_tot+1UL];

    // this is given by constructor, no memory allocation necessary
    void *tmp_prt_properties[ParticleFields::Nfields];

    // stuff that happens during construction
    void compute_prt_indices ();
    void sort_prt_indices ();
    void reorder_prt_properties ();
    void compute_offsets ();

    class Geometry
    {
        // these two helper functions bring the geometric setup into a canonical form
        static void mod_translations (const typename GroupFields::coord_t grp_coord[GroupFields::dims[0]],
                                      typename GroupFields::coord_t cub_coord[GroupFields::dims[0]],
                                      typename GroupFields::coord_t periodicity);
        static void mod_reflections (typename GroupFields::coord_t cub_coord[GroupFields::dims[0]]);
    public :
        // assumes that the parameters passed have units such that the cube sidelength is unity
        // The cub_coord array will be modified!
        static bool sph_cub_intersect (const typename GroupFields::coord_t grp_coord[GroupFields::dims[0]],
                                       typename GroupFields::coord_t cub_coord[GroupFields::dims[0]],
                                       typename GroupFields::coord_t grp_R,
                                       typename GroupFields::coord_t periodicity);
        Geometry () = delete;
    };

public :
    Sorting (size_t Nprt_, typename ParticleFields::coord_t Bsize_, void **tmp_prt_properties_);
    Sorting () = delete;
    ~Sorting ();

    // store the sorted properties here (instance must allocate memory for this!)
    // user can access these
    void *tmp_prt_properties_sorted[ParticleFields::Nfields];

    // user can use this function to find the indices of all particles that may
    // belong to a given group
    std::vector<std::pair<size_t, size_t>> prt_idx_ranges
        (const typename GroupFields::coord_t grp_coord[GroupFields::dims[0]],
         const typename GroupFields::coord_t R) const;
};// }}}

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
void
Workspace<GroupFields,ParticleFields>::Sorting::compute_prt_indices ()
{// {{{
    prt_indices.reserve(Nprt);
    typename ParticleFields::coord_t *prt_coord
        = (typename ParticleFields::coord_t *)tmp_prt_properties[0];

    #define GRID(x, dir) ((size_t)(x[dir] / acell))

    for (size_t prt_idx=0; prt_idx != Nprt; ++prt_idx, prt_coord+=ParticleFields::dims[0])
        prt_indices.emplace_back(prt_idx, Ncells_side * Ncells_side * GRID(prt_coord, 0)
                                          +             Ncells_side * GRID(prt_coord, 1)
                                          +                           GRID(prt_coord, 2));

    #undef GRID
}// }}}

template<typename GroupFields, typename ParticleFields>
void
Workspace<GroupFields,ParticleFields>::Sorting::sort_prt_indices ()
{// {{{
    std::sort(prt_indices.begin(), prt_indices.end(),
              [](std::pair<size_t,size_t> a,
                 std::pair<size_t,size_t> b)
              { return a.second < b.second; } );
}// }}}

template<typename GroupFields, typename ParticleFields>
void
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
void
Workspace<GroupFields,ParticleFields>::Sorting::compute_offsets ()
{// {{{
    for (size_t ii=0; ii != Ncells_tot; ++ii)
    {
        offsets[ii] = Nprt;

        // find the index where the last cell started,
        // taking into account possibly empty cells
        // Note that the loop will not be executed at ii=0
        // and be careful with unsigned integer arithmetic
        size_t jj = 0UL;
        for (size_t kk=ii-1UL; kk < Ncells_tot; --kk)
            if (offsets[kk] < Nprt)
            {
                jj = offsets[kk];
                break;
            }

        for (; jj != Nprt; ++jj)
            if (prt_indices[jj].second == ii)
            {
                offsets[ii] = jj;
                break;
            }
    }

    offsets[Ncells_tot] = Nprt;
}// }}}

template<typename GroupFields, typename ParticleFields>
std::vector<std::pair<size_t, size_t>>
Workspace<GroupFields,ParticleFields>::Sorting::prt_idx_ranges
    (const typename GroupFields::coord_t grp_coord[GroupFields::dims[0]],
     typename GroupFields::coord_t R) const
{// {{{
    std::vector<std::pair<size_t, size_t>> out;

    typename GroupFields::coord_t grp_coord_normalized[GroupFields::dims[0]];
    for (size_t ii=0; ii != GroupFields::dims[0]; ++ii)
        grp_coord_normalized[ii] = grp_coord[ii] / acell;

    const typename GroupFields::coord_t R_normalized = R / acell;

    // loop over all cells -- there's a more efficient way of doing this,
    //                        by some preselection based on R,
    //                        but this is fast enough
    for (size_t ii=0; ii != Ncells_tot; ++ii)
    {
        typename GroupFields::coord_t cub_coord[]
            = { ii / (Ncells_side*Ncells_side),
                (ii/Ncells_side) % Ncells_side,
                ii % Ncells_side };

        if (// check if group has overlap with this cell
            Geometry::sph_cub_intersect(grp_coord_normalized, cub_coord,
                                        R_normalized,
                                        (typename GroupFields::coord_t)Bsize)
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

template<typename GroupFields, typename ParticleFields>
inline bool
Workspace<GroupFields,ParticleFields>::Sorting::Geometry::sph_cub_intersect
    (const typename GroupFields::coord_t grp_coord[GroupFields::dims[0]],
     typename GroupFields::coord_t cub_coord[GroupFields::dims[0]],
     typename GroupFields::coord_t grp_R,
     typename GroupFields::coord_t periodicity)
{// {{{
    mod_translations(grp_coord, cub_coord, periodicity);
    mod_reflections(cub_coord);

    return std::hypot(std::max((typename GroupFields::coord_t)0.0, cub_coord[0]),
                      std::max((typename GroupFields::coord_t)0.0, cub_coord[1]),
                      std::max((typename GroupFields::coord_t)0.0, cub_coord[2]) ) < grp_R;
}// }}}

template<typename GroupFields, typename ParticleFields>
inline void
Workspace<GroupFields,ParticleFields>::Sorting::Geometry::mod_translations
    (const typename GroupFields::coord_t grp_coord[GroupFields::dims[0]],
     typename GroupFields::coord_t cub_coord[GroupFields::dims[0]],
     typename GroupFields::coord_t periodicity)
{// {{{
    for (size_t ii=0; ii != GroupFields::dims[0]; ++ii)
    {
        cub_coord[ii] -= grp_coord[ii];
        if (cub_coord[ii] > 0.5 * periodicity)
            cub_coord[ii] -= periodicity;
        else if (cub_coord[ii] < -0.5 * periodicity)
            cub_coord[ii] += periodicity;
    }
}// }}}

template<typename GroupFields, typename ParticleFields>
inline void
Workspace<GroupFields,ParticleFields>::Sorting::Geometry::mod_reflections
    (typename GroupFields::coord_t cub_coord[GroupFields::dims[0]])
{// {{{
    for (size_t ii=0; ii != GroupFields::dims[0]; ++ii)
        if (cub_coord[ii] < -0.5)
            cub_coord[ii] = - ( cub_coord[ii] + 1.0 );
}// }}}

#endif // WORKSPACE_SORTING_HPP
