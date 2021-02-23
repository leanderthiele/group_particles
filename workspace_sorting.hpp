#ifndef WORKSPACE_SORTING_HPP
#define WORKSPACE_SORTING_HPP

#include <vector>
#include <utility>
#include <tuple>
#include <cstdlib>
#include <algorithm>
#include <cstring>
#include <cassert>
#include <cmath>

#include "fields.hpp"
#include "workspace.hpp"
#include "geom_utils.hpp"
#include "timing.hpp"

// TODO
// put parallel execution policy into std::sort
// (depending on preprocessor flag)

template<typename AFields>
class Workspace<AFields>::Sorting
{// {{{
    coord_t Bsize;
    const size_t Nprt;

    static constexpr const size_t Ncells_side = 8UL;
    static constexpr const size_t Ncells_tot  = Ncells_side * Ncells_side * Ncells_side;
    coord_t acell;
    
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
        static void mod_translations (const coord_t grp_coord[3], coord_t cub_coord[3]);
        static void mod_reflections (coord_t cub_coord[3]);
    public :
        // this function does not take the periodic boundary conditions into account
        static bool sph_cub_intersect (const coord_t grp_coord[3],
                                       coord_t cub_coord[3],
                                       coord_t grp_Rsq);
        Geometry () = delete;
    };

public :
    Sorting (size_t Nprt_, coord_t Bsize_, void **tmp_prt_properties_);
    Sorting () = delete;
    ~Sorting ();

    // store the sorted properties here (instance must allocate memory for this!)
    // user can access these
    void *tmp_prt_properties_sorted[AFields::ParticleFields::Nfields];

    // user can use this function to find the indices of all particles that may
    // belong to a given group
    std::vector<std::tuple<size_t, size_t, std::array<int,3>>> prt_idx_ranges
        (const coord_t grp_coord[3],
         const coord_t Rsq) const;

    // user can use this function to find the indices of all particles that may
    // belong to a given group
    // same functionality as above, but hopefully correct implementation
    std::vector<std::tuple<size_t, size_t, std::array<int,3>>> prt_idx_ranges
        (const coord_t grp_coord[3],
         const coord_t R, const coord_t Rsq) const;
};// }}}

// ----- Implementation -----

template<typename AFields>
Workspace<AFields>::Sorting::Sorting (size_t Nprt_,
                                      coord_t Bsize_,
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
        tmp_prt_properties_sorted[ii] = std::malloc(Nprt * AFields::ParticleFields::strides_fcoord[ii]);
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
    auto *prt_coord = (coord_t *)tmp_prt_properties[0];

    #define GRID(x, dir) ((size_t)(x[dir] / acell))

    for (size_t prt_idx=0; prt_idx != Nprt;
         ++prt_idx, prt_coord += 3)
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

        for (size_t prt_idx=0; prt_idx != Nprt; ++prt_idx, dest += AFields::ParticleFields::strides_fcoord[ii])
            std::memcpy(dest, (char *)(tmp_prt_properties[ii])
                              + prt_indices[prt_idx].first * AFields::ParticleFields::strides_fcoord[ii],
                        AFields::ParticleFields::strides_fcoord[ii]);
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
std::vector<std::tuple<size_t, size_t, std::array<int,3>>>
Workspace<AFields>::Sorting::prt_idx_ranges
    (const coord_t grp_coord[3], coord_t R, coord_t Rsq) const
{// {{{
    std::vector<std::tuple<size_t, size_t, std::array<int,3>>> out;

    coord_t grp_coord_normalized[3];
    for (size_t ii=0; ii != 3; ++ii)
        grp_coord_normalized[ii] = grp_coord[ii] / acell;

    const coord_t R_normalized = R / acell;
    const coord_t Rsq_normalized = Rsq / (acell*acell);

    for (int64_t  xx  = (int64_t)(grp_coord_normalized[0]-R_normalized) - 1L;
                  xx <= (int64_t)(grp_coord_normalized[0]+R_normalized);
                ++xx)
    {
        int64_t idx_x = Ncells_side * Ncells_side
                        * ((Ncells_side+xx%Ncells_side) % Ncells_side);

        for (int64_t  yy  = (int64_t)(grp_coord_normalized[1]-R_normalized) - 1L;
                      yy <= (int64_t)(grp_coord_normalized[1]+R_normalized);
                    ++yy)
        {
            int64_t idx_y = idx_x + Ncells_side 
                                    * ((Ncells_side+yy%Ncells_side) % Ncells_side);

            for (int64_t  zz  = (int64_t)(grp_coord_normalized[2]-R_normalized) - 1L;
                          zz <= (int64_t)(grp_coord_normalized[2]+R_normalized);
                        ++zz)
            {
                // this is the index in the flattened array of cells
                int64_t ii = idx_y + ((Ncells_side+zz%Ncells_side) % Ncells_side);
                
                coord_t cub_coord[] = { (coord_t)xx, (coord_t)yy, (coord_t)zz };

                if (Geometry::sph_cub_intersect(grp_coord_normalized, cub_coord, Rsq_normalized)
                    && offsets[ii] < Nprt)
                {
                    #define PER(x) ((x>=Ncells_side) ? 1 : (x<0) ? -1 : 0)
                    std::array<int,3> periodic_to_add { PER(xx), PER(yy), PER(zz) };
                    #undef PER

                    std::tuple<size_t, size_t, std::array<int,3>> idx_range
                        { offsets[ii], Nprt, periodic_to_add };

                    // find the last one -- taking possibly empty cells that follow into account!
                    for (size_t jj=(size_t)ii+1UL; jj != Ncells_tot; ++jj)
                        if (offsets[jj] < Nprt)
                        {
                            std::get<1>(idx_range) = offsets[jj];
                            break;
                        }
                    
                    out.push_back(idx_range);
                }
            }
        }
    }

    return out;
}// }}}

template<typename AFields>
inline bool
Workspace<AFields>::Sorting::Geometry::sph_cub_intersect
    (const coord_t grp_coord[3],
     coord_t cub_coord[3],
     coord_t grp_Rsq)
{
    mod_translations(grp_coord, cub_coord);
    mod_reflections(cub_coord);

    return GeomUtils::hypotsq(std::max((coord_t)0.0, cub_coord[0]),
                              std::max((coord_t)0.0, cub_coord[1]),
                              std::max((coord_t)0.0, cub_coord[2])
                             ) < grp_Rsq;
}

// no periodic boudary conditions!!!
template<typename AFields>
inline void
Workspace<AFields>::Sorting::Geometry::mod_translations
    (const coord_t grp_coord[3], coord_t cub_coord[3])
{// {{{
    for (size_t ii=0; ii != 3; ++ii)
        cub_coord[ii] -= grp_coord[ii];
}// }}}

template<typename AFields>
inline void
Workspace<AFields>::Sorting::Geometry::mod_reflections
    (coord_t cub_coord[3])
{// {{{
    for (size_t ii=0; ii != 3; ++ii)
        if (cub_coord[ii] < (coord_t)-0.5)
            cub_coord[ii] = - ( cub_coord[ii] + (coord_t)1.0 );
}// }}}

#endif // WORKSPACE_SORTING_HPP
