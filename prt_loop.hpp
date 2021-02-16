#ifndef PRT_LOOP_HPP
#define PRT_LOOP_HPP

#include <cmath>
#include <cstdio>

#include "callback.hpp"
#include "fields.hpp"
#include "read_fields.hpp"
#include "workspace.hpp"
#include "workspace_memory.hpp"
#include "workspace_sorting.hpp"
#include "geom_utils.hpp"
#include "timing.hpp"

template<typename AFields>
void
Workspace<AFields>::prt_loop ()
{// {{{
    #ifndef NDEBUG
    std::fprintf(stderr, "Started Workspace::prt_loop ...\n");
    #endif // NDEBUG

    // the file name for the current chunk will be written here
    std::string fname;

    // loop until the callback function returns false
    for (size_t chunk_idx=0; callback.prt_chunk(chunk_idx, fname); ++chunk_idx)
    {
        #ifndef NDEBUG
        TIME_PT(t1);
        #endif // NDEBUG

        // open the hdf5 file
        auto fptr = std::make_shared<H5::H5File>(fname, H5F_ACC_RDONLY);

        // read metadata
        size_t Nprt_this_file;
        float Bsize_this_file;
        callback.read_prt_meta(chunk_idx, fptr, Bsize_this_file, Nprt_this_file);

        if (chunk_idx==0)
            Bsize = Bsize_this_file;
        else
            assert(std::fabs(Bsize/Bsize_this_file - 1.0F) < 1e-5F);

        if (!Nprt_this_file) continue;

        // allocate storage
        #ifndef NDEBUG
        TIME_PT(t2);
        #endif // NDEBUG
        realloc_tmp_storage<typename AFields::ParticleFields>(Nprt_this_file, tmp_prt_properties);
        #ifndef NDEBUG
        TIME_MSG(t2, "prt_loop memory allocation for particle chunk data");
        #endif // NDEBUG

        // read the file data
        #ifndef NDEBUG
        TIME_PT(t3);
        #endif // NDEBUG
        read_fields<AFields, typename AFields::ParticleFields>(callback, fptr, Nprt_this_file,
                                                               tmp_prt_properties);
        #ifndef NDEBUG
        TIME_MSG(t3, "prt_loop read_fields for particle chunk data");
        #endif // NDEBUG

        // run the loop
        #ifndef NAIVE
        #   ifndef NDEBUG
        TIME_PT(t4);
        #   endif // NDEBUG
        prt_loop_sorted(Nprt_this_file);
        #   ifndef NDEBUG
        TIME_MSG(t4, "prt_loop->prt_loop_sorted");
        #   endif // NDEBUG
        #else // NAIVE
        #   ifndef NDEBUG
        TIME_PT(t4);
        #   endif // NDEBUG
        #   warning "Compiling with the naive particle loop instead of the (much faster) sorted one."
        prt_loop_naive(Nprt_this_file);
        #   ifndef NDEBUG
        TIME_MSG(t4, "prt_loop->prt_loop_naive");
        #   endif // NDEBUG
        #endif // NAIVE

        #ifndef NDEBUG
        TIME_MSG(t1, "chunk %lu in Workspace::prt_loop", chunk_idx+1UL);
        #endif // NDEBUG

        #ifndef NDEBUG
        std::fprintf(stderr, "In Workspace::prt_loop : did %lu chunks.\n", chunk_idx+1UL);
        #endif // NDEBUG
    }// for chunk_idx

    // save memory by shrinking the temporary particle storage
    realloc_tmp_storage<typename AFields::ParticleFields>(1, tmp_prt_properties);
}// }}}

template<typename AFields>
void
Workspace<AFields>::prt_loop_naive (size_t Nprt_this_file)
{// {{{
    typename Callback<AFields>::PrtProperties prt (tmp_prt_properties);

    // loop over particles
    for (size_t prt_idx=0; prt_idx != Nprt_this_file; ++prt_idx, prt.advance())
    {
        typename Callback<AFields>::GrpProperties grp (grp_properties);

        // now loop over groups
        for (size_t grp_idx=0; grp_idx != Ngrp; ++grp_idx, grp.advance())
            // do stuff
            prt_loop_inner(grp_idx, grp, prt);
    }// for prt_idx
}// }}}

template<typename AFields>
void
Workspace<AFields>::prt_loop_sorted (size_t Nprt_this_file)
{// {{{
    // create a Sorting instance, constructing it will perform the main work
    // associated with this object
    #ifndef NDEBUG
    TIME_PT(t1);
    #endif // NDEBUG

    Sorting prt_sort (Nprt_this_file, Bsize, tmp_prt_properties);

    #ifndef NDEBUG
    TIME_MSG(t1, "initialization of Sorting instance (Nprt=%lu)", Nprt_this_file);
    #endif // NDEBUG

    typename Callback<AFields>::GrpProperties grp (grp_properties);

    // loop over groups
    for (size_t grp_idx=0; grp_idx != Ngrp; ++grp_idx, grp.advance())
    {
        // compute which cells have intersection with this group
        std::vector<std::pair<size_t,size_t>> prt_idx_ranges
            = prt_sort.prt_idx_ranges(grp.coord(), grp_radii_sq[grp_idx]);

        // no particles in the vicinity of this group
        if (prt_idx_ranges.empty()) continue;
        
        // loop over cells
        for (auto &prt_idx_range : prt_idx_ranges)
        {
            typename Callback<AFields>::PrtProperties prt (prt_sort.tmp_prt_properties_sorted,
                                                           prt_idx_range.first);
            // loop over particles
            for (size_t prt_idx=prt_idx_range.first;
                        prt_idx != prt_idx_range.second;
                        ++prt_idx, prt.advance())
                prt_loop_inner(grp_idx, grp, prt);
        }// for prt_idx_range
    }// for grp_idx
}// }}}

template<typename AFields>
inline void
Workspace<AFields>::prt_loop_inner
    (size_t grp_idx,
     const typename Callback<AFields>::GrpProperties &grp,
     const typename Callback<AFields>::PrtProperties &prt)
{// {{{
    // figure out the distance between group and particle,
    float Rsq = GeomUtils::periodic_hypotsq<float>(grp.coord(), prt.coord(), Bsize);

    // check if this particle belongs to the group
    if (Rsq > grp_radii_sq[grp_idx]
        || !callback.prt_select(grp_idx, grp, prt, Rsq))
        return;

    // particle belongs to group: do the user-defined thing with it
    callback.prt_action(grp_idx, grp, prt, Rsq);
}// }}}

#endif // PRT_LOOP_HPP
