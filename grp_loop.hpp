#ifndef GRP_LOOP_HPP
#define GRP_LOOP_HPP

#include <memory>
#include <cstring>
#include <cstdio>

#include "H5Cpp.h"

#include "workspace.hpp"
#include "workspace_memory.hpp"
#include "hdf5_utils.hpp"
#include "callback.hpp"

template<typename AFields>
void
Workspace<AFields>::grp_loop ()
{
    #ifndef NDEBUG
    std::fprintf(stderr, "Started Workspace::grp_loop ...\n");
    #endif // NDEBUG

    // the file name for the current chunk will be written here
    std::string fname;

    // loop until the callback function returns false
    for (size_t chunk_idx=0; callback.grp_chunk(chunk_idx, fname); ++chunk_idx)
    {
        // open the hdf5 file
        auto fptr = std::make_shared<H5::H5File>(fname, H5F_ACC_RDONLY);

        // read metadata
        size_t Ngrp_this_file;
        callback.read_grp_meta(chunk_idx, fptr, Ngrp_this_file);

        if (!Ngrp_this_file) continue;

        // allocate storage
        realloc_tmp_storage<typename AFields::GroupFields>(Ngrp_this_file, tmp_grp_properties);

        // read the file data
        hdf5Utils::read_fields<AFields, typename AFields::GroupFields>(callback, fptr, Ngrp_this_file, tmp_grp_properties);

        // file not needed anymore
        fptr->close();

        // convert coordinates to global type
        AFields::GroupFields::convert_coords(Ngrp_this_file, tmp_grp_properties[0]);

        typename Callback<AFields>::GrpProperties grp (tmp_grp_properties);

        // now loop over groups to see which ones belong into permanent storage
        for (size_t grp_idx=0; grp_idx != Ngrp_this_file; ++grp_idx, grp.advance())
            if (callback.grp_select(grp))
            {
                realloc_grp_storage_if_necessary ();

                // let the user do some stuff
                callback.grp_action(grp);
                
                // compute group radius
                grp_radii[Ngrp]    = callback.grp_radius(grp);
                grp_radii_sq[Ngrp] = grp_radii[Ngrp] * grp_radii[Ngrp];

                // copy properties into permanent storage
                for (size_t ii=0; ii != AFields::GroupFields::Nfields; ++ii)
                    std::memcpy((char *)(grp_properties[ii]) + Ngrp * AFields::GroupFields::strides_fcoord[ii],
                                grp[ii], AFields::GroupFields::strides_fcoord[ii]);
                
                // advance the counter
                ++Ngrp;
            }

        #ifndef NDEBUG
        std::fprintf(stderr, "In Workspace::grp_loop : did %lu chunks.\n", chunk_idx+1UL);
        #endif // NDEBUG
    }// for chunk_idx

    // save memory by shrinking the temporary buffers
    realloc_tmp_storage<typename AFields::GroupFields>(1, tmp_grp_properties);

    // save memory by reallocating the perhaps too large buffers
    shrink_grp_storage();

    #ifndef NDEBUG
    std::fprintf(stderr, "Ended Workspace::grp_loop, %lu groups loaded.\n", Ngrp);
    #endif // NDEBUG
}


#endif // GRP_LOOP_HPP
