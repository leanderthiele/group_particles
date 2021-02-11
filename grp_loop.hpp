#ifndef GRP_LOOP_HPP
#define GRP_LOOP_HPP

#include <memory>
#include <cstring>
#include <cstdio>

#include "H5Cpp.h"

#include "workspace.hpp"
#include "workspace_memory.hpp"
#include "read_fields.hpp"
#include "callback.hpp"

template<typename GroupFields, typename ParticleFields>
void
Workspace<GroupFields,ParticleFields>::grp_loop (Callback &callback)
{
    #ifndef NDEBUG
    std::fprintf(stderr, "Started Workspace::grp_loop ...\n");
    #endif // NDEBUG

    // the file name for the current chunk will be written here
    std::string fname;

    // the number of groups in this file
    size_t Ngrp_this_file;

    // a temporary array to store pointers to group properties
    void *this_grp_properties[GroupFields::Nfields];

    // loop until the callback function returns false
    for (size_t chunk_idx=0; callback.grp_chunk(chunk_idx, fname); ++chunk_idx)
    {
        // open the hdf5 file
        auto fptr = std::make_shared<H5::H5File>(fname, H5F_ACC_RDONLY);

        // read metadata
        callback.read_grp_meta(chunk_idx, fptr, Ngrp_this_file);

        if (!Ngrp_this_file) continue;

        // allocate storage
        realloc_tmp_storage<GroupFields>(Ngrp_this_file, tmp_grp_properties);

        // read the file data
        read_fields<FieldTypes::GrpFld, GroupFields>(callback, fptr, Ngrp_this_file, tmp_grp_properties);

        // now loop over groups to see which ones belong into permanent storage
        for (size_t grp_idx=0; grp_idx != Ngrp_this_file; ++grp_idx)
        {
            for (size_t ii=0; ii != GroupFields::Nfields; ++ii)
                // we make use of the fact that a char is one byte wide
                this_grp_properties[ii] = (char *)(tmp_grp_properties[ii])
                                          + grp_idx * GroupFields::strides[ii];
            
            if (callback.grp_select(this_grp_properties))
            {
                realloc_grp_storage_if_necessary ();

                // let the user do some stuff
                callback.grp_action(this_grp_properties);
                
                // compute group radius
                grp_radii[Ngrp] = callback.grp_radius(tmp_grp_properties);

                // copy properties into permanent storage
                for (size_t ii=0; ii != GroupFields::Nfields; ++ii)
                    std::memcpy((char *)(grp_properties[ii]) + Ngrp * GroupFields::strides[ii],
                                tmp_grp_properties[ii],
                                GroupFields::strides[ii]);
                
                // advance the counter
                ++Ngrp;
            }
        }// for grp_idx

        #ifndef NDEBUG
        std::fprintf(stderr, "In Workspace::grp_loop : did %lu chunks.\n", chunk_idx+1UL);
        #endif // NDEBUG
    }// for chunk_idx

    // save memory by shrinking the temporary buffers
    realloc_tmp_storage<GroupFields>(1, tmp_grp_properties);

    // save memory by reallocating the perhaps too large buffers
    shrink_grp_storage();

    #ifndef NDEBUG
    std::fprintf(stderr, "Ended Workspace::grp_loop, %lu groups loaded.\n", Ngrp);
    #endif // NDEBUG
}


#endif // GRP_LOOP_HPP
