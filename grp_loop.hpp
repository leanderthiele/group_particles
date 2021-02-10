#ifndef GRP_LOOP_HPP
#define GRP_LOOP_HPP

#include <memory>
#include <cstring>

#include "H5Cpp.h"

#include "workspace.hpp"
#include "workspace_memory.hpp"
#include "read_fields.hpp"
#include "callback.hpp"

template<typename GroupFields, typename ParticleFields>
void
Workspace::grp_loop (Callback &callback)
{
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
        realloc_tmp_grp_storage(Ngrp_this_file);

        // read the file data
        read_fields<FieldTypes::GrpFld, GroupFields>(callback, fptr, Ngrp, tmp_grp_properties);

        // now loop over groups to see which ones belong into permanent storage
        for (size_t grp_idx=0; grp_idx != Ngrp_this_file; ++grp_idx)
        {
            for (size_t ii=0; ii != GroupFields::Nfields; ++ii)
                this_grp_properties[ii] = tmp_grp_properties[ii]
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
                    std::memcpy(grp_properties[ii] + Ngrp * GroupFields::strides[ii],
                                tmp_grp_properties[ii],
                                GroupFields::strides[ii]);
                
                // advance the counter
                ++Ngrp;
            }
        }
    }
}


#endif // GRP_LOOP_HPP
