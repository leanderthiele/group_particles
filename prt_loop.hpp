#ifndef PRT_LOOP_HPP
#define PRT_LOOP_HPP

#include <cmath>

#include "callback.hpp"
#include "fields.hpp"
#include "read_fields.hpp"
#include "workspace.hpp"
#include "workspace_memory.hpp"

template<typename GroupFields, typename ParticleFields>
void
Workspace<GroupFields,ParticleFields>::prt_loop (Callback &callback)
{
    // the file name for the current chunk will be written here
    std::string fname;

    // the number of particles in this file
    size_t Nprt_this_file;
    float Bsize_this_file;

    // temporary arrays to store pointers to properties
    void *this_prt_properties[ParticleFields::Nfields];
    void *this_grp_properties[GroupFields::Nfields];

    // for convenience
    float *this_prt_coords;
    float *this_grp_coords;

    // loop until the callback function returns false
    for (size_t chunk_idx=0; callback.prt_chunk(chunk_idx, fname); ++chunk_idx)
    {
        // open the hdf5 file
        auto fptr = std::make_shared<H5::H5File>(fname, H5F_ACC_RDONLY);

        // read metadata
        callback.read_prt_meta(chunk_idx, fptr, Bsize_this_file, Nprt_this_file);

        if (chunk_idx==0)
            Bsize == Bsize_this_file;
        else
            assert(std::fabs(Bsize/Bsize_this_file - 1.0F) < 1e-5F);

        if (!Nprt_this_file) continue;

        // allocate storage
        realloc_tmp_storage<ParticleFields>(Nprt_this_file, tmp_prt_properties);

        // read the file data
        read_fields<FieldTypes::PrtFld, ParticleFields>(callback, fptr, Nprt_this_file, tmp_prt_properties);

        // now loop over particles (TODO : if we have Paco's sorting thing implemented, reverse loops)
        for (size_t prt_idx=0; prt_idx != Nprt_this_file; ++prt_idx)
        {
            for (size_t ii=0; ii != ParticleFields::Nfields; ++ii)
                this_prt_properties[ii] = tmp_prt_properties[ii]
                                          + prt_idx * ParticleFields::strides[ii];
            this_prt_coords = tmp_prt_properties[0]
                              + prt_idx * ParticleFields::strides[0];

            // now loop over groups
            for (size_t grp_idx=0; grp_idx != Ngrp; ++grp_idx)
            {
                for (size_t ii=0; ii != GroupFields::Nfields; ++ii)
                    this_grp_properties[ii] = grp_properties[ii]
                                              + grp_idx * GroupFields::strides[ii];
                this_grp_coords = grp_properties[0]
                                  + grp_idx * GroupFields::strides[0];

                // figure out the distance between group and particle,
                // remembering periodic boundary conditions
                #define PERIODIC (dist)                   \
                    (dist > 0.5F * Bsize) ? dist-Bsize    \
                    : (dist < -0.5F * Bsize) ? dist+Bsize \
                    : dist
                #define DIST (dir)                              \
                    this_prt_coords[dir] - this_grp_coords[dir]
                float R = std::hypot(PERIODIC(DIST(0)),
                                     PERIODIC(DIST(1)),
                                     PERIODIC(DIST(2)) );
                #undef PERIODIC
                #undef DIST

                // check if this particle belongs to the group
                if (R > grp_radii[grp_idx]
                    || !callback.prt_select(this_grp_properties,
                                            this_prt_properties))
                    continue;

                // particle belongs to group: do the user-defined thing with it
                callback.prt_action(grp_idx, this_grp_properties, this_prt_properties);
            }// for grp_idx
        }// for prt_idx
    }// for chunk_idx

    // save memory by shrinking the temporary particle storage
    realloc_tmp_storage<ParticleFields>(1, tmp_prt_properties);
}

#endif // PRT_LOOP_HPP
