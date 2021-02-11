#ifndef CALLBACK_HPP
#define CALLBACK_HPP

#include <string>
#include <memory>

#include "H5Cpp.h"

// base class for all user-defined call-back functions
// User should inherit from this class and override the
// member functions
struct Callback
{
    // these functions write the chunk file name corresponding
    // to chunk_idx into return argument and return true.
    // If no chunk corresponds to chunk_idx, false is returned
    virtual bool grp_chunk (size_t chunk_idx, std::string &fname) const = 0;
    virtual bool prt_chunk (size_t chunk_idx, std::string &fname) const = 0;

    // Where to find the group fields in the hdf5 file
    // (for Illustris: "Group/")
    virtual std::string grp_name () const = 0;

    // Where to find the particle fields in the hdf5 file
    // (for Illustris: "PartType<>/")
    virtual std::string prt_name () const = 0;

    // allows the user to read metadata for each group chunk
    // (thus this function is not marked as constant)
    // This function is required to write the number of groups in this chunk into the
    // return value.
    virtual void read_grp_meta (size_t chunk_idx, std::shared_ptr<H5::H5File> fptr,
                                size_t &Ngroups) = 0;

    // allows the user to read metadata for each particle chunk
    // (thus this function is not marked as constant)
    // This function is required to write the box size and number of particles in this
    // chunk into the return values.
    virtual void read_prt_meta (size_t chunk_idx, std::shared_ptr<H5::H5File> fptr,
                                float &Bsize, size_t &Npart) = 0;

    // returns whether the group described by the argument should be
    // considered
    virtual bool grp_select (void **grp_properties) const = 0;

    // this function will be called for all groups for which grp_select
    // returns true.
    // It will be called in order of groups encountered
    // (thus, it is NOT required to be thread safe!)
    // It is not marked as const since it probably stores some data in the child class.
    virtual void grp_action (void **grp_properties) = 0;

    // returns the group radius. Only particles falling within the radius
    // should be considered
    virtual float grp_radius (void **grp_properties) const = 0;

    // returns whether a particle should be considered for the halo passed.
    // Only particles falling within grp_radius will be passed here.
    virtual bool prt_select (void **grp_properties, void **prt_properties, float R) const = 0;

    // this function will be called for all particles for which prt_select
    // returns true.
    // It will not be called in order of particles encountered
    // and should thus be thread safe.
    // grp_idx corresponds to the order in which groups were passed to grp_action.
    // R is the radius from the group center (as defined by the user)
    // It is not marked as const since it probably stores some data in the child class.
    virtual void prt_action (size_t grp_idx, void **grp_properties, void **prt_properties, float R) = 0;
};

#endif // CALLBACK_HPP
