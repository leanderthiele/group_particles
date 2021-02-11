#ifndef WORKSPACE_HPP
#define WORKSPACE_HPP

#include "callback.hpp"
#include "fields.hpp"

template<typename GroupFields, typename ParticleFields>
class Workspace
{
    // where to find the user's functions
    Callback &callback;

    // metadata we need to store permanently
    float Bsize;

    // data we need to store permanently
    // (acoording to user-defined selection and radius calculation)
    size_t Ngrp = 0UL;
    size_t alloced_grp = 0UL;
    void *grp_properties[GroupFields::Nfields];
    float *grp_radii;

    // temporary buffers
    void *tmp_grp_properties[GroupFields::Nfields];
    void *tmp_prt_properties[ParticleFields::Nfields];

    void realloc_grp_storage (size_t new_size);

    template<typename Fields>
    void realloc_tmp_storage (size_t new_size, void **buf);

    // we can also use this function for the initial malloc
    void realloc_grp_storage_if_necessary ();

    void shrink_grp_storage ();

    // sets the pointers in dest to the locations in src indexed
    // by idx
    // T is one of GroupFields, ParticleFields
    template<typename T>
    void collect_properties (void **dest, void **src, size_t idx);

    // everything we need to sort particles
    class Sorting;

    // --- helper functions for the loops ---
    
    // the inner action, invariant under how we do the loops
    void prt_loop_inner (size_t grp_idx,
                         void **this_grp_properties,
                         void **this_prt_properties);
    
    // computes the distance between a particle and a group,
    // taking into account periodic boundary conditions
    float prt_grp_dist (typename GroupFields::coord_t *grp_coord,
                        typename ParticleFields::coord_t *prt_coord);

    // the simple loop over all particles
    void prt_loop_naive (size_t Nprt_this_file);

    // the more sophisticated loop grouping particles into cells
    // and considering only a subset for each group
    void prt_loop_sorted (size_t Nprt_this_file);

public :
    void grp_loop ();
    void prt_loop ();

    Workspace (Callback &callback_);

    ~Workspace ();

    // Note : we do not implement copy and move constructors,
    //        despite explicit memory management.
    //        The caller should simply not use these.
};

#endif // WORKSPACE_HPP
