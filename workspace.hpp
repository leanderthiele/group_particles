#ifndef WORKSPACE_HPP
#define WORKSPACE_HPP

#include "callback.hpp"
#include "fields.hpp"

template<typename AFields>
class Workspace
{
    // where to find the user's functions
    Callback<AFields> &callback;

    // metadata we need to store permanently
    float Bsize;

    // data we need to store permanently
    // (acoording to user-defined selection and radius calculation)
    size_t Ngrp = 0UL;
    size_t alloced_grp = 0UL;
    void *grp_properties[AFields::GroupFields::Nfields];
    float *grp_radii;

    // temporary buffers
    void *tmp_grp_properties[AFields::GroupFields::Nfields];
    void *tmp_prt_properties[AFields::ParticleFields::Nfields];

    void realloc_grp_storage (size_t new_size);

    // T is one of GroupFields, ParticleFields
    template<typename T>
    void realloc_tmp_storage (size_t new_size, void **buf);

    // we can also use this function for the initial malloc
    void realloc_grp_storage_if_necessary ();

    void shrink_grp_storage ();

    // everything we need to sort particles
    class Sorting;

    // --- helper functions for the loops ---
    
    // the inner action, invariant under how we do the loops
    void prt_loop_inner (size_t grp_idx,
                         const typename Callback<AFields>::GrpProperties &grp,
                         const typename Callback<AFields>::PrtProperties &prt);
    
    // computes the distance between a particle and a group,
    // taking into account periodic boundary conditions
    float prt_grp_dist (typename AFields::GroupFields::coord_t *grp_coord,
                        typename AFields::ParticleFields::coord_t *prt_coord);

    // the simple loop over all particles
    void prt_loop_naive (size_t Nprt_this_file);

    // the more sophisticated loop grouping particles into cells
    // and considering only a subset for each group
    void prt_loop_sorted (size_t Nprt_this_file);

public :
    void grp_loop ();
    void prt_loop ();

    Workspace (Callback<AFields> &callback_);

    ~Workspace ();

    // Note : we do not implement copy and move constructors,
    //        despite explicit memory management.
    //        The caller should simply not use these.
};

#endif // WORKSPACE_HPP
